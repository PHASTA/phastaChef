#include "pcAdapter.h"
#include <MeshSimAdapt.h>
#include <SimUtil.h>
#include <SimPartitionedMesh.h>
#include <SimDiscrete.h>
#include "apfSIM.h"
#include "gmi_sim.h"
#include <PCU.h>
#include <cassert>

extern void MSA_setBLSnapping(pMSAdapt, int onoff);

namespace pc {

  apf::Field* convertField(apf::Mesh* m,
    const char* inFieldname,
    const char* outFieldname) {
    apf::Field* inf = m->findField(inFieldname);
    assert(inf);
    int size = apf::countComponents(inf);
    apf::Field* outf = m->findField(outFieldname);
    if (outf)
      apf::destroyField(outf);
    outf = apf::createPackedField(m, outFieldname, size);
    double* inVal = new double[size];
    double* outVal = new double[size];
    apf::MeshEntity* vtx;
    apf::MeshIterator* it = m->begin(0);
    while ((vtx = m->iterate(it))) {
      apf::getComponents(inf, vtx, 0, inVal);
      for (int i = 0; i < size; i++){
        outVal[i] = inVal[i];
      }
      apf::setComponents(outf,vtx, 0, outVal);
    }
    m->end(it);
    apf::destroyField(inf);
    return outf;
  }

  int getSimFields(apf::Mesh2*& m, int simFlag, pField* sim_flds) {
    int num_flds = 0;
    if (m->findField("solution")) {
      num_flds += 3;
      sim_flds[0] = apf::getSIMField(chef::extractField(m,"solution","pressure",1,apf::SCALAR,simFlag));
      sim_flds[1] = apf::getSIMField(chef::extractField(m,"solution","velocity",2,apf::VECTOR,simFlag));
      sim_flds[2] = apf::getSIMField(chef::extractField(m,"solution","temperature",5,apf::SCALAR,simFlag));
      apf::destroyField(m->findField("solution"));
    }

    if (m->findField("time derivative of solution")) {
      num_flds += 3;
      sim_flds[3] = apf::getSIMField(chef::extractField(m,"time derivative of solution","der_pressure",1,apf::SCALAR,simFlag));
      sim_flds[4] = apf::getSIMField(chef::extractField(m,"time derivative of solution","der_velocity",2,apf::VECTOR,simFlag));
      sim_flds[5] = apf::getSIMField(chef::extractField(m,"time derivative of solution","der_temperature",5,apf::SCALAR,simFlag));
      apf::destroyField(m->findField("time derivative of solution"));
    }

    if (m->findField("mesh_vel")) {
      num_flds += 1;
      sim_flds[6] = apf::getSIMField(chef::extractField(m,"mesh_vel","mesh_vel_sim",1,apf::VECTOR,simFlag));
      apf::destroyField(m->findField("mesh_vel"));
    }
    return num_flds;
  }

  /* unpacked solution into serveral fields,
     put these field explicitly into pPList */
  pPList getSimFieldList(ph::Input& in, apf::Mesh2*& m){
    pField* sim_flds = new pField[7]; // Hardcoding
    int num_flds = getSimFields(m, in.simmetrixMesh, sim_flds);
    assert(num_flds == PList_size(sim_fld_lst));
    pPList sim_fld_lst = PList_new();
    for (int i = 0; i < num_flds; i++) {
      PList_append(sim_fld_lst, sim_flds[i]);
    }
    return sim_fld_lst;
  }

  void transferSimFields(apf::Mesh2*& m) {
    if (m->findField("pressure")) // assume we had solution before
      chef::combineField(m,"solution","pressure","velocity","temperature");
    if (m->findField("der_pressure")) // assume we had time derivative of solution before
      chef::combineField(m,"time derivative of solution","der_pressure","der_velocity","der_temperature");
    if (m->findField("mesh_vel_sim"))
      convertField(m, "mesh_vel_sim", "mesh_vel");
  }

  void runMeshAdapter(ph::Input& in, apf::Mesh2*& m, apf::Field*& orgSF, int step) {
    if (m->findField("material_type"))
      apf::destroyField(m->findField("material_type"));
    if (m->findField("meshQ"))
      apf::destroyField(m->findField("meshQ"));
    in.writeGeomBCFiles = 1; //write GeomBC file for visualization

    /* use the size field of the mesh before mesh motion */
    apf::Field* szFld = orgSF;

    if(in.simmetrixMesh == 1) {
      apf::MeshSIM* sim_m = dynamic_cast<apf::MeshSIM*>(m);
      pParMesh sim_pm = sim_m->getMesh();
      pMesh pm = PM_mesh(sim_pm,0);

      // declaration
      VIter vIter;
      pVertex meshVertex;

      /* create the Simmetrix adapter */
      if(!PCU_Comm_Self())
        printf("Start mesh adapt\n");
      pMSAdapt adapter = MSA_new(sim_pm, 1);
      MSA_setAdaptBL(adapter, 1);
      MSA_setExposedBLBehavior(adapter,BL_DisallowExposed);
//      MSA_setNoMigration(adapter,1); // hack; since split/cut mesh is not supported with adaptation
//      MSA_setBLSnapping(adapter, 0); // currently needed for parametric model
      MSA_setBLMinLayerAspectRatio(adapter, 0.0); // needed in parallel

      /* use size field before mesh motion */
//      double* inVal = new double[apf::countComponents(szFld)];
      if(!PCU_Comm_Self())
        printf("Start mesh adapt of setting size field\n");
      vIter = M_vertexIter(pm);
      while((meshVertex = VIter_next(vIter))){
//        apf::getComponents(szFld, reinterpret_cast<apf::MeshEntity*>(meshVertex), 0, inVal);
//        MSA_setVertexSize(adapter, meshVertex, inVal[0]); // use the prescribed size field
        MSA_scaleVertexSize(adapter, meshVertex, 1.0); // use the size field of the mesh before mesh motion
      }
      VIter_delete(vIter);

      /* set fields to be mapped */
      if (in.solutionMigration) {
        pPList sim_fld_lst = getSimFieldList(in, m);
        MSA_setMapFields(adapter, sim_fld_lst);
        PList_delete(sim_fld_lst);
      }

      /* run the adapter */
      if(!PCU_Comm_Self())
        printf("do real mesh adapt\n");
      pProgress progress = Progress_new();
      MSA_adapt(adapter, progress);
      MSA_delete(adapter);

      // write mesh
      if(!PCU_Comm_Self())
        printf("write mesh for mesh adapt\n");
      PM_write(sim_pm, "adapted_mesh.sms", progress);
      Progress_delete(progress);

      /* transfer data back to apf */
      if (in.solutionMigration) {
        transferSimFields(m);
      }
    }
    else {
      assert(szFld);
      apf::synchronize(szFld);
      apf::synchronize(m->getCoordinateField());
      /* do SCOREC mesh adaptation */
      chef::adapt(m,szFld,in);
    }
    m->verify();
    chef::balance(in,m);
  }

}
