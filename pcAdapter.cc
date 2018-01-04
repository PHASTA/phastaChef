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

  /* unpacked solution into serveral fields,
     put these field explicitly into pPList */
  pPList getSimFieldList(ph::Input& in, apf::Mesh2*& m){
    int num_flds = 0;
    pField* sim_flds = new pField[7]; // Hardcoding
    pPList sim_fld_lst = PList_new();
    if (m->findField("solution")) {
      num_flds += 3;
      sim_flds[0] = apf::getSIMField(chef::extractField(m,"solution","pressure",1,apf::SCALAR,in.simmetrixMesh));
      sim_flds[1] = apf::getSIMField(chef::extractField(m,"solution","velocity",2,apf::VECTOR,in.simmetrixMesh));
      sim_flds[2] = apf::getSIMField(chef::extractField(m,"solution","temperature",5,apf::SCALAR,in.simmetrixMesh));
      for (int i = 0; i < 3; i++)
        PList_append(sim_fld_lst, sim_flds[i]);
    }

    if (m->findField("time derivative of solution")) {
      num_flds += 3;
      sim_flds[3] = apf::getSIMField(chef::extractField(m,"time derivative of solution","der_pressure",1,apf::SCALAR,in.simmetrixMesh));
      sim_flds[4] = apf::getSIMField(chef::extractField(m,"time derivative of solution","der_velocity",2,apf::VECTOR,in.simmetrixMesh));
      sim_flds[5] = apf::getSIMField(chef::extractField(m,"time derivative of solution","der_temperature",5,apf::SCALAR,in.simmetrixMesh));
      for (int i = 3; i < 6; i++)
        PList_append(sim_fld_lst, sim_flds[i]);
    }

    if (m->findField("mesh_vel")) {
      num_flds += 1;
      sim_flds[6] = apf::getSIMField(chef::extractField(m,"mesh_vel","mesh_vel_sim",1,apf::VECTOR,in.simmetrixMesh));
      PList_append(sim_fld_lst, sim_flds[6]);
    }

    assert(num_flds == PList_size(sim_fld_lst));
    return sim_fld_lst;
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
      if(!PCU_Comm_Self())
        printf("Start mesh adapt of setting size field\n");
      vIter = M_vertexIter(pm);
      while((meshVertex = VIter_next(vIter))){
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
        if (m->findField("pressure")) // assume we had solution before
          chef::combineField(m,"solution","pressure","velocity","temperature");
        if (m->findField("der_pressure")) // assume we had time derivative of solution before
          chef::combineField(m,"time derivative of solution","der_pressure","der_velocity","der_temperature");
        if (m->findField("mesh_vel_sim"))
          convertField(m, "mesh_vel_sim", "mesh_vel");
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
