#include "pcAdapter.h"
#include <MeshSimAdapt.h>
#include <SimUtil.h>
#include <SimPartitionedMesh.h>
#include <cassert>

namespace pc {
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
      printf("Start mesh adapt\n");
      pMSAdapt adapter = MSA_new(sim_pm, 1);
      MSA_setAdaptBL(adapter, 0);
      MSA_setExposedBLBehavior(adapter,BL_DisallowExposed);

      /* use size field before mesh motion */
      printf("Start mesh adapt of setting size field\n");
      vIter = M_vertexIter(pm);
      while((meshVertex = VIter_next(vIter))){
        MSA_scaleVertexSize(adapter, meshVertex, 1.0); // use the size field of the mesh before mesh motion
//        MSA_scaleVertexSize(adapter, meshVertex, 2.0); // double size field
//        MSA_scaleVertexSize(adapter, meshVertex, 0.5); // half size field
      }
      VIter_delete(vIter);

      /* copy the size field from APF to the Simmetrix adapter */
//      apf::MeshEntity* v;
//      apf::MeshIterator* it = m->begin(0);
//      while ((v = m->iterate(it))) {
//        double size = apf::getScalar(szFld, v, 0);
//        MSA_setVertexSize(adapter, (pVertex) v, size);
//      }
//      m->end(it);
//      apf::destroyField(szFld);

      /* unpacked solution into serveral fields */
      apf::Field* pre_field = chef::extractField(m,"solution","pressure",1,apf::SCALAR,in.simmetrixMesh);
      apf::Field* vel_field = chef::extractField(m,"solution","velocity",2,apf::VECTOR,in.simmetrixMesh);
      apf::Field* tem_field = chef::extractField(m,"solution","temperature",5,apf::SCALAR,in.simmetrixMesh);

      /* put these field explicitly into pPList */
      int num_flds = 3; // for now
      pField* sim_flds = new pField[num_flds];
      sim_flds[0] = apf::getSIMField(pre_field);
      sim_flds[1] = apf::getSIMField(vel_field);
      sim_flds[2] = apf::getSIMField(tem_field);
      pPList sim_fld_lst = PList_new();
      for (int i = 0; i < num_flds; i++)
        PList_append(sim_fld_lst, sim_flds[i]);
      assert(num_flds == PList_size(sim_fld_lst));

      /* set fields to be mapped */
      MSA_setMapFields(adapter, sim_fld_lst);
      PList_delete(sim_fld_lst);

//      PM_write(sim_pm, "before_adapt.sms", NULL);

      /* run the adapter */
      printf("do real mesh adapt\n");
      pProgress progress = Progress_new();
      MSA_adapt(adapter, progress);
      MSA_delete(adapter);

      // write model and mesh
      printf("write mesh for mesh adapt\n");
      M_write(pm, "after_adapt.sms", 0, progress);
      Progress_delete(progress);
      /* transfer data back to apf */
      apf::Field* solution = chef::combineField(m,"solution","pressure","velocity","temperature");
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
