#include "pcUpdateMesh.h"
#include <SimPartitionedMesh.h>
#include "SimAdvMeshing.h"
#include "SimModel.h"
#include "SimUtil.h"
#include "SimParasolidKrnl.h"
#include "SimMeshMove.h"
#include "SimMeshTools.h"
#include "MeshSim.h"
#include "MeshSimAdapt.h"
#include "apfSIM.h"
#include "gmi_sim.h"
#include <cassert>

namespace pc {
  bool updateAPFCoord(apf::Mesh2* m) {
    apf::Field* f = m->findField("motion_coords");
    assert(f);
    double* vals = new double[apf::countComponents(f)];
    assert(apf::countComponents(f) == 3);
    apf::MeshEntity* vtx;
    apf::Vector3 points;
    apf::MeshIterator* itr = m->begin(0);
    while( (vtx = m->iterate(itr)) ) {
      apf::getComponents(f, vtx, 0, vals);
      for ( int i = 0; i < 3; i++ )  points[i] = vals[i];
      m->setPoint(vtx, 0, points);
    }
    m->end(itr);
    delete [] vals;
    return true;
  }

  bool updateSIMCoord(apf::Mesh2* m, int step, int caseId) {
    MS_init();
    SimAdvMeshing_start();
    SimParasolid_start(1);
    SimMeshTools_start();
    pProgress progress = Progress_new();

    apf::MeshSIM* apf_msim = dynamic_cast<apf::MeshSIM*>(m);
    pParMesh ppm = apf_msim->getMesh();
    pMesh pm = PM_mesh(ppm,0);

    M_write(pm, "before_overwrite.sms", 0, progress);

    gmi_model* gmiModel = apf_msim->getModel();
    pGModel model = gmi_export_sim(gmiModel);
//    pNativeModel nmodel = ParasolidNM_createFromFile("model_nat.x_t", 0);
//    pGModel model = GM_load("model.smd", nmodel, progress);

    apf::Field* f = m->findField("motion_coords");
    assert(f);
    double* vals = new double[apf::countComponents(f)];
    assert(apf::countComponents(f) == 3);

    // case specific variables
    const double trans[3] = {2e-4, 0.0, 0.0};
    double disp, sfct, rang;
    if (caseId == 1) {
      disp = 2e-4 * (int)(step / 2 - 1);
      printf("current step is %d; disp = %f\n",step,disp);
      sfct = 0.9025;
      rang = 2.0;
    }
    else if (caseId == 2) {
      disp = 4e-4;
      sfct = 0.9025;
      rang = 2.0;
    }
    else {
      printf("wrong case id\n");
      assert(0);
    }

    // prepare variables
    int gas = 1;        // hardcoding
    int grain11 = 1339; // hardcoding
    int grain12 = 1036; // hardcoding
    int grain13 = 733;  // hardcoding
    int grain21 = 1440; // hardcoding
    int grain22 = 1137; // hardcoding
    int grain23 = 834;  // hardcoding
    int grain31 = 1238; // hardcoding
    int grain32 = 935;  // hardcoding
    int grain33 = 632;  // hardcoding
    int faceC1  = 54;   // hardcoding
    int faceC2  = 46;   // hardcoding
    int faceC3  = 41;   // hardcoding
    const double rotax[3] = {0.0, 0.0, 1.0};
    const double rotax2[3] = {0.0, 0.0, -1.0};
    const double rotpt11[3] = {0.5e-3 + disp,   1.125e-3, 0.0};
    const double rotpt12[3] = {2.625e-3 + disp, 1.125e-3, 0.0};
    const double rotpt13[3] = {4.75e-3 + disp,  1.125e-3, 0.0};
    const double rotpt21[3] = {0.5e-3 + disp,   0.0,      0.0};
    const double rotpt22[3] = {2.625e-3 + disp, 0.0,      0.0};
    const double rotpt23[3] = {4.75e-3 + disp,  0.0,      0.0};
    const double rotpt31[3] = {0.5e-3 + disp,  -1.125e-3, 0.0};
    const double rotpt32[3] = {2.625e-3 + disp,-1.125e-3, 0.0};
    const double rotpt33[3] = {4.75e-3 + disp, -1.125e-3, 0.0};

    // declaration
    GRIter grIter;
    pGRegion modelRegion;
    GFIter gfIter;
    pGFace modelFace;
    VIter vIter;
    pVertex meshVertex;
    double newpt[3];
    double newpar[2];
    double xyz[3];

    // start mesh mover
    printf("Start mesh mover\n");
    pMeshMover mmover = MeshMover_new(pm, 0);

    // mesh motion of moving body
    grIter = GM_regionIter(model);
    while((modelRegion = GRIter_next(grIter))){
      if (GEN_tag(modelRegion) == grain11) {
        printf("Start mesh mover on grain11\n");
        MeshMover_setTransform(mmover, modelRegion, trans, rotax, rotpt11, rang, sfct);
      }
      if (GEN_tag(modelRegion) == grain12) {
        printf("Start mesh mover on grain12\n");
        MeshMover_setTransform(mmover, modelRegion, trans, rotax, rotpt12, rang, sfct);
      }
      if (GEN_tag(modelRegion) == grain13) {
        printf("Start mesh mover on grain13\n");
        MeshMover_setTransform(mmover, modelRegion, trans, rotax, rotpt13, rang, sfct);
      }
      if (GEN_tag(modelRegion) == grain21) {
        printf("Start mesh mover on grain21\n");
        MeshMover_setTransform(mmover, modelRegion, trans, rotax, rotpt21, 0.0, sfct);
      }
      if (GEN_tag(modelRegion) == grain22) {
        printf("Start mesh mover on grain22\n");
        MeshMover_setTransform(mmover, modelRegion, trans, rotax, rotpt22, 0.0, sfct);
      }
      if (GEN_tag(modelRegion) == grain23) {
        printf("Start mesh mover on grain23\n");
        MeshMover_setTransform(mmover, modelRegion, trans, rotax, rotpt23, 0.0, sfct);
      }
      if (GEN_tag(modelRegion) == grain31) {
        printf("Start mesh mover on grain31\n");
        MeshMover_setTransform(mmover, modelRegion, trans, rotax2, rotpt31, rang, sfct);
      }
      if (GEN_tag(modelRegion) == grain32) {
        printf("Start mesh mover on grain32\n");
        MeshMover_setTransform(mmover, modelRegion, trans, rotax2, rotpt32, rang, sfct);
      }
      if (GEN_tag(modelRegion) == grain33) {
        printf("Start mesh mover on grain33\n");
        MeshMover_setTransform(mmover, modelRegion, trans, rotax2, rotpt33, rang, sfct);
      }
    }
    GRIter_delete(grIter);

    // mesh motion of vertices on surfaces
    printf("Start mesh mover on surfaces and edges\n");
    gfIter = GM_faceIter(model);
    while((modelFace = GFIter_next(gfIter))){
      if(GEN_tag(modelFace) == faceC1 ||
         GEN_tag(modelFace) == faceC2 ||
         GEN_tag(modelFace) == faceC3){
        vIter = M_classifiedVertexIter(pm, modelFace, 1);
        while((meshVertex = VIter_next(vIter))){
          V_coord(meshVertex, xyz);
          apf::MeshEntity* vtx = reinterpret_cast<apf::MeshEntity*>(meshVertex);
          apf::getComponents(f, vtx, 0, vals);
          const double disp[3] = {vals[0]-xyz[0], vals[1]-xyz[1], vals[2]-xyz[2]};
          V_movedParamPoint(meshVertex,disp,newpar,newpt);
          MeshMover_setSurfaceMove(mmover,meshVertex,newpar,newpt);
        }
        VIter_delete(vIter);
      }
    }
    GFIter_delete(gfIter);

    // mesh motion of vertices in region
    printf("Start mesh mover in region\n");
    grIter = GM_regionIter(model);
    while((modelRegion = GRIter_next(grIter))){
      if(GEN_tag(modelRegion) == gas){
        vIter = M_classifiedVertexIter(pm, modelRegion, 0);
        while((meshVertex = VIter_next(vIter))){
          if (EN_isBLEntity(meshVertex)) continue;
          apf::MeshEntity* vtx = reinterpret_cast<apf::MeshEntity*>(meshVertex);
          apf::getComponents(f, vtx, 0, vals);
//          V_coord(meshVertex, xyz);
//          printf("old position: (%f,%f,%f)\n",xyz[0],xyz[1],xyz[2]);
//          printf("new position: (%f,%f,%f)\n",vals[0],vals[1],vals[2]);
          const double newloc[3] = {vals[0], vals[1], vals[2]};
          MeshMover_setVolumeMove(mmover,meshVertex,newloc);
        }
        VIter_delete(vIter);
      }
    }
    GRIter_delete(grIter);

    // do real work
   printf("do real mesh mover\n");
   assert(MeshMover_run(mmover, progress));
   MeshMover_delete(mmover);

    // write model and mesh
    printf("write model for mesh mover\n");
    GM_write(model, "after_mover.smd", 0, progress);
    printf("write mesh for mesh mover\n");
    M_write(pm, "after_mover.sms", 0, progress);

    Progress_delete(progress);
    SimMeshTools_stop();
    SimParasolid_stop(1);
    SimAdvMeshing_stop();
    MS_exit();
    return true;
  }

  void updateMeshCoord(ph::Input& in, apf::Mesh2* m, int step, int caseId) {
    bool done = false;
    if (in.simmetrixMesh)
      done = updateSIMCoord(m,step,caseId);
    else
      done = updateAPFCoord(m);
    assert(done);
  }

}
