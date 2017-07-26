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
#include <list>

namespace pc {
  struct movingBodyMotion {
    int tag;
    double trans[3];
    double rotaxis[3];
    double rotpt[3];
    double rotang;
    double scale;
  };

  struct meshMotion {
    std::list<movingBodyMotion> movingBodyMotions;
    std::list<int> surfaceTags;
    std::list<int> edgeTags;
    std::list<int> regionTags;
  };

  /* ideally, this comes from some input file */
  meshMotion configureMotion(int caseId, int step) {
    double disp, sfct, rang;
    sfct = 0.9025;
    rang = 2.0;
    if (caseId == 1) {
      disp = 2e-4 * (int)(step / 2 - 1);
      printf("current step is %d; disp = %f\n",step,disp);
    }
    else if (caseId == 2) {
      disp = 4e-4;
    }
    else {
      printf("wrong case id\n");
      assert(0);
    }
    meshMotion mm;
    movingBodyMotion mbm;
    mm.surfaceTags.push_back(54);
    mm.surfaceTags.push_back(46);
    mm.surfaceTags.push_back(41);
    mm.regionTags.push_back(1);
// grain11
    mbm.tag = 1339;
    mbm.trans[0] = 2e-4;
    mbm.rotaxis[2] = 1.0;
    mbm.rotpt[0] = 0.5e-3 + disp;
    mbm.rotpt[1] = 1.125e-3;
    mbm.rotang = rang;
    mbm.scale = sfct;
    mm.movingBodyMotions.push_back(mbm);
// grain12
    mbm.tag = 1036;
    mbm.trans[0] = 2e-4;
    mbm.rotaxis[2] = 1.0;
    mbm.rotpt[0] = 2.625e-3 + disp;
    mbm.rotpt[1] = 1.125e-3;
    mbm.rotang = rang;
    mbm.scale = sfct;
    mm.movingBodyMotions.push_back(mbm);
// grain13
    mbm.tag = 733;
    mbm.trans[0] = 2e-4;
    mbm.rotaxis[2] = 1.0;
    mbm.rotpt[0] = 4.75e-3 + disp;
    mbm.rotpt[1] = 1.125e-3;
    mbm.rotang = rang;
    mbm.scale = sfct;
    mm.movingBodyMotions.push_back(mbm);
// grain21
    mbm.tag = 1440;
    mbm.trans[0] = 2e-4;
    mbm.rotaxis[2] = 1.0;
    mbm.rotpt[0] = 0.5e-3 + disp;
    mbm.rotang = 0.0;
    mbm.scale = sfct;
    mm.movingBodyMotions.push_back(mbm);
// grain22
    mbm.tag = 1137;
    mbm.trans[0] = 2e-4;
    mbm.rotaxis[2] = 1.0;
    mbm.rotpt[0] = 2.625e-3 + disp;
    mbm.rotang = 0.0;
    mbm.scale = sfct;
    mm.movingBodyMotions.push_back(mbm);
// grain23
    mbm.tag = 834;
    mbm.trans[0] = 2e-4;
    mbm.rotaxis[2] = 1.0;
    mbm.rotpt[0] = 4.75e-3 + disp;
    mbm.rotang = 0.0;
    mbm.scale = sfct;
    mm.movingBodyMotions.push_back(mbm);
// grain31
    mbm.tag = 1238;
    mbm.trans[0] = 2e-4;
    mbm.rotaxis[2] = -1.0;
    mbm.rotpt[0] = 0.5e-3 + disp;
    mbm.rotpt[1] = -1.125e-3;
    mbm.rotang = rang;
    mbm.scale = sfct;
    mm.movingBodyMotions.push_back(mbm);
// grain32
    mbm.tag = 935;
    mbm.trans[0] = 2e-4;
    mbm.rotaxis[2] = -1.0;
    mbm.rotpt[0] = 2.625e-3 + disp;
    mbm.rotpt[1] = -1.125e-3;
    mbm.rotang = rang;
    mbm.scale = sfct;
    mm.movingBodyMotions.push_back(mbm);
// grain33
    mbm.tag = 632;
    mbm.trans[0] = 2e-4;
    mbm.rotaxis[2] = -1.0;
    mbm.rotpt[0] = 4.75e-3 + disp;
    mbm.rotpt[1] = -1.125e-3;
    mbm.rotang = rang;
    mbm.scale = sfct;
    mm.movingBodyMotions.push_back(mbm);
// return
    return mm;
  }

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

  bool updateSIMDiscreteCoord(apf::Mesh2* m) {
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

    apf::Field* f = m->findField("motion_coords");
    assert(f);
    double* vals = new double[apf::countComponents(f)];
    assert(apf::countComponents(f) == 3);

    // declaration
    pGRegion modelRegion;
    VIter vIter;
    pVertex meshVertex;
    double newpt[3];

    // start mesh mover
    printf("Start mesh mover\n");
    pMeshMover mmover = MeshMover_new(pm, 0);

    // mesh motion of vertices in region
    vIter = M_vertexIter(pm);
    while((meshVertex = VIter_next(vIter))){
      apf::MeshEntity* vtx = reinterpret_cast<apf::MeshEntity*>(meshVertex);
      apf::getComponents(f, vtx, 0, vals);
      newpt[0] = vals[0];
      newpt[1] = vals[1];
      newpt[2] = vals[2];
      pPList closureRegion = GEN_regions(EN_whatIn(meshVertex));
      assert(PList_size(closureRegion));
      modelRegion = (pGRegion) PList_item(closureRegion, 0);
      printf("set move on (discrete) model: region %d; coord=(%12.16e,%12.16e,%12.16e)\n", GEN_tag(modelRegion),vals[0],vals[1],vals[2]);
      MeshMover_setDiscreteDeformMove(mmover,modelRegion,meshVertex,newpt);
      PList_delete(closureRegion);
    }
    VIter_delete(vIter);

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

    apf::Field* f = m->findField("motion_coords");
    assert(f);
    double* vals = new double[apf::countComponents(f)];
    assert(apf::countComponents(f) == 3);

    // declaration
    pGRegion modelRegion;
    pGFace modelFace;
    VIter vIter;
    pVertex meshVertex;
    double newpt[3];
    double newpar[2];
    double xyz[3];

    // start mesh mover
    printf("Start mesh mover\n");
    pMeshMover mmover = MeshMover_new(pm, 0);

    // configure mesh motion
    meshMotion mm = configureMotion(caseId, step);

    // mesh motion of moving body
    for (std::list<movingBodyMotion>::iterator mit = mm.movingBodyMotions.begin(); mit != mm.movingBodyMotions.end(); ++mit) {
      assert(modelRegion = (pGRegion) GM_entityByTag(model, 3, mit->tag));
      printf("set moving body: region %d\n",mit->tag);
      MeshMover_setTransform(mmover, modelRegion, mit->trans, mit->rotaxis, mit->rotpt, mit->rotang, mit->scale);
    }

    // mesh motion of vertices on surfaces
    std::list<int>::iterator lit;
    for (lit = mm.surfaceTags.begin(); lit != mm.surfaceTags.end(); ++lit) {
      assert(modelFace = (pGFace) GM_entityByTag(model, 2, *lit));
      printf("set move on surface: face %d\n",*lit);
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

    // mesh motion of vertices in region
    for (lit = mm.regionTags.begin(); lit != mm.regionTags.end(); ++lit) {
      assert(modelRegion = (pGRegion) GM_entityByTag(model, 3, *lit));
      printf("set move on region: region %d\n",*lit);
      vIter = M_classifiedVertexIter(pm, modelRegion, 0);
      while((meshVertex = VIter_next(vIter))){
        apf::MeshEntity* vtx = reinterpret_cast<apf::MeshEntity*>(meshVertex);
        apf::getComponents(f, vtx, 0, vals);
        const double newloc[3] = {vals[0], vals[1], vals[2]};
        MeshMover_setVolumeMove(mmover,meshVertex,newloc);
      }
      VIter_delete(vIter);
    }

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
    if (in.simmetrixMesh) {
      /* assumption: parametric model always needs a caseId
         when caseId = 0, it is supposed to be discrete model */
      if (!caseId)
        done = updateSIMDiscreteCoord(m);
      else
        done = updateSIMCoord(m,step,caseId);
    }
    else {
      done = updateAPFCoord(m);
    }
    assert(done);
  }

}
