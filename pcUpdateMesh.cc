#include "pcUpdateMesh.h"
#include "pcAdapter.h"
#include "pcWriteFiles.h"
#include <SimPartitionedMesh.h>
#include "SimAdvMeshing.h"
#include "SimModel.h"
#include "SimUtil.h"
#include "SimParasolidKrnl.h"
#include "SimField.h"
#include "SimMeshMove.h"
#include "SimMeshTools.h"
#include "SimDiscrete.h"
#include "MeshSim.h"
#include "MeshSimAdapt.h"
#include "apfSIM.h"
#include "gmi_sim.h"
#include <PCU.h>
#include <cassert>
#include <list>
#include <cstring>
#include <cstdlib>

extern void MSA_setBLSnapping(pMSAdapt, int onoff);

namespace pc {
  struct movingBodyMotion {
    movingBodyMotion(int t = 0, double r = 0.0, double s = 0.0)
    {
      tag = t;
      memset(trans, 0.0, sizeof trans);
      memset(rotaxis, 0.0, sizeof rotaxis);
      memset(rotpt, 0.0, sizeof rotpt);
      rotang = r;
      scale = s;
    }
    int tag;
    double trans[3];
    double rotaxis[3];
    double rotpt[3];
    double rotang;
    double scale;
    void set_trans(double x, double y, double z){
      trans[0] = x; trans[1] = y; trans[2] = z;
    }
    void set_rotaxis(double x, double y, double z){
      rotaxis[0] = x; rotaxis[1] = y; rotaxis[2] = z;
    }
    void set_rotpt(double x, double y, double z){
      rotpt[0] = x; rotpt[1] = y; rotpt[2] = z;
    }
  };

  struct meshMotion {
    std::list<movingBodyMotion> movingBodyMotions;
    std::list<int> surfaceTags;
    std::list<int> edgeTags;
    std::list<int> regionTags;
  };

  /* ideally, this comes from some input file */
  meshMotion configureMotion(int caseId, int step) {
    double offset, disp, rang, sfct;
    offset = 0.0;
    disp = 2e-4 + 2e-5;
    rang = 2.0;
    sfct = 0.9;
    if (caseId == 1) {
      offset = 2e-4 * (double)(step - 1) + 2e-5;
      printf("current step is %d; offset = %f\n",step,offset);
    }
    else if (caseId == 2) {
      offset = 4e-4;
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
    mbm = movingBodyMotion(1339, rang, sfct);
    mbm.set_trans(disp, 0.0, 0.0);
    mbm.set_rotaxis(0.0, 0.0, 1.0);
    mbm.set_rotpt(0.5e-3 + offset, 1.125e-3, 0.0);
    mm.movingBodyMotions.push_back(mbm);
// grain12
    mbm = movingBodyMotion(1036, rang, sfct);
    mbm.set_trans(disp, 0.0, 0.0);
    mbm.set_rotaxis(0.0, 0.0, 1.0);
    mbm.set_rotpt(2.625e-3 + offset, 1.125e-3, 0.0);
    mm.movingBodyMotions.push_back(mbm);
// grain13
    mbm = movingBodyMotion(733, rang, sfct);
    mbm.set_trans(disp, 0.0, 0.0);
    mbm.set_rotaxis(0.0, 0.0, 1.0);
    mbm.set_rotpt(4.75e-3 + offset, 1.125e-3, 0.0);
    mm.movingBodyMotions.push_back(mbm);
// grain21
    mbm = movingBodyMotion(1440, 0.0, sfct);
    mbm.set_trans(disp, 0.0, 0.0);
    mbm.set_rotaxis(0.0, 0.0, 1.0);
    mbm.set_rotpt(0.5e-3 + offset, 0.0, 0.0);
    mm.movingBodyMotions.push_back(mbm);
// grain22
    mbm = movingBodyMotion(1137, 0.0, sfct);
    mbm.set_trans(disp, 0.0, 0.0);
    mbm.set_rotaxis(0.0, 0.0, 1.0);
    mbm.set_rotpt(2.625e-3 + offset, 0.0, 0.0);
    mm.movingBodyMotions.push_back(mbm);
// grain23
    mbm = movingBodyMotion(834, 0.0, sfct);
    mbm.set_trans(disp, 0.0, 0.0);
    mbm.set_rotaxis(0.0, 0.0, 1.0);
    mbm.set_rotpt(4.75e-3 + offset, 0.0, 0.0);
    mm.movingBodyMotions.push_back(mbm);
// grain31
    mbm = movingBodyMotion(1238, rang, sfct);
    mbm.set_trans(disp, 0.0, 0.0);
    mbm.set_rotaxis(0.0, 0.0, -1.0);
    mbm.set_rotpt(0.5e-3 + offset, -1.125e-3, 0.0);
    mm.movingBodyMotions.push_back(mbm);
// grain32
    mbm = movingBodyMotion(935, rang, sfct);
    mbm.set_trans(disp, 0.0, 0.0);
    mbm.set_rotaxis(0.0, 0.0, -1.0);
    mbm.set_rotpt(2.625e-3 + offset, -1.125e-3, 0.0);
    mm.movingBodyMotions.push_back(mbm);
// grain33
    mbm = movingBodyMotion(632, rang, sfct);
    mbm.set_trans(disp, 0.0, 0.0);
    mbm.set_rotaxis(0.0, 0.0, -1.0);
    mbm.set_rotpt(4.75e-3 + offset, -1.125e-3, 0.0);
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

  void addImproverInMover(pMeshMover mmover, pPList sim_fld_lst) {
    // mesh improver
    if(!PCU_Comm_Self())
      printf("Add mesh improver\n");
    pVolumeMeshImprover vmi = MeshMover_createImprover(mmover);
    VolumeMeshImprover_setModifyBL(vmi, 1);
    VolumeMeshImprover_setShapeMetric(vmi, ShapeMetricType_VolLenRatio, 0.3);

    // set field to be mapped
    if (sim_fld_lst)
      VolumeMeshImprover_setMapFields(vmi, sim_fld_lst);
  }

  void addAdapterInMover(pMeshMover mmover,  pPList sim_fld_lst, apf::Mesh2*& m) {
    // mesh adapter
    if(!PCU_Comm_Self())
      printf("Add mesh adapter\n");
    pMSAdapt msa = MeshMover_createAdapter(mmover);
    MSA_setAdaptBL(msa, 1);
    MSA_setExposedBLBehavior(msa,BL_DisallowExposed);
    MSA_setBLSnapping(msa, 0); // currently needed for parametric model
    MSA_setBLMinLayerAspectRatio(msa, 0.0); // needed in parallel

    // use current size field
    apf::MeshEntity* v;
    apf::MeshIterator* vit = m->begin(0);
    while ((v = m->iterate(vit))) {
      pVertex meshVertex = reinterpret_cast<pVertex>(v);
      MSA_scaleVertexSize(msa, meshVertex, 1.0); // use the size field of the mesh before mesh motion
    }
    m->end(vit);

    // set field to be mapped
    if (sim_fld_lst)
      MSA_setMapFields(msa, sim_fld_lst);
  }

// temporarily used to write serial moved mesh and model
  bool updateAndWriteSIMDiscreteModel(apf::Mesh2* m) {
    Sim_logOn("updateAndWriteSIMDiscreteModel.log");

    pProgress progress = Progress_new();
    Progress_setDefaultCallback(progress);

    apf::MeshSIM* apf_msim = dynamic_cast<apf::MeshSIM*>(m);
    gmi_model* gmiModel = apf_msim->getModel();
    pGModel model = gmi_export_sim(gmiModel);
    pParMesh ppm = apf_msim->getMesh();

    PM_write(ppm, "before_mover_parallel.sms", progress);

    // declaration
    pGRegion modelRegion;
    VIter vIter;
    pVertex meshVertex;
    double disp[3];
    double xyz[3];
    double newpt[3];

    pMesh pm = M_createFromParMesh(ppm,3,progress);
//    M_release(ppm);
    if (pm)
      M_write(pm, "before_mover_serial.sms", 0, progress);

    FILE* sFile;
    int id;
if (pm) {
    // open file for writing displacements
    sFile = fopen ("allrank_id_disp.dat", "w");

    apf::Field* f = m->findField("motion_coords");
    assert(f);
    double* vals = new double[apf::countComponents(f)];
    assert(apf::countComponents(f) == 3);

    // write id and displacement to file
    id = 0;
    vIter = M_vertexIter(pm);
    while((meshVertex = VIter_next(vIter))){
      apf::MeshEntity* vtx = reinterpret_cast<apf::MeshEntity*>(meshVertex);
      apf::getComponents(f, vtx, 0, vals);
      V_coord(meshVertex, xyz);
      disp[0] = vals[0] - xyz[0];
      disp[1] = vals[1] - xyz[1];
      disp[2] = vals[2] - xyz[2];

      fprintf(sFile, "%010d ", id);
      for(int ii = 0; ii < 3; ii++) {
        fprintf(sFile, " %f", disp[ii]);
      }
      fprintf(sFile, "\n");

      id++;
    }
    VIter_delete(vIter);

    // close file
    fclose (sFile);

    // write serial mesh and model
    printf("write discrete model and serial mesh\n");
    GM_write(M_model(pm), "discreteModel_serial.smd", 0, progress);
    M_write(pm, "mesh_serial.sms", 0, progress);

    // load serial mesh and displacement
    pDiscreteModel dmodel = (pDiscreteModel) GM_load("discreteModel_serial.smd", 0, progress);
    pMesh mesh = M_load("mesh_serial.sms", dmodel, progress);
    sFile = fopen("allrank_id_disp.dat", "r");

    printf("start mesh mover on the serial mesh\n");
    pMeshMover mmover = MeshMover_new(mesh, 0);

    // mesh motion of vertices in region
    int counter = 0;
    vIter = M_vertexIter(mesh);
    while((meshVertex = VIter_next(vIter))){

      fscanf(sFile, "%010d", &id);
      assert(counter == id);
      for (int i = 0; i < 3; i++)
        fscanf(sFile, "%lf", &disp[i]);
      counter++;

      V_coord(meshVertex, xyz);
      newpt[0] = xyz[0] + disp[0];
      newpt[1] = xyz[1] + disp[1];
      newpt[2] = xyz[2] + disp[2];

      pPList closureRegion = GEN_regions(EN_whatIn(meshVertex));
      assert(PList_size(closureRegion));
      modelRegion = (pGRegion) PList_item(closureRegion, 0);
      assert(GEN_isDiscreteEntity(modelRegion));
      MeshMover_setDiscreteDeformMove(mmover,modelRegion,meshVertex,newpt);
      PList_delete(closureRegion);
    }
    VIter_delete(vIter);

    // do real work
    printf("do real mesh mover\n");
    assert(MeshMover_run(mmover, progress));
    MeshMover_delete(mmover);

    // write model and mesh
    printf("write updated discrete model and serial mesh\n");
    GM_write(M_model(pm), "updated_model.smd", 0, progress);
    M_write(pm, "after_mover_serial.sms", 0, progress);
    exit(0);
}

    Progress_delete(progress);
    return true;
  }



// temporarily used to write displacement field
  bool updateAndWriteSIMDiscrete(apf::Mesh2* m) {
    Sim_logOn("updateAndWriteSIMDiscrete.log");

    pProgress progress = Progress_new();
    Progress_setDefaultCallback(progress);

    apf::MeshSIM* apf_msim = dynamic_cast<apf::MeshSIM*>(m);
    gmi_model* gmiModel = apf_msim->getModel();
    pGModel model = gmi_export_sim(gmiModel);
    pParMesh ppm = apf_msim->getMesh();

    // declaration
    pGRegion modelRegion;
    VIter vIter;
    pVertex meshVertex;
    double disp[3];
    double xyz[3];

/*
    pMesh pm = M_createFromParMesh(ppm,3,progress);
    M_release(ppm);
*/
/*
    // migrate mesh to part 0
    pGEntMeshMigrator gmig = GEntMeshMigrator_new(ppm, 3);
    GRIter grIter = GM_regionIter(model);
    while((modelRegion = GRIter_next(grIter))){
      int gid = PMU_gid(0, 0);
      GEntMeshMigrator_add(gmig, modelRegion, gid);
    }
    GRIter_delete(grIter);
    GEntMeshMigrator_run(gmig, progress);
    GEntMeshMigrator_delete(gmig);
*/
    pMesh pm = PM_mesh(ppm,0);

    if(!PCU_Comm_Self())
      printf("write mesh just after migration: before_mover.sms\n");
    PM_write(ppm, "before_mover.sms", progress);

    apf::Field* f = m->findField("motion_coords");
    assert(f);
    double* vals = new double[apf::countComponents(f)];
    assert(apf::countComponents(f) == 3);

    // start transfer to pField
    if(!PCU_Comm_Self())
      printf("Start transfer to pField\n");
    pPolyField pf = PolyField_new(1, 0);
    pField dispFd = Field_new(ppm, 3, "disp", "displacement", ShpLagrange, 1, 1, 1, pf);
    Field_apply(dispFd, 3, progress);
    pDofGroup dof;

    // mesh motion of vertices in region
    vIter = M_vertexIter(pm);
    while((meshVertex = VIter_next(vIter))){
      apf::MeshEntity* vtx = reinterpret_cast<apf::MeshEntity*>(meshVertex);
      apf::getComponents(f, vtx, 0, vals);
      V_coord(meshVertex, xyz);
      disp[0] = vals[0] - xyz[0];
      disp[1] = vals[1] - xyz[1];
      disp[2] = vals[2] - xyz[2];

      for(int eindex = 0; (dof = Field_entDof(dispFd,meshVertex,eindex)); eindex++) {
        int dofs_per_node = DofGroup_numComp(dof);
        assert(dofs_per_node == 3);
        for(int ii = 0; ii < dofs_per_node; ii++)
          DofGroup_setValue(dof,ii,0,disp[ii]);
      }
    }
    VIter_delete(vIter);

    // do real work
    if(!PCU_Comm_Self())
      printf("write sim field\n");
    Field_write(dispFd, "dispFd.fld", 0, NULL, progress);

    // used to write discrete model at a certain time step
//    pDiscreteModel dmodel;

/*
    dmodel = DM_createFromMesh(pm, 0, progress);
    // define the Discrete model
    DM_findEdgesByFaceNormals(dmodel, 0, progress);
    DM_eliminateDanglingEdges(dmodel, progress);
    assert(!DM_completeTopology(dmodel, progress));
*/
/*
    if(!PCU_Comm_Self()) {
      dmodel = DM_createFromModel(model, pm);
      assert(!DM_completeTopology(dmodel, progress));
      printf("write extracted discrete model\n");
      GM_write(dmodel, "extracted_model.smd", 0, progress);
    }
*/
    // write model and mesh
    if(!PCU_Comm_Self())
      printf("write mesh after creating discrete model: after_mover.sms\n");
    GM_write(model, "updated_model.smd", 0, progress);
    PM_write(ppm, "after_mover.sms", progress);

    Progress_delete(progress);
    return true;
  }

  bool updateSIMDiscreteCoord(ph::Input& in, apf::Mesh2* m, int cooperation) {
    Sim_logOn("updateSIMDiscreteCoord.log");

    pProgress progress = Progress_new();
    Progress_setDefaultCallback(progress);

    apf::MeshSIM* apf_msim = dynamic_cast<apf::MeshSIM*>(m);
    pParMesh ppm = apf_msim->getMesh();
    pMesh pm = PM_mesh(ppm,0);

    PM_write(ppm, "before_mover.sms", progress);

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
    if(!PCU_Comm_Self())
      printf("Start mesh mover\n");
    pMeshMover mmover = MeshMover_new(ppm, 0);

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
      assert(GEN_isDiscreteEntity(modelRegion));
      MeshMover_setDiscreteDeformMove(mmover,modelRegion,meshVertex,newpt);
      PList_delete(closureRegion);
    }
    VIter_delete(vIter);

    // add mesh improver and solution transfer
    if (cooperation) {
      pPList sim_fld_lst;
      if (in.solutionMigration)
        sim_fld_lst = getSimFieldList(in, m);
      addImproverInMover(mmover, sim_fld_lst);
      addAdapterInMover(mmover, sim_fld_lst, m);
    }

    // do real work
    if(!PCU_Comm_Self())
      printf("do real mesh mover\n");
    assert(MeshMover_run(mmover, progress));
    MeshMover_delete(mmover);

    // transfer sim fields to apf fields
    if (cooperation) {
      pc::transferSimFields(m);
    }

    // write model and mesh
    if(!PCU_Comm_Self())
      printf("write discrete model and mesh for mesh mover\n");
    GM_write(model, "updated_model.smd", 0, progress);
    PM_write(ppm, "after_mover.sms", progress);

    Progress_delete(progress);
    return true;
  }

  bool updateSIMCoord(ph::Input& in, apf::Mesh2* m, int step, int caseId, int cooperation) {
    Sim_logOn("updateSIMCoord.log");

    pProgress progress = Progress_new();
    Progress_setDefaultCallback(progress);

    apf::MeshSIM* apf_msim = dynamic_cast<apf::MeshSIM*>(m);
    pParMesh ppm = apf_msim->getMesh();
    pMesh pm = PM_mesh(ppm,0);

    PM_write(ppm, "before_mover.sms", progress);

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
    if(!PCU_Comm_Self())
      printf("Start mesh mover\n");
    pMeshMover mmover = MeshMover_new(ppm, 0);

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

    // add mesh improver and solution transfer
    if (cooperation) {
      pPList sim_fld_lst;
      if (in.solutionMigration)
        sim_fld_lst = getSimFieldList(in, m);
      addImproverInMover(mmover, sim_fld_lst);
      addAdapterInMover(mmover, sim_fld_lst, m);
    }

    // do real work
    if(!PCU_Comm_Self())
      printf("do real mesh mover\n");
    assert(MeshMover_run(mmover, progress));
    MeshMover_delete(mmover);

    // transfer sim fields to apf fields
    if (cooperation) {
      pc::transferSimFields(m);
    }

    // write model and mesh
    if(!PCU_Comm_Self())
      printf("write model and mesh for mesh mover\n");
    GM_write(model, "updated_model.smd", 0, progress);
    PM_write(ppm, "after_mover.sms", progress);

    Progress_delete(progress);
    return true;
  }

  void runMeshMover(ph::Input& in, apf::Mesh2* m, int step, int caseId, int cooperation) {
    bool done = false;
    if (in.simmetrixMesh) {
      /* assumption: parametric model always needs a caseId
         when caseId = 0, it is supposed to be discrete model */
      if (caseId == 0)
        done = updateSIMDiscreteCoord(in, m, cooperation);
      else if (caseId == 10) // hack!
        done = updateAndWriteSIMDiscrete(m); // hack!
      else if (caseId == 20) // hack!
        done = updateAndWriteSIMDiscreteModel(m); // hack!
      else
        done = updateSIMCoord(in, m, step, caseId, cooperation);
    }
    else {
      done = updateAPFCoord(m);
    }
    assert(done);
  }

  void updateMesh(ph::Input& in, apf::Mesh2* m, apf::Field* szFld, int step, int caseId, int cooperation) {
    if (in.simmetrixMesh && cooperation) {
      pc::runMeshMover(in,m,step,caseId,1);
      m->verify();
      pc::writeSequence(m,step,"updated_mesh_");
    }
    else {
      pc::runMeshMover(in,m,step,caseId);
      m->verify();
      pc::writeSequence(m,step,"after_mover_");

      pc::runMeshAdapter(in,m,szFld,step);
      pc::writeSequence(m,step,"after_adapter_");
    }
  }

}
