#include "pcUpdateMesh.h"
#include "pcAdapter.h"
#include "pcSmooth.h"
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
#include <phastaChef.h>
#include <string.h>
#include <cassert>
#include <cstdio>

extern void MSA_setBLSnapping(pMSAdapt, int onoff);

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

  void addImproverInMover(pMeshMover& mmover, pPList sim_fld_lst) {
    // mesh improver
    if(!PCU_Comm_Self())
      printf("Add mesh improver attributes\n");
    pVolumeMeshImprover vmi = MeshMover_createImprover(mmover);
    VolumeMeshImprover_setModifyBL(vmi, 1);
    VolumeMeshImprover_setShapeMetric(vmi, ShapeMetricType_VolLenRatio, 0.3); // should be user-defined

    // set field to be mapped
    if (PList_size(sim_fld_lst))
      VolumeMeshImprover_setMapFields(vmi, sim_fld_lst);
  }

  void addAdapterInMover(pMeshMover& mmover,  pPList sim_fld_lst, apf::Mesh2*& m) {
    apf::Field* sizes = m->findField("sizes");
    apf::Vector3 v_mag = apf::Vector3(0.0,0.0,0.0);

    // mesh adapter
    if(!PCU_Comm_Self())
      printf("Add mesh adapter attributes\n");
    pMSAdapt msa = MeshMover_createAdapter(mmover);
    MSA_setAdaptBL(msa, 1);
    MSA_setExposedBLBehavior(msa,BL_DisallowExposed);
    MSA_setBLSnapping(msa, 0); // currently needed for parametric model
    MSA_setBLMinLayerAspectRatio(msa, 0.0); // needed in parallel

    // use current size field
    apf::MeshEntity* v;
    apf::MeshIterator* vit = m->begin(0);
    while ((v = m->iterate(vit))) {
      apf::getVector(sizes,v,0,v_mag);
      pVertex meshVertex = reinterpret_cast<pVertex>(v);
      MSA_setVertexSize(msa, meshVertex, v_mag[0]);
    }
    m->end(vit);

    // set field to be mapped
    if (PList_size(sim_fld_lst))
      MSA_setMapFields(msa, sim_fld_lst);

    // destroy mesh size field
    apf::destroyField(m->findField("sizes"));
  }

// temporarily used to write serial moved mesh and model
// it also writes coordinates to file
  bool updateAndWriteSIMDiscreteCoord(apf::Mesh2* m) {
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
    delete [] vals;

    // close file
    fclose (sFile);

    // write serial mesh and model
    if(!PCU_Comm_Self())
      printf("write discrete model and serial mesh\n");
    GM_write(M_model(pm), "discreteModel_serial.smd", 0, progress);
    M_write(pm, "mesh_serial.sms", 0, progress);
/*
    // load serial mesh and displacement
    pDiscreteModel dmodel = (pDiscreteModel) GM_load("discreteModel_serial.smd", 0, progress);
    pMesh mesh = M_load("mesh_serial.sms", dmodel, progress);
    sFile = fopen("allrank_id_disp.dat", "r");

    if(!PCU_Comm_Self())
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

    // close file
    fclose (sFile);

    // do real work
    if(!PCU_Comm_Self())
      printf("do real mesh mover\n");
    int isRunMover = MeshMover_run(mmover, progress);
    assert(isRunMover);
    MeshMover_delete(mmover);

    // write model and mesh
    if(!PCU_Comm_Self())
      printf("write updated discrete model and serial mesh\n");
    GM_write(M_model(pm), "updated_model.smd", 0, progress);
    M_write(pm, "after_mover_serial.sms", 0, progress);
    exit(0);
*/
}
    PCU_Barrier();
    Progress_delete(progress);
    return true;
  }

// temporarily used to write displacement field
  bool updateAndWriteSIMDiscreteField(apf::Mesh2* m) {
    Sim_logOn("updateAndWriteSIMDiscreteField.log");

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
      printf("write mesh: before_mover.sms\n");
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
    delete [] vals;

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

// hardcoding {
  void prescribe_proj_mesh_size(
    pGModel model, pMesh pm, apf::Mesh2* m,
    apf::Field* sizes, double proj_disp
  ) {
    // find the range of the refinement zone around the projectile
    pGFace p_left = (pGFace)GM_entityByTag(model, 2, 1);  // tail face
    pGFace p_rigt = (pGFace)GM_entityByTag(model, 2, 29); // head face
    double x_min = 1.6;
    double x_max = 0.0;

    double xyz[3];
    VIter vIter;
    pVertex meshVertex;
    vIter = M_classifiedVertexIter(pm, p_left, 1);
    while((meshVertex = VIter_next(vIter))){
      V_coord(meshVertex, xyz);
      if (xyz[0] < x_min) x_min = xyz[0];
    }
    VIter_delete(vIter);

    vIter = M_classifiedVertexIter(pm, p_rigt, 1);
    while((meshVertex = VIter_next(vIter))){
      V_coord(meshVertex, xyz);
      if (xyz[0] > x_max) x_max = xyz[0];
    }
    VIter_delete(vIter);

    PCU_Min_Doubles(&x_min, 1);
    PCU_Max_Doubles(&x_max, 1);

    x_min = x_min - 0.02 + proj_disp;
    x_max = x_max + 0.11 + proj_disp;

    printf("x_min and x_max = %f and %f\n",x_min,x_max);

    // set refinement zone 1: around the projectile size = D/40
    double ref_r = 0.048 + 0.0001;
    apf::Vector3 v_mag = apf::Vector3(0.0,0.0,0.0);
    apf::MeshEntity* v;
    apf::MeshIterator* vit = m->begin(0);
    while ((v = m->iterate(vit))) {
      pVertex meshVertex = reinterpret_cast<pVertex>(v);
      double xyz[3];
      V_coord(meshVertex, xyz);
      if (xyz[1]*xyz[1]+xyz[2]*xyz[2] <= ref_r*ref_r) {
        if (xyz[0] >= 0.0 && xyz[0] < x_min) {
          v_mag = apf::Vector3(0.012,0.012,0.012);
          apf::setVector(sizes,v,0,v_mag);
        }
        else if (xyz[0] >= x_min && xyz[0] <= x_max) {
          v_mag = apf::Vector3(0.003,0.003,0.003);
          apf::setVector(sizes,v,0,v_mag);
        }
        else if (xyz[0] > x_max && xyz[0] <= 1.6) {
          v_mag = apf::Vector3(0.012,0.012,0.012);
          apf::setVector(sizes,v,0,v_mag);
        }
      }
    }
    m->end(vit);
  }
// hardcoding }

// check if a model entity is (on) a rigid body
  int isOnRigidBody(pGModel model, pGEntity modelEnt, std::vector<ph::rigidBodyMotion> rbms) {
    for(unsigned id = 0; id < rbms.size(); id++)
      if(GEN_inClosure(GM_entityByTag(model, 3, rbms[id].tag), modelEnt)) return (int)id;
    // not find
    return -1;
  }

// auto detect non-rigid body model entities
  bool updateSIMCoordAuto(ph::Input& in, apf::Mesh2* m, int cooperation) {
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
    apf::NewArray<double> vals(apf::countComponents(f));
    assert(apf::countComponents(f) == 3);

    // declaration
    pGRegion modelRegion;
    pGFace modelFace;
    pGEdge modelEdge;
    pGVertex modelVertex;
    VIter vIter;
    pVertex meshVertex;
    double newpt[3];
    double newpar[2];
    double xyz[3];

    // start mesh mover
    if(!PCU_Comm_Self())
      printf("Start mesh mover\n");
    pMeshMover mmover = MeshMover_new(ppm, 0);

    std::vector<ph::rigidBodyMotion> rbms;
    core_get_rbms(rbms);
    // loop over model regions
    GRIter grIter = GM_regionIter(model);
    while((modelRegion=GRIter_next(grIter))){
      int id = isOnRigidBody(model, modelRegion, rbms);
      if(id >= 0) {
        assert(!GEN_isDiscreteEntity(modelRegion)); // should be parametric geometry
        printf("set rigid body motion: region %d; disp = (%e,%e,%e)\n",
                GEN_tag(modelRegion), rbms[id].trans[0], rbms[id].trans[1], rbms[id].trans[2]);
        MeshMover_setTransform(mmover, modelRegion, rbms[id].trans, rbms[id].rotaxis,
                                   rbms[id].rotpt, rbms[id].rotang, rbms[id].scale);
      }
      else {
        if (!GEN_isDiscreteEntity(modelRegion)) { // parametric
          vIter = M_classifiedVertexIter(pm, modelRegion, 0);
          while((meshVertex = VIter_next(vIter))){
            apf::MeshEntity* vtx = reinterpret_cast<apf::MeshEntity*>(meshVertex);
            apf::getComponents(f, vtx, 0, &vals[0]);
            const double newloc[3] = {vals[0], vals[1], vals[2]};
            MeshMover_setVolumeMove(mmover,meshVertex,newloc);
          }
          VIter_delete(vIter);
        }
        else { // discrete
          vIter = M_classifiedVertexIter(pm, modelRegion, 1);
          while((meshVertex = VIter_next(vIter))){
            apf::MeshEntity* vtx = reinterpret_cast<apf::MeshEntity*>(meshVertex);
            apf::getComponents(f, vtx, 0, &vals[0]);
            const double newloc[3] = {vals[0], vals[1], vals[2]};
            MeshMover_setDiscreteDeformMove(mmover,modelRegion,meshVertex,newloc);
          }
          VIter_delete(vIter);
        }
      }
    }
    GRIter_delete(grIter);

    // loop over model surfaces
    GFIter gfIter = GM_faceIter(model);
    while((modelFace=GFIter_next(gfIter))){
      int id = isOnRigidBody(model, modelFace, rbms);
      if(id >= 0) {
        assert(!GEN_isDiscreteEntity(modelFace)); // should be parametric geometry
        continue;
      }
      else {
        if (!GEN_isDiscreteEntity(modelFace)) { // parametric
          vIter = M_classifiedVertexIter(pm, modelFace, 0);
          while((meshVertex = VIter_next(vIter))){
            V_coord(meshVertex, xyz);
            apf::MeshEntity* vtx = reinterpret_cast<apf::MeshEntity*>(meshVertex);
            apf::getComponents(f, vtx, 0, &vals[0]);
            const double disp[3] = {vals[0]-xyz[0], vals[1]-xyz[1], vals[2]-xyz[2]};
            V_movedParamPoint(meshVertex,disp,newpar,newpt);
            MeshMover_setSurfaceMove(mmover,meshVertex,newpar,newpt);
          }
          VIter_delete(vIter);
        }
      }
    }
    GFIter_delete(gfIter);

    // loop over model edges
    GEIter geIter = GM_edgeIter(model);
    while((modelEdge=GEIter_next(geIter))){
      int id = isOnRigidBody(model, modelEdge, rbms);
      if(id >= 0) {
        assert(!GEN_isDiscreteEntity(modelEdge)); // should be parametric geometry
        continue;
      }
      else {
        if (!GEN_isDiscreteEntity(modelEdge)) { // parametric
          vIter = M_classifiedVertexIter(pm, modelEdge, 0);
          while((meshVertex = VIter_next(vIter))){
            V_coord(meshVertex, xyz);
            apf::MeshEntity* vtx = reinterpret_cast<apf::MeshEntity*>(meshVertex);
            apf::getComponents(f, vtx, 0, &vals[0]);
            const double disp[3] = {vals[0]-xyz[0], vals[1]-xyz[1], vals[2]-xyz[2]};
            V_movedParamPoint(meshVertex,disp,newpar,newpt);
            MeshMover_setSurfaceMove(mmover,meshVertex,newpar,newpt);
          }
          VIter_delete(vIter);
        }
      }
    }
    GEIter_delete(geIter);

    // loop over model vertices
    GVIter gvIter = GM_vertexIter(model);
    while((modelVertex=GVIter_next(gvIter))){
      int id = isOnRigidBody(model, modelVertex, rbms);
      if(id >= 0) {
        assert(!GEN_isDiscreteEntity(modelVertex)); // should be parametric geometry
        continue;
      }
      else {
        if (!GEN_isDiscreteEntity(modelVertex)) { // parametric
          vIter = M_classifiedVertexIter(pm, modelVertex, 0);
          while((meshVertex = VIter_next(vIter))){
            V_coord(meshVertex, xyz);
            apf::MeshEntity* vtx = reinterpret_cast<apf::MeshEntity*>(meshVertex);
            apf::getComponents(f, vtx, 0, &vals[0]);
            const double disp[3] = {vals[0]-xyz[0], vals[1]-xyz[1], vals[2]-xyz[2]};
            assert(sqrt(disp[0]*disp[0] + disp[1]*disp[1] + disp[2]*disp[2]) < 1e-10); // threshold 1e-10
          }
          VIter_delete(vIter);
        }
      }
    }
    GVIter_delete(gvIter);

    // add mesh improver and solution transfer
    pPList sim_fld_lst = PList_new();
    PList_clear(sim_fld_lst);
    if (cooperation) {
      if (in.solutionMigration)
        sim_fld_lst = getSimFieldList(in, m);
// create a field to store mesh size
      if(m->findField("sizes")) apf::destroyField(m->findField("sizes"));
      if(m->findField("frames")) apf::destroyField(m->findField("frames"));
      apf::Field* sizes  = apf::createSIMFieldOn(m, "sizes", apf::VECTOR);
      apf::Field* frames = apf::createSIMFieldOn(m, "frames", apf::MATRIX);
      ph::attachSIMSizeField(m, sizes, frames);

// prescribe mesh size field for the projectile case
// this is hardcoded, please comment out this call for other usage
//      prescribe_proj_mesh_size(model, pm, m, sizes, in.rbParamData[0]);

// add mesh smooth/gradation function here
//      addSmootherInMover(m);

      addImproverInMover(mmover, sim_fld_lst);
      addAdapterInMover(mmover, sim_fld_lst, m);
    }
    PList_clear(sim_fld_lst);
    PList_delete(sim_fld_lst);

    // do real work
    if(!PCU_Comm_Self())
      printf("do real mesh mover\n");
    int isRunMover = MeshMover_run(mmover, progress);
    assert(isRunMover);
    MeshMover_delete(mmover);

    // transfer sim fields to apf fields
    if (cooperation) {
      pc::transferSimFields(m);
    }

    // set rigid body total disp to be zero
    for (size_t i_rbpd = 0; (int)i_rbpd < 3*in.nRigidBody; i_rbpd++)
      in.rbParamData[i_rbpd] = 0.0;

    // write model and mesh
    if(!PCU_Comm_Self())
      printf("write model and mesh for mesh mover\n");
    GM_write(model, "updated_model.smd", 0, progress);
    PM_write(ppm, "after_mover.sms", progress);

    Progress_delete(progress);
    return true;
  }



  void runMeshMover(ph::Input& in, apf::Mesh2* m, int step, int cooperation) {
    bool done = false;
    if (in.simmetrixMesh) {
      done = updateSIMCoordAuto(in, m, cooperation);
    }
    else {
      done = updateAPFCoord(m);
    }
    assert(done);
  }

  void updateMesh(ph::Input& in, apf::Mesh2* m, apf::Field* szFld, int step, int cooperation) {
    if (in.simmetrixMesh && cooperation) {
      pc::runMeshMover(in,m,step,cooperation);
      m->verify();
      pc::writeSequence(m,step,"updated_mesh_");
    }
    else {
      pc::runMeshMover(in,m,step);
      m->verify();
      pc::writeSequence(m,step,"after_mover_");

      pc::runMeshAdapter(in,m,szFld,step);
      pc::writeSequence(m,step,"after_adapter_");
    }
  }

}
