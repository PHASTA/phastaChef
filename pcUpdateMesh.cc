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

  bool updateAPFCoord(ph::Input& in, apf::Mesh2* m) {
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
    pc::writeSequence(m, in.timeStepNumber, "pvtu_mesh_");
    return true;
  }

  void addImproverInMover(pMeshMover& mmover, pPList sim_fld_lst) {
    // mesh improver
    if(!PCU_Comm_Self())
      printf("Add mesh improver attributes\n");
    pVolumeMeshImprover vmi = MeshMover_createImprover(mmover);
    pc::setupSimImprover(vmi, sim_fld_lst);
  }

  void addAdapterInMover(pMeshMover& mmover,  pPList& sim_fld_lst, ph::Input& in, apf::Mesh2*& m) {
    // mesh adapter
    if(!PCU_Comm_Self())
      printf("Add mesh adapter attributes\n");
    pMSAdapt msa = MeshMover_createAdapter(mmover);
    pc::setupSimAdapter(msa, in, m, sim_fld_lst);
  }

  void balanceEqualWeights(pParMesh pmesh, pProgress progress) {
    // get total number of processors
    int totalNumProcs = PMU_size();
    // get total number of parts
    int totalNumParts = PM_totalNumParts(pmesh);

    // we assume one part per processor
    if (totalNumProcs != totalNumParts) {
      if( !PCU_Comm_Self() )
        fprintf(stderr, "Error: N of procs %d not equal to N of parts %d\n",
                totalNumProcs, totalNumParts);
    }

    // start load balance
    if(!PCU_Comm_Self())
      printf("Start load balance\n");
    pPartitionOpts pOpts = PartitionOpts_new();
    // Set total no. of parts
    PartitionOpts_setTotalNumParts(pOpts, totalNumParts);
    // Sets processes to be equally weighted
    PartitionOpts_setProcWtEqual(pOpts);
    PM_partition(pmesh, pOpts, progress);     // Do the partitioning
    PartitionOpts_delete(pOpts);              // Done with options
    // print out elements of each part
    pMesh mesh = PM_mesh(pmesh,0);
    int numElmOnPart = M_numRegions(mesh);
    long numTolElm = PCU_Add_Long(numElmOnPart);
    if(!PCU_Comm_Self())
      printf("Total No. of Elm: %d\n", numTolElm);
  }


// temporarily used to write serial moved mesh and model
// it also writes coordinates to file
  bool updateAndWriteSIMDiscreteCoord(apf::Mesh2* m) {
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
    writeSIMModel(model, in.timeStepNumber, "sim_model_");
    writeSIMMesh(ppm, in.timeStepNumber, "sim_mesh_");

    Progress_delete(progress);
    return true;
  }

// hardcoding {
  bool isInCircularSector(double y, double z) {
    double PI = 3.14159265;
    double degree = atan2(z,y) * 180.0 / PI + 180.0;
    if ( (degree >= 0.0        && degree <= 3.3)
      || (degree >= -3.3+45.0  && degree <= 3.3+45.0)
      || (degree >= -3.3+90.0  && degree <= 3.3+90.0)
      || (degree >= -3.3+135.0 && degree <= 3.3+135.0)
      || (degree >= -3.3+180.0 && degree <= 3.3+180.0)
      || (degree >= -3.3+225.0 && degree <= 3.3+225.0)
      || (degree >= -3.3+270.0 && degree <= 3.3+270.0)
      || (degree >= -3.3+315.0 && degree <= 3.3+315.0)
      || (degree >= -3.3+360.0 && degree <= 360.0) ) {
      return true;
    }
    return false;
  }

  void prescribe_proj_mesh_size(
    apf::Mesh2* m, apf::Field* sizes, double proj_disp)
  {
    // hardcoded parameters
    int tailFace = 110; // 2mm case
    int headFace = 145; // 2mm case
//    int tailFace = 473; // touch case
//    int headFace = 553; // touch case
    double tail = 2.0;
    double head = 0.0;
    double t_min = 0.0;
    double t_max = 2.0;
    double fgap_zone_left_dis_tail = 0.002;
    double fgap_zone_rigt_dis_tail = 0.028;
    double step_fine_rig1_dis_tail = 0.435;
    double step_fine_rig2_dis_tail = 0.581;
    double step_mode_rig1_dis_tail = 0.411;
    double step_mode_rig2_dis_tail = 0.597;
    double step_cors_rig1_dis_tail = 0.363;
    double step_cors_rig2_dis_tail = 0.629;
    double fine_zone_left_dis_tail = 0.008;
    double fine_zone_rigt_dis_tail = 0.052;
    double mode_zone_left_dis_tail = 0.020;
    double mode_zone_rigt_dis_head = 0.110;
    double fgap_r = 0.052;
    double zone_r = 0.044;
    double tube_r = 0.060;
    double ref_tol = 0.0001;
    tube_r = tube_r + ref_tol;
    apf::Vector3 fgap_size = apf::Vector3(0.0005,0.0005,0.0005);
    apf::Vector3 stfn_size = apf::Vector3(0.0005,0.0005,0.0005);
    apf::Vector3 stmd_size = apf::Vector3(0.001,0.001,0.001);
    apf::Vector3 stcr_size = apf::Vector3(0.002,0.002,0.002);
    apf::Vector3 fine_size = apf::Vector3(0.002,0.002,0.002);
    apf::Vector3 mode_size = apf::Vector3(0.004,0.004,0.004);
    apf::Vector3 cors_size = apf::Vector3(0.008,0.008,0.008);

    // load model and mesh
    apf::MeshSIM* apf_msim = dynamic_cast<apf::MeshSIM*>(m);
    pParMesh ppm = apf_msim->getMesh();
    pMesh pm = PM_mesh(ppm,0);

    gmi_model* gmiModel = apf_msim->getModel();
    pGModel model = gmi_export_sim(gmiModel);

    // find the range of the refinement zone around the projectile
    pGFace p_left = (pGFace)GM_entityByTag(model, 2, tailFace); // tail face
    pGFace p_rigt = (pGFace)GM_entityByTag(model, 2, headFace); // head face

    double xyz[3];
    VIter vIter;
    pVertex meshVertex;
    vIter = M_classifiedVertexIter(pm, p_left, 1);
    while((meshVertex = VIter_next(vIter))){
      V_coord(meshVertex, xyz);
      if (xyz[0] < tail) tail = xyz[0];
    }
    VIter_delete(vIter);

    vIter = M_classifiedVertexIter(pm, p_rigt, 1);
    while((meshVertex = VIter_next(vIter))){
      V_coord(meshVertex, xyz);
      if (xyz[0] > head) head = xyz[0];
    }
    VIter_delete(vIter);

    PCU_Min_Doubles(&tail, 1);
    PCU_Max_Doubles(&head, 1);

    double fgap_zone_min = proj_disp + tail - fgap_zone_left_dis_tail;
    double fgap_zone_max = proj_disp + tail + fgap_zone_rigt_dis_tail;
    double step_fine_min = proj_disp + tail + step_fine_rig1_dis_tail;
    double step_fine_max = proj_disp + tail + step_fine_rig2_dis_tail;
    double step_mode_min = proj_disp + tail + step_mode_rig1_dis_tail;
    double step_mode_max = proj_disp + tail + step_mode_rig2_dis_tail;
    double step_cors_min = proj_disp + tail + step_cors_rig1_dis_tail;
    double step_cors_max = proj_disp + tail + step_cors_rig2_dis_tail;
    double fine_zone_min = proj_disp + tail - fine_zone_left_dis_tail;
    double fine_zone_max = proj_disp + tail + fine_zone_rigt_dis_tail;
    double mode_zone_min = proj_disp + tail - mode_zone_left_dis_tail;
    double mode_zone_max = proj_disp + head + mode_zone_rigt_dis_head;

    if(!PCU_Comm_Self())
      printf("tail and head and proj_disp = %f and %f and %f\n",tail,head,proj_disp);

    apf::Vector3 cur_size = apf::Vector3(0.0,0.0,0.0);
    apf::MeshEntity* v;
    apf::MeshIterator* vit = m->begin(0);
    while ((v = m->iterate(vit))) {
      pVertex meshVertex = reinterpret_cast<pVertex>(v);
      double xyz[3];
      V_coord(meshVertex, xyz);
      if (xyz[1]*xyz[1]+xyz[2]*xyz[2] <= zone_r*zone_r) {
        if (xyz[0] >= t_min && xyz[0] < mode_zone_min) {
          apf::setVector(sizes,v,0,cors_size);
        }
        else if (xyz[0] >= mode_zone_min && xyz[0] < fine_zone_min) {
          apf::setVector(sizes,v,0,mode_size);
        }
        else if (xyz[0] >= fine_zone_min && xyz[0] <= fine_zone_max) {
          apf::setVector(sizes,v,0,fine_size);
        }
        else if (xyz[0] > fine_zone_max && xyz[0] < mode_zone_max) {
          apf::setVector(sizes,v,0,mode_size);
        }
        else if (xyz[0] >= mode_zone_max && xyz[0] < t_max) {
          apf::setVector(sizes,v,0,cors_size);
        }
      }
      else if (xyz[1]*xyz[1]+xyz[2]*xyz[2] >  zone_r*zone_r
            && xyz[1]*xyz[1]+xyz[2]*xyz[2] <= tube_r*tube_r ) {
        if (xyz[1]*xyz[1]+xyz[2]*xyz[2] >  fgap_r*fgap_r
         && xyz[0] >= fgap_zone_min && xyz[0] <= fgap_zone_max
         && isInCircularSector(xyz[1],xyz[2]) ) {
          apf::setVector(sizes,v,0,fgap_size);
        }
        else if (xyz[0] >= step_fine_min && xyz[0] <= step_fine_max) {
          apf::setVector(sizes,v,0,stfn_size);
        }
        else if (xyz[0] >= step_mode_min && xyz[0] <= step_mode_max) {
          apf::setVector(sizes,v,0,stmd_size);
        }
        else if (xyz[0] >= step_cors_min && xyz[0] <= step_cors_max) {
          apf::setVector(sizes,v,0,stcr_size);
        }
        else if (xyz[0] >= t_min && xyz[0] < fine_zone_min) {
          apf::setVector(sizes,v,0,mode_size);
        }
        else if (xyz[0] >= fine_zone_min && xyz[0] <= fine_zone_max) {
          apf::setVector(sizes,v,0,fine_size);
        }
        else if (xyz[0] > fine_zone_max && xyz[0] < t_max) {
          apf::setVector(sizes,v,0,mode_size);
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
    if (in.writeSimLog)
      Sim_logOn("updateSIMCoord.log");

    pProgress progress = Progress_new();
    Progress_setDefaultCallback(progress);

    apf::MeshSIM* apf_msim = dynamic_cast<apf::MeshSIM*>(m);
    pParMesh ppm = apf_msim->getMesh();
    pMesh pm = PM_mesh(ppm,0);

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
    if (in.nRigidBody > 0) {
      core_get_rbms(rbms);
    }
    else {
      rbms.clear();
    }
    // loop over model regions
    GRIter grIter = GM_regionIter(model);
    while((modelRegion=GRIter_next(grIter))){
      int id = isOnRigidBody(model, modelRegion, rbms);
      if(id >= 0) {
        assert(!GEN_isDiscreteEntity(modelRegion)); // should be parametric geometry
/*
        printf("set rigid body motion: region %d; disp = (%e,%e,%e); rotaxis = (%e,%e,%e); rotpt = (%e,%e,%e); rotang = %f; scale = %f\n",
                GEN_tag(modelRegion), rbms[id].trans[0], rbms[id].trans[1], rbms[id].trans[2],
                rbms[id].rotaxis[0], rbms[id].rotaxis[1], rbms[id].rotaxis[2],
                rbms[id].rotpt[0], rbms[id].rotpt[1], rbms[id].rotpt[2],
                rbms[id].rotang, rbms[id].scale);
*/
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
      addAdapterInMover(mmover, sim_fld_lst, in, m);
      addImproverInMover(mmover, sim_fld_lst);
    }

    // do real work
    if(!PCU_Comm_Self())
      printf("do real mesh mover\n");
    int isRunMover = MeshMover_run(mmover, progress);
    assert(isRunMover);
    MeshMover_delete(mmover);

    PList_clear(sim_fld_lst);
    PList_delete(sim_fld_lst);

    if (cooperation) {
      // load balance
      balanceEqualWeights(ppm, progress);

      // transfer sim fields to apf fields
      if (in.solutionMigration)
        pc::transferSimFields(m);
    }

    // set rigid body total disp to be zero
    for (size_t i_nrbs = 0; (int)i_nrbs < in.nRigidBody; i_nrbs++) {
      for (size_t i_rbpd = 0; (int)i_rbpd < 3; i_rbpd++)
        in.rbParamData[i_rbpd+i_nrbs*(size_t)in.nRBParam] = 0.0;
      for (size_t i_rbpd = 12; (int)i_rbpd < 13; i_rbpd++)
        in.rbParamData[i_rbpd+i_nrbs*(size_t)in.nRBParam] = 0.0;
    }

    // write model and mesh
    if(!PCU_Comm_Self())
      printf("write model and mesh after mesh modification\n");
    writeSIMModel(model, in.timeStepNumber, "sim_model_");
    if (cooperation)
      writeSIMMesh(ppm, in.timeStepNumber, "sim_mesh_");
    else
      writeSIMMesh(ppm, in.timeStepNumber, "sim_moved_mesh_");

    Progress_delete(progress);
    return true;
  }



  void runMeshMover(ph::Input& in, apf::Mesh2* m, int step, int cooperation) {
    bool done = false;
    if (in.simmetrixMesh) {
      done = updateSIMCoordAuto(in, m, cooperation);
    }
    else {
      done = updateAPFCoord(in, m);
    }
    assert(done);
  }

  void updateMesh(ph::Input& in, apf::Mesh2* m, apf::Field* szFld, int step, int cooperation) {
    if (in.simmetrixMesh && cooperation) {
      pc::runMeshMover(in,m,step,cooperation);
      m->verify();
    }
    else {
      pc::runMeshMover(in,m,step);
      m->verify();
      pc::runMeshAdapter(in,m,szFld,step);
      m->verify();
    }
  }

}
