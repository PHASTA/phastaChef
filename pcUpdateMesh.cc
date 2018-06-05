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
#include <phastaChef.h>
#include <string.h>
#include <cassert>
#include <cstdio>

extern void MSA_setBLSnapping(pMSAdapt, int onoff);

namespace pc {
/*
  int getValueByString(char* header, char* token, double* val) {
    if (token == NULL || token[0] == '\n') return 0;
    int counter = 0;
    if (strcmp(header, token) == 0)
      printf("found token: %s\n", token); //debugging
    else
      return 0;
    while(token != NULL && token[0] != '\n') {
      token = strtok (NULL, " \n");
      if (token != NULL && token[0] != '\n') {
        val[counter] = atof(token);
        counter++;
      }
    }
    return counter;
  }

  bool configureMotion(int step, meshMotion& mm) {
    FILE* pFile = fopen("meshMotion.inp", "r");
    if(!pFile) {
      printf("cannot find file: meshMotion.inp\n");
      return false;
    }
    char target[1024] = "\0";
    char* token;
    double* val = new double[100];
    double* pdVal = new double[100];
    int i = 0;
    int nt = 0;
    int npd = 0;
    int caseId = 0;
    while( fgets(target, 1024, pFile) && !feof(pFile) ) {
      if (target[0] == '#') continue; // ignore comment line
      token = strtok ( target, ":\n" );

      nt = getValueByString("Parametric deforming", token, val);
      if(nt > 0) npd = nt;
      for(i = 0; i < nt; i++)
        pdVal[i] = val[i];

      nt = getValueByString("Time depended parametric deforming case ID", token, val);
      if(nt)
        caseId = (int)val[0];

      nt = getValueByString("Parametric surrounding regions", token, val);
      for(i = 0; i < nt; i++)
        mm.parSurRegions.push_back((int)val[i]);

      nt = getValueByString("Parametric surrounding faces", token, val);
      for(i = 0; i < nt; i++)
        mm.parSurFaces.push_back((int)val[i]);

      nt = getValueByString("Parametric surrounding edges", token, val);
      for(i = 0; i < nt; i++)
        mm.parSurEdges.push_back((int)val[i]);

      nt = getValueByString("Discrete deforming", token, val);
      for(i = 0; i < nt; i++)
        mm.disDefRegions.push_back((int)val[i]);

      nt = getValueByString("Discrete surrounding", token, val);
      for(i = 0; i < nt; i++)
        mm.disSurRegions.push_back((int)val[i]);
    }

    if (caseId > 0) {
      mm = pc::getTDMeshMotion(caseId, step);
    }
    else {
      if (npd) {
        rewind(pFile);
        while( fgets(target, 1024, pFile) && !feof(pFile) ) {
          if (target[0] == '#') continue; // ignore comment line
          token = strtok ( target, ":\n" );
          char buf[50] = "\0";
          for(i = 0; i < npd; i++){
            sprintf(buf, "Parametric deforming tag %d", (int)pdVal[i]);
            nt = getValueByString(buf, token, val);
            if (nt > 0) {
              assert(nt == 11);
              for(int ii = 0; ii < nt; ii++) {
                printf("val[%d] = %f\n",ii,val[ii]);
              }
              rigidBodyMotion rbm = rigidBodyMotion((int)pdVal[i], val[9], val[10]);
              rbm.set_trans(val[0], val[1], val[2]);
              rbm.set_rotaxis(val[3], val[4], val[5]);
              rbm.set_rotpt(val[6], val[7], val[8]);
              mm.rigidBodyMotions.push_back(rbm);
            }
          }
        }
      }
      return true;
    }
  }
*/

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
    VolumeMeshImprover_setShapeMetric(vmi, ShapeMetricType_VolLenRatio, 0.3);

    // set field to be mapped
    if (PList_size(sim_fld_lst))
      VolumeMeshImprover_setMapFields(vmi, sim_fld_lst);
  }

  void addAdapterInMover(pMeshMover& mmover,  pPList sim_fld_lst, apf::Mesh2*& m) {
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
      pVertex meshVertex = reinterpret_cast<pVertex>(v);
      MSA_scaleVertexSize(msa, meshVertex, 1.0); // use the size field of the mesh before mesh motion
    }
    m->end(vit);

    // set field to be mapped
    if (PList_size(sim_fld_lst))
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
    delete [] vals;

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

    // close file
    fclose (sFile);

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



/*
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

    // start mesh mover
    if(!PCU_Comm_Self())
      printf("Start mesh mover\n");
    pMeshMover mmover = MeshMover_new(ppm, 0);

    // mesh motion of vertices in region
    vIter = M_vertexIter(pm);
    while((meshVertex = VIter_next(vIter))){
      apf::MeshEntity* vtx = reinterpret_cast<apf::MeshEntity*>(meshVertex);
      apf::getComponents(f, vtx, 0, vals);
      const double newloc[3] = {vals[0], vals[1], vals[2]};
      pPList closureRegion = GEN_regions(EN_whatIn(meshVertex));
      assert(PList_size(closureRegion));
      modelRegion = (pGRegion) PList_item(closureRegion, 0);
      assert(GEN_isDiscreteEntity(modelRegion));
      MeshMover_setDiscreteDeformMove(mmover,modelRegion,meshVertex,newloc);
      PList_delete(closureRegion);
    }
    VIter_delete(vIter);
    delete [] vals;

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

  bool updateSIMCoord(ph::Input& in, apf::Mesh2* m, int cooperation, meshMotion& mm) {
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
    pGEdge modelEdge;
    VIter vIter;
    pVertex meshVertex;
    double newpt[3];
    double newpar[2];
    double xyz[3];

    // start mesh mover
    if(!PCU_Comm_Self())
      printf("Start mesh mover\n");
    pMeshMover mmover = MeshMover_new(ppm, 0);

    // mesh motion of parametric deforming model region
    for (std::list<rigidBodyMotion>::iterator mit = mm.rigidBodyMotions.begin(); mit != mm.rigidBodyMotions.end(); ++mit) {
      assert(modelRegion = (pGRegion) GM_entityByTag(model, 3, mit->tag));
      assert(!GEN_isDiscreteEntity(modelRegion));
      printf("set rigid body motion: region %d\n",mit->tag);
      MeshMover_setTransform(mmover, modelRegion, mit->trans, mit->rotaxis, mit->rotpt, mit->rotang, mit->scale);
    }

    std::list<int>::iterator lit;
    // mesh motion of parametric surrounding model region
    for (lit = mm.parSurRegions.begin(); lit != mm.parSurRegions.end(); ++lit) {
      assert(modelRegion = (pGRegion) GM_entityByTag(model, 3, *lit));
      assert(!GEN_isDiscreteEntity(modelRegion));
      printf("set move on parametric surrounding region %d\n",*lit);
      vIter = M_classifiedVertexIter(pm, modelRegion, 0);
      while((meshVertex = VIter_next(vIter))){
        apf::MeshEntity* vtx = reinterpret_cast<apf::MeshEntity*>(meshVertex);
        apf::getComponents(f, vtx, 0, vals);
        const double newloc[3] = {vals[0], vals[1], vals[2]};
        MeshMover_setVolumeMove(mmover,meshVertex,newloc);
      }
      VIter_delete(vIter);
    }

    // mesh motion of parametric surrounding model face
    for (lit = mm.parSurFaces.begin(); lit != mm.parSurFaces.end(); ++lit) {
      assert(modelFace = (pGFace) GM_entityByTag(model, 2, *lit));
      assert(!GEN_isDiscreteEntity(modelFace));
      printf("set move on parametric face %d\n",*lit);
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

    // mesh motion of parametric surrounding model edge
    for (lit = mm.parSurEdges.begin(); lit != mm.parSurEdges.end(); ++lit) {
      assert(modelEdge = (pGEdge) GM_entityByTag(model, 1, *lit));
      assert(!GEN_isDiscreteEntity(modelEdge));
      printf("set move on parametric edge %d\n",*lit);
      vIter = M_classifiedVertexIter(pm, modelEdge, 1);
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

    // mesh motion of discrete deforming model region
    for (lit = mm.disDefRegions.begin(); lit != mm.disDefRegions.end(); ++lit) {
      assert(modelRegion = (pGRegion) GM_entityByTag(model, 3, *lit));
      assert(GEN_isDiscreteEntity(modelRegion));
      printf("set move on discrete deforming region %d\n",*lit);
      vIter = M_classifiedVertexIter(pm, modelRegion, 1);
      while((meshVertex = VIter_next(vIter))){
        apf::MeshEntity* vtx = reinterpret_cast<apf::MeshEntity*>(meshVertex);
        apf::getComponents(f, vtx, 0, vals);
        const double newloc[3] = {vals[0], vals[1], vals[2]};
        MeshMover_setDiscreteDeformMove(mmover,modelRegion,meshVertex,newloc);
      }
      VIter_delete(vIter);
    }

    // mesh motion of discrete surrounding model region
    for (lit = mm.disSurRegions.begin(); lit != mm.disSurRegions.end(); ++lit) {
      assert(modelRegion = (pGRegion) GM_entityByTag(model, 3, *lit));
      assert(GEN_isDiscreteEntity(modelRegion));
      printf("set move on discrete surrounding region %d\n",*lit);
      vIter = M_classifiedVertexIter(pm, modelRegion, 0);
      while((meshVertex = VIter_next(vIter))){
        apf::MeshEntity* vtx = reinterpret_cast<apf::MeshEntity*>(meshVertex);
        apf::getComponents(f, vtx, 0, vals);
        const double newloc[3] = {vals[0], vals[1], vals[2]};
        MeshMover_setVolumeMove(mmover,meshVertex,newloc);
      }
      VIter_delete(vIter);
    }

    // delete vals
    delete [] vals;

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
*/

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
          printf("set move on parametric region %d\n", GEN_tag(modelRegion));
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
          printf("set move on discrete region %d\n", GEN_tag(modelRegion));
          vIter = M_classifiedVertexIter(pm, modelRegion, 0);
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
        printf("skip rigid body motion: surface %d\n", GEN_tag(modelFace));
        continue;
      }
      else {
        if (!GEN_isDiscreteEntity(modelFace)) { // parametric
          printf("set move on parametric face %d\n", GEN_tag(modelFace));
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
        else { // discrete
          printf("set move on discrete face %d\n", GEN_tag(modelFace));
          vIter = M_classifiedVertexIter(pm, modelFace, 0);
          while((meshVertex = VIter_next(vIter))){
            apf::MeshEntity* vtx = reinterpret_cast<apf::MeshEntity*>(meshVertex);
            apf::getComponents(f, vtx, 0, &vals[0]);
            const double newloc[3] = {vals[0], vals[1], vals[2]};
            pPList closureRegion = GEN_regions(modelFace);
            assert(PList_size(closureRegion));
            modelRegion = (pGRegion) PList_item(closureRegion, 0);
            assert(GEN_isDiscreteEntity(modelRegion));
            MeshMover_setDiscreteDeformMove(mmover,modelRegion,meshVertex,newloc);
            PList_delete(closureRegion);
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
        printf("skip rigid body motion: edge %d\n", GEN_tag(modelEdge));
        continue;
      }
      else {
        if (!GEN_isDiscreteEntity(modelEdge)) { // parametric
          printf("set move on parametric edge %d\n", GEN_tag(modelEdge));
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
        else { // discrete
          printf("set move on discrete edge %d\n", GEN_tag(modelEdge));
          vIter = M_classifiedVertexIter(pm, modelEdge, 0);
          while((meshVertex = VIter_next(vIter))){
            apf::MeshEntity* vtx = reinterpret_cast<apf::MeshEntity*>(meshVertex);
            apf::getComponents(f, vtx, 0, &vals[0]);
            const double newloc[3] = {vals[0], vals[1], vals[2]};
            pPList closureRegion = GEN_regions(modelEdge);
            assert(PList_size(closureRegion));
            modelRegion = (pGRegion) PList_item(closureRegion, 0);
            assert(GEN_isDiscreteEntity(modelRegion));
            MeshMover_setDiscreteDeformMove(mmover,modelRegion,meshVertex,newloc);
            PList_delete(closureRegion);
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
        printf("skip rigid body motion: vertex %d\n", GEN_tag(modelVertex));
        continue;
      }
      else {
        if (!GEN_isDiscreteEntity(modelVertex)) { // parametric
          printf("set move on parametric vertex %d\n", GEN_tag(modelVertex));
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
        else { // discrete
          printf("set move on discrete vertex %d\n", GEN_tag(modelVertex));
          vIter = M_classifiedVertexIter(pm, modelRegion, 0);
          while((meshVertex = VIter_next(vIter))){
            apf::MeshEntity* vtx = reinterpret_cast<apf::MeshEntity*>(meshVertex);
            apf::getComponents(f, vtx, 0, &vals[0]);
            const double newloc[3] = {vals[0], vals[1], vals[2]};
            pPList closureRegion = GEN_regions(modelVertex);
            assert(PList_size(closureRegion));
            modelRegion = (pGRegion) PList_item(closureRegion, 0);
            assert(GEN_isDiscreteEntity(modelRegion));
            MeshMover_setDiscreteDeformMove(mmover,modelRegion,meshVertex,newloc);
            PList_delete(closureRegion);
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
      addImproverInMover(mmover, sim_fld_lst);
      addAdapterInMover(mmover, sim_fld_lst, m);
    }
    PList_clear(sim_fld_lst);
    PList_delete(sim_fld_lst);

    // do real work
    if(!PCU_Comm_Self())
      printf("do real mesh mover\n");
    assert(MeshMover_run(mmover, progress));
    MeshMover_delete(mmover);
//    MSA_delete(msa); // fix not unique id bug; but not proper
//    VolumeMeshImprover_delete(vmi); // // fix not unique id bug; but not proper

    // transfer sim fields to apf fields
    if (cooperation) {
      pc::transferSimFields(m);
    }

    // set rigid body total disp to be zero
    for (size_t i_rbpd = 0; i_rbpd < 3*in.nRigidBody; i_rbpd++)
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
//      meshMotion mm;
//      if (configureMotion(step, mm))
//        done = updateSIMCoord(in, m, cooperation, mm);
//      else
//        done = updateSIMDiscreteCoord(in, m, cooperation);
        done = updateSIMCoordAuto(in, m, cooperation);
//      else if (caseId == 10) // hack!
//        done = updateAndWriteSIMDiscrete(m); // hack!
//      else if (caseId == 20) // hack!
//        done = updateAndWriteSIMDiscreteModel(m); // hack!
    }
    else {
      done = updateAPFCoord(m);
    }
    assert(done);
  }

  void updateMesh(ph::Input& in, apf::Mesh2* m, apf::Field* szFld, int step, int cooperation) {
    if (in.simmetrixMesh && cooperation) {
      pc::runMeshMover(in,m,step,1);
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
