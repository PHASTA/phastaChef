#include "chefPhasta.h"
#include "sam.h"
#include "samSz.h"
#include <PCU.h>
#include <pcu_util.h>
#include <pcu_io.h>
#include <chef.h>
#include <spr.h>
#include <phasta.h>
#include "phIO.h"
#include <ph.h>
#include <phBC.h>
#include <phstream.h>
#include <sam.h>
#include <phastaChef.h>
#include <apfMDS.h>
#include <stdlib.h>
#include <unistd.h>

#include <apfSIM.h>
#include <gmi_sim.h>
#include <SimPartitionedMesh.h>
#include <MeshSimAdapt.h>
#include <SimField.h>
#include <SimAdvMeshing.h>
#include "SimMeshTools.h"
#include "SimParasolidKrnl.h"
#include <SimDiscrete.h>

#include "pcWriteFiles.h"
#include "pcUpdateMesh.h"
#include "pcAdapter.h"

#include <cstring>
#include <cassert>

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }

  /* may not need this */
  static FILE* openfile_read(ph::Input&, const char* path) {
    FILE* f = NULL;
    f = pcu_group_open(path, false);
    return f;
  }

  /* may not need this */
  static FILE* openstream_read(ph::Input& in, const char* path) {
    std::string fname(path);
    std::string restartStr("restart");
    FILE* f = NULL;
    if( fname.find(restartStr) != std::string::npos )
      f = openRStreamRead(in.rs);
    else {
      fprintf(stderr,
        "ERROR %s type of stream %s is unknown... exiting\n",
        __func__, fname.c_str());
      exit(1);
    }
    return f;
  }

  void setupChef(ph::Input& ctrl, int step) {
    PCU_ALWAYS_ASSERT(step > 0);
    //don't split or tetrahedronize
    ctrl.splitFactor = 1;
    ctrl.tetrahedronize = 0;
    ctrl.solutionMigration = 1;
    ctrl.adaptFlag = 0;
    ctrl.writeGeomBCFiles = 1;
  }

  static apf::Mesh2* loadMesh(gmi_model*& g, ph::Input& in) {
    const char* meshfile = in.meshFileName.c_str();
    if (ph::mesh_has_ext(meshfile, "sms")) {
      if (in.simmetrixMesh == 0) {
        if (PCU_Comm_Self()==0)
          fprintf(stderr, "oops, turn on flag: simmetrixMesh\n");
        in.simmetrixMesh = 1;
        in.filterMatches = 0; //not support
      }
    }
    return ph::loadMesh(g, meshfile);
  }

  void projectAndAttachFields(apf::Mesh2*& src_m, apf::Mesh2*& dst_m) {
    pProgress progress = Progress_new();
    Progress_setDefaultCallback(progress);

    // record num of parts
    int numParts = PCU_Comm_Peers();

    // load src model and mesh
    apf::MeshSIM* src_apf_msim = dynamic_cast<apf::MeshSIM*>(src_m);
    pParMesh src_ppm = src_apf_msim->getMesh();
    gmi_model* src_gmiModel = src_apf_msim->getModel();
    pGModel src_model = gmi_export_sim(src_gmiModel);
    PM_write(src_ppm, "src_mesh.sms", progress);

    pGRegion modelRegion;
    GRIter grIter;
    pGEntMeshMigrator gmig;
    // migrate all elements to rank 0
    gmig = GEntMeshMigrator_new(src_ppm, 3);
    grIter = GM_regionIter(src_model);
    while((modelRegion = GRIter_next(grIter))){
      int gid = PMU_gid(0, 0);
      GEntMeshMigrator_add(gmig, modelRegion, gid);
    }
    GRIter_delete(grIter);
    GEntMeshMigrator_run(gmig, progress);
    GEntMeshMigrator_delete(gmig);
    pMesh src_mesh = PM_mesh(src_ppm,0);

    // load dest model and mesh
    apf::MeshSIM* dst_apf_msim = dynamic_cast<apf::MeshSIM*>(dst_m);
    pParMesh dst_ppm = dst_apf_msim->getMesh();
    gmi_model* dst_gmiModel = dst_apf_msim->getModel();
    pGModel dst_model = gmi_export_sim(dst_gmiModel);

    // migrate all elements to rank 0
    gmig = GEntMeshMigrator_new(dst_ppm, 3);
    grIter = GM_regionIter(dst_model);
    while((modelRegion = GRIter_next(grIter))){
      int gid = PMU_gid(0, 0);
      GEntMeshMigrator_add(gmig, modelRegion, gid);
    }
    GRIter_delete(grIter);
    GEntMeshMigrator_run(gmig, progress);
    GEntMeshMigrator_delete(gmig);
    pMesh dst_mesh = PM_mesh(dst_ppm,0);

    // create mesh data for processing flag
    pMeshDataId mdid = MD_newMeshDataId("done_flag");

    // load fields
    phSolver::Input inp("solver.inp", "input.config");
    int num_flds = pc::getNumOfMappedFields(src_m);
    pField* src_flds = new pField[num_flds];
    pc::removeOtherFields(src_m, inp);
    int num_chck = pc::getSimFields(src_m, 1, src_flds, inp);
    PCU_ALWAYS_ASSERT(num_chck == num_flds);

    // create fields on destination mesh
    int valueType = 0;
    for(int i = 0; i < num_flds; i++) {
      if (Field_numComp(src_flds[i]) == 1)
        valueType = apf::SCALAR;
      else if (Field_numComp(src_flds[i]) == 3)
        valueType = apf::VECTOR;
      else if (Field_numComp(src_flds[i]) == 9)
        valueType = apf::MATRIX;
      else {
        printf("error: number of components is not correct!\n");
        assert(0);
      }
      if(dst_m->findField(Field_name(src_flds[i])))
        apf::destroyField(dst_m->findField(Field_name(src_flds[i])));
      if(!PCU_Comm_Self())
        printf("create a sim field %s on destination mesh\n", Field_name(src_flds[i]));
      apf::Field* rf = apf::createSIMFieldOn(dst_m, Field_name(src_flds[i]), valueType);
    }

    PCU_Barrier();

    // loop over model regions
    double xyz[3];
    pGRegion src_modelRegion;
    pGRegion dst_modelRegion;
    pRegion dst_meshRegion;
    pEdge dst_meshEdge;
    pVertex dst_meshVertex;
    pPList dvList = PList_new();
    pPList deList = PList_new();
    grIter = GM_regionIter(src_model);
    while((src_modelRegion = GRIter_next(grIter))){
    // add model region to domain
      pGDomain gdom = GDomain_new();
      GDomain_addModelEntity(gdom, src_modelRegion, 1);

    // start mesh region finder
      pMeshRegionFinder mrf = MeshRegionFinder_new(src_mesh, 1.0, gdom, progress);

    // we assume that the "same" model region has the same tag
      dst_modelRegion = (pGRegion) GM_entityByTag(dst_model, 3, GEN_tag(src_modelRegion));

    // loop over destination mesh regions
      RIter rIter = M_classifiedRegionIter(dst_mesh, dst_modelRegion);
      while((dst_meshRegion = RIter_next(rIter))){
        int useFirFlag = 1;
        double min = 0.0;
        deList = R_edges(dst_meshRegion, 1);
        void *eiter = 0;
        while((dst_meshEdge = (pEdge)PList_next(deList, &eiter))){
          double el = E_length(dst_meshEdge);
          if (useFirFlag) {
            min = el;
            useFirFlag = 0;
          }
          else {
            if (el < min)
              min = el;
          }
        }
        double tol = 0.5 * min; // hardcoding

    // loop over downward vertices
        dvList = R_vertices(dst_meshRegion, 1);
        void *viter = 0; // must initialize to 0
        while((dst_meshVertex = (pVertex)PList_next(dvList, &viter))){
          if (EN_getDataInt(dst_meshVertex, mdid, NULL))
            continue;
          V_coord(dst_meshVertex, xyz);
          const double loc[3] = {xyz[0], xyz[1], xyz[2]};

          double params[3] = {0.0, 0.0, 0.0};
          double cent[3] = {0.0, 0.0, 0.0};
          double dist[1] = {0.0};
          pRegion foundMR = MeshRegionFinder_find(mrf, loc, params, dist);
//          if(foundMR) {
//            EN_centroid(foundMR, cent);
//            printf("\nrank: %d; we found: (%f, %f, %f)\n", PCU_Comm_Self(), xyz[0], xyz[1], xyz[2]);
//            printf("  center of source region: (%f, %f, %f)\n", cent[0], cent[1], cent[2]);
//            printf("  parametric of point: (%f, %f, %f)\n", params[0], params[1], params[2]);
//            printf("  distance to domain: %12.16e; (tol = %12.16e)\n", dist[0], tol);
//          }
//          if(foundMR && (dist[0] < tol)) {
          if(foundMR) {
            // loop over fields
            for(int i = 0; i < num_flds; i++) {
            // get value from source mesh
              int numOfComp = Field_numComp(src_flds[i]);
              double* inVal = new double[numOfComp];
              pInterp intp = Field_entInterp(src_flds[i], foundMR);
              Interp_deriv0(intp, params, 0, inVal);

            // set value on destination mesh
              apf::Field* dst_fld = dst_m->findField(Field_name(src_flds[i]));
              double* outVal = new double[numOfComp];
              for (int j = 0; j < numOfComp; j++){
                outVal[j] = inVal[j];
              }
              apf::MeshEntity* vtx = reinterpret_cast<apf::MeshEntity*> (dst_meshVertex);
              apf::setComponents(dst_fld, vtx, 0, outVal);
            }
          // mark this destination vertex
            EN_attachDataInt(dst_meshVertex, mdid, 1);
          }
          else {
            if (foundMR == 0)
              printf("cannot find mesh region by point (%f,%f,%f)\n",loc[0],loc[1],loc[2]);
            else
              printf("dist %f is larger than tol %f\n",dist[0],tol);
          }
        }
        PList_delete(dvList);
      }
      RIter_delete(rIter);
      MeshRegionFinder_delete(mrf);
      GDomain_delete(gdom);
    }
    GRIter_delete(grIter);

    // delete mesh data for processing flag
    MD_deleteMeshDataId(mdid);

    // partition the dst mesh
    if(!PCU_Comm_Self())
      printf("start to partition the dst mesh\n");
    pPartitionOpts pOpts = PartitionOpts_new();
    PartitionOpts_setTotalNumParts(pOpts, numParts);
    PM_partition(dst_ppm, pOpts, progress);
    PartitionOpts_delete(pOpts);
    PM_write(dst_ppm, "dst_mesh.sms", progress);

    // transfer fields
    PCU_Barrier();
    printf("rank %d transfer sim fields to apf fields\n", PCU_Comm_Self());
    pc::transferSimFields(dst_m);

    Progress_delete(progress);
  }

} //end namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_Protect();
  if( argc != 3 ) {
    if(!PCU_Comm_Self())
      fprintf(stderr, "Usage: %s <dst_model.smd> <dst_mesh.sms> \n",argv[0]);
    exit(EXIT_FAILURE);
  }
  const char* attribFilename = argv[1];
  const char* meshFilename = argv[2];

  rstream rs = makeRStream();
  rstream dst_rs = makeRStream();
  grstream grs = makeGRStream();
  grstream dst_grs = makeGRStream();
  ph::Input ctrl;
  ph::Input dst_ctrl;
  ctrl.load("adapt.inp");
  dst_ctrl.load("adapt.inp");
  chefPhasta::initModelers(ctrl.writeSimLog);
  dst_ctrl.attributeFileName = attribFilename;
  dst_ctrl.meshFileName = meshFilename;
  dst_ctrl.solutionMigration = 0;

  /* load the model and mesh */
  gmi_model* dst_g = 0;
  apf::Mesh2* dst_m = 0;
  chef::cook(dst_g, dst_m, dst_ctrl, dst_grs); // used to load model and mesh
  dst_m->verify();

  gmi_model* g = 0;
  apf::Mesh2* m = 0;
  int step = ctrl.timeStepNumber;
  setupChef(ctrl, step);
  chef::cook(g, m, ctrl, grs); // used to load model and mesh
  m->verify();

  ctrl.rs = rs;
  dst_ctrl.rs = dst_rs;
  clearGRStream(grs);
  clearGRStream(dst_grs);

  /* transfer fields to mesh */
  chef::readAndAttachFields(ctrl,m);

  /* update model and write new model */
//  pc::runMeshMover(ctrl,m,step);
  pc::updateAPFCoord(ctrl,m);
  m->verify();

  /* project solution to new mesh */
  projectAndAttachFields(m, dst_m);
  printf("rank %d done with projectAndAttachFields\n", PCU_Comm_Self());

  /* write geombc and restart */
  dst_ctrl.solutionMigration = 1;
  chef::preprocess(dst_m, dst_ctrl, dst_grs);
  clearRStream(rs);
  clearRStream(dst_rs);

  destroyGRStream(grs);
  destroyGRStream(dst_grs);
  destroyRStream(rs);
  destroyRStream(dst_rs);
  freeMesh(m);
  chefPhasta::finalizeModelers(ctrl.writeSimLog);
  PCU_Comm_Free();
  MPI_Finalize();
}
