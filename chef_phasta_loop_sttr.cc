#include "chefPhasta.h"
#include "sam.h"
#include "samSz.h"
#include <PCU.h>
#include <chef.h>
#include <spr.h>
#include <phasta.h>
#include "phIO.h"
#include <phstream.h>
#include <sam.h>
#include <phastaChef.h>
#include <apfMDS.h>
#include <stdlib.h>
#include <assert.h>
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
#include <cstring>

#include "pcWriteFiles.h"
#include "pcUpdateMesh.h"
#include "pcAdapter.h"

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }

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
    //don't split or tetrahedronize
    ctrl.splitFactor = 1;
    ctrl.tetrahedronize = 0;
    ctrl.timeStepNumber = step;
    ctrl.solutionMigration = 1;
    if(step>0) {
      if(!PCU_Comm_Self()) {
        fprintf(stderr, "STATUS error based adapt %d\n", step);
        fprintf(stderr, "STATUS ctrl.attributeFileName %s step %d\n",
            ctrl.attributeFileName.c_str(), step);
      }
      ctrl.adaptStrategy = 7;
      ctrl.adaptFlag = 1;
      ctrl.writeGeomBCFiles = 1;
    }
  }

} //end namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_Protect();
  if( argc != 3 ) {
    if(!PCU_Comm_Self())
      fprintf(stderr, "Usage: %s <maxTimeStep> <motion case id>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  int maxStep = atoi(argv[1]);
  int caseId  = atoi(argv[2]);
  chefPhasta::initModelers();
  rstream rs = makeRStream();
  grstream grs = makeGRStream();
  ph::Input ctrl;
  ctrl.load("adapt.inp");
  /* load the model and mesh */
  gmi_model* g = 0;
  apf::Mesh2* m = 0;
  chef::cook(g, m, ctrl, rs, grs);
  /* load input file for solver */
  phSolver::Input inp("solver.inp", "input.config");
  int step = 0; int phtStep = 0; int seq  = 0;
  pc::writeSequence(m,seq,"test_"); seq++;
  do {
    m->verify();
    pass_info_to_phasta(m, ctrl);
    /* take the initial mesh as size field */
    apf::Field* szFld = samSz::isoSize(m);
    step = phasta(inp,grs,rs);
    ctrl.rs = rs;
    clearGRStream(grs);
    if(!PCU_Comm_Self())
      fprintf(stderr, "STATUS ran to step %d\n", step);
    if( step >= maxStep )
      break;
    setupChef(ctrl,step);
    chef::readAndAttachFields(ctrl,m);
    pc::updateMeshCoord(ctrl,m,step,caseId);
    m->verify();

    pc::writePHTfiles(phtStep, step-phtStep, PCU_Comm_Peers()); phtStep = step;
    pc::writeSequence(m,seq,"test_"); seq++;
    /* do mesh adaptation */
    pc::runMeshAdapter(ctrl,m,szFld,step);
    pc::writeSequence(m,seq,"test_"); seq++;

    chef::preprocess(m,ctrl,grs);
    clearRStream(rs);
  } while( step < maxStep );
  destroyGRStream(grs);
  destroyRStream(rs);
  freeMesh(m);
  chefPhasta::finalizeModelers();
  PCU_Comm_Free();
  MPI_Finalize();
}
