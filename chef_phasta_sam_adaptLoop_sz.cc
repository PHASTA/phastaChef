#include "chefPhasta.h"
#include <PCU.h>
#include <pumi_version.h>
#include <chef.h>
#include <phasta.h>
#include <phstream.h>
#include <apfMDS.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

/** \file chef_phasta_sam_adaptLoop_sz.cc
    \brief Example in-memory driver for adaptive loops
    \remark Runs Chef and then PHASTA until the user-specified maximum
            PHASTA time step is reached.  Size fields to drive adaptation
            are defined using SAM (from
            <a href=https://github.com/SCOREC/core>core</a>)
            by reading the PHASTA "errors" field.
*/

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }

  void setupChef(ph::Input& ctrl, int step) {
    //don't split or tetrahedronize
    ctrl.printIOtime = 1; //report time spent streaming
    ctrl.splitFactor = 1;
    ctrl.tetrahedronize = 0;
    ctrl.timeStepNumber = step;
    ctrl.solutionMigration = 1;
    if(step>1) {
      if(!PCU_Comm_Self()) {
        fprintf(stderr, "STATUS error based adapt %d\n", step);
        fprintf(stderr, "STATUS ctrl.attributeFileName %s step %d\n",
            ctrl.attributeFileName.c_str(), step);
      }
      ctrl.adaptStrategy = 1; //error field adapt
      ctrl.adaptFlag = 1;
    }
  }
}
int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_Protect();
  if( argc != 3 ) {
    if(!PCU_Comm_Self())
      fprintf(stderr, "Usage: %s <maxTimeStep> <chef input config>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  if( !PCU_Comm_Self() )
    printf("PUMI Git hash %s\n", pumi_version());
  int maxStep = atoi(argv[1]);
  const char* chefinp = argv[2];

  double t0 = PCU_Time();

  chefPhasta::initModelers();
  grstream grs = makeGRStream();
  ph::Input ctrl;
  ctrl.load(chefinp);
  /* load the model and mesh */
  gmi_model* mdl = gmi_load(ctrl.modelFileName.c_str());
  apf::Mesh2* m = NULL;
  /* read restart files (and split if requested) */
  chef::cook(mdl,m,ctrl,grs);
  assert(m);
  rstream rs = makeRStream();
  phSolver::Input inp("solver.inp", "input.config");
  int step = 0;
  do {
    double cycleStart = PCU_Time();
    step = phasta(inp,grs,rs);
    clearGRStream(grs);
    if(!PCU_Comm_Self())
      fprintf(stderr, "STATUS ran to step %d\n", step);
    if( step >= maxStep )
      break;
    setupChef(ctrl,step);
    chef::cook(mdl,m,ctrl,rs,grs);
    clearRStream(rs);
    if(!PCU_Comm_Self())
      fprintf(stderr, "STATUS cycle time %f seconds\n", PCU_Time()-cycleStart);
  } while( step < maxStep );
  destroyGRStream(grs);
  destroyRStream(rs);
  freeMesh(m);
  chefPhasta::finalizeModelers();
  if(!PCU_Comm_Self())
    fprintf(stderr, "STATUS done %f seconds\n", PCU_Time()-t0);
  PCU_Comm_Free();
  MPI_Finalize();
}
