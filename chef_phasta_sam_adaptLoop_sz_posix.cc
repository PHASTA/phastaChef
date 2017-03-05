#include "chefPhasta.h"
#include <PCU.h>
#include <chef.h>
#include <phasta.h>
#include <phstream.h>
#include <sam.h>
#include <apfMDS.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

/** \file chef_phasta_sam_adaptLoop_sz_posix.cc
    \brief Example POSIX file-based driver for adaptive loops
    \remark Runs Chef and then PHASTA until the user-specified maximum
            PHASTA time step is reached.  Size fields to drive adaptation
            are defined using SAM (from
            <a href=https://github.com/SCOREC/core>core</a>)
            by reading the PHASTA "errors" field.
            This example also demonstrates the use of the fine grained
            chef.h and phasta.h APIs.
*/

namespace {
  void pwd() {
    if(!PCU_Comm_Self()) {
      char path[4096];
      getcwd(path,4096);
      fprintf(stderr, "STATUS changed to %s dir\n", path);
    }
  }

  void mychdir(int step) {
    std::stringstream path;
    path << step;
    string s = path.str();
    const int fail = chdir(s.c_str());
    if(fail) {
      fprintf(stderr, "ERROR failed to change to %d dir... exiting\n", step);
      exit(1);
    }
    pwd();
  }

  std::string makeMeshName(int step) {
    std::stringstream meshname;
    meshname  << "bz2:" << "t" << step << "p" << PCU_Comm_Peers() << "/";
    return meshname.str();
  }

  std::string makeRestartName() {
    std::stringstream restartname;
    restartname << PCU_Comm_Peers() << "-procs_case/restart";
    return restartname.str();
  }

  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }

  static FILE* openfile_read(ph::Input&, const char* path) {
    return fopen(path, "r");
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
    ctrl.restartFileName = makeRestartName();
    ctrl.outMeshFileName = makeMeshName(step);
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
  int maxStep = atoi(argv[1]);
  const char* chefinp = argv[2];

  chefPhasta::initModelers();
  /* read chef config */
  ph::Input ctrl;
  ctrl.load(chefinp);
  /* load the model and mesh */
  gmi_model* mdl = gmi_load(ctrl.modelFileName.c_str());
  apf::Mesh2* m = NULL;
  /* read restart files (and split if requested) */
  ctrl.outMeshFileName = makeMeshName(ctrl.timeStepNumber);
  chef::cook(mdl,m,ctrl);
  freeMesh(m); m = NULL;
  phSolver::Input inp("solver.inp", "input.config");
  int step = ctrl.timeStepNumber;
  do {
    pwd();
    /* The next Chef run needs to load the mesh from 
     * the solve that is about to run.*/
    ctrl.meshFileName = makeMeshName(step);
    step = phasta(inp);
    if(!PCU_Comm_Self())
      fprintf(stderr, "STATUS ran to step %d\n", step);
    if( step >= maxStep )
      break;
    setupChef(ctrl,step);
    m = apf::loadMdsMesh(ctrl.modelFileName.c_str(),
                         ctrl.meshFileName.c_str());
    chef::cook(mdl,m,ctrl);
    freeMesh(m); m = NULL;
    mychdir(step);
  } while( step < maxStep );
  chefPhasta::finalizeModelers();
  if(!PCU_Comm_Self())
    fprintf(stderr, "STATUS done");
  PCU_Comm_Free();
  MPI_Finalize();
}
