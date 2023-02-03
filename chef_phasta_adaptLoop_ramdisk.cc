#include <PCU.h>
#include <chef.h>
#include <phasta.h>
#include <phstream.h>
#include <iostream>
#include <sstream>
#include "chefPhasta.h"
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>

/** \file chef_phasta_adaptLoop_ramdisk.cc
    \brief Example file-based driver for adaptive loops using a ramdisk
    \remark Runs Chef and then PHASTA until the user-specified maximum
            PHASTA time step is reached.
*/

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }

  void mychdir(const char* path) {
    const int fail = chdir(path);
    if(fail) {
      fprintf(stderr, "ERROR failed to change to %s dir... exiting\n", path);
      exit(1);
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
  }

  std::string makeMeshName(int step) {
    std::stringstream meshname;
    meshname  << "bz2:" << "t" << step << "p" << PCU_Comm_Peers() << "_.smb";
    return meshname.str();
  }

  std::string makeRestartName() {
    std::stringstream restartname;
    restartname << PCU_Comm_Peers() << "-procs_case/restart";
    return restartname.str();
  }

  void setupChef(ph::Input& ctrl, int step) {
    //don't split or tetrahedronize
    ctrl.splitFactor = 1;
    ctrl.tetrahedronize = 0;
    ctrl.timeStepNumber = step;
    ctrl.solutionMigration = 1;
    if(step>1) {
      ctrl.adaptStrategy = 2; //error field adapt
      ctrl.adaptErrorThreshold = 1e-2;
      ctrl.adaptErrorFieldName = "errors";
      ctrl.adaptErrorFieldIndex = 0;
      ctrl.adaptFlag = 1;
    }
    ctrl.outMeshFileName = makeMeshName(step);
    ctrl.restartFileName = makeRestartName();
  }

}
int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_Protect();
  if( argc != 6 ) {
    if(!PCU_Comm_Self())
      fprintf(stderr, "Usage: %s <maxTimeStep> <ramdisk path> <solver.inp> <input.config> <adapt.inp>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  int maxStep = atoi(argv[1]);
  const char* ramdisk = argv[2];
  const char* solverinp = argv[3];
  const char* inputcfg = argv[4];
  const char* adaptinp = argv[5];
  gmi_model* g = 0;
  apf::Mesh2* m = 0;

  mychdir(ramdisk);
  ph::Input ctrl;
  ctrl.load(adaptinp);
  ctrl.outMeshFileName = makeMeshName(0);
  chefPhasta::initModelers();
  chef::cook(g,m,ctrl);
  freeMesh(m); m = NULL;
  phSolver::Input inp(solverinp,inputcfg);
  int step = 0;
  do {
    ctrl.meshFileName = makeMeshName(step);
    step = phasta(inp);
    assert(step >= 0);
    if(!PCU_Comm_Self())
      fprintf(stderr, "CAKE ran to step %d\n", step);
    setupChef(ctrl,step);
    chef::cook(g,m,ctrl);
    freeMesh(m); m = NULL;
    mychdir(step);
  } while( step < maxStep );
  chefPhasta::finalizeModelers();
  PCU_Comm_Free();
  MPI_Finalize();
}
