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

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
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
    meshname  << "bz2:" << "t" << step << "p" << PCU_Comm_Peers() << "/";
    return meshname.str();
  }

  std::string makeRestartName() {
    std::stringstream restartname;
    restartname << PCU_Comm_Peers() << "-procs_case/restart";
    return restartname.str();
  }

  std::string prefixCwd(std::string name) {
    char cwd[4096] = "\0";
    getcwd(cwd,4096);
    std::stringstream s;
    s << cwd << "/" << name;
    return s.str();
  }

  void setupChef(ph::Input& ctrl, int step) {
    //don't split or tetrahedronize
    ctrl.splitFactor = 1;
    ctrl.recursivePtn = 0;
    ctrl.tetrahedronize = 0;
    ctrl.timeStepNumber = step;
    ctrl.solutionMigration = 1;
    if(step>1) {
      if(!PCU_Comm_Self()) {
        fprintf(stderr, "CAKE error based adapt %d\n", step);
        fprintf(stderr, "CAKE ctrl.attributeFileName %s step %d\n",
            ctrl.attributeFileName.c_str(), step);
      }
      ctrl.adaptStrategy = 1; //error field adapt
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
  if( argc != 2 ) {
    if(!PCU_Comm_Self())
      fprintf(stderr, "Usage: %s <maxTimeStep>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  int maxStep = atoi(argv[1]);
  chefPhasta::initModelers();
  gmi_model* g = 0;
  apf::Mesh2* m = 0;
  ph::Input ctrl;
  ctrl.load("adapt.inp");
  ctrl.modelFileName = prefixCwd(ctrl.modelFileName);
  ctrl.attributeFileName = prefixCwd(ctrl.attributeFileName);
  ctrl.outMeshFileName = makeMeshName(0);
  chef::cook(g,m,ctrl);
  freeMesh(m); m = NULL;
  phSolver::Input inp("solver.inp", "input.config");
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
