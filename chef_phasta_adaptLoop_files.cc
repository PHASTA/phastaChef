#include <PCU.h>
#include <chef.h>
#include <phasta.h>
#include <phstream.h>
#include <iostream>
#include <sstream>
#include "chefPhasta.h"
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
    //std::stringstream meshname;
    //meshname  << "bz2:t" << step << "p" << PCU_Comm_Peers() << "/";
    //ctrl.outMeshFileName = meshname.str();
  }
}
int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  fprintf(stdout, "----you made it----\n");
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
  chef::cook(g,m,ctrl);
  phSolver::Input inp("solver.inp", "input.config");
  int step = 0;
  do {
    step = phasta(inp);
    if(!PCU_Comm_Self())
      fprintf(stderr, "CAKE ran to step %d\n", step);
    setupChef(ctrl,step);
    chef::cook(g,m,ctrl);
    mychdir(step);
  } while( step < maxStep );
  freeMesh(m);
  chefPhasta::finalizeModelers();
  PCU_Comm_Free();
  MPI_Finalize();
}
