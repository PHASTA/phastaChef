#include <PCU.h>
#include <chef.h>
#include <phasta.h>
#include <phstream.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include "chefPhasta.h"

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }
}

void setupChef(ph::Input& ctrl, int step) {
  //don't split or tetrahedronize
  ctrl.splitFactor = 1;
  ctrl.recursivePtn = 0;
  ctrl.tetrahedronize = 0;
  ctrl.timeStepNumber = step;
  ctrl.solutionMigration = 1;
  if(!PCU_Comm_Self()) {
    fprintf(stderr, "CAKE enabling UR at step %d\n", step);
    fprintf(stderr, "CAKE ctrl.attributeFileName %s step %d\n",
        ctrl.attributeFileName.c_str(), step);
  }
  if(4==step) {
    ctrl.adaptStrategy = 7; //UR
    ctrl.adaptFlag = 1;
  } else {
    ctrl.adaptFlag = 0;
  }
  std::stringstream meshname;
  meshname  << "bz2:t" << step << "p" << PCU_Comm_Peers() << "/";
  ctrl.outMeshFileName = meshname.str();
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
  grstream grs = makeGRStream();
  ph::Input ctrl;
  ctrl.load("adapt.inp");
  chef::cook(g,m,ctrl,grs);
  rstream rs = makeRStream();
  phSolver::Input inp("solver.inp", "input.config");
  int step = 0;
  do { 
    step = phasta(inp,grs,rs);
    clearGRStream(grs);
    if(!PCU_Comm_Self())
      fprintf(stderr, "CAKE ran to step %d\n", step);
    setupChef(ctrl,step);
    chef::cook(g,m,ctrl,rs,grs);
    clearRStream(rs);
  } while( step < maxStep );
  destroyGRStream(grs);
  destroyRStream(rs);
  freeMesh(m);
  chefPhasta::finalizeModelers();
  PCU_Comm_Free();
  MPI_Finalize();
}
