#include <PCU.h>
#include <chef.h>
#include <phasta.h>
#include <phstream.h>
#include <iostream>
#include <sstream>
#include "chefPhasta.h"

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }
}

void setupChefParma(ph::Input& ctrl, int step) {
  //don't split or tetrahedronize
  ctrl.splitFactor = 1;
  ctrl.recursivePtn = 0;
  ctrl.tetrahedronize = 0;
  ctrl.timeStepNumber = step;
  ctrl.solutionMigration = 1;
  std::stringstream meshname;
  meshname  << "bz2:t" << step << "p" << PCU_Comm_Peers() << "/";
  ctrl.outMeshFileName = meshname.str();
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_Protect();
  chefPhasta::initModelers();
  gmi_model* g = 0;
  apf::Mesh2* m = 0;
  grstream grs = makeGRStream();
  ph::Input ctrl;
  ctrl.load("adapt.inp");
  chef::cook(g,m,ctrl,grs);
  rstream rs = makeRStream();
  phSolver::Input inp("solver.inp", "input.config");
  int step = phasta(inp,grs,rs);
  if(!PCU_Comm_Self())
    fprintf(stderr, "CAKE ran to step %d\n", step);
  setupChefParma(ctrl,step);
  chef::cook(g,m,ctrl,rs);
  destroyGRStream(grs);
  destroyRStream(rs);
  freeMesh(m);
  chefPhasta::finalizeModelers();
  PCU_Comm_Free();
  MPI_Finalize();
}
