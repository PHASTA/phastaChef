#include <PCU.h>
#include <chef.h>
#include <phasta.h>
#include <phstream.h>
#include "chefPhasta.h"

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_Protect();
  chefPhasta::initModelers();
  gmi_model* g = 0;
  apf::Mesh2* m = 0;
  GRStream* grs = makeGRStream();
  ph::Input ctrl;
  ctrl.load("adapt.inp");
  chef::cook(g,m,ctrl,grs);
  phSolver::Input inp("solver.inp", "input.config");
  int step = phasta(inp,grs);
  if(!PCU_Comm_Self())
    fprintf(stderr, "CAKE ran to step %d\n", step);
  destroyGRStream(grs);
  freeMesh(m);
  chefPhasta::finalizeModelers();
  PCU_Comm_Free();
  MPI_Finalize();
}
