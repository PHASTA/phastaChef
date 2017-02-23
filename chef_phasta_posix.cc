#include <PCU.h>
#include <chef.h>
#include <phasta.h>
#include "chefPhasta.h"

/** \file chef_phasta_posix.cc
    \brief Example file-based driver for running Chef then PHASTA
    \remark Runs Chef and then PHASTA using POSIX files to transfer
            mesh and field information.
*/

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
  chef::cook(g,m);
  phasta(argc,argv);
  freeMesh(m);
  chefPhasta::finalizeModelers();
  PCU_Comm_Free();
  MPI_Finalize();
}
