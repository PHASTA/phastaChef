#include "phasta.h"
#include "chef.h"

int main(int argc, char** argv) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  assert(provided == MPI_THREAD_MULTIPLE);
  PCU_Comm_Init();
  PCU_Protect();
  chef(argc,argv);
  phasta();
  PCU_Comm_Free();
  MPI_Finalize();
}




