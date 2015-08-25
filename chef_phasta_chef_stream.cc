#include <PCU.h>
#include <chef.h>
#include <gmi_mesh.h>
#include <phasta.h>
#include <phstream.h>

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
  gmi_register_mesh();
  gmi_model* g = 0;
  apf::Mesh2* m = 0;
  grstream grs = makeGRStream();
  chef::cook(g,m,"adapt.inp",grs);
  rstream rs = makeRStream();
  phasta(argc,argv,grs,rs);
  chef::cook(g,m,"adaptParma.inp",rs);
  destroyGRStream(grs);
  destroyRStream(rs);
  freeMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
