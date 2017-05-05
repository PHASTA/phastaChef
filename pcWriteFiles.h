#ifndef PC_WRITEFILES_H
#define PC_WRITEFILES_H

#include <apf.h>
#include <apfMesh2.h>

namespace pc {
  void writeSequence(apf::Mesh2* m, int step, const char* filename);

  void writePHTfiles(int step, int nstep, int nproc);
}

#endif
