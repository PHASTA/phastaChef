#ifndef PC_WRITEFILES_H
#define PC_WRITEFILES_H

#include <phasta.h>
#include <apf.h>
#include <apfMesh2.h>
#include <chef.h>
#include <SimPartitionedMesh.h>
#include "SimModel.h"

namespace pc {
  void writeSequence(apf::Mesh2* m, int step, const char* filename);

  void writeSIMModel (pGModel model, int step, const char* filename);

  void writeSIMMesh (pParMesh mesh, int step, const char* filename);

  void writePHTfiles (int old_step, int cur_step, phSolver::Input& inp);
}

#endif
