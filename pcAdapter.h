#ifndef PC_ADAPTER_H
#define PC_ADAPTER_H

#include <apf.h>
#include <apfMesh2.h>
#include <apfSIM.h>
#include <apfMDS.h>
#include <chef.h>

namespace pc {
  void runMeshAdapter(ph::Input& in, apf::Mesh2*& m, apf::Field*& orgSF, int step);
}

#endif
