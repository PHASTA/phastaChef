#ifndef PC_ERROR_H
#define PC_ERROR_H

#include "pcWriteFiles.h"
#include <SimField.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfSIM.h>
#include <apfMDS.h>
#include <chef.h>
#include <phasta.h>

namespace pc {
  double getShortestEdgeLength(apf::Mesh* m, apf::MeshEntity* elm);

  void attachVMSSizeField(apf::Mesh2*& m, ph::Input& in, phSolver::Input& inp);
}

#endif
