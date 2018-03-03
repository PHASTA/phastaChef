#ifndef PC_UPDATEMESH_H
#define PC_UPDATEMESH_H

#include <apf.h>
#include <apfMesh2.h>
#include <apfSIM.h>
#include <apfMDS.h>
#include <chef.h>

namespace pc {
  void runMeshMover(ph::Input& in, apf::Mesh2* m, int step, int caseId);

  void updateMesh(ph::Input& in, apf::Mesh2* m, apf::Field* szFld, int step, int caseId, int cooperation = 1);
}

#endif
