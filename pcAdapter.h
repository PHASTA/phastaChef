#ifndef PC_ADAPTER_H
#define PC_ADAPTER_H

#include <SimField.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfSIM.h>
#include <apfMDS.h>
#include <chef.h>

namespace pc {
  int getSimFields(apf::Mesh2*& m, int simFlag, pField* sim_flds);

  pPList getSimFieldList(ph::Input& in, apf::Mesh2*& m);

  void transferSimFields(apf::Mesh2*& m);

  void runMeshAdapter(ph::Input& in, apf::Mesh2*& m, apf::Field*& orgSF, int step);
}

#endif
