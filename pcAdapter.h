#ifndef PC_ADAPTER_H
#define PC_ADAPTER_H

#include <SimField.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfSIM.h>
#include <apfMDS.h>
#include <chef.h>
#include <phasta.h>

namespace pc {
  int getNumOfMappedFields(phSolver::Input& inp);

  void removeOtherFields(apf::Mesh2*& m, phSolver::Input& inp);

  int getSimFields(apf::Mesh2*& m, int simFlag, pField* sim_flds, phSolver::Input& inp);

  pPList getSimFieldList(ph::Input& in, apf::Mesh2*& m);

  void transferSimFields(apf::Mesh2*& m);

  void runMeshAdapter(ph::Input& in, apf::Mesh2*& m, apf::Field*& orgSF, int step);
}

#endif
