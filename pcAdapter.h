#ifndef PC_ADAPTER_H
#define PC_ADAPTER_H

#include "pcWriteFiles.h"
#include <SimField.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfSIM.h>
#include <apfMDS.h>
#include <chef.h>
#include <phasta.h>
#include <MeshSimAdapt.h>

namespace pc {

  void attachMeshSizeField(apf::Mesh2*& m, ph::Input& in, phSolver::Input& inp);

  int getNumOfMappedFields(phSolver::Input& inp);

  void removeOtherFields(apf::Mesh2*& m, phSolver::Input& inp);

  int getSimFields(apf::Mesh2*& m, int simFlag, pField* sim_flds, phSolver::Input& inp);

  pPList getSimFieldList(ph::Input& in, apf::Mesh2*& m);

  void attachMinSizeFlagField(apf::Mesh2*& m, ph::Input& in);

  void transferSimFields(apf::Mesh2*& m);

  void attachCurrentSizeField(apf::Mesh2*& m);

  void setupSimImprover(pVolumeMeshImprover vmi, pPList sim_fld_lst);

  void setupSimAdapter(pMSAdapt adapter, ph::Input& in, apf::Mesh2*& m, pPList& sim_fld_lst);

  void runMeshAdapter(ph::Input& in, apf::Mesh2*& m, apf::Field*& orgSF, int step);
}

#endif
