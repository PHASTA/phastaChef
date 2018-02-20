#ifndef CHEF_PHASTA_H 
#define CHEF_PHASTA_H

#include <gmi_mesh.h>
#ifdef GMI_SIM_FOUND
#include <gmi_sim.h>
#include <SimUtil.h>
#include <SimPartitionedMesh.h>
#include <SimField.h>
#include <SimAdvMeshing.h>
#include <SimDiscrete.h>
#include "SimMeshTools.h"
#include "SimParasolidKrnl.h"
#include <MeshSim.h>
#include <MeshSimAdapt.h>
namespace chefPhasta {
  void initModelers() {
    MS_init();
    SimModel_start();
    Sim_readLicenseFile(0);
    SimPartitionedMesh_start(0, 0);
    Sim_logOn("phastaChef.log");
    SimParasolid_start(1);
    SimField_start();
    SimAdvMeshing_start();
    SimMeshTools_start();
    SimDiscrete_start(0);
    gmi_sim_start();
    gmi_register_sim();
    gmi_register_mesh();
  }
  void finalizeModelers() {
    SimDiscrete_stop(0);
    SimMeshTools_stop();
    SimAdvMeshing_stop();
    gmi_sim_stop();
    SimField_stop();
    Sim_logOff();
    SimPartitionedMesh_stop();
    SimParasolid_stop(1);
    Sim_unregisterAllKeys();
    SimModel_stop();
    MS_exit();
  }
}
#else
namespace chefPhasta {
  void initModelers() {
    gmi_register_mesh();
  }
  void finalizeModelers() {
  }
}
#endif

#endif
