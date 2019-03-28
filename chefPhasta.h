#ifndef CHEF_PHASTA_H 
#define CHEF_PHASTA_H

#include <PCU.h>
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
  void initModelers(int simLog = 0) {
    if (simLog)
      Sim_logOn("phastaChef.log");
    MS_init();
    SimModel_start();
    Sim_readLicenseFile(0);
    SimPartitionedMesh_start(0, 0);
    SimParasolid_start(1);
    SimField_start();
    SimAdvMeshing_start();
    SimMeshTools_start();
    SimDiscrete_start(0);
    gmi_sim_start();
    gmi_register_sim();
    gmi_register_mesh();
    if(!PCU_Comm_Self())
      printf("SimModSuite Version: %s\n",Sim_buildID());
  }
  void finalizeModelers(int simLog = 0) {
    SimDiscrete_stop(0);
    SimMeshTools_stop();
    SimAdvMeshing_stop();
    gmi_sim_stop();
    SimField_stop();
    SimPartitionedMesh_stop();
    SimParasolid_stop(1);
    Sim_unregisterAllKeys();
    SimModel_stop();
    MS_exit();
    if (simLog)
      Sim_logOff();
  }
}
#else
namespace chefPhasta {
  void initModelers(int simLog = 0) {
    (void) simLog;
    gmi_register_mesh();
  }
  void finalizeModelers(int simLog = 0) {
    (void) simLog;
  }
}
#endif

#endif
