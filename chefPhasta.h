#ifndef CHEF_PHASTA_H 
#define CHEF_PHASTA_H

#include <gmi_mesh.h>
#ifdef GMI_SIM_FOUND
#include <gmi_sim.h>
#include <SimUtil.h>
namespace chefPhasta {
  void initModelers() {
    Sim_readLicenseFile(0);
    gmi_sim_start();
    gmi_register_sim();
    gmi_register_mesh();
  }
  void finalizeModelers() {
    gmi_sim_stop();
    Sim_unregisterAllKeys();
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
