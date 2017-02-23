#ifndef CHEF_PHASTA_H 
#define CHEF_PHASTA_H

/** \file chefPhasta.h
    \brief The Chef interface */

#include <gmi_mesh.h>
#ifdef GMI_SIM_FOUND
#include <gmi_sim.h>
#include <SimUtil.h>
namespace chefPhasta {
  /* \brief initialize the Simmetrix interfaces to geometric modelers */
  void initModelers() {
    Sim_readLicenseFile(0);
    gmi_sim_start();
    gmi_register_sim();
    gmi_register_mesh();
  }
  /* \brief finalize the Simmetrix interfaces to geometric modelers */
  void finalizeModelers() {
    gmi_sim_stop();
    Sim_unregisterAllKeys();
  }
}
#else
namespace chefPhasta {
  /* \brief initialize the SCOREC interfaces to geometric modelers */
  void initModelers() {
    gmi_register_mesh();
  }
  /* \brief empty place holder function */
  void finalizeModelers() {
  }
}
#endif

#endif
