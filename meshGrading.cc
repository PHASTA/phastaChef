#include "chefPhasta.h"
#include "sam.h"
#include "samSz.h"
#include <PCU.h>
#include <pcu_util.h>
#include <pcu_io.h>
#include <chef.h>
#include <spr.h>
#include <phasta.h>
#include "phIO.h"
#include <ph.h>
#include <phBC.h>
#include <phstream.h>
#include <sam.h>
#include <phastaChef.h>
#include <apfMDS.h>
#include <stdlib.h>
#include <unistd.h>
#include "lionPrint.h"

#include <apfSIM.h>
#include <gmi_sim.h>
#include <SimPartitionedMesh.h>
#include <MeshSimAdapt.h>
#include <SimField.h>
#include <SimAdvMeshing.h>
#include "SimMeshTools.h"
#include "SimParasolidKrnl.h"
#include <SimDiscrete.h>

#include "pcSmooth.h"

#include <cstring>
#include <cassert>

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }
} //end namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_Protect();
  lion_set_verbosity(1);
  if( argc != 5 ) {
    if(!PCU_Comm_Self())
      fprintf(stderr, "Usage: %s <model.smd> <mesh.sms> <graded_mesh.sms> <gradation factor>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  const char* modelFilename = argv[1];
  const char* meshFilename = argv[2];
  const char* outputFilename = argv[3];
  double gradingFactor = atof(argv[4]);
  chefPhasta::initModelers();
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  // load mesh and model
  pGModel model = GM_load(modelFilename, NULL, progress);
  pParMesh pmesh = PM_load(meshFilename, model, progress);
  apf::Mesh2* m = apf::createMesh(pmesh);

  // apply mesh gradation
  pc::meshGradation(m, gradingFactor);

  // write out mesh
  PM_write(pmesh, outputFilename, progress);

  freeMesh(m);
  Progress_delete(progress);
  chefPhasta::finalizeModelers();
  PCU_Comm_Free();
  MPI_Finalize();
}
