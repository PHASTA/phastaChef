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
#include <phiotimer.h>
#include <phstream.h>
#include <sam.h>
#include <phastaChef.h>
#include <apfMDS.h>
#include <stdlib.h>
#include <unistd.h>

#include <apfSIM.h>
#include <gmi_sim.h>
#include <SimPartitionedMesh.h>
#include <MeshSimAdapt.h>
#include <SimField.h>
#include <SimAdvMeshing.h>
#include "SimMeshTools.h"
#include "SimParasolidKrnl.h"
#include <SimDiscrete.h>
#include <cstring>

#include "pcWriteFiles.h"
#include "pcUpdateMesh.h"
#include "pcAdapter.h"

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }

  /* may not need this */
  static FILE* openfile_read(ph::Input&, const char* path) {
    FILE* f = NULL;
    PHASTAIO_OPENTIME(f = pcu_group_open(path, false);)
    return f;
  }

  /* may not need this */
  static FILE* openstream_read(ph::Input& in, const char* path) {
    std::string fname(path);
    std::string restartStr("restart");
    FILE* f = NULL;
    if( fname.find(restartStr) != std::string::npos )
      PHASTAIO_OPENTIME(f = openRStreamRead(in.rs);)
    else {
      fprintf(stderr,
        "ERROR %s type of stream %s is unknown... exiting\n",
        __func__, fname.c_str());
      exit(1);
    }
    return f;
  }

  void setupChef(ph::Input& ctrl, int step) {
    PCU_ALWAYS_ASSERT(step > 0);
    //don't split or tetrahedronize
    ctrl.splitFactor = 1;
    ctrl.tetrahedronize = 0;
    ctrl.solutionMigration = 1;
    ctrl.adaptFlag = 0;
    ctrl.writeGeomBCFiles = 1;
  }

  static apf::Mesh2* loadMesh(gmi_model*& g, ph::Input& in) {
    const char* meshfile = in.meshFileName.c_str();
    if (ph::mesh_has_ext(meshfile, "sms")) {
      if (in.simmetrixMesh == 0) {
        if (PCU_Comm_Self()==0)
          fprintf(stderr, "oops, turn on flag: simmetrixMesh\n");
        in.simmetrixMesh = 1;
        in.filterMatches = 0; //not support
      }
    }
    return ph::loadMesh(g, meshfile);
  }

} //end namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_Protect();
  if( argc != 2 ) {
    if(!PCU_Comm_Self())
      fprintf(stderr, "Usage: %s <motion case id>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  int caseId  = atoi(argv[1]);
  chefPhasta::initModelers();
  rstream rs = makeRStream();
  grstream grs = makeGRStream();
  ph::Input ctrl;
  ctrl.load("adapt.inp");

  /* load the model and mesh */
  gmi_model* g = 0;
  apf::Mesh2* m = 0;
  int step = ctrl.timeStepNumber;
  setupChef(ctrl, step);
  chef::cook(g, m, ctrl, grs);
  m->verify();

  /* take the initial mesh as size field */
  apf::Field* szFld = samSz::isoSize(m);
  ctrl.rs = rs;
  clearGRStream(grs);

  /* transfer fields to mesh */
  chef::readAndAttachFields(ctrl,m);

  /* update model and write new model */
  pc::updateMeshCoord(ctrl,m,step,caseId);
  m->verify();

  /* do mesh adaptation and write new mesh */
  pc::runMeshAdapter(ctrl,m,szFld,step);

  /* write geombc and restart */
  chef::preprocess(m,ctrl,grs);
  clearRStream(rs);

  destroyGRStream(grs);
  destroyRStream(rs);
  freeMesh(m);
  chefPhasta::finalizeModelers();
  PCU_Comm_Free();
  MPI_Finalize();
}
