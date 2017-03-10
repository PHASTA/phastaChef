#include "chefPhasta.h"
#include <PCU.h>
#include <pumi_version.h>
#include <chef.h>
#include <phasta.h>
#include <phstream.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <unistd.h>

/** \file chef_phasta_sam_adaptLoop_sz_ramdisk.cc
    \brief Example POSIX file-based driver using a ramdisk for adaptive loops
    \remark Runs Chef and then PHASTA until the user-specified maximum
            PHASTA time step is reached.  Size fields to drive adaptation
            are defined using SAM (from
            <a href=https://github.com/SCOREC/core>core</a>) in Chef
            by reading the PHASTA "errors" field.
*/

namespace {
  void pwd() {
    if(!PCU_Comm_Self()) {
      char path[4096];
      getcwd(path,4096);
      fprintf(stderr, "STATUS changed to %s dir\n", path);
    }
  }

  void mychdir(const char* path) {
    const int fail = chdir(path);
    if(fail) {
      fprintf(stderr, "ERROR failed to change to %s dir... exiting\n", path);
      exit(1);
    }
  }

  void mychdir(int step) {
    std::stringstream path;
    path << step;
    string s = path.str();
    const int fail = chdir(s.c_str());
    if(fail) {
      fprintf(stderr, "ERROR failed to change to %d dir... exiting\n", step);
      exit(1);
    }
    pwd();
  }

  std::string makeMeshName(int step) {
    std::stringstream meshname;
    meshname  << "bz2:" << "t" << step << "p" << PCU_Comm_Peers() << "_.smb";
    return meshname.str();
  }

  std::string makeRestartName() {
    std::stringstream restartname;
    restartname << PCU_Comm_Peers() << "-procs_case/restart";
    return restartname.str();
  }

  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }


  void setupChef(ph::Input& ctrl, int step) {
    //don't split or tetrahedronize
    ctrl.printIOtime = 1; //report time spent streaming
    ctrl.splitFactor = 1;
    ctrl.tetrahedronize = 0;
    ctrl.timeStepNumber = step;
    ctrl.solutionMigration = 1;
    ctrl.adaptStrategy = 1; //error field adapt
    ctrl.adaptFlag = 1;
    ctrl.restartFileName = makeRestartName();
    ctrl.outMeshFileName = makeMeshName(step);
    if(!PCU_Comm_Self()) {
      fprintf(stderr, "STATUS error based adapt %d\n", step);
      fprintf(stderr, "STATUS ctrl.attributeFileName %s step %d\n",
          ctrl.attributeFileName.c_str(), step);
    }
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_Protect();
  if( argc != 6 ) {
    if(!PCU_Comm_Self())
      fprintf(stderr, "Usage: %s <maxTimeStep> <ramdisk path> <solver.inp> <input.config> <adapt.inp>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  if( !PCU_Comm_Self() )
    printf("PUMI Git hash %s\n", pumi_version());
  int maxStep = atoi(argv[1]);
  const char* ramdisk = argv[2];
  const char* solverinp = argv[3];
  const char* inputcfg = argv[4];
  const char* chefinp = argv[5];

  double t0 = PCU_Time();

  mychdir(ramdisk);

  chefPhasta::initModelers();
  /* read chef config */
  ph::Input ctrl;
  ctrl.load(chefinp);
  /* read restart files (and split if requested) */
  ctrl.outMeshFileName = makeMeshName(ctrl.timeStepNumber);
  gmi_model* g = NULL;
  apf::Mesh2* m = NULL;
  chef::cook(g,m,ctrl);
  freeMesh(m); m = NULL;
  phSolver::Input inp(solverinp,inputcfg);
  int step = ctrl.timeStepNumber;
  do {
    double cycleStart = PCU_Time();
    pwd();
    /* The next Chef run needs to load the mesh from 
     * the solve that is about to run.*/
    ctrl.meshFileName = makeMeshName(step);
    step = phasta(inp);
    if(!PCU_Comm_Self())
      fprintf(stderr, "STATUS ran to step %d\n", step);
    if( step >= maxStep )
      break;
    setupChef(ctrl,step);
    chef::cook(g,m,ctrl);
    freeMesh(m); m = NULL;
    if(ctrl.adaptFlag) mychdir(step);
    if(!PCU_Comm_Self())
      fprintf(stderr, "cycle time %f seconds\n", PCU_Time()-cycleStart);
  } while( step < maxStep );
  chefPhasta::finalizeModelers();
  if(!PCU_Comm_Self())
    fprintf(stderr, "STATUS done %f seconds\n", PCU_Time()-t0);
  PCU_Comm_Free();
  MPI_Finalize();
}
