#include "chefPhasta.h"
#include <PCU.h>
#include <chef.h>
#include <phasta.h>
#include <phastaIO.h>
#include <phstream.h>
#include <sam.h>
#include <apfMDS.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }

  apf::Field* getField(apf::Mesh* m) {
    /* Hijacking threshold and converting it to same as
     * chefs runFromGivenSize in core:phasta/phAdapt.cc
     * 6th entry (idx 5) is an isotropic size.
     * copied from Chef RunFromGivenSize */
    const unsigned idx = 5;
    const char* fldName = "errors";
    return sam::specifiedIso(m,fldName,idx);
  }

  static FILE* openfile_read(ph::Input&, const char* path) {
    return fopen(path, "r");
  }

  static FILE* openstream_read(ph::Input& in, const char* path) {
    std::string fname(path);
    std::string restartStr("restart");
    FILE* f = NULL;
    if( fname.find(restartStr) != std::string::npos )
      f = openRStreamRead(in.rs);
    else {
      fprintf(stderr,
        "ERROR %s type of stream %s is unknown... exiting\n",
        __func__, fname.c_str());
      exit(1);
    }
    return f;
  }

  void setupChef(ph::Input& ctrl, int step) {
    //don't split or tetrahedronize
    ctrl.splitFactor = 1;
    ctrl.tetrahedronize = 0;
    ctrl.timeStepNumber = step;
    ctrl.solutionMigration = 1;
    if(step>1) {
      if(!PCU_Comm_Self()) {
        fprintf(stderr, "STATUS error based adapt %d\n", step);
        fprintf(stderr, "STATUS ctrl.attributeFileName %s step %d\n",
            ctrl.attributeFileName.c_str(), step);
      }
      ctrl.adaptStrategy = 1; //error field adapt
      ctrl.adaptFlag = 1;
    }
  }
}
int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_Protect();
  if( argc != 2 ) {
    if(!PCU_Comm_Self())
      fprintf(stderr, "Usage: %s <maxTimeStep>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  int maxStep = atoi(argv[1]);
  chefPhasta::initModelers();
  grstream grs = makeGRStream();
  ph::Input ctrl;
  ctrl.load("samAdaptLoop.inp");
  /* setup file reading */
  ctrl.openfile_read = openfile_read;
  /* load the model and mesh */
  gmi_model* mdl = gmi_load(ctrl.modelFileName.c_str());
  apf::Mesh2* m = apf::loadMdsMesh(mdl,ctrl.meshFileName.c_str());
  /* read restart files (and split if requested) */
  chef::cook(mdl,m,ctrl,grs);
  rstream rs = makeRStream();
  /* setup stream reading */
  ctrl.openfile_read = openstream_read;
  ctrl.rs = rs;
  phSolver::Input inp("solver.inp", "input.config");
  int step = 0;
  do {
    step = phasta(inp,grs,rs);
    clearGRStream(grs);
    if(!PCU_Comm_Self())
      fprintf(stderr, "STATUS ran to step %d\n", step);
    if( step >= maxStep )
      break;
    if( step >= 325 ) {
      int istop=1;
     }
      
    setupChef(ctrl,step);
    chef::readAndAttachFields(ctrl,m);
    apf::Field* szFld = getField(m);
    assert(szFld);
    chef::adapt(m,szFld);
    apf::destroyField(szFld);
    chef::balanceAndReorder(ctrl,m);
    chef::preprocess(m,ctrl,grs);
    clearRStream(rs);
    FILE* f = openGRStream(grs,"restart");
    int params[128];
    readHeader(f,"solution",params,3,"binary");
    int size = params[0]*params[1];
    double* data = (double*) calloc(size,sizeof(double));
    readDataBlock(f,data,size,"double","binary");
    free(data);
    fclose(f); 
  } while( step < maxStep );
  destroyGRStream(grs);
  destroyRStream(rs);
  freeMesh(m);
  chefPhasta::finalizeModelers();
  PCU_Comm_Free();
  MPI_Finalize();
  fprintf(stderr, "STATUS done");
}
