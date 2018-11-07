#include "chefPhasta.h"
#include <PCU.h>
#include <lionPrint.h>
#include <pumi_version.h>
#include <chef.h>
#include <phasta.h>
#include <phstream.h>

#include <sstream>
#include <iostream>
#include <fstream>
//#include <cstring>

#include <apfMDS.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

/** \file chef_phasta_sam_adaptLoop_sz.cc
    \brief Example in-memory driver for adaptive loops
    \remark Runs Chef and then PHASTA until the user-specified maximum
            PHASTA time step is reached.  Size fields to drive adaptation
            are defined using SAM (from
            <a href=https://github.com/SCOREC/core>core</a>)
            by reading the PHASTA "errors" field.
*/

namespace {
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
    if(step>1) {
      if(!PCU_Comm_Self()) {
        fprintf(stdout, "STATUS error based adapt %d\n", step);
        fprintf(stdout, "STATUS ctrl.attributeFileName %s step %d\n",
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
  lion_set_verbosity(1);
  if( argc != 3 ) {
    if(!PCU_Comm_Self())
      fprintf(stdout, "Usage: %s <nAdaptCycle> <chef input config>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  if( !PCU_Comm_Self() )
    printf("PUMI Git hash %s\n", pumi_version());
  int nAdaptCycle = atoi(argv[1]);
  const char* chefinp = argv[2];

  double t0 = PCU_Time();

  chefPhasta::initModelers();
  grstream grs = makeGRStream();
  ph::Input ctrl;
  ctrl.load(chefinp);
  /* load the model and mesh */
  gmi_model* mdl = gmi_load(ctrl.modelFileName.c_str());
  apf::Mesh2* m = NULL;

// override what was in inp file for starting mesh and restart so that chain jobs continue
  double lastTime = PCU_Time();
  double newTime; 
  int Rstep,Mstep;
  int step = ctrl.timeStepNumber;
  std::string numAdaptMesh="numAdaptMesh.dat"; // this will need to be written after a successful adapt but assumes
                                               // we have at least one restart that we adapted to.
  if(step !=0) {
    std::string numRestart="numRestart.dat";   // if you link to <numprocs>-procs_case/numstart.dat no need to write this
    std::ifstream numR(numRestart.c_str());
    numR >> Rstep;
    step=Rstep;  // otherwise the code uses non updated step from inp file
    numR.close();

    std::ifstream numM(numAdaptMesh.c_str());
    numM >> Mstep;
    numM.close();
  
    ctrl.timeStepNumber=Rstep;
    std::stringstream ss;
    ss << "bz2:" << Mstep << "/mdsMesh_bz2/";
    ctrl.meshFileName = ss.str();
    std::string printme=ss.str();
    if(!PCU_Comm_Self()) fprintf(stdout,"meshFileName to read is %s \n", ctrl.meshFileName.c_str());
    if(nAdaptCycle == 1) {
      ctrl.outMeshFileName="bz2:mdsMesh_bz2/"; // this will get the step number inserted within Chef
      ctrl.writeGeomBCFiles= 1; // Assuming we will need this viz and above for adapt 
    }
  }

  newTime = PCU_Time();
  if(!PCU_Comm_Self())
      fprintf(stdout, "prep for pre-cook %f seconds\n", newTime-lastTime);
  lastTime=newTime;

  int iAdaptCycle=0;
  /* read restart files (and split and adapt if requested) */
  chef::cook(mdl,m,ctrl,grs);

  newTime = PCU_Time();
  if(!PCU_Comm_Self())
      fprintf(stdout, "First cook %f seconds\n", newTime-lastTime);
  lastTime=newTime;

  if(!PCU_Comm_Self()) {   // only rank 0 needs to write this as otherwise file system will hate you (and go slow)
    if(ctrl.adaptFlag ==1) { // there was an adapt at start so write the step number to a file in case it is one cycle run
      std::ofstream numM(numAdaptMesh.c_str());
      numM << step;
      numM.close();
    }
  }
  assert(m);
  rstream rs = makeRStream();

  newTime = PCU_Time();
  if(!PCU_Comm_Self())
      fprintf(stdout, "First Stream Prep %f seconds\n", newTime-lastTime);
  lastTime=newTime;

  phSolver::Input inp("solver.inp", "input.config");

  newTime = PCU_Time();
  if(!PCU_Comm_Self())
      fprintf(stdout, "First inp %f seconds\n", newTime-lastTime);
  lastTime=newTime;

  do {
    double cycleStart = PCU_Time();
    step = phasta(inp,grs,rs);

  newTime = PCU_Time();
  if(!PCU_Comm_Self())
      fprintf(stdout, " phasta Run  %f seconds\n", newTime-lastTime);
  lastTime=newTime;

    iAdaptCycle++;
    clearGRStream(grs);
    if(!PCU_Comm_Self())
      fprintf(stdout, "STATUS ran to step %d\n", step);
    if( iAdaptCycle >= nAdaptCycle )
      break;
    if(nAdaptCycle - iAdaptCycle == 1) { // this will be the last adapt so turn on mesh write flags to enable chain
      ctrl.outMeshFileName="bz2:mdsMesh_bz2/"; // this will get the step number inserted within Chef
      ctrl.writeGeomBCFiles= 1; // Assuming we will need this viz and above for adapt 
    }

  newTime = PCU_Time();
  if(!PCU_Comm_Self())
      fprintf(stdout, " Prep before setupChef  %f seconds\n", newTime-lastTime);
  lastTime=newTime;

    setupChef(ctrl,step);

  newTime = PCU_Time();
  if(!PCU_Comm_Self())
      fprintf(stdout, " setupChef  %f seconds\n", newTime-lastTime);
  lastTime=newTime;

    chef::cook(mdl,m,ctrl,rs,grs);

  newTime = PCU_Time();
  if(!PCU_Comm_Self())
      fprintf(stdout, " loop cook  %f seconds\n", newTime-lastTime);
  lastTime=newTime;


// write the step number to a file
  if(!PCU_Comm_Self()) {   // only rank 0 needs to write this as otherwise file system will hate you (and go slow)
    std::ofstream numM(numAdaptMesh.c_str());
    numM << step;
    numM.close();
  }
    clearRStream(rs);
    if(!PCU_Comm_Self())
      fprintf(stdout, "STATUS cycle time %f seconds\n", PCU_Time()-cycleStart);

  newTime = PCU_Time();
  if(!PCU_Comm_Self())
      fprintf(stdout, " Stream end of while  %f seconds\n", newTime-lastTime);
  lastTime=newTime;

  } while( iAdaptCycle < nAdaptCycle ); // could rewrite this as a for loop but ok for now to leave as is
  destroyGRStream(grs);
  destroyRStream(rs);
  freeMesh(m);
  chefPhasta::finalizeModelers();
  if(!PCU_Comm_Self())
    fprintf(stdout, "STATUS done %f seconds\n", PCU_Time()-t0);
  PCU_Comm_Free();
  MPI_Finalize();
}
