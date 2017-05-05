#include "chefPhasta.h"
#include "sam.h"
#include "samSz.h"
#include <PCU.h>
#include <chef.h>
#include <spr.h>
#include <phasta.h>
#include "phIO.h"
#include <phstream.h>
#include <sam.h>
#include <apfMDS.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

#include <apfSIM.h>
#include <gmi_sim.h>
#include <SimPartitionedMesh.h>
#include <MeshSimAdapt.h>
#include <SimField.h>
#include <cstring>

#include "pcWriteFiles.h"
#include "pcUpdateMesh.h"
#include "pcAdapter.h"

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }

  apf::Field* refineProjSF(apf::Mesh2* m, apf::Field* orgSF, int step) {
    /* define geom of projectile */
    double orgCenter = 0.435;
    double acc = 60000.0;
    double dt = 0.00001;
    double curCenter = orgCenter + 0.5*acc*(step-1.0)*(step-1.0)*dt*dt;

    /* define imaginary cylinder and refine size */
    double box[] = {0.85, 0.0, 0.0, 0.855, 0.065};
    double ref = 0.006;
    double cor = 0.02;

    /* define size field based on current center */
    apf::Field* newSz = apf::createFieldOn(m,"refineProjSF",apf::SCALAR);
    apf::Vector3 points;
    double h = 0.0;
    double f = 0.0;
    apf::MeshEntity* vtx;
    apf::MeshIterator* itr = m->begin(0);
    while( (vtx = m->iterate(itr)) ) {
      m->getPoint(vtx, 0, points);
      if ( fabs(points[0]- box[0]) < box[3] &&
           sqrt((points[1]-box[1])*(points[1]-box[1]) +
                (points[2]-box[2])*(points[2]-box[2])) < box[4]){
        if ( points[0]- curCenter >= -0.235 && points[0]- curCenter <= 0.0 ) // near proj rear part
        {
          h = ref;
        }
        else if( points[0]- curCenter <= 0.235 && points[0]- curCenter > 0.0 ) // near proj front part
        {
          h = ref+0.002;
        }
        else if( points[0]- curCenter < -0.235 ) // rear part
        {
          f = ((curCenter-0.235)-points[0])/(curCenter-0.235+0.005);
          h = ref * (1-f) + cor * f;
        }
        else if( points[0]- curCenter >  0.235 )  // front part
        {
          f = (points[0]-(curCenter+0.235))/(1.705-(curCenter+0.235));
          h = (ref+0.004) * (1-f) + cor * f;
        }
        else
          printf("surprise! we should not fall into here\n");
      }
      else {
        h = apf::getScalar(orgSF,vtx,0);
        if (h < cor) h = cor;
      }
      apf::setScalar(newSz,vtx,0,h);
    }
    m->end(itr);
    apf::destroyField(orgSF);
    return newSz;
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
    if(step>0) {
      if(!PCU_Comm_Self()) {
        fprintf(stderr, "STATUS error based adapt %d\n", step);
        fprintf(stderr, "STATUS ctrl.attributeFileName %s step %d\n",
            ctrl.attributeFileName.c_str(), step);
      }
      ctrl.adaptStrategy = 1; //error field adapt
      ctrl.adaptFlag = 1;
      ctrl.writeGeomBCFiles = 0;
    }
  }

  int isMeshqGood(apf::Mesh* m, double crtn) {
    apf::Field* meshq = m->findField("meshQ");
    if (!meshq) {
      if (!PCU_Comm_Self())
        fprintf(stderr, "Not find meshQ field.\n");
      assert(meshq);
    }
    int meshGood = 1;
    apf::MeshEntity* elm;
    apf::MeshIterator* itr = m->begin(m->getDimension());
    while( (elm = m->iterate(itr)) ) {
      if (apf::getScalar(meshq, elm, 0) < crtn) {
        meshGood = 0;
        break;
      }
    }
    m->end(itr);
    PCU_Barrier();
    if (PCU_Min_Int(meshGood) && !PCU_Comm_Self())
      printf("Mesh is Good. No need for adaptation!\n");
    return PCU_Min_Int(meshGood);
  }

} //end namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_Protect();
//debugging for simmetrix mesh
  Sim_readLicenseFile(0);
  SimPartitionedMesh_start(0, 0);
  Sim_logOn("loopDriver.log");
  SimUtil_start();
  SimModel_start();
  SimField_start();
  gmi_sim_start();
  gmi_register_sim();
  MS_init();
//end debugging
  if( argc != 2 ) {
    if(!PCU_Comm_Self())
      fprintf(stderr, "Usage: %s <maxTimeStep>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  int maxStep = atoi(argv[1]);
  chefPhasta::initModelers();
  rstream rs = makeRStream();
  grstream grs = makeGRStream();
  ph::Input ctrl;
  ctrl.load("samAdaptLoop.inp");
  /* load the model and mesh */
  gmi_register_mesh();
  gmi_model* g = 0;
  apf::Mesh2* m = 0;
  chef::cook(g, m, ctrl, rs, grs);
  /* load input file for solver */
  phSolver::Input inp("solver.inp", "input.config");
  int step = 0; int phtStep = 0; int seq  = 0;
  pc::writeSequence(m,seq,"test_"); seq++;
  do {
    m->verify();
    /* take the initial mesh as size field */
    apf::Field* szFld = samSz::isoSize(m);
    step = phasta(inp,grs,rs);
    ctrl.rs = rs;
    clearGRStream(grs);
    if(!PCU_Comm_Self())
      fprintf(stderr, "STATUS ran to step %d\n", step);
    if( step >= maxStep )
      break;
    setupChef(ctrl,step);
    chef::readAndAttachFields(ctrl,m);
    pc::updateMeshCoord(ctrl,m);
//    int doAdaptation = !isMeshqGood(m, ctrl.meshqCrtn);
    int doAdaptation = 1; // make it run anyway

    m->verify();
    if ( doAdaptation ) {
      pc::writePHTfiles(phtStep, step-phtStep, PCU_Comm_Peers()); phtStep = step;
      pc::writeSequence(m,seq,"test_"); seq++;
      /* do mesh adaptation */
      pc::runMeshAdapter(ctrl,m,szFld,step);
      pc::writeSequence(m,seq,"test_"); seq++;
    }
    chef::preprocess(m,ctrl,grs);
    clearRStream(rs);
  } while( step < maxStep );
  destroyGRStream(grs);
  destroyRStream(rs);
  freeMesh(m);
  chefPhasta::finalizeModelers();
//debugging for simmetrix mesh
  gmi_sim_stop();
  SimField_stop();
  SimModel_stop();
  SimUtil_stop();
  Sim_logOff();
  SimPartitionedMesh_stop();
  Sim_unregisterAllKeys();
//end debugging
  PCU_Comm_Free();
  MPI_Finalize();
}
