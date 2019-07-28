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

#include "pcWriteFiles.h"
#include "pcUpdateMesh.h"
#include "pcAdapter.h"

#include <cstring>
#include <cassert>

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }

  /* may not need this */
  static FILE* openfile_read(ph::Input&, const char* path) {
    FILE* f = NULL;
    f = pcu_group_open(path, false);
    return f;
  }

  /* may not need this */
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

  /* the following size field is hardcoded
     for the 6-grain problem */
  apf::Field* getPreSz(apf::Mesh* m) {
    double isoSize;
    double size[1];
    double anisosize[3][3];
    double center[3] = {0.0, 0.0, 0.0};
    double r = 0.0;
    double R = 0.456;
    apf::Field* szFld = createFieldOn(m, "isoSize", apf::SCALAR);;
    apf::MeshEntity* v;
    apf::MeshIterator* vit = m->begin(0);
    while ((v = m->iterate(vit))) {
      apf::Vector3 p;
      m->getPoint(v, 0, p);
      pVertex meshVertex = reinterpret_cast<pVertex>(v);
      pGEntity clsModel = EN_whatIn(meshVertex);
      int dim = GEN_type(clsModel);
      int tag = GEN_tag(clsModel);
      if ( EN_isBLEntity(meshVertex) ||
           dim == Gvertex ||
           dim == Gedge ||
           dim == Gface ||
          (dim == Gregion && tag == 40 ) ) {
        int szType = V_size(meshVertex, size, anisosize);
        assert(szType == 1);
        if (size[0] > 0.5) {
          printf("size is too large: %f on model (%d, %d)\n", size[0], dim, tag);
          size[0] = size[0] / 2.0;
        }
        isoSize = size[0];
      }
      else {
        assert(dim == Gregion);
        switch(tag) {
          case 110:
            center[0] = -2.6;
            center[1] = -0.8;
            break;
          case 163:
            center[0] =  2.6;
            center[1] = -0.8;
            break;
          case 223:
            center[0] =  0.0;
            center[1] = -0.8;
            break;
          case 290:
            center[0] = -2.6;
            center[1] =  0.8;
            break;
          case 364:
            center[0] =  2.6;
            center[1] =  0.8;
            break;
          case 445:
            center[0] =  0.0;
            center[1] =  0.8;
            break;
          default:
            printf("not quite as the plan. tag = %d\n",tag);
            assert(0);
        }
        if( p[0] <= (center[0] - 0.5) ) {
          r = sqrt( (p[0]-(center[0]-0.5))*(p[0]-(center[0]-0.5)) +
                    (p[1]-center[1])*(p[1]-center[1]) +
                    (p[2]-center[2])*(p[2]-center[2]) );
        }
        else if ( p[0] >= (center[0] + 0.5) ) {
          r = sqrt( (p[0]-(center[0]+0.5))*(p[0]-(center[0]+0.5)) +
                    (p[1]-center[1])*(p[1]-center[1]) +
                    (p[2]-center[2])*(p[2]-center[2]) );
        }
        else {
          r = sqrt( (p[1]-center[1])*(p[1]-center[1]) +
                    (p[2]-center[2])*(p[2]-center[2]) );
        }
        isoSize = 0.2 - 0.15 * r / R;
      }
      apf::setScalar(szFld,v,0,isoSize);
    }
    m->end(vit);
    return szFld;
  }

} //end namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_Protect();
  lion_set_verbosity(1);
  if( argc != 2 ) {
    if(!PCU_Comm_Self())
      fprintf(stderr, "Usage: %s <mode id>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  int modeId  = atoi(argv[1]);
  rstream rs = makeRStream();
  grstream grs = makeGRStream();
  ph::Input ctrl;
  ctrl.load("adapt.inp");
  chefPhasta::initModelers(ctrl.writeSimLog);

  /* load the model and mesh */
  gmi_model* g = 0;
  apf::Mesh2* m = 0;
  int step = ctrl.timeStepNumber;
  setupChef(ctrl, step);
  chef::cook(g, m, ctrl, grs);
  m->verify();

  /* take the initial mesh as size field */
  apf::Field* szFld = samSz::isoSize(m);

  /* use the prescribed size field */
//  apf::Field* szFld = getPreSz(m);

  ctrl.rs = rs;
  clearGRStream(grs);

  /* transfer fields to mesh */
  chef::readAndAttachFields(ctrl,m);

  /* update model and write new model */
  if(modeId == 0) {
    pc::updateMesh(ctrl,m,szFld,step);
    /* write geombc and restart files */
    ctrl.writeRestartFiles = 1;
    chef::preprocess(m,ctrl,grs);
  }
  else if(modeId == 1) {
    pc::updateAndWriteSIMDiscreteCoord(m);
  }
  else if(modeId == 2) {
    pc::updateAndWriteSIMDiscreteField(m);
  }

  clearRStream(rs);
  destroyGRStream(grs);
  destroyRStream(rs);
  freeMesh(m);
  chefPhasta::finalizeModelers(ctrl.writeSimLog);
  PCU_Comm_Free();
  MPI_Finalize();
}
