#ifndef PC_UPDATEMESH_H
#define PC_UPDATEMESH_H

#include <apf.h>
#include <apfMesh2.h>
#include <apfSIM.h>
#include <apfMDS.h>
#include <chef.h>
#include <list>
#include <cstring>
#include <cstdlib>

namespace pc {

  struct rigidBodyMotion {
    rigidBodyMotion(int t = 0, double r = 0.0, double s = 0.0)
    {
      tag = t;
      memset(trans, 0.0, sizeof trans);
      memset(rotaxis, 0.0, sizeof rotaxis);
      memset(rotpt, 0.0, sizeof rotpt);
      rotang = r;
      scale = s;
    }
    int tag;
    double trans[3];
    double rotaxis[3];
    double rotpt[3];
    double rotang;
    double scale;
    void set_trans(double x, double y, double z){
      trans[0] = x; trans[1] = y; trans[2] = z;
    }
    void set_rotaxis(double x, double y, double z){
      rotaxis[0] = x; rotaxis[1] = y; rotaxis[2] = z;
    }
    void set_rotpt(double x, double y, double z){
      rotpt[0] = x; rotpt[1] = y; rotpt[2] = z;
    }
  };

  struct meshMotion {
    int caseId;  // may not be used
    std::list<rigidBodyMotion> rigidBodyMotions;
    std::list<int> disDefRegions;
    std::list<int> disSurRegions;
    std::list<int> parSurRegions;
    std::list<int> parSurFaces;
    std::list<int> parSurEdges;
  };

  meshMotion getTDMeshMotion(int caseId, int step);

  bool updateAPFCoord(ph::Input& in, apf::Mesh2* m);

  bool updateAndWriteSIMDiscreteCoord(apf::Mesh2* m);

  bool updateAndWriteSIMDiscreteField(apf::Mesh2* m);

  void runMeshMover(ph::Input& in, apf::Mesh2* m, int step, int cooperation = 0);

  void updateMesh(ph::Input& in, apf::Mesh2* m, apf::Field* szFld, int step, int cooperation = 1);

  void balanceEqualWeights(pParMesh pmesh, pProgress progress);

// hardcoding {
  void prescribe_proj_mesh_size(
    apf::Mesh2* m, apf::Field* sizes, double proj_disp);
// hardcoding }

}

#endif
