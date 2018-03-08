#include "pcUpdateMesh.h"
#include <PCU.h>
#include <cassert>

namespace pc {
  /* ideally, this comes from some input file */
  meshMotion getTDMeshMotion(int caseId, int step) {
    double offset, disp, rang, sfct;
    offset = 0.0;
    disp = 2e-4 + 2e-5;
    rang = 2.0;
    sfct = 0.9;
    if (caseId == 1) {
      offset = 2e-4 * (double)(step - 1) + 2e-5;
      printf("current step is %d; offset = %f\n",step,offset);
    }
    else if (caseId == 2) {
      offset = 4e-4;
    }
    else {
      printf("wrong case id\n");
      assert(0);
    }
    meshMotion mm;
    rigidBodyMotion mbm;
    mm.parSurFaces.push_back(54);
    mm.parSurFaces.push_back(46);
    mm.parSurFaces.push_back(41);
    mm.parSurRegions.push_back(1);
// grain11
    mbm = rigidBodyMotion(1339, rang, sfct);
    mbm.set_trans(disp, 0.0, 0.0);
    mbm.set_rotaxis(0.0, 0.0, 1.0);
    mbm.set_rotpt(0.5e-3 + offset, 1.125e-3, 0.0);
    mm.rigidBodyMotions.push_back(mbm);
// grain12
    mbm = rigidBodyMotion(1036, rang, sfct);
    mbm.set_trans(disp, 0.0, 0.0);
    mbm.set_rotaxis(0.0, 0.0, 1.0);
    mbm.set_rotpt(2.625e-3 + offset, 1.125e-3, 0.0);
    mm.rigidBodyMotions.push_back(mbm);
// grain13
    mbm = rigidBodyMotion(733, rang, sfct);
    mbm.set_trans(disp, 0.0, 0.0);
    mbm.set_rotaxis(0.0, 0.0, 1.0);
    mbm.set_rotpt(4.75e-3 + offset, 1.125e-3, 0.0);
    mm.rigidBodyMotions.push_back(mbm);
// grain21
    mbm = rigidBodyMotion(1440, 0.0, sfct);
    mbm.set_trans(disp, 0.0, 0.0);
    mbm.set_rotaxis(0.0, 0.0, 1.0);
    mbm.set_rotpt(0.5e-3 + offset, 0.0, 0.0);
    mm.rigidBodyMotions.push_back(mbm);
// grain22
    mbm = rigidBodyMotion(1137, 0.0, sfct);
    mbm.set_trans(disp, 0.0, 0.0);
    mbm.set_rotaxis(0.0, 0.0, 1.0);
    mbm.set_rotpt(2.625e-3 + offset, 0.0, 0.0);
    mm.rigidBodyMotions.push_back(mbm);
// grain23
    mbm = rigidBodyMotion(834, 0.0, sfct);
    mbm.set_trans(disp, 0.0, 0.0);
    mbm.set_rotaxis(0.0, 0.0, 1.0);
    mbm.set_rotpt(4.75e-3 + offset, 0.0, 0.0);
    mm.rigidBodyMotions.push_back(mbm);
// grain31
    mbm = rigidBodyMotion(1238, rang, sfct);
    mbm.set_trans(disp, 0.0, 0.0);
    mbm.set_rotaxis(0.0, 0.0, -1.0);
    mbm.set_rotpt(0.5e-3 + offset, -1.125e-3, 0.0);
    mm.rigidBodyMotions.push_back(mbm);
// grain32
    mbm = rigidBodyMotion(935, rang, sfct);
    mbm.set_trans(disp, 0.0, 0.0);
    mbm.set_rotaxis(0.0, 0.0, -1.0);
    mbm.set_rotpt(2.625e-3 + offset, -1.125e-3, 0.0);
    mm.rigidBodyMotions.push_back(mbm);
// grain33
    mbm = rigidBodyMotion(632, rang, sfct);
    mbm.set_trans(disp, 0.0, 0.0);
    mbm.set_rotaxis(0.0, 0.0, -1.0);
    mbm.set_rotpt(4.75e-3 + offset, -1.125e-3, 0.0);
    mm.rigidBodyMotions.push_back(mbm);
// return
    return mm;
  }

}
