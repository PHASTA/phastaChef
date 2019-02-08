#include "pcError.h"
#include "pcUpdateMesh.h"
#include "pcSmooth.h"
#include "pcWriteFiles.h"
#include <MeshSimAdapt.h>
#include <SimUtil.h>
#include <SimPartitionedMesh.h>
#include <SimDiscrete.h>
#include "SimMeshMove.h"
#include "SimMeshTools.h"
#include "apfSIM.h"
#include "gmi_sim.h"
#include <PCU.h>
#include <cassert>
#include <phastaChef.h>

namespace pc {
  getCellMeshSize() {
    //read phasta cell-based field errorH1
    //loop over elements
    //get shortest height of this element
    //get error
    //get desired error
    //get new size
    //set new size
  }

  getNodalMeshSize() {
    //read cell-based size
    //loop over vertices
    //loop over adjacent elements
    //get weighted size and weight
    //get size of this vertex
    //set new size
  }

  getVMSMeshSize() {
    getCellMeshSize();
    getNodalMeshSize();
  }
}
