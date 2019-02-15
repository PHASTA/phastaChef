#include "pcError.h"
#include <MeshSimAdapt.h>
#include <SimUtil.h>
#include <SimPartitionedMesh.h>
#include "SimMeshMove.h"
#include "SimMeshTools.h"
#include "apfSIM.h"
#include "gmi_sim.h"
#include "apfShape.h"
#include <PCU.h>
#include <cassert>
#include <phastaChef.h>

namespace pc {
  void attachVMSSizeFieldH1(apf::Mesh2*& m, ph::Input& in) {
    //read phasta cell-based field errorH1
    apf::Field* err = m->findField("errorH1");
    //get nodal-based mesh size field
    apf::Field* sizes = m->findField("sizes");
    //create a field to store cell-based mesh size
    int nsd = m->getDimension();
    apf::Field* cell_size = apf::createField(m, "cell_size", apf::SCALAR, apf::getConstant(nsd));

    //get desired error
    double desr_err[3];
    desr_err[0] = in.simAdaptDesiredErrorMass;
    desr_err[1] = in.simAdaptDesiredErrorMomt;
    desr_err[2] = in.simAdaptDesiredErrorEnrg;

    //loop over elements
    apf::NewArray<double> curr_err(apf::countComponents(err));
    apf::MeshElement *me;
    apf::MeshEntity* elm;
    apf::MeshIterator* it = m->begin(nsd);
    while ((elm = m->iterate(it))) {
      double h_old = 0.0;
      double h_new = 0.0;
      me = apf::createMeshElement(m, elm);
      //get shortest height of this element
      if (nsd == 2)
        h_old = sqrt(apf::measure(me) * 4 / sqrt(3));
      else
        h_old = apf::computeShortestHeightInTet(m,elm);
      //get error
      apf::getComponents(err, elm, 0, &curr_err[0]);
      //get new size
      //currently, we only focus on the momemtum error // debugging
      double factor = 0.0;
      if (desr_err[1] / curr_err[1] > 100.0)
        factor = 100.0;
      else
        factor = desr_err[1] / curr_err[1];
      h_new = h_old * pow(factor, 2.0/(2.0*(1.0)+nsd));
      //set new size
      apf::setScalar(cell_size, elm, 0, h_new);
      apf::destroyMeshElement(me);
    }
    m->end(it);

    //loop over vertices
    apf::MeshEntity* vtx;
    it = m->begin(0);
    while ((vtx = m->iterate(it))) {
      apf::Adjacent adj_elm;
      m->getAdjacent(vtx, m->getDimension(), adj_elm);
      double weightedSize = 0.0;
      double totalError = 0.0;
      //loop over adjacent elements
      for (std::size_t i = 0; i < adj_elm.getSize(); ++i) {
        //get weighted size and weight
        apf::getComponents(err, adj_elm[i], 0, &curr_err[0]);
        double curr_size = apf::getScalar(cell_size, adj_elm[i], 0);
        //currently, we only focus on the momemtum error // debugging
//        weightedSize += apf::getScalar(cell_size,adj_elm[i],0)*curr_err[1];
        weightedSize += curr_size*curr_err[1];
        totalError += curr_err[1];
      }
      //get size of this vertex
      weightedSize = weightedSize / totalError;
      //set new size
      apf::Vector3 v_mag;
      v_mag[0] = weightedSize;
      v_mag[1] = weightedSize;
      v_mag[2] = weightedSize;
      apf::setVector(sizes, vtx, 0, v_mag);
    }
    m->end(it);

    //delete cell-based error and mesh size
    apf::destroyField(err);
    apf::destroyField(cell_size);
  }

  void attachVMSSizeField(apf::Mesh2*& m, ph::Input& in) {
    attachVMSSizeFieldH1(m, in);
  }
}
