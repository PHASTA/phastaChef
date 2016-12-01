#include "customSF.h"

namespace {

  apf::Field* getSprSF(apf::Mesh2* m) {
    const int order = 2;
    double adaptRatio = 0.1;
//    apf::Field* szFld;
    apf::Field* temperature = chef::extractField(m,"solution","temperature",5,apf::SCALAR);
    assert(temperature);
    apf::Field* eps = spr::getGradIPField(temperature, "eps", order);
//    apf::Field* test = apf::createFieldOn(m, "TEST", apf::SCALAR);
//    apf::zeroField(test);
//    apf::Field* eps_star = spr::recoverField(eps);
    apf::writeVtkFiles("test_eps",m);
    apf::destroyField(temperature);
    apf::Field* szFld = spr::getSPRSizeField(eps,adaptRatio);
    apf::destroyField(eps);
    return szFld;
  }

  apf::Field* multipleSF(apf::Mesh* m, apf::Field* sf, double factor) {
    apf::Field* sz = createFieldOn(m, "multipliedSize", apf::SCALAR);
    apf::MeshEntity* vtx;
    apf::MeshIterator* itr = m->begin(0);
    while( (vtx = m->iterate(itr)) ) {
      double h = apf::getScalar(sf,vtx,0);
      apf::setScalar(sz,vtx,0,h*factor);
    }
    m->end(itr);
    return sz; 
  }

  /* cylinder is defined as [center x, center y, center z, 
                             half of x length, radius] */
  apf::Field* cylRefineSF(apf::Mesh2* m, double* box, double ref, double cor) {
    apf::Field* newSz = apf::createFieldOn(m,"cylRefineIso",apf::SCALAR);
    apf::Vector3 points;
    double h = 0.0;
    assert(box[3] > 0 && box[4] > 0);
    apf::MeshEntity* vtx;
    apf::MeshIterator* itr = m->begin(0);
    while( (vtx = m->iterate(itr)) ) {
      m->getPoint(vtx, 0, points);
      /* assume the cylinder normal is x direction */
      double xDis = fabs(points[0]-box[0]);
      double rDis = sqrt((points[1]-box[1])*(points[1]-box[1]) +  
                         (points[2]-box[2])*(points[2]-box[2]));
      if ( xDis < box[3] && rDis < box[4] ) {
        h = ref;
      }
      else if (xDis < 3*box[3] && rDis < 3*box[4]) {
        double xFtr = 0.0;
        double rFtr = 0.0;
        /* set x factor */
        if (xDis < box[3]) 
          xFtr = 1.0;
        else 
          xFtr = (3.0*box[3]-xDis)/2.0/box[3];
        /* set radius factor */
        if (rDis < box[4]) 
          rFtr = 1.0;
        else 
          rFtr = (3.0*box[4]-rDis)/2.0/box[4];
        h = ref*xFtr*rFtr + cor*(1.0-xFtr*rFtr);
      }
      else {
        h = cor;
      }
      apf::setScalar(newSz,vtx,0,h);
    }
    m->end(itr);
    return newSz;
  }

} //end namespace


 
