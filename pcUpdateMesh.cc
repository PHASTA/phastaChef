#include "pcUpdateMesh.h"
#include <SimPartitionedMesh.h>
#include <cassert>

namespace pc {
  bool updateAPFCoord(apf::Mesh2* m) {
    apf::Field* f = m->findField("motion_coords");
    assert(f);
    double* vals = new double[apf::countComponents(f)];
    assert(apf::countComponents(f) == 3);
    apf::MeshEntity* vtx;
    apf::Vector3 points;
    apf::MeshIterator* itr = m->begin(0);
    while( (vtx = m->iterate(itr)) ) {
      apf::getComponents(f, vtx, 0, vals);
      for ( int i = 0; i < 3; i++ )  points[i] = vals[i];
      m->setPoint(vtx, 0, points);
    }
    m->end(itr);
    delete [] vals;
    return true;
  }

  bool updateSIMCoord(apf::Mesh2* m) {
    apf::MeshSIM* apf_msim = dynamic_cast<apf::MeshSIM*>(m);
    pParMesh ppm = apf_msim->getMesh();
    pMesh pm = PM_mesh(ppm,0);

    M_write(pm, "before_overwrite.sms", 0, NULL);

    apf::Field* f = m->findField("motion_coords");
	assert(f);
	double* vals = new double[apf::countComponents(f)];
	assert(apf::countComponents(f) == 3);

    fprintf(stderr, "not finished yet\n");
    assert(0);

    apf::MeshEntity* vtx;
	apf::MeshIterator* itr = m->begin(0);
	while( (vtx = m->iterate(itr)) ) {
	  apf::getComponents(f, vtx, 0, vals);
	  pVertex simVtx = reinterpret_cast<pVertex>(vtx);
//	  if(simVtx is classified on moving body ) continue;
      switch(V_whatInType(simVtx)){
	    case Gedge :
		  break;
		case Gface :
		  break;
		case Gregion :
		  break;
	  }
	}
	m->end(itr);
  }

  void updateMeshCoord(ph::Input& in, apf::Mesh2* m) {
    bool done = false;
    if (in.simmetrixMesh)
      done = updateSIMCoord(m);
    else
      done = updateAPFCoord(m);
    assert(done);
  }

}
