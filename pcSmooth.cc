#include "pcSmooth.h"
#include <MeshSimAdapt.h>
#include <SimUtil.h>
#include <SimPartitionedMesh.h>
#include <SimDiscrete.h>
#include "apfSIM.h"
#include "gmi_sim.h"
#include <PCU.h>
#include <cassert>
#include <queue>
#include <iostream>
#include <algorithm>
#include <ph.h>
#include <phastaChef.h>

/* This part of code was originally written by Alvin Zhang */

using namespace std;

namespace pc {

bool isInCylinder(apf::MeshEntity* en) {
  double x_min = 0.0;
  double x_max = 2.0;
  double r_max = 0.065;

  apf::Adjacent enAdjVert;
  m->getAdjacent(en, 0, enAdjVert);
  apf::Vector3 xyz = apf::Vector3(0.0,0.0,0.0);

  for (size_t i=0; i < enAdjVert.getSize(); ++i) {
    m->getPoint(enAdjVert[i],0,xyz);
    double xyz_r = sqrt(xyz[1]*xyz[1] + xyz[2]*xyz[2]);
    if (xyz[0] > x_min && xyz[0] < x_max && xyz_r < r_max)
      return true;
  }
  return false;
}

int gradeSizeModify(apf::Mesh* m, double gradingFactor,
    double size[2], apf::Adjacent edgAdjVert,
    apf::Adjacent vertAdjEdg,
    std::queue<apf::MeshEntity*> &markedEdges,
    apf::MeshTag* isMarked,
    int idxFlag)
{
  //Determine a switching scheme depending on which vertex needs a modification
  int idx1,idx2;
  if(idxFlag == 0){
    idx1=0;
    idx2=1;
  }
  else{
    idx1=1;
    idx2 = 0;
  }

  int marker[3] = {0,1,0};
  double marginVal = 0.01;
  int needsParallel=0;
  apf::Vector3 v_mag = apf::Vector3(0.0,0.0,0.0);

  apf::Field* sizes = m->findField("sizes");

  if(size[idx1]>(gradingFactor*size[idx2])*(1+marginVal))
  {
    if(m->isOwned(edgAdjVert[idx1]))
    {
      size[idx1] = gradingFactor*size[idx2];
      v_mag = apf::Vector3(size[idx1], size[idx1], size[idx1]);
      apf::setVector(sizes,edgAdjVert[idx1],0,v_mag);
      m->getAdjacent(edgAdjVert[idx1], 1, vertAdjEdg);
      for (std::size_t i=0; i<vertAdjEdg.getSize();++i){
        m->getIntTag(vertAdjEdg[i],isMarked,&marker[2]);
        //if edge is not already marked
        if(!marker[2]){
          if (isInCylinder(vertAdjEdg[i])) { // for projectile case only
            m->setIntTag(vertAdjEdg[i],isMarked,&marker[1]);
            markedEdges.push(vertAdjEdg[i]);
          }
        }
      }
    } //end isOwned
    else
    { //Pack information to owning processor
      needsParallel=1;
      apf::Copies remotes;
      m->getRemotes(edgAdjVert[idx1],remotes);
      double newSize = gradingFactor*size[idx2];
      int owningPart=m->getOwner(edgAdjVert[idx1]);
      PCU_COMM_PACK(owningPart, remotes[owningPart]);
      PCU_COMM_PACK(owningPart,newSize);
    }
  }

  return needsParallel;
}

void markEdgesInitial(apf::Mesh* m, std::queue<apf::MeshEntity*> &markedEdges,double gradingFactor)
{
  //marker structure for 0) not marked 1) marked 2)storage
  int marker[3] = {0,1,0};

  double size[2];
  apf::Vector3 v_mag = apf::Vector3(0.0,0.0,0.0);
  apf::MeshTag* isMarked = m->findTag("isMarked");
  apf::Field* sizes = m->findField("sizes");
  apf::Adjacent edgAdjVert;
  apf::MeshEntity* edge;
  apf::MeshIterator* it = m->begin(1);
  while((edge=m->iterate(it))){
    m->getAdjacent(edge, 0, edgAdjVert);

    // for projectile case only
    if (isInCylinder(edge)) {

    for (std::size_t i=0; i < edgAdjVert.getSize(); ++i){
      apf::getVector(sizes,edgAdjVert[i],0,v_mag);
      size[i] = v_mag[0];
    }
    if( (size[0] > gradingFactor*size[1]) || (size[1] > gradingFactor*size[0]) ){
      //add edge to a queue
      markedEdges.push(edge);
      //tag edge to indicate that it is part of queue
      m->setIntTag(edge,isMarked,&marker[1]);
    }

    } // end isInCylinder
    else{
      m->setIntTag(edge,isMarked,&marker[0]);
    }
  }
  m->end(it);
}

int serialGradation(apf::Mesh* m, std::queue<apf::MeshEntity*> &markedEdges,double gradingFactor)
{
  double size[2];
  apf::Vector3 v_mag = apf::Vector3(0.0,0.0,0.0);
  //marker structure for 0) not marked 1) marked 2)storage
  int marker[3] = {0,1,0};
  apf::MeshTag* isMarked = m->findTag("isMarked");
  apf::Field* sizes = m->findField("sizes");
  apf::Adjacent edgAdjVert;
  apf::Adjacent vertAdjEdg;
  apf::MeshEntity* edge;
  apf::MeshIterator* it = m->begin(1);
  int needsParallel=0;

  //perform serial gradation while packing necessary info for parallel
  while(!markedEdges.empty()){
    edge = markedEdges.front();
    m->getAdjacent(edge, 0, edgAdjVert);
    for (std::size_t i=0; i < edgAdjVert.getSize(); ++i){
      apf::getVector(sizes,edgAdjVert[i],0,v_mag);
      size[i] = v_mag[0];
    }

    needsParallel+=gradeSizeModify(m, gradingFactor, size, edgAdjVert,
      vertAdjEdg, markedEdges, isMarked, 0);
    needsParallel+=gradeSizeModify(m, gradingFactor, size, edgAdjVert,
      vertAdjEdg, markedEdges, isMarked, 1);

    m->setIntTag(edge,isMarked,&marker[0]);
    markedEdges.pop();
  }
  return needsParallel;
}

void addSmoother(apf::Mesh2* m, double gradingFactor) {
  meshGradation(m, gradingFactor);
}

void meshGradation(apf::Mesh2* m, double gradingFactor)
{
  if(!PCU_Comm_Self())
    std::cout<<"Starting grading\n";
  apf::MeshEntity* edge;
  apf::Adjacent edgAdjVert;
  apf::Adjacent vertAdjEdg;
  double size[2];
  std::queue<apf::MeshEntity*> markedEdges;
  apf::MeshTag* isMarked = m->createIntTag("isMarked",1);
  apf::Field* sizes = m->findField("sizes");

  //marker structure for 0) not marked 1) marked 2)storage
  int marker[3] = {0,1,0};

  apf::MeshIterator* it;
  markEdgesInitial(m,markedEdges,gradingFactor);

  int needsParallel=1;
  int nCount=1;
  while(needsParallel)
  {
    PCU_Comm_Begin();
    needsParallel = serialGradation(m,markedEdges,gradingFactor);

    PCU_Add_Ints(&needsParallel,1);
    if(!PCU_Comm_Self())
      std::cerr<<"Sending size info for gradation"<<std::endl;
    PCU_Comm_Send();

    apf::MeshEntity* ent;
    double receivedSize;
    double currentSize;
    double newSize;
    apf::Vector3 v_mag = apf::Vector3(0.0,0.0,0.0);

    //Need a container to get all entitites that need to be updated on remotes
    std::queue<apf::MeshEntity*> updateRemoteVertices;

    apf::Copies remotes;
    //owning copies are receiving
    while(PCU_Comm_Receive())
    {
      PCU_COMM_UNPACK(ent);
      PCU_COMM_UNPACK(receivedSize);

      if(!m->isOwned(ent)){
        std::cout<<"THERE WAS AN ERROR"<<std::endl;
        std::exit(1);
      }

      apf::getVector(sizes,ent,0,v_mag);
      currentSize = v_mag[0];
      newSize = std::min(receivedSize,currentSize);
      v_mag = apf::Vector3(newSize, newSize, newSize);
      apf::setVector(sizes,ent,0,v_mag);

      //add adjacent edges into Q
      m->getAdjacent(ent, 1, vertAdjEdg);
      for (std::size_t i=0; i<vertAdjEdg.getSize();++i)
      {
        edge = vertAdjEdg[i];
        m->getIntTag(vertAdjEdg[i],isMarked,&marker[2]);
        if(!marker[2])
        {
          if (isInCylinder(edge)) { // for projectile case only
            markedEdges.push(edge);
            //tag edge to indicate that it is part of queue
            m->setIntTag(edge,isMarked,&marker[1]);
          }
        }
      }
      updateRemoteVertices.push(ent);
    }

    PCU_Comm_Begin();

    while(!updateRemoteVertices.empty())
    {
      ent = updateRemoteVertices.front();
      //get remote copies and send updated mesh sizes
      m->getRemotes(ent,remotes);
      for(apf::Copies::iterator iter=remotes.begin(); iter!=remotes.end();++iter)
      {
        PCU_COMM_PACK(iter->first, iter->second);
      }
      updateRemoteVertices.pop();
    }

    PCU_Comm_Send();
    //while remote copies are receiving
    while(PCU_Comm_Receive())
    {
      //unpack
      PCU_COMM_UNPACK(ent);
      //PCU_COMM_UNPACK(receivedSize);
      assert(!m->isOwned(ent));

      if(m->isOwned(ent)){
        std::cout<<"Problem occurred\n";
        std::exit(1);
      }

      //add adjacent edges into Q
      m->getAdjacent(ent, 1, vertAdjEdg);
      for (std::size_t i=0; i<vertAdjEdg.getSize();++i)
      {
        edge = vertAdjEdg[i];
        m->getIntTag(vertAdjEdg[i],isMarked,&marker[2]);
        if(!marker[2])
        {
          if (isInCylinder(edge)) { // for projectile case only
            markedEdges.push(edge);
            //tag edge to indicate that it is part of queue
            m->setIntTag(edge,isMarked,&marker[1]);
          }
        }
      }
    }
    apf::synchronize(sizes);

  } //end outer while

  //Cleanup of edge marker field
  it = m->begin(1);
  while((edge=m->iterate(it))){
    m->removeTag(edge,isMarked);
  }
  m->end(it);
  m->destroyTag(isMarked);

  if(!PCU_Comm_Self())
    std::cout<<"Completed grading\n";
}

} // end namespace pc
