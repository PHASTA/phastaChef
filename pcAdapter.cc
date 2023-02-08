#include "pcAdapter.h"
#include "pcError.h"
#include "pcUpdateMesh.h"
#include "pcSmooth.h"
#include "pcWriteFiles.h"
#include <SimUtil.h>
#include <SimPartitionedMesh.h>
#include <SimDiscrete.h>
#include "SimMeshMove.h"
#include "SimMeshTools.h"
#include <SimAdvMeshing.h>
#include "apfSIM.h"
#include "gmi_sim.h"
#include <PCU.h>
#include <cassert>
#include <phastaChef.h>
#include <maStats.h>
#include <apfShape.h>
#include <math.h>
#include <fstream>
#include <limits>
#include <iomanip>
#include <set>
#include <map>
#include <spr.h>

extern void MSA_setBLSnapping(pMSAdapt, int onoff);

namespace pc {

     std::vector<std::string> ParseDeliminatedString(std::string& FullLine, std::string& delim){
          std::vector<std::string> output;
          size_t pos=0;
          std::string token;

          while ((pos = FullLine.find(delim)) != std::string::npos) {
               token = FullLine.substr(0, pos);
               output.push_back(token);
               FullLine.erase(0, pos + delim.length());
          }
          output.push_back(token);


          return output;          
     }

     std::set<apf::MeshEntity*> TraceSurf(std::map<apf::MeshEntity*, std::pair<std::set<apf::MeshEntity*>,int> >& LocGraph,
                                        std::set<apf::MeshEntity*>& StartSet){
          
          std::set<apf::MeshEntity*> output=StartSet;
          for (auto itr1=StartSet.begin(); itr1!=StartSet.end(); itr1++){
               for (auto nei=LocGraph[*itr1].first.begin(); nei!=LocGraph[*itr1].first.end(); nei++){
                    output.insert(*nei); 
               }     
          }
     return output;
     }

     void PropagateID1(std::map<apf::MeshEntity*, std::pair<set<apf::MeshEntity*>, int> >& LocGraph,
                         std::set<apf::MeshEntity*>& InputSet,
                         int ID){
          if(InputSet.size()==0) {return;}
          std::set<apf::MeshEntity*> NextSet;
          for(auto itr=InputSet.begin(); itr!=InputSet.end(); itr++){
               //set ID
               LocGraph[*itr].second=ID;
               //check neighbours
               for(auto nei=LocGraph[*itr].first.begin(); nei!=LocGraph[*itr].first.end(); nei++ ){
                    if (LocGraph[*nei].second!=0){//skip set IDs
                         NextSet.insert(*nei);     
                    }
               }
          }
          PropagateID1(LocGraph,NextSet,ID);
          return;
     }


     void PropagateID(std::set<apf::MeshEntity*>& Start, std::map<apf::MeshEntity*, std::pair<set<apf::MeshEntity*>, int> >& LocGraph, int ID){
         int check=0; std::set<apf::MeshEntity*> Fill=Start;
         while (check!=Fill.size()){
              check=Fill.size();
              Fill=TraceSurf(LocGraph, Fill);
          }
          for (auto itr=Fill.begin(); itr!=Fill.end(); itr++){
               LocGraph[*itr].second=ID;     
          }
     }


  std::set<apf::MeshEntity*> GetNeighs(std::set<apf::MeshEntity*>& Input, apf::Mesh2* m){
     std::set<apf::MeshEntity*> Output = Input;
     for (auto itr1=Input.begin(); itr1!=Input.end(); itr1++){
          apf::Adjacent Adja; apf::MeshEntity* tmp;
          getBridgeAdjacent(m, *itr1, 2,3, Adja); //get adjacent elements;
          for(size_t i=0; i<Adja.getSize(); i++){
               tmp=Adja[i];
               Output.insert(tmp);
          }
     }
     return Output;
  } 

  std::set<apf::MeshEntity*> GetN_Neighs(apf::MeshEntity* elm, apf::Mesh2* m, int n){
     int check=0; std::set<apf::MeshEntity*> Fill={elm};
     
     while (check !=n){
          Fill=GetNeighs(Fill, m);
          check++;
     }
     return Fill;
  }


  bool vertexIsInCylinder(apf::MeshEntity* v) {
    double x_min = 0.0;
    double x_max = 2.01;
    double r_max = 0.07;//outer radius

    apf::Vector3 xyz = apf::Vector3(0.0,0.0,0.0);
    m->getPoint(v,0,xyz);
    double xyz_r = sqrt(xyz[1]*xyz[1] + xyz[2]*xyz[2]);
    if (xyz[0] > x_min && xyz[0] < x_max && xyz_r < r_max)
      return true;
    return false;
  }


  apf::Field* convertField(apf::Mesh* m,
    const char* inFieldname,
    const char* outFieldname) {
    apf::Field* inf = m->findField(inFieldname);
    assert(inf);
    int size = apf::countComponents(inf);
    apf::Field* outf = m->findField(outFieldname);
    if (outf)
      apf::destroyField(outf);
    outf = apf::createPackedField(m, outFieldname, size);
    apf::NewArray<double> inVal(size);
    apf::NewArray<double> outVal(size);
    apf::MeshEntity* vtx;
    apf::MeshIterator* it = m->begin(0);
    while ((vtx = m->iterate(it))) {
      apf::getComponents(inf, vtx, 0, &inVal[0]);
      for (int i = 0; i < size; i++){
        outVal[i] = inVal[i];
      }
      apf::setComponents(outf,vtx, 0, &outVal[0]);
    }
    m->end(it);
    apf::destroyField(inf);
    return outf;
  }

  apf::Field* convertVtxFieldToElm(apf::Mesh* m,
                  const char* vFieldname,
                  const char* eFieldname) {
    apf::Field* vf = m->findField(vFieldname);
    assert(vf);
    int size = apf::countComponents(vf);
    apf::Field* ef = m->findField(eFieldname);
    if (ef) apf::destroyField(ef);
    ef = apf::createPackedField(m, eFieldname, size, apf::getConstant(m->getDimension()));
    apf::NewArray<double> vVal(size);
    apf::NewArray<double> eVal(size);
    apf::MeshEntity* elm;
    apf::MeshIterator* it = m->begin(m->getDimension());
    while ((elm = m->iterate(it))) {
      for (int i = 0; i < size; i++) eVal[i] = 0.0;
      apf::Downward vtx;
      int nbv = m->getDownward(elm, 0, vtx);
      for (int j = 0; j < nbv; j++){
        apf::getComponents(vf, vtx[j], 0, &vVal[0]);
        for (int i = 0; i < size; i++){
          eVal[i] += vVal[i]/(double)nbv;
        }
      }
      apf::setComponents(ef, elm, 0, &eVal[0]);
    }
    m->end(it);
    apf::destroyField(vf);
    return ef;
  }

  void attachMeshSizeField(apf::Mesh2*& m, ph::Input& in, phSolver::Input& inp) {
    /* create a field to store mesh size */
    if(m->findField("sizes")) apf::destroyField(m->findField("sizes"));
    apf::Field* sizes = apf::createSIMFieldOn(m, "sizes", apf::VECTOR);
    /* switch between VMS error mesh size and initial mesh size */
    if((string)inp.GetValue("Error Estimation Option") != "False") {
      pc::attachVMSSizeField(m, in, inp);
    }
    else {
      if(m->findField("frames")) apf::destroyField(m->findField("frames"));
      apf::Field* frames = apf::createSIMFieldOn(m, "frames", apf::MATRIX);
      ph::attachSIMSizeField(m, sizes, frames);
    }
// prescribe mesh size field for the projectile case
// this is hardcoded, please comment out this call for other usage
//    pc::prescribe_proj_mesh_size(m, sizes, in.rbParamData[0]);
  }

  int getNumOfMappedFields(apf::Mesh2*& m) {
    /* initially, we have 7 fields: pressure, velocity, temperature,
       time der of pressure, time der of velocity, time der of temperature,
       ,mesh velocity and 1 optional field: time resource bound factor field */
    int numOfMappedFields;
    if (m->findField("ctcn_elm")) numOfMappedFields = 8;
    else numOfMappedFields = 7;
    return numOfMappedFields;
  }

  /* remove all fields except for solution, time
     derivative of solution, mesh velocity and keep
     certain field if corresponding option is on */
  void removeOtherFields(apf::Mesh2*& m, phSolver::Input& inp) {
    int index = 0;
    int numOfPackFields = 4;
    while (m->countFields() > numOfPackFields) {
      apf::Field* f = m->getField(index);
      if ( f == m->findField("solution") ||
           f == m->findField("time derivative of solution") ||
           f == m->findField("mesh_vel") ||
           f == m->findField("ctcn_elm") ) {
        index++;
        continue;
      }
      m->removeField(f);
      apf::destroyField(f);
    }
    m->verify();
  }

  int getSimFields(apf::Mesh2*& m, int simFlag, pField* sim_flds, phSolver::Input& inp) {
    int num_flds = 0;
    if (m->findField("solution")) {
      num_flds += 3;
      sim_flds[0] = apf::getSIMField(chef::extractField(m,"solution","pressure",1,apf::SCALAR,simFlag));
      sim_flds[1] = apf::getSIMField(chef::extractField(m,"solution","velocity",2,apf::VECTOR,simFlag));
      sim_flds[2] = apf::getSIMField(chef::extractField(m,"solution","temperature",5,apf::SCALAR,simFlag));
      apf::destroyField(m->findField("solution"));
    }

    if (m->findField("time derivative of solution")) {
      num_flds += 3;
      sim_flds[3] = apf::getSIMField(chef::extractField(m,"time derivative of solution","der_pressure",1,apf::SCALAR,simFlag));
      sim_flds[4] = apf::getSIMField(chef::extractField(m,"time derivative of solution","der_velocity",2,apf::VECTOR,simFlag));
      sim_flds[5] = apf::getSIMField(chef::extractField(m,"time derivative of solution","der_temperature",5,apf::SCALAR,simFlag));
      apf::destroyField(m->findField("time derivative of solution"));
    }

    if (m->findField("mesh_vel")) {
      num_flds += 1;
      sim_flds[6] = apf::getSIMField(chef::extractField(m,"mesh_vel","mesh_vel_sim",1,apf::VECTOR,simFlag));
      apf::destroyField(m->findField("mesh_vel"));
    }

    if (m->findField("ctcn_elm")) {
      num_flds += 1;
      sim_flds[7] = apf::getSIMField(chef::extractField(m,"ctcn_elm","ctcn_elm_sim",1,apf::SCALAR,simFlag));
      apf::destroyField(m->findField("ctcn_elm"));
    }

    return num_flds;
  }

  /* unpacked solution into serveral fields,
     put these field explicitly into pPList */
  pPList getSimFieldList(ph::Input& in, apf::Mesh2*& m){
    /* load input file for solver */
    phSolver::Input inp("solver.inp", "input.config");
    int num_flds = getNumOfMappedFields(m);
    removeOtherFields(m,inp);
    pField* sim_flds = new pField[num_flds];
    getSimFields(m, in.simmetrixMesh, sim_flds, inp);
    pPList sim_fld_lst = PList_new();
    for (int i = 0; i < num_flds; i++) {
      PList_append(sim_fld_lst, sim_flds[i]);
    }
    assert(num_flds == PList_size(sim_fld_lst));
    delete [] sim_flds;
    return sim_fld_lst;
  }

  void transferSimFields(apf::Mesh2*& m) {
    if (m->findField("pressure")) // assume we had solution before
      chef::combineField(m,"solution","pressure","velocity","temperature");
    if (m->findField("der_pressure")) // assume we had time derivative of solution before
      chef::combineField(m,"time derivative of solution","der_pressure","der_velocity","der_temperature");
    if (m->findField("mesh_vel_sim"))
      convertField(m, "mesh_vel_sim", "mesh_vel");
    if (m->findField("ctcn_elm_sim"))
      convertVtxFieldToElm(m, "ctcn_elm_sim", "err_tri_f");
    // destroy mesh size field
    if(m->findField("sizes"))  apf::destroyField(m->findField("sizes"));
    if(m->findField("frames")) apf::destroyField(m->findField("frames"));
  }

  void attachCurrentSizeField(apf::Mesh2*& m) {
    int  nsd = m->getDimension();
    if(m->findField("cur_size")) apf::destroyField(m->findField("cur_size"));
    apf::Field* cur_size = apf::createField(m, "cur_size", apf::SCALAR, apf::getConstant(nsd));
    // loop over non-BL elements
    apf::MeshEntity* e;
    apf::MeshIterator* eit = m->begin(nsd);
    while ((e = m->iterate(eit))) {
      pRegion meshRegion = reinterpret_cast<pRegion>(e);
      if (EN_isBLEntity(meshRegion)) continue;
      // set mesh size field
      double h = 0.0;
      if (m->getType(e) == apf::Mesh::TET)
        h = apf::computeShortestHeightInTet(m,e) * sqrt(3.0);
      else
        h = pc::getShortestEdgeLength(m,e);
      apf::setScalar(cur_size, e, 0, h);
    }
    m->end(eit);

    // get sim model
    apf::MeshSIM* sim_m = dynamic_cast<apf::MeshSIM*>(m);
    pParMesh sim_pm = sim_m->getMesh();
    pMesh pm = PM_mesh(sim_pm,0);

    gmi_model* gmiModel = sim_m->getModel();
    pGModel model = gmi_export_sim(gmiModel);

    // loop over model faces
    pGFace modelFace;
    pFace meshFace;
    pFace blFace;
    pRegion blRegion;
    FIter fIter;
    pEntity seed;
    pPList growthRegions = PList_new();
    pPList growthFaces = PList_new();
    GFIter gfIter = GM_faceIter(model);
    while((modelFace=GFIter_next(gfIter))){
      // loop over mesh faces on model face
      fIter = M_classifiedFaceIter(pm, modelFace, 1);
      while((meshFace = FIter_next(fIter))){
        apf::MeshEntity* apf_f = reinterpret_cast<apf::MeshEntity*>(meshFace);
        // check if BL base
        if (BL_isBaseEntity(meshFace, modelFace)) {
          // loop over BL regions and layers
          for(int faceSide = 0; faceSide < 2; faceSide++){
            int hasSeed = BL_stackSeedEntity(meshFace, modelFace, faceSide, NULL, &seed);
            if (hasSeed) {
              BL_growthRegionsAndLayerFaces((pRegion)seed, growthRegions, growthFaces, Layer_Entity);
              if (PList_size(growthRegions) >= (PList_size(growthFaces)-1)*3) { // tet
                for(int i = 0; i < PList_size(growthFaces); i++) {
                  blFace = (pFace)PList_item(growthFaces,i);
                  apf::MeshEntity* apf_f = reinterpret_cast<apf::MeshEntity*>(blFace);
                  double h = apf::computeShortestHeightInTri(m,apf_f) * sqrt(2.0);
                  for(int j = 0; j < 3; j++) {
                    if (i*3+j == PList_size(growthRegions)) break;
                    blRegion = (pRegion)PList_item(growthRegions,i*3+j);
                    // set mesh size field
                    apf::MeshEntity* apf_r = reinterpret_cast<apf::MeshEntity*>(blRegion);
                    apf::setScalar(cur_size, apf_r, 0, h);
                  }
                }
              }
              else if (PList_size(growthRegions) >= (PList_size(growthFaces)-1)) { // wedge
                for(int i = 0; i < PList_size(growthFaces); i++) {
                  if (i == PList_size(growthRegions)) break;
                  blFace = (pFace)PList_item(growthFaces,i);
                  apf::MeshEntity* apf_f = reinterpret_cast<apf::MeshEntity*>(blFace);
                  double h = apf::computeShortestHeightInTri(m,apf_f) * sqrt(2.0);
                  blRegion = (pRegion)PList_item(growthRegions,i);
                  // set mesh size field
                  apf::MeshEntity* apf_r = reinterpret_cast<apf::MeshEntity*>(blRegion);
                  apf::setScalar(cur_size, apf_r, 0, h);
                }
              }
            }
            else if (hasSeed < 0) {
              printf("not support blending BL mesh or miss some info!\n");
              exit(0);
            }
          }
        }
      }
      FIter_delete(fIter);
    }
    GFIter_delete(gfIter);
  }

  double estimateAdaptedMeshElements(apf::Mesh2*& m, apf::Field* sizes) {
    attachCurrentSizeField(m);
    apf::Field* cur_size = m->findField("cur_size");
    assert(cur_size);

    double estElm = 0.0;
    apf::Vector3 v_mag  = apf::Vector3(0.0, 0.0, 0.0);
    int num_dims = m->getDimension();
    assert(num_dims == 3); // only work for 3D mesh
    apf::Vector3 xi = apf::Vector3(0.25, 0.25, 0);
    apf::MeshEntity* en;
    apf::MeshIterator* eit = m->begin(num_dims);
    while ((en = m->iterate(eit))) {
      apf::MeshElement* elm = apf::createMeshElement(m,en);
      apf::Element* fd_elm = apf::createElement(sizes,elm);
      apf::getVector(fd_elm,xi,v_mag);
      double h_old = apf::getScalar(cur_size,en,0);
      if(EN_isBLEntity(reinterpret_cast<pEntity>(en))) {
        estElm = estElm + (h_old/v_mag[0])*(h_old/v_mag[0]);
      }
      else {
        estElm = estElm + (h_old/v_mag[0])*(h_old/v_mag[0])*(h_old/v_mag[0]);
      }
    }
    m->end(eit);

    apf::destroyField(cur_size);

    double estTolElm = PCU_Add_Double(estElm);
    return estTolElm;
  }

  void initializeCtCn(apf::Mesh2*& m) {
    apf::Field* ctcn = apf::createSIMFieldOn(m, "ctcn_elm", apf::SCALAR);
    apf::MeshEntity* v;
    apf::MeshIterator* vit = m->begin(0);
    while ((v = m->iterate(vit))) {
      apf::setScalar(ctcn,v,0,1.0);
    }
    m->end(vit);
  }

  void applyMaxSizeBound(apf::Mesh2*& m, apf::Field* sizes, ph::Input& in) {
    //int barrelTag = 69;  // 2mm case
    //int domainTag = 113; // 2mm case
    //int projctTag = 108; // 2mm case
    //int domainTag = 1; //hollow barrel case
    //int a_barrelTag = 285; //back of barrel
    //int b_barrelTag = 377; //front of barrel


    // int domainTag = 171; //diamond (unfitted)
    //  int domainTag = 195; //capsule (unfitted)
    // int domainTag = 152; //capsule (m2_fit)
    // int domainTag = 153; //cpasule (m_fit_long)
    int domainTag = 92; //bullet stationary cube
    // int domainTag = 254;
    //int a_barrelTag = 185;
    //int b_barrelTag = 1;
    //int domainTag = 128;//2d dw case
    
    //apf::ModelEntity* bme = m->findModelEntity(3, barrelTag);
    
    apf::ModelEntity* dme = m->findModelEntity(3, domainTag); //for barrel problems
    
    //apf::ModelEntity* pme = m->findModelEntity(3, projctTag);

    //apf::ModelEntity* a_bme = m->findModelEntity(3,a_barrelTag); // for hb problems
    //apf::ModelEntity* b_bme = m->findModelEntity(3,b_barrelTag);//

    apf::Vector3 v_mag = apf::Vector3(0.0,0.0,0.0);
    //double dmeSizeUpperBound = 0.1; //upper bound for farfield

    apf::MeshEntity* v;
    apf::MeshIterator* vit = m->begin(0);
    while ((v = m->iterate(vit))) {
      apf::getVector(sizes,v,0,v_mag);
      for (int i = 0; i < 3; i++) {
        apf::ModelEntity* me = m->toModel(v);
        //if ((m->isInClosureOf(me, b_bme) || m->isInClosureOf(me, a_bme) ) && v_mag[i] > 0.016) v_mag[i] = 0.016;
        if (m->isInClosureOf(me, dme) && v_mag[i] > in.simSizeUpperBound) v_mag[i] = in.simSizeUpperBound;
        //if (m->isInClosureOf(me, pme) && v_mag[i] > 0.004) v_mag[i] = 0.004;
        /*if(vertexIsInCylinder(v) && v_mag[i] > in.simSizeUpperBound)
          v_mag[i] = in.simSizeUpperBound; //set in adapt.inp*/
      }
      apf::setVector(sizes,v,0,v_mag);
    }
    m->end(vit);
  }

  double applyMaxNumberElement(apf::Mesh2*& m, apf::Field* sizes, ph::Input& in)  {
    /* scale mesh if number of elements exceeds threshold */
    double N_est = estimateAdaptedMeshElements(m, sizes);
    double cn = N_est / (double)in.simMaxAdaptMeshElements;
    cn = (cn>1.0)?cbrt(cn):1.0;
    if(!PCU_Comm_Self())
      printf("Estimated No. of Elm: %f and c_N = %f\n", N_est, cn);
    apf::Field* sol = m->findField("solution");
    apf::Field* ctcn = m->findField("ctcn_elm");
    assert(sol);
    assert(ctcn);
    apf::Vector3 v_mag = apf::Vector3(0.0,0.0,0.0);
    apf::MeshEntity* v;
    apf::MeshIterator* vit = m->begin(0);
    while ((v = m->iterate(vit))) {
      apf::getVector(sizes,v,0,v_mag);
      for (int i = 0; i < 3; i++)
        v_mag[i] = v_mag[i] * cn;
      apf::setVector(sizes,v,0,v_mag);

      double f = apf::getScalar(ctcn,v,0);
      f = f * cn;
      apf::setScalar(ctcn,v,0,f);
    }
    m->end(vit);
    return cn;
  }

  void applyMaxTimeResource(apf::Mesh2*& m, apf::Field* sizes,
                            ph::Input& in, phSolver::Input& inp) {
    apf::Field* sol = m->findField("solution");
    apf::Field* ctcn = m->findField("ctcn_elm");
    assert(sol);
    assert(ctcn);
    apf::Vector3 v_mag = apf::Vector3(0.0,0.0,0.0);
    apf::NewArray<double> s(in.ensa_dof);
    double maxCt = 1.0;
    double minCtH = 1.0e16;
    apf::MeshEntity* v;
    apf::MeshIterator* vit = m->begin(0);
    while ((v = m->iterate(vit))) {
      apf::getComponents(sol, v, 0, &s[0]);
      double u = sqrt(s[1]*s[1]+s[2]*s[2]+s[3]*s[3]);
      double c = sqrt(1.4*8.3145*s[4]/0.029); // ideal air assumed here
      double t = inp.GetValue("Time Step Size");
      double h_min = (u+c)*t/in.simCFLUpperBound;
      if (h_min < in.simSizeLowerBound) h_min = in.simSizeLowerBound;
      apf::getVector(sizes,v,0,v_mag);
      double f = apf::getScalar(ctcn,v,0);
      for (int i = 0; i < 3; i++) {
        if(v_mag[i] < h_min) {
          if(h_min/(v_mag[i]) > maxCt) maxCt = h_min/(v_mag[i]);
          apf::setScalar(ctcn,v,0,h_min/v_mag[i]*f);
          if(h_min < minCtH) minCtH = h_min;
          v_mag[i] = h_min;
        }
      }
      apf::setVector(sizes,v,0,v_mag);
    }
    m->end(vit);

    double maxCtAll  = PCU_Max_Double(maxCt);
    double minCtHAll = PCU_Min_Double(minCtH);
    if (!PCU_Comm_Self())
      printf("max time resource bound factor and min reached size: %f and %f\n",maxCtAll,minCtHAll);
  }


  void syncMeshSize(apf::Mesh2*& m, apf::Field* sizes) {
    PCU_Comm_Begin();
    apf::Copies remotes;
    apf::Vector3 v_mag = apf::Vector3(0.0,0.0,0.0);
    apf::MeshEntity* v;
    apf::MeshIterator* vit = m->begin(0);
    while ((v = m->iterate(vit))) {
      apf::getVector(sizes,v,0,v_mag);
      if(m->isShared(v)) {
        m->getRemotes(v, remotes);
        APF_ITERATE(apf::Copies, remotes, rit) {
          PCU_COMM_PACK(rit->first, rit->second);
          PCU_Comm_Pack(rit->first, &(v_mag[0]), sizeof(double));
        }
      }
    }
    m->end(vit);

    PCU_Comm_Send();
    while (PCU_Comm_Receive()) {
      apf::MeshEntity* rv;
      PCU_COMM_UNPACK(rv);
      double rv_mag;
      PCU_Comm_Unpack(&(rv_mag), sizeof(double));
      apf::getVector(sizes,v,0,v_mag);
      if(rv_mag < v_mag[0]) { // smaller wins
        v_mag = apf::Vector3(rv_mag,rv_mag,rv_mag);
        apf::setVector(sizes,rv,0,v_mag);
      }
    }
  }

  void setupSimImprover(pVolumeMeshImprover vmi, pPList sim_fld_lst) {
    VolumeMeshImprover_setModifyBL(vmi, 1);
    VolumeMeshImprover_setShapeMetric(vmi, ShapeMetricType_VolLenRatio, 0.3);
    VolumeMeshImprover_setSmoothType(vmi, 1); // 0:Laplacian-based; 1:Gradient-based

    /* set fields to be mapped */
    if (PList_size(sim_fld_lst))
      VolumeMeshImprover_setMapFields(vmi, sim_fld_lst);
  }

  void setupSimAdapter(pMSAdapt adapter, ph::Input& in, apf::Mesh2*& m, pPList& sim_fld_lst) {
    MSA_setAdaptBL(adapter, 1);
    MSA_setExposedBLBehavior(adapter,BL_DisallowExposed);
    MSA_setBLSnapping(adapter, 0); // currently needed for parametric model
    MSA_setBLMinLayerAspectRatio(adapter, 0.0); // needed in parallel
    MSA_setSizeGradation(adapter, 1, 0.0);
	

    /* attach mesh size field */
    phSolver::Input inp("solver.inp", "input.config");
    attachMeshSizeField(m, in, inp);
    apf::Field* sizes = m->findField("sizes");
    assert(sizes);
	
    /* initial ctcn field */
    pc::initializeCtCn(m);

    /* apply upper bound */
    // pc::applyMaxSizeBound(m, sizes, in);

    /* apply max number of element */
    // double cn = pc::applyMaxNumberElement(m, sizes, in);

    /* scale mesh if reach time resource bound */
    // pc::applyMaxTimeResource(m, sizes, in, inp);

    /* apply upper bound */
    pc::applyMaxSizeBound(m, sizes, in);

    apf::Field* position = m->findField("motion_coords");
    apf::Vector3 pos;
    apf::Field* PG_avg = m->findField("PG_avg");
    apf::Vector3 PG;
    double aspect_ratio = 8;
    apf::Field* shock_param = m->findField("Shock Param");

    apf::Vector3 v_mag = apf::Vector3(0.0,0.0,0.0);
    apf::MeshEntity* v;
    apf::MeshIterator* vit = m->begin(0);
    while ((v = m->iterate(vit))) {
      /*
      Game plan: at each vertex, check if an adjacent element contains a shock
      if yes - perform anisotropic refinement by PG direction and VMS size for thickness
              tangent direction size defined by aspect ratio
      else - perform isotropic refinement by the VMS size
      */
      apf::Adjacent adj_elm;
      m->getAdjacent(v, m->getDimension(), adj_elm);
      bool shock = false;
      for (std::size_t i = 0; i < adj_elm.getSize(); ++i) {
        double param = apf::getScalar(shock_param, adj_elm[i], 0);  
        if(param > 0.5) shock = true;
      }
 
      apf::getVector(sizes,v,0,v_mag);
      pVertex meshVertex = reinterpret_cast<pVertex>(v);
      bool aniso = in.anisotropicShockAdaptation;

      if(!EN_isOnPartBdry(meshVertex)){
        if(aniso && shock){
          // anisotropic refinement based on PG direction
          double aspect_ratio = in.anisotropicShockAR;

          apf::getComponents(PG_avg,v,0,&PG[0]);
          double PG_mag = sqrt(PG[0]*PG[0]+PG[1]*PG[1]+PG[2]*PG[2]);
          double n1[3] = {PG[0]/PG_mag, PG[1]/PG_mag, PG[2]/PG_mag};

          double e[3] = {1,0,0}; 
          if(n1[1] < n1[0]){
            if(n1[2] < n1[1]){
              e[0] = 0;
              e[2] = 1;
            }
            else{
              e[0] = 0;
              e[1] = 1;
            } 
          }
          else if (n1[2] < n1[0]){
            e[0] = 0;
            e[2] = 1;
          }

          double t1[3] = {n1[1]*e[2]-e[1]*n1[2], -(n1[0]*e[2]-e[0]*n1[2]), n1[0]*e[1]-e[0]*n1[1]};
          double t1_mag = sqrt(t1[0]*t1[0]+t1[1]*t1[1]+t1[2]*t1[2]);
          t1[0] /= t1_mag; t1[1] /= t1_mag; t1[2] /= t1_mag;
          double t2[3] = {n1[1]*t1[2]-t1[1]*n1[2], -(n1[0]*t1[2]-t1[0]*n1[2]), n1[0]*t1[1]-t1[0]*n1[1]};
          double t2_mag = sqrt(t2[0]*t2[0]+t2[1]*t2[1]+t2[2]*t2[2]);
          t2[0] /= t2_mag; t2[1] /= t2_mag; t2[2] /= t2_mag;

          // set vector mags for element size in those directions
          // n1[0] *= v_mag[0]; n1[1] *= v_mag[0]; n1[2] *= v_mag[0];
          // t1[0] *= v_mag[0]*aspect_ratio; t1[1] *= v_mag[0]*aspect_ratio; t1[2] *= v_mag[0]*aspect_ratio;
          // t2[0] *= v_mag[0]*aspect_ratio; t2[1] *= v_mag[0]*aspect_ratio; t2[2] *= v_mag[0]*aspect_ratio;
          // set vector mags for element size in those directions
          // n1[0] *= v_mag[0]/aspect_ratio; n1[1] *= v_mag[0]/aspect_ratio; n1[2] *= v_mag[0]/aspect_ratio;
          // t1[0] *= v_mag[0]; t1[1] *= v_mag[0]; t1[2] *= v_mag[0];
          // t2[0] *= v_mag[0]; t2[1] *= v_mag[0]; t2[2] *= v_mag[0];

          double n_val = v_mag[0];
          double t_val = v_mag[0]*in.anisotropicShockAR;


          // double n_val = 0.025/8;
          // double t_val = 0.025;
          n1[0] *= n_val; n1[1] *= n_val; n1[2] *= n_val;
          t1[0] *= t_val; t1[1] *= t_val; t1[2] *= t_val;
          t2[0] *= t_val; t2[1] *= t_val; t2[2] *= t_val;


          t1_mag = sqrt(t1[0]*t1[0]+t1[1]*t1[1]+t1[2]*t1[2]);
          t2_mag = sqrt(t2[0]*t2[0]+t2[1]*t2[1]+t2[2]*t2[2]);
          double n1_mag = sqrt(n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2]);

          double anisoSize[3][3] = {{n1[0],n1[1],n1[2]},{t1[0],t1[1],t1[2]},{t2[0],t2[1],t2[2]}}; 
          // if(PCU_Comm_Self()==0){
          //   std::cout << "n1:\t" <<  n1[0] << "\t" << n1[1] << "\t" << n1[2] << "\tmag: " << n1_mag << std::endl;
          //   std::cout << "t1:\t" <<  t1[0] << "\t" << t1[1] << "\t" << t1[2] << "\tmag: " << t1_mag << std::endl;
          //   std::cout << "t2:\t" <<  t2[0] << "\t" << t2[1] << "\t" << t2[2] << "\tmag: " << t2_mag << std::endl << std::endl;
          // }

          MSA_setAnisoVertexSize(adapter, meshVertex, anisoSize);
        }
        else {
          // isotropic refinement based on VMS error
          MSA_setVertexSize(adapter, meshVertex, v_mag[0]);
        }
      }
      else if(!aniso){
        MSA_setVertexSize(adapter,meshVertex,v_mag[0]);
      }

      // MS_setMaxAnisoRatio(pACase cs, double ratio);

      // apf::getComponents(position,v,0,&pos[0]);
      // double h_max = 0.0875/10;
      // // double h_min = 0.0875/96;
      // double h_min = h_max / 32;

      // double L1 = (-0.0875+-0.02)/2;
      // double delta = 0.0015;
      // double theta = 45*3.14159265/180;
      // double tan_theta = 1;


      // apf::getVector(sizes,v,0,v_mag);
      // v_mag[0] = h_min;
      // double anisoSize[3][3] = {{v_mag[0]*32,0,0},{0,v_mag[0]/sqrt(2),-v_mag[0]/sqrt(2)},{0,v_mag[0]*32/sqrt(2),v_mag[0]*32/sqrt(2)}};
      // v_mag[0] = h_max;
      // // pVertex meshVertex = reinterpret_cast<pVertex>(v);
      // if (pos[1] < L1+delta +pos[2]/tan_theta  && pos[1] > L1-delta + pos[2]/tan_theta){
      //   MSA_setAnisoVertexSize(adapter, meshVertex, anisoSize);
      // }
      // else{
      //   MSA_setVertexSize(adapter, meshVertex, v_mag[0]);
      // }
    }
    m->end(vit);

     /* add mesh smooth/gradation function here */
    pc::addSmoother(m, in.gradingFactor);

    /* sync mesh size over partitions */
    // pc::syncMeshSize(m, sizes);

    /* use current size field */
    if(!PCU_Comm_Self())
      printf("Start mesh adapt of setting size field\n");


    /* write error and mesh size */
    pc::writeSequence(m, in.timeStepNumber, "error_mesh_size_");

    /* set fields to be mapped */
    PList_clear(sim_fld_lst);
    if (in.solutionMigration) {
      sim_fld_lst = getSimFieldList(in, m);
      MSA_setMapFields(adapter, sim_fld_lst);
    }
  }

  void runMeshAdapter(ph::Input& in, apf::Mesh2*& m, apf::Field*& orgSF, int step) {
    /* use the size field of the mesh before mesh motion */
    apf::Field* szFld = orgSF;

    if(in.simmetrixMesh == 1) {
      if (in.writeSimLog)
        Sim_logOn("sim_mesh_adaptation.log");
      pProgress progress = Progress_new();
      Progress_setDefaultCallback(progress);

      apf::MeshSIM* sim_m = dynamic_cast<apf::MeshSIM*>(m);
      pParMesh sim_pm = sim_m->getMesh();
      pMesh pm = PM_mesh(sim_pm,0);

      gmi_model* gmiModel = sim_m->getModel();
      pGModel model = gmi_export_sim(gmiModel);

      // declaration
      VIter vIter;
      pVertex meshVertex;

      pProgress progress_tmp = Progress_new();
      Progress_setDefaultCallback(progress_tmp);
      apf::MeshSIM* apf_msim = dynamic_cast<apf::MeshSIM*>(m);
      pParMesh ppm = apf_msim->getMesh();
  
      
      /*Begin shock detection code: Written by Isaac Tam, 2020*/ 
      
      PM_write(ppm, "ShockInd.sms",progress_tmp);
 
      apf::Field* S_Ind = m->findField("Shock_Ind");
      apf::Field* P_Filt = m->findField("P_Filt");
      apf::Field* vms_err = m->findField("VMS_error");\

      apf::Field* sol = m->findField("solution");
      apf::Field* td_sol = m->findField("time derivative of solution");
     

      //mesh resampling hacks
      bool mesh_resample = false; // resample mesh from solution transfer in paraview
      if (mesh_resample && step == 1){ //just the initial step
          
          //read in the resampled data
          std::ifstream in_stream("outpoints.csv");
          if (!in_stream.good()){
               std::cerr<< "can't find the resampled data"<<std::endl;
               exit(1);
          }
          std::string debugstring1; in_stream >> debugstring1; //remove header
          std::string Lin_Parse; std::string delimin=",";
          
          std::map < int, std::vector<double> >  resample_map; // vert_ID, solution vec

          while(in_stream >> Lin_Parse){
               std::vector <std::string > working_vec = ParseDeliminatedString(Lin_Parse, delimin);
               // populate map
               if(PCU_Comm_Self() == std::stoi(working_vec[5] )){
                   std::vector< double> dropin;
                   dropin.push_back(std::stod(working_vec[0])); dropin.push_back(std::stod(working_vec[1]));
                   dropin.push_back(std::stod(working_vec[2])); dropin.push_back(std::stod(working_vec[3]));
                   dropin.push_back(std::stod(working_vec[4]));

                   resample_map[std::stoi(working_vec[6])] = dropin;

                   //if (!PCU_Comm_Self()) std::cout << working_vec[5] <<" " <<working_vec[6]<< " " << dropin[0] <<std::endl;
               }
          }
          

          // solution can be changed here
          apf::MeshEntity* vert_tmp;
          apf::MeshIterator* v_iter = m->begin(0);
          while ((vert_tmp = m->iterate(v_iter))) {
               // grab solution          
               apf::NewArray<double> solu_tmp(apf::countComponents(sol));
               apf::getComponents(sol, vert_tmp, 0, &solu_tmp[0]);
               pEntity ent1 = reinterpret_cast<pEntity>(vert_tmp);
               int v_ID = EN_id(ent1);

               std::vector<double> insert_vec= resample_map[v_ID];

               solu_tmp[0]=insert_vec[0]; solu_tmp[1]=insert_vec[1];
               solu_tmp[2]=insert_vec[2]; solu_tmp[3]=insert_vec[3];
               solu_tmp[4]=insert_vec[4];

               apf::setComponents(sol, vert_tmp, 0, &solu_tmp[0]);
          }
        pc::writeSequence(m,in.timeStepNumber, "resampled_prior_");
        // end mesh resampling hacks
      }

      /* Get smooth pressure gradient field, chef level shock detector*/
     
      apf::Field* PG_avg= apf::createFieldOn(m,"PG_avg",apf::VECTOR);
      apf::Field* shk_det= apf::createFieldOn(m,"shk_det",apf::SCALAR);
      apf::Field* num_elms= apf::createFieldOn(m,"num_elms",apf::SCALAR);
 
      apf::NewArray<double> holder(apf::countComponents(P_Filt));
      apf::Vector3 PG_add = apf::Vector3(0.0,0.0,0.0);

      //parallel implementation for pressure gradient field - Steven Spreizer, 2023
      PCU_Comm_Begin();
 
      apf::MeshEntity* v_tmp;
      apf::MeshIterator* v_itr = m->begin(0);
      while ((v_tmp = m->iterate(v_itr))) {
        //get adjacent nodes
        PG_add[0]=0.0; PG_add[1]=0.0; PG_add[2]=0.0;
        apf::Adjacent Adja;
        m->getAdjacent(v_tmp,3,Adja);
        int num_elm = Adja.getSize(); //number of elements for current sum
          
        //loop for contributions from neighboring elements
        for (size_t i=0; i<num_elm; i++){
          apf::getComponents(P_Filt, Adja[i], 0, &holder[0]);
          PG_add[0]=PG_add[0]+holder[0];
          PG_add[1]=PG_add[1]+holder[1];
          PG_add[2]=PG_add[2]+holder[2];
        }

        //set value
        apf::setVector(PG_avg, v_tmp, 0, PG_add); 
        apf::setScalar(num_elms,v_tmp,0, Adja.getSize());

        //pack up data for communications
        if(!m->isOwned(v_tmp)){
          apf::Copies remotes;
          m->getRemotes(v_tmp,remotes);
          int owningPart = m->getOwner(v_tmp);

          PCU_COMM_PACK(owningPart,remotes[owningPart]); //send entity
          PCU_COMM_PACK(owningPart,num_elm); //send int 
          PCU_COMM_PACK(owningPart,PG_add); //send apf::Vector3
        }
      }
      m->end(v_itr);

      PCU_Comm_Send();

      apf::MeshEntity* ent;
      int received_num_elm;
      apf::Vector3 received_PG;
      apf::Vector3 current_PG;
      int current_num_elm;

      while(PCU_Comm_Receive()){
        //receive comms 
        PCU_COMM_UNPACK(ent);
        PCU_COMM_UNPACK(received_num_elm);
        PCU_COMM_UNPACK(received_PG);

        if(!m->isOwned(ent)){
          std::cout << "ERROR: Data sent to non-owner entity" << std::endl;
          std::exit(1);
        }

        apf::getComponents(PG_avg,ent,0,&current_PG[0]);
        current_num_elm = apf::getScalar(num_elms,ent,0);
        int total_elms = current_num_elm + received_num_elm;
        //compute accumulation from remotes
        current_PG[0] = current_PG[0] + received_PG[0];
        current_PG[1] = current_PG[1] + received_PG[1];
        current_PG[2] = current_PG[2] + received_PG[2];

        //update current fields (on owner)
        apf::setComponents(PG_avg, ent, 0, &current_PG[0]); 
        apf::setScalar(num_elms,ent,0, total_elms);
      } //end receive comms

      //perform averaging
      int total_elms;
      v_itr = m->begin(0);
      while((v_tmp = m->iterate(v_itr))){
        if(m->isOwned(v_tmp)){
          total_elms = apf::getScalar(num_elms,v_tmp,0);
          apf::getComponents(PG_avg,v_tmp,0,&current_PG[0]);
          current_PG[0] /= (double)total_elms;
          current_PG[1] /= (double)total_elms;
          current_PG[2] /= (double)total_elms;
          apf::setComponents(PG_avg, v_tmp, 0, &current_PG[0]); 
        }
      }
      m->end(v_itr);

      //synchronize data on all processes
      apf::synchronize(PG_avg);
      apf::synchronize(num_elms);

      /*
      Verification code for PG parallel implementation
      - send from all remotes back to owner and make sure values are consistent
      */
      bool verify = false;
      if(verify){
        PCU_Comm_Begin();
        v_itr = m->begin(0);
        while((v_tmp = m->iterate(v_itr))){
          //pack messages on all non-owners
          if(!m->isOwned(v_tmp)){
            apf::getComponents(PG_avg,v_tmp,0,&current_PG[0]);
            current_num_elm = apf::getScalar(num_elms,v_tmp,0);

            apf::Copies remotes;
            m->getRemotes(v_tmp,remotes);
            int owningPart = m->getOwner(v_tmp);

            PCU_COMM_PACK(owningPart,remotes[owningPart]); //send entity
            PCU_COMM_PACK(owningPart,current_num_elm); //send int 
            PCU_COMM_PACK(owningPart,current_PG); //send apf::Vector3
          }
        }
        m->end(v_itr);

        PCU_Comm_Send();

        while(PCU_Comm_Receive()){
          //receive comms 
          PCU_COMM_UNPACK(ent);
          PCU_COMM_UNPACK(received_num_elm);
          PCU_COMM_UNPACK(received_PG);

          if(!m->isOwned(ent)){
            std::cout << "ERROR: Data sent to non-owner entity" << std::endl;
            std::exit(1);
          }

          apf::getComponents(PG_avg,ent,0,&current_PG[0]);
          current_num_elm = apf::getScalar(num_elms,ent,0);

          assert(received_num_elm == current_num_elm);
          assert(abs(current_PG[0] - received_PG[0]) < 1e-12);
          assert(abs(current_PG[1] - received_PG[1]) < 1e-12);
          assert(abs(current_PG[2] - received_PG[2]) < 1e-12);
        }
      } //end if verify
      /*
      End pressure gradient parallel comms
      */

      v_itr = m->begin(0);
      while((v_tmp = m->iterate(v_itr))){
        double loc_det=0.0;
        apf::NewArray<double> sol_tmp(apf::countComponents(sol));
        apf::NewArray<double> td_sol_tmp(apf::countComponents(td_sol));

        apf::getComponents(sol, v_tmp, 0, &sol_tmp[0]);
        apf::getComponents(td_sol, v_tmp, 0, &td_sol_tmp[0]);

        loc_det= sol_tmp[1]*PG_add[0] + sol_tmp[2]*PG_add[1] + sol_tmp[3]*PG_add[2];//term 2
        loc_det= loc_det + td_sol_tmp[0];//term 1
        loc_det= loc_det/ ( sqrt(1.4*287*sol_tmp[4])  );//speed of sound
        loc_det= loc_det/ ( sqrt(PG_add[0]*PG_add[0] + PG_add[1]*PG_add[1] + PG_add[2]*PG_add[2]) );
        // P Grad Mag
        
        apf::setScalar(shk_det, v_tmp, 0, loc_det);
            
      }
      m->end(v_itr);

      
      /* Loop through elements and output IDs that contain a shock after filtering */
      
      int nsd = 3;
      apf::MeshEntity* elm;

      
      apf::NewArray<double> Shock_Ind(apf::countComponents(S_Ind));
      apf::NewArray<double> P_filter(apf::countComponents(P_Filt));
      apf::NewArray<double> VMS_err(apf::countComponents(vms_err));

      //default values for filtering
      double P_thres_max   = 10000000000000000000000.0; // accept the highest values
      double P_thres_min   = 0.0; // accept lowest values
      //P_thres_min = 3529981.18288112; //debug val for DW case @ step 488
      double VMS_thres_max = 10000000000000000000000.0;
      double VMS_thres_min = 0.0;
      
      //read from "Shock.inp"
      std::ifstream in_str("Shock.inp");
      if (!in_str.good()){
          if(!PCU_Comm_Self()){
               std::cout << "Can't open Shock.inp to read. Using defaults.\n" <<std::endl;
          }
      }else{
      std::string parse; 
          while(in_str >> parse){
                 if (parse == "P_thres_max"){
                     in_str >> P_thres_max;
                } else if (parse == "P_thres_min"){
                     in_str >> P_thres_min;
                } else if (parse == "VMS_thres_max"){
                     in_str >> VMS_thres_max;
                } else if (parse == "VMS_thres_min"){
                      in_str >> VMS_thres_min;
                }
          }
      }
      
      //Setup for shock parameter recording
      std::vector<int> ShkIDs;
      std::vector<apf::Vector3 > ShkLocs;
      std::vector<double> ShkFilt;
      
      apf::Field* Shock_Param = apf::createField(m, "Shock Param", apf::SCALAR, apf::getConstant(nsd));
      apf::Field* Shock_IDs   = apf::createField(m, "Shock ID", apf::SCALAR, apf::getConstant(nsd));


      apf::MeshIterator* it = m->begin(nsd);
      while ((elm = m->iterate(it))) {
          apf::getComponents(S_Ind, elm, 0, &Shock_Ind[0]);//phasta det calc
          apf::getComponents(P_Filt, elm, 0, &P_filter[0]);
          apf::getComponents(vms_err, elm, 0, &VMS_err[0]);

          //get momentum vms err
          double moment_err = 0.0;
          moment_err = sqrt(VMS_err[1]*VMS_err[1]
                           +VMS_err[2]*VMS_err[2]
                           +VMS_err[3]*VMS_err[3]);

          //actual filtering
          if (P_thres_max > P_filter[3] && P_filter[3] > P_thres_min){//pressure filter, between max and min
           if ( VMS_thres_max > moment_err && moment_err > VMS_thres_min) {//similar VMS_filter
               
               apf::Adjacent Adja;
               m->getAdjacent(elm,0,Adja);
               //find iso surface elms
               bool loc_gt = false; bool loc_lt=false; double loc_check=0.0;
               for (size_t i=0; i<Adja.getSize(); i++){
                    apf::getComponents(shk_det, Adja[i], 0, &loc_check);
                    if (loc_check >= 1.0){
                         loc_gt=true;
                    }
                    if (loc_check <= 1.0){
                         loc_lt=true;
                    }

               }

               //if (Shock_Ind[0] >1 && Shock_Ind[1] <1){//iso surface elements from phasta
               if (loc_gt && loc_lt){ // iso surface elements from chef          
                         
                         apf::Vector3 x3p = apf::Vector3(0.0,0.0,0.0);
                         m->getPoint(Adja[0],0,x3p);//for geom filter
                         double x3p_r= sqrt(x3p[1]*x3p[1] + x3p[2]*x3p[2]);
                         
                         if(true){// !(x3p[0] < 2.01 && x3p_r < 0.07)){ // geom filtering
                              pEntity ent = reinterpret_cast<pEntity>(elm);
                              int rID = EN_id(ent);
                              
                              apf::Vector3 CellCent = apf::getLinearCentroid(m, elm);
                              ShkLocs.push_back(CellCent);

                              ShkIDs.push_back(rID);
                              ShkFilt.push_back(P_filter[3]);
                              apf::setScalar(Shock_Param, elm, 0, 1.0);
                              apf::setScalar(Shock_IDs, elm, 0, rID);
                         }
               }
             }
          }
      }
      m->end(it);
      
      std::string ShkFile ="ShockElms-"+ std::to_string(PCU_Comm_Self()) + ".txt";
      std::ofstream ShockOut(ShkFile);
      std::setprecision(std::numeric_limits<double>::digits10+1);
      if(!ShockOut.good()){
          std::cerr << "Can't open "<< ShkFile <<" to write.\n" <<std::endl;
          exit(1);
      }
      double tmpO1; double tmpO2; double tmpO3;
      for (unsigned int i=0; i<ShkIDs.size();i++){
          ShockOut <<ShkIDs[i] << ",";
          tmpO1=ShkLocs[i][0]; tmpO2=ShkLocs[i][1]; tmpO3=ShkLocs[i][2];
          ShockOut << tmpO1 << "," << tmpO2 <<"," <<tmpO3 << "," << ShkFilt[i] <<std::endl;
      }
      

     /* Shock System Identification process  */

      apf::Field* SurfID= apf::createField(m,"Surf ID", apf::SCALAR, apf::getConstant(nsd)); 

      double work_val=0.0; std::vector<apf::MeshEntity*> LocShkElms;

      it=m->begin(nsd);
      while(elm=m->iterate(it)){ //loop over elms
          apf::setScalar(SurfID,elm,0,0.0); //set surfID to zero
          apf::getComponents(Shock_Param, elm, 0, &work_val);
          if (work_val>0.0){
               LocShkElms.push_back(elm);     
          }

      }
      m->end(it);
      
      double all_s_elms=LocShkElms.size();
      PCU_Add_Doubles(&all_s_elms,1);

      if(!PCU_Comm_Self()){
          std::cout<< "Shock Elms count: "  <<all_s_elms<<std::endl;     
      }
     

      /*
      // LocGraph[elm]= pair (Neighbors set, surfID) 
      std::map<apf::MeshEntity*, std::pair<std::set<apf::MeshEntity*>,int> > LocGraph;

      int n_layers=5; int maxcount=0;
      for (int i=0; i<LocShkElms.size();i++){//populate the graph

          std::set<apf::MeshEntity*> RawSet = GetN_Neighs(LocShkElms[i],m,n_layers);
          LocGraph[LocShkElms[i]].first=RawSet;
          LocGraph[LocShkElms[i]].second=0;

          for (std::set<apf::MeshEntity*>::iterator i1=RawSet.begin(); i1!=RawSet.end(); i1++){
                    double work_val=0.0; apf::getComponents(Shock_Param,*i1,0,&work_val);
                    if ( work_val > 0.0){
                         LocGraph[LocShkElms[i]].first.insert(*i1);   
                    }
          }
          if (maxcount < LocGraph[LocShkElms[i]].first.size()){ maxcount = LocGraph[LocShkElms[i]].first.size();}
      }


      if(!PCU_Comm_Self()){
          std::cout << "most neighbours: " << maxcount<<std::endl;
          apf::MeshIterator* it1 = m->begin(2);//loop over faces
          apf::MeshEntity* fac; int count=0; int counter=0;
          std::set<apf::MeshEntity*> B_Elms1, B_Elms2; //elms and faces
          while(fac=m->iterate(it1)){
               apf::Adjacent Adja;
               m->getAdjacent(fac,3,Adja);
               if (m->isShared(fac)){//count boundary faces, used to vertices
                    count++;     
               }
               for (size_t i=0; i<Adja.getSize(); i++){
                    if(m->isShared(fac)){
                         B_Elms1.insert(Adja[i]);
                         B_Elms2.insert(fac);
                         apf::setScalar(SurfID,Adja[i],0,-7.0);
                    }
                       
               } 
          }
      
          std::cout <<"Boundary faces | elms | faces: "<< count<<" " << B_Elms1.size()<<" "<< B_Elms2.size() <<std::endl;
          m->end(it1);
      }
      pc::writeSequence(m,in.timeStepNumber, "shock_field_");

      
          Solve the connected components problem.
          loc_suf_ID=1
          loop over shock elms
               if it has non-zero loc_surf_ID, continue to next shock elm
               else, this one belongs to loc_surf_ID
                    propagate the ID out
               loc_suf_ID++;

       
      
      bool ProcessSurfs=false;
      if (ProcessSurfs){
      int cur_surf=1;
      for ( int i=0; i<LocShkElms.size(); i++){
          if (LocGraph[LocShkElms[i]].second!=0){ continue;}
          else{
               PropagateID(LocGraph[LocShkElms[i]].first,LocGraph,cur_surf);//recursively go through neighs, assign surf ID
               cur_surf++;//increment surfID
          }
      }
      double tmp2;
      for (int i=0; i<LocShkElms.size();i++){
          tmp2=LocGraph[LocShkElms[i]].second;
          apf::setScalar(SurfID,LocShkElms[i],0,tmp2);
      }
      }
      */
      
      // create the Simmetrix adapter *
      if(!PCU_Comm_Self())
        printf("Start mesh adapt\n");
      pMSAdapt adapter = MSA_new(sim_pm, 1);
      pPList sim_fld_lst = PList_new();
      setupSimAdapter(adapter, in, m, sim_fld_lst);

      // run the adapter *
      if(!PCU_Comm_Self())
        printf("do real mesh adapt\n");
      MSA_adapt(adapter, progress);
      MSA_delete(adapter);

      // create Simmetrix improver *
      pVolumeMeshImprover vmi = VolumeMeshImprover_new(sim_pm);
      setupSimImprover(vmi, sim_fld_lst);

      // run the improver *
      VolumeMeshImprover_execute(vmi, progress);
      VolumeMeshImprover_delete(vmi);

      PList_clear(sim_fld_lst);
      PList_delete(sim_fld_lst);

      // load balance *
      pc::balanceEqualWeights(sim_pm, progress);


      // write mesh *
      if(!PCU_Comm_Self())
        printf("write mesh after mesh adaptation\n");
      writeSIMMesh(sim_pm, in.timeStepNumber, "sim_mesh_");
      
      
      Progress_delete(progress);

      // transfer data back to apf *
      if (in.solutionMigration)
        transferSimFields(m); 
    }
    else {
      assert(szFld);
      apf::synchronize(szFld);
      apf::synchronize(m->getCoordinateField());
      /* do SCOREC mesh adaptation */
      chef::adapt(m,szFld,in);
      chef::balance(in,m);
    }
    m->verify();
  }

}
