#include "pcAdapter.h"
#include "pcShock.h"
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

// // for shock surface ID
// #include <armadillo>
// #include <unordered_set>
#include <vector>

// typedef unordered_set<int>::iterator nbd_iter;
// typedef unordered_set<int>& set_ref;

extern void MSA_setBLSnapping(pMSAdapt, int onoff);

namespace pc {

  /*
  Begin porting of Isaac's shock surface ID code
  */
  // struct sElm{
  //   int ElmID;
  //   int SurfID = 0;

  //   double x_pos; 
  //   double y_pos; 
  //   double z_pos; //xyz coords
  //   double svd_val = 0.0; // planarity indicator

  //   std::unordered_set<int> nbs; //neighbours
  // };
  

  // void ReadData(std::vector<sElm> &AllElms, int proc, int cum_elements)
  // {
  //   std::string file = "ShockElms-" + std::to_string(proc) + ".txt";
  //   std::ifstream in_str(file);
  //   if (!in_str)
  //   {
  //     cerr << "Could not open " << file << " to read\n";
  //     exit(1);
  //   }

  //   string parse;
  //   cout << "Reading " << proc << endl;
  //   string delimiter = ",";
  //   int id_iter = cum_elements;

  //   while (in_str >> parse)
  //   {
  //     vector<string> work_vec;

  //     string s = parse;
  //     size_t pos = 0;
  //     string token;
  //     while ((pos = s.find(delimiter)) != std::string::npos)
  //     {
  //       token = s.substr(0, pos);
  //       work_vec.push_back(token);
  //       s.erase(0, pos + delimiter.length());
  //     }
  //     work_vec.push_back(s); // work vec has the x,y,z coords

  //     sElm drop;
  //     drop.x_pos = stod(work_vec[1]);
  //     drop.y_pos = stod(work_vec[2]);
  //     drop.z_pos = stod(work_vec[3]);

  //     drop.ElmID = id_iter;
  //     id_iter++;

  //     AllElms.push_back(drop);
  //   }

  //   cout << "Done reading " << proc << endl;
  //   return;
  // }

  // void GetNbs(std::vector<sElm>& Elms, double spacing){
  //   //get neighbors = elements within a set distance between them
  //   float cur2 = 0.0;
  //   int comp1 = 0;

  //   std::cout << "Elements size: " << Elms.size() << std::endl;
  //   for (int i = 0; i < Elms.size(); i++)
  //   {
  //     for (int j = i+1; j < Elms.size(); j++)
  //     { // account for symmetry?
  //       if (i == j)
  //       {
  //         continue;
  //       }
  //       double dist = (Elms[i].x_pos - Elms[j].x_pos) * (Elms[i].x_pos - Elms[j].x_pos) 
  //                   + (Elms[i].y_pos - Elms[j].y_pos) * (Elms[i].y_pos - Elms[j].y_pos) 
  //                   + (Elms[i].z_pos - Elms[j].z_pos) * (Elms[i].z_pos - Elms[j].z_pos);
  //       dist = sqrt(dist); // absolute dist between elms
  //       if (dist < spacing)
  //       {
  //         Elms[i].nbs.insert(Elms[j].ElmID);
  //         Elms[j].nbs.insert(Elms[i].ElmID);
  //       }
  //     }
  //     // cur2 = 10 * i / Elms.size();
  //     // cur2 = std::ceil(cur2);
  //     // if (comp1 < cur2)
  //     // {
  //     //   comp1++;
  //     //   std::cout << 10 * cur2 << "%" << std::endl;
  //     // }
  //   }

  //   std::cout << "Done finding neighbors" << std::endl;
  //   return;
  // }

  // void get_n_layers(int n, set_ref cur_set, set_ref work_set,set_ref ret_set, std::vector<sElm>& Elms){
  //   if (n == 0)
  //     return; // recursive stopping condition

  //   for (nbd_iter i = cur_set.begin(); i != cur_set.end(); i++)
  //   {
  //     // insert all of cur_set into ret, since we ignore it in work
  //     if (ret_set.find(*i) == ret_set.end())
  //       ;
  //     ret_set.insert(*i);
  //     for (nbd_iter j = Elms[*i].nbs.begin(); j != Elms[*i].nbs.end(); j++)
  //     {
  //       // insert cur's nbs into work and ret
  //       if (ret_set.find(*j) == ret_set.end())
  //       {
  //         work_set.insert(*j);
  //         ret_set.insert(*j);
  //       }
  //     }
  //   }

  //   n--; // decrement layer counter
  //   cur_set = work_set;
  //   work_set.clear(); // generally set up for the recursion
  //   get_n_layers(n, cur_set, work_set, ret_set, Elms);

  //   return;
  // }
  
  // void CalcPlanarity(std::vector<sElm>& Elms){
  //   // std::cout << "getting planarity" << std::endl;

  //   float cur3 = 0.0;
  //   int percent_done = 0;
  //   int comp2 = 0;
  //   for (int i = 0; i < Elms.size(); i++)
  //   {
  //     percent_done = ceil(100*i/Elms.size());
  //     if(percent_done%10==0){
  //       std::cout << "Planarity calc " << percent_done << "\% completed" << std::endl;
  //     }


  //     // cur3=100*i/Elms.size();
  //     // cur3=ceil(cur3);
  //     // if (i % 100 == 0)
  //     // {
  //     //   // comp2++;
  //     //   cout << i << " so far" << endl;
  //     // }

  //     std::unordered_set<int> tmp1;
  //     std::unordered_set<int> tmp2;
  //     std::unordered_set<int> friends;
  //     tmp1.insert(Elms[i].ElmID);
  //     get_n_layers(7, tmp1, tmp2, friends, Elms); // friends includes itself

  //     arma::mat locMeme(3, friends.size());
  //     int iter = 0; // for the matrix insertion memes

  //     double x_c = 0.0; // Elms[i].x_pos;
  //     double y_c = 0.0; // Elms[i].y_pos;
  //     double z_c = 0.0; // Elms[i].z_pos;

  //     for (nbd_iter j = friends.begin(); j != friends.end(); j++)
  //     { // populate matrix for svd
  //       locMeme(0, iter) = Elms[*j].x_pos - x_c;
  //       locMeme(1, iter) = Elms[*j].y_pos - y_c;
  //       locMeme(2, iter) = Elms[*j].z_pos - z_c;
  //       iter++;
  //     }

  //     for (int j = 0; j < friends.size(); j++)
  //     { // comment loop out to remove centering
  //       x_c += locMeme(0, j);
  //       y_c += locMeme(1, j);
  //       z_c += locMeme(2, j);
  //     }
  //     x_c = x_c / friends.size();
  //     y_c = y_c / friends.size();
  //     z_c = z_c / friends.size(); // this too
  //     for (int k = 0; k < friends.size(); k++)
  //     {
  //       locMeme(0, k) = locMeme(0, k) - x_c;
  //       locMeme(1, k) = locMeme(1, k) - y_c;
  //       locMeme(2, k) = locMeme(2, k) - z_c;
  //     }

  //     arma::vec s = arma::svd(locMeme);
  //     double pln = 0.0; // default to as planar as possible
  //     if (friends.size() > 3)
  //     { // need more than 3 points for a guess
  //       pln = (s(2) * s(2)) / ((s(0) * s(0)) + (s(1) * s(1)) + (s(2) * s(2)));
  //     }
  //     Elms[i].svd_val = pln; // svd_val is the normalized 3rd singular value
  //     // cout << friends.size() << " svd val: " <<pln<<endl;
  //   }
  //   std::cout << "got planarity" << std::endl;

  //   return;
  // }

  // std::unordered_set<int> outliner(std::unordered_set<int> &start, std::vector<sElm> &All_Elms, int cur_ID)
  // {
  //   // int start_size = 0;
  //   // //count number of neighbors of each element in set and sum
  //   // for (nbd_iter i = start.begin(); i != start.end(); i++)
  //   // {
  //   //   start_size += All_Elms[*i].nbs.size();
  //   // }
  //   std::unordered_set<int> ret;
  //   // ret.rehash(start_size); //start with the max size, save rehash time

  //   //iterate over neighboring elements
  //   for (nbd_iter i = start.begin(); i != start.end(); i++)
  //   { 
  //     //grab all of current element's neighbors and iterate
  //     std::unordered_set<int> work = All_Elms[*i].nbs;
  //     for (nbd_iter j = work.begin(); j != work.end(); j++)
  //     { 
  //       //search current set of neighbors and exclude if contained
  //       if (start.find(*j) == start.end())
  //       { 
  //         //check that element is unnumbered
  //         if (All_Elms[*j].SurfID == 0)
  //         { 
  //           //insert element to numbering set
  //           ret.insert(*j);
  //         }
  //       }
  //     }
  //   }
  //   // a set full of the unnumbered neighbours
  //   return ret;
  // }

  // void SurfTracer(std::unordered_set<int> &cur_set, std::vector<sElm> &All_Elms, int cur_ID)
  // {
  //   //loop neighbors and set surface ID to current ID
  //   for (nbd_iter i = cur_set.begin(); i != cur_set.end(); i++)
  //   {
  //     All_Elms[*i].SurfID = cur_ID; // set the ID
  //   }
  //   //fill set with next group of elements to add to current shock system 
  //   std::unordered_set<int> next_set = outliner(cur_set, All_Elms, cur_ID);
  //   if (next_set.size() > 0)
  //   { 
  //     // recurse if set is not empty
  //     SurfTracer(next_set, All_Elms, cur_ID);
  //   }
  //   return;
  // }

  // std::vector<int> SurfaceSorter(std::vector<sElm> &All_Elms)
  // {
  //   int SurfID_iter = 1;
  //   //loop over all elements
  //   for (int i = 0; i < All_Elms.size(); i++)
  //   {
  //     //if non surface ID has been assigned
  //     if (All_Elms[i].SurfID == 0)
  //     {
  //       //create set for all elements in shock system
  //       std::unordered_set<int> begin;
  //       begin.insert(i); // initialize the function w/ a size 1 set

  //       SurfTracer(begin, All_Elms, SurfID_iter);
  //       SurfID_iter++;
  //     }
  //   }

  //   //vector of number of elements in each shock system
  //   std::vector<int> ret;
  //   ret.resize(SurfID_iter);
  //   for (int i = 0; i < All_Elms.size(); i++)
  //   {
  //     ret[All_Elms[i].SurfID]++;
  //   }

  //   return ret;
  // }

  // void WriteData(std::vector<sElm> &Elms){
  //   std::ofstream out_str("processed_shock_elements.txt");
  //   if (!out_str.good())
  //   {
  //     cerr << "Can't open " << "processed_shock_elements.txt" << " to write." << endl;
  //     exit(1);
  //   }

  //   out_str << setprecision(std::numeric_limits<double>::digits10 + 1);
  //   // out_str << "\"vtkOriginalPointIds\",\"Points:0\",\"Points:1\",\"Points:2\",\"Surface_ID\",\"Planarity_Ind\" " << std::endl;

  //   for (int i = 0; i < Elms.size(); i++)
  //   {
  //     out_str << i << ",";
  //     out_str << Elms[i].x_pos << ",";
  //     out_str << Elms[i].y_pos << ",";
  //     out_str << Elms[i].z_pos << ",";
  //     out_str << Elms[i].SurfID << ",";
  //     out_str << Elms[i].svd_val << ",";
  //     out_str << "dummy" << std::endl;
  //   }

  //   return;
  // }

  // void assemble_pts(arma::mat &main, std::vector<arma::vec> &drop)
  // {

  //   for (int i = 0; i < drop.size(); i++)
  //   {
  //     main = arma::join_rows(main, drop[i]);
  //   }
  //   main.shed_col(0);
  // }

  // void shock_id_main(char* argv[]){
  //   std::vector<sElm> All_Elms;
  //   ReadData(All_Elms, argv); // intake data

  //   GetNbs(All_Elms, argv); // Establish connectivity data

  //   CalcPlanarity(All_Elms); // step 3 away and get the svd val

  //   WriteData(All_Elms, argv);
  // }

  /*
  End shock ID code
  */

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
        // if (m->isInClosureOf(me, dme) && v_mag[i] > in.simSizeUpperBound) v_mag[i] = in.simSizeUpperBound;
        if (v_mag[i] > in.simSizeUpperBound) v_mag[i] = in.simSizeUpperBound;
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

  void syncMeshSize(apf::Mesh2*& m){
    apf::Field* sizes = m->findField("sizes");
    assert(sizes);

    if(!PCU_Comm_Self())
      std::cerr << "Begin sync mesh size" << std::endl;
    apf::synchronize(sizes);
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

    // double grad_rate = 1.0;
    // MSA_setSizeGradation(adapter, 0, grad_rate);
	

    /* attach mesh size field */
    phSolver::Input inp("solver.inp", "input.config");
    attachMeshSizeField(m, in, inp);
    apf::Field* sizes = m->findField("sizes");
    assert(sizes);
	
    /* initial ctcn field */
    pc::initializeCtCn(m);

    /* apply upper bound */
    pc::applyMaxSizeBound(m, sizes, in);

    /* attach current size field for verification */
    pc::attachCurrentSizeField(m);

    /* apply max number of element */
    // double cn = pc::applyMaxNumberElement(m, sizes, in);

    /* scale mesh if reach time resource bound */
    // pc::applyMaxTimeResource(m, sizes, in, inp);

    /* apply upper bound */
    // pc::applyMaxSizeBound(m, sizes, in);

    /* add mesh smooth/gradation function here */
    // pc::addSmoother(m, in.gradingFactor);

    pc::syncMeshSize(m);

    apf::Field* position = m->findField("motion_coords");
    apf::Vector3 pos;
    apf::Field* PG_avg = m->findField("PG_avg");
    apf::Vector3 PG;
    apf::Field* shock_param = m->findField("Shock Param");
    apf::Field* surf_id = m->findField("Shock_ID");
    apf::Field* shock_vert = m->findField("shock_vert");
    apf::Field* plan_field = m->findField("plan_vert");
    apf::Field* shock_line = m->findField("shock_line");
    apf::Field* shock_line_marker = m->findField("shock_line_marker");

    apf::Field* aniso_size = apf::createFieldOn(m,"aniso_size",apf::MATRIX);
    apf::Field* shk_id = apf::createField(m,"shk_id",apf::SCALAR,apf::getConstant(3));

    apf::Matrix3x3 an_size;

    apf::Vector3 v_mag = apf::Vector3(0.0,0.0,0.0);
    apf::MeshEntity* v;
    apf::MeshIterator* vit = m->begin(0);
    while ((v = m->iterate(vit))) {
      //set isotropic/anisotropic size field at vertices

      double shock = apf::getScalar(shock_vert,v,0);
      double planarity = apf::getScalar(plan_field,v,0);
      double scale_factor = in.sizeScaleFactor;
 
      apf::getVector(sizes,v,0,v_mag);
      pVertex meshVertex = reinterpret_cast<pVertex>(v);
      bool aniso = in.anisotropicShockAdaptation;
      double aspect_ratio = in.anisotropicShockAR;

      apf::Vector3 pt;
      m->getPoint(v,0,pt);
      double iso_marker = apf::getScalar(shock_line_marker,v,0);

      if(aniso && shock && iso_marker < 0.4){
        //this section is unverified
        //vertices which get 2 short component, 1 long (intersecting shocks)
        apf::getComponents(PG_avg,v,0,&PG[0]);
        apf::Vector3 n = PG.normalize();
        apf::Vector3 t1;
        apf::getVector(shock_line,v,0,t1);
        apf::Vector3 t2 = apf::cross(n,t1).normalize();

        n = n * v_mag[0] / aspect_ratio * scale_factor;
        t1 = t1 * v_mag[0] * scale_factor;
        t2 = t2 * v_mag[0] / aspect_ratio * scale_factor;

        double anisoSize[3][3] = {{n[0],n[1],n[2]},{t1[0],t1[1],t1[2]},{t2[0],t2[1],t2[2]}}; 

        MSA_setAnisoVertexSize(adapter, meshVertex, anisoSize);
      }
      else if(aniso && shock){
        // vertices which get 1 short component, 2 long
        // rewrite to use apf::Vector calculations as opposed to by hand
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

        double n_val = v_mag[0]/aspect_ratio*scale_factor;
        double t_val = v_mag[0]*scale_factor;
        
        n1[0] *= n_val; n1[1] *= n_val; n1[2] *= n_val;
        t1[0] *= t_val; t1[1] *= t_val; t1[2] *= t_val;
        t2[0] *= t_val; t2[1] *= t_val; t2[2] *= t_val;


        t1_mag = sqrt(t1[0]*t1[0]+t1[1]*t1[1]+t1[2]*t1[2]);
        t2_mag = sqrt(t2[0]*t2[0]+t2[1]*t2[1]+t2[2]*t2[2]);
        double n1_mag = sqrt(n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2]);

        apf::Vector3 Vecs[3];
        Vecs[0] = apf::Vector3(n1[0],n1[1],n1[2]);
        Vecs[1] = apf::Vector3(t1[0],t1[1],t1[2]);
        Vecs[2] = apf::Vector3(t2[0],t2[1],t2[2]);
        an_size = apf::Matrix3x3(Vecs);

        apf::setMatrix(aniso_size,v,0,an_size);

        double anisoSize[3][3] = {{n1[0],n1[1],n1[2]},{t1[0],t1[1],t1[2]},{t2[0],t2[1],t2[2]}}; 
        v_mag[0] = n_val; v_mag[1] = t_val; v_mag[2] = t_val;
        MSA_setAnisoVertexSize(adapter, meshVertex, anisoSize);
      }
      else {
        // isotropic refinement based on VMS error
        an_size = apf::Matrix3x3(v_mag[0],0,0,0,v_mag[0],0,0,0,v_mag[0]);
        apf::setMatrix(aniso_size,v,0,an_size);
        MSA_setVertexSize(adapter, meshVertex, v_mag[0]*scale_factor);
      }
    } //end iterate over vertices
    m->end(vit);

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
      // bool mesh_resample = false; // resample mesh from solution transfer in paraview
      int mesh_resample = in.manualMeshResample;
      if (mesh_resample && step == 1){ // just the initial step
        std::cout << "Enter manual solution migration" << std::endl;

        // read in the resampled data
        std::ifstream in_stream("outpoints.csv");
        if (!in_stream.good())
        {
          std::cerr << "can't find the resampled data" << std::endl;
          exit(1);
        }
        std::string debugstring1;
        in_stream >> debugstring1; // remove header
        std::string Lin_Parse;
        std::string delimin = ",";

        // std::map<int, std::vector<double>> resample_map; // vert_ID, solution vec
        std::vector<std::vector<double>> resample_map;

        std::cout << "Begin loop file" << std::endl;
        while (in_stream >> Lin_Parse)
        {
          std::vector<std::string> working_vec = ParseDeliminatedString(Lin_Parse, delimin);
          // populate map

          if (PCU_Comm_Self() == std::stoi(working_vec[5]))
          {
            std::vector<double> dropin;
            dropin.push_back(std::stod(working_vec[0]));
            dropin.push_back(std::stod(working_vec[1]));
            dropin.push_back(std::stod(working_vec[2]));
            dropin.push_back(std::stod(working_vec[3]));
            dropin.push_back(std::stod(working_vec[4]));

            // resample_map[std::stoi(working_vec[6])] = dropin;

            resample_map.push_back(dropin);
            // if (!PCU_Comm_Self()) std::cout << working_vec[5] <<" " <<working_vec[6]<< " " << dropin[0] <<std::endl;
          }
        }
        std::cout << "Begin iterate mesh to change solution" << std::endl;

        // solution can be changed here
        apf::MeshEntity *vert_tmp;
        apf::MeshIterator *v_iter = m->begin(0);
        while ((vert_tmp = m->iterate(v_iter)))
        {
          // grab solution
          apf::NewArray<double> solu_tmp(apf::countComponents(sol));
          apf::getComponents(sol, vert_tmp, 0, &solu_tmp[0]);
          pEntity ent1 = reinterpret_cast<pEntity>(vert_tmp);
          int v_ID = EN_id(ent1);

          std::vector<double> insert_vec = resample_map[v_ID];

          // std::cout <<  insert_vec[0] << std::endl;
          // std::cout <<  insert_vec[1] << std::endl;
          // std::cout <<  insert_vec[2] << std::endl;
          // std::cout <<  insert_vec[3] << std::endl;
          // std::cout <<  insert_vec[4] << std::endl;

          solu_tmp[0] = insert_vec[0];
          solu_tmp[1] = insert_vec[1];
          solu_tmp[2] = insert_vec[2];
          solu_tmp[3] = insert_vec[3];
          solu_tmp[4] = insert_vec[4];

          apf::setComponents(sol, vert_tmp, 0, &solu_tmp[0]);
        }
        pc::writeSequence(m, in.timeStepNumber, "resampled_prior_");
      }// end mesh resampling hacks

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
      } //end PG accumulation and comm packaging iteration
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
      } //end averaging iteration
      m->end(v_itr);

      //synchronize data on all processes
      apf::synchronize(PG_avg);
      apf::destroyField(num_elms);    

      //iterate over vertices and calculate normal mach number (Lovely & Haimes), store in shk_det
      apf::Vector3 PG_tmp = apf::Vector3(0.0,0.0,0.0);
      v_itr = m->begin(0);
      while((v_tmp = m->iterate(v_itr))){
        double loc_det=0.0;
        apf::NewArray<double> sol_tmp(apf::countComponents(sol));
        apf::NewArray<double> td_sol_tmp(apf::countComponents(td_sol));

        apf::getComponents(sol, v_tmp, 0, &sol_tmp[0]);
        apf::getComponents(td_sol, v_tmp, 0, &td_sol_tmp[0]);
        apf::getComponents(PG_avg, v_tmp, 0, &PG_tmp[0]);

        loc_det= sol_tmp[1]*PG_tmp[0] + sol_tmp[2]*PG_tmp[1] + sol_tmp[3]*PG_tmp[2];//term 2
        loc_det= loc_det + td_sol_tmp[0];//term 1
        loc_det= loc_det/ ( sqrt(1.4*287*sol_tmp[4])  );//speed of sound
        loc_det= loc_det/ ( sqrt(PG_tmp[0]*PG_tmp[0] + PG_tmp[1]*PG_tmp[1] + PG_tmp[2]*PG_tmp[2]) );
        // P Grad Mag
        
        apf::setScalar(shk_det, v_tmp, 0, loc_det);
      } //end iterate setting normal mach number
      m->end(v_itr);

      /* Loop through elements and output IDs that contain a shock after filtering */      
      int nsd = 3;
      apf::MeshEntity* elm;

      apf::NewArray<double> Shock_Ind(apf::countComponents(S_Ind));
      apf::NewArray<double> P_filter(apf::countComponents(P_Filt));
      apf::NewArray<double> VMS_err(apf::countComponents(vms_err));

      //default values for filtering
      //pressure gradient threshold
      double P_thres_max   = 10000000000000000000000.0; // accept the highest values
      double P_thres_min   = 0.0; // accept lowest values
      //VMS error threshold
      double VMS_thres_max = 10000000000000000000000.0;
      double VMS_thres_min = 0.0;
      
      //read from "Shock.inp"
      std::ifstream in_str("Shock.inp");
      if (!in_str.good()){
        if (!PCU_Comm_Self()){
          std::cout << "Can't open Shock.inp to read. Using defaults.\n" 
                    << std::endl;
        }
      }
      else{
        std::string parse;
        while (in_str >> parse){
          if (parse == "P_thres_max"){
            in_str >> P_thres_max;
          }
          else if (parse == "P_thres_min"){
            in_str >> P_thres_min;
          }
          else if (parse == "VMS_thres_max"){
            in_str >> VMS_thres_max;
          }
          else if (parse == "VMS_thres_min"){
            in_str >> VMS_thres_min;
          }
        }
      }

      //Setup for shock parameter recording
      std::vector<int> ShkIDs;
      std::vector<apf::Vector3 > ShkLocs;
      std::vector<double> ShkFilt;
      
      //1 or 0 indicator of shock in element
      apf::Field* Shock_Param = apf::createField(m, "Shock_Param", apf::SCALAR, apf::getConstant(nsd));
      //shock system ID 
      apf::Field* Shock_IDs   = apf::createField(m, "Shock_ID", apf::SCALAR, apf::getConstant(nsd));
      //planarity indicator
      apf::Field* plan_val = apf::createField(m, "planarity", apf::SCALAR, apf::getConstant(nsd));
      
      //iterate over elements to set shock parameter
      int num_shock_elms = 0;
      apf::MeshIterator* it = m->begin(nsd);
      while ((elm = m->iterate(it))) {
        apf::getComponents(S_Ind, elm, 0, &Shock_Ind[0]);//phasta det calc
        apf::getComponents(P_Filt, elm, 0, &P_filter[0]);
        apf::getComponents(vms_err, elm, 0, &VMS_err[0]);

        apf::setScalar(Shock_Param,elm,0,0.0);
        apf::setScalar(Shock_IDs,elm,0,0.0);
        apf::setScalar(plan_val,elm,0,-1.0);

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
            //loop adjacent vertices to see if element contains a normal 
            //mach number of 1 => element contains shock
            for (size_t i=0; i<Adja.getSize(); i++){
              apf::getComponents(shk_det, Adja[i], 0, &loc_check);
              if (loc_check >= 1.0){
                loc_gt=true;
              }
              if (loc_check <= 1.0){
                loc_lt=true;
              }
            }

            if (loc_gt && loc_lt){ // iso surface elements from chef    
              /* option for geometric filtering as well:
              if desired add if condtiional to only look at elements in filter          
              apf::Vector3 x3p = apf::Vector3(0.0,0.0,0.0);
              m->getPoint(Adja[0],0,x3p);//for geom filter
              double x3p_r= sqrt(x3p[1]*x3p[1] + x3p[2]*x3p[2]);
              !(x3p[0] < 2.01 && x3p_r < 0.07)){ // geom filtering
              */
              pEntity ent = reinterpret_cast<pEntity>(elm);
              int rID = EN_id(ent);
              
              apf::Vector3 CellCent = apf::getLinearCentroid(m, elm);
              ShkLocs.push_back(CellCent);

              ShkIDs.push_back(rID);
              ShkFilt.push_back(P_filter[3]);
              apf::setScalar(Shock_Param, elm, 0, 1.0);
              // apf::setScalar(Shock_IDs, elm, 0, rID);
              apf::setScalar(Shock_IDs, elm, 0, 0);
              ++num_shock_elms;
            }
          }
        }
      } //end iterate over elements to set shock parameter
      m->end(it);

      //option to output shock information to file
      bool shock_to_file = true;
      if(shock_to_file){
        std::string ShkFile ="ShockElms-"+ std::to_string(PCU_Comm_Self()) + ".txt";
        std::ofstream ShockOut(ShkFile);
        ShockOut << std::setprecision(std::numeric_limits<double>::digits10 + 1);
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
      }

      PCU_Add_Int(num_shock_elms);
      if(!PCU_Comm_Self()){
        std::cout<< "Shock Elms count: " << num_shock_elms <<std::endl;     
      }

      // if(true){
      if(in.shockIDSegment){
        std::cerr << "Entering shock system ID and segmentation code " << std::endl;
        /* Shock system ID and segmentation using Armadillo */
        if(!PCU_Comm_Self()){ //proceed in serial 
          std::vector<sElm> all_shock_elements;
          double spacing = 0.1; //closeness parameter
          
          for(int i = 0; i < PCU_Comm_Peers();i++){
            //loop over each file and add to working vector
            ReadData(all_shock_elements,i,all_shock_elements.size());
          }
          GetNbs(all_shock_elements,spacing);

          std::vector<int> SurfSizes;
          std::cout << "Beginning SurfaceSorter" << std::endl;
          SurfSizes = SurfaceSorter(all_shock_elements);
          std::cout << "SurfaceSorter complete" << std::endl;
 
          int threshold = 10; //threshold num elements to consider shock system as just noise
          for (int i = 1; i < SurfSizes.size(); i++)
          { // check surface sizes
            if(SurfSizes[i] < threshold){
              for(int j = 0; j < all_shock_elements.size(); j++){
                if(all_shock_elements[j].SurfID == i){
                  all_shock_elements[j].SurfID = 0;
                }
              }
            }
            else{
              std::cout << i << " " << SurfSizes[i] << std::endl;
            }
          }
        
          //segmentation procedure
          CalcPlanarity(all_shock_elements);
          WriteData(all_shock_elements);
        }

        //add code to now read back in file, compare element centroid position with positions on current process
        //if element is from that process, read the shock surface ID and write to the field data

        PCU_Barrier();

        std::string file = "processed_shock_elements.txt";
        in_str = std::ifstream(file);
        if(!in_str){
          std::cerr << "Could not open " << file << " to read." << std::endl;
        }

        std::string parse;
        string delimiter = ",";
        int id_iter = 0;

        if(!PCU_Comm_Self()){
          std::cout << "Looping shock elements to add ID to field data" << std::endl;
        }

        while (in_str >> parse)
        {
          std::vector<std::string> work_vec;

          std::string s = parse;
          size_t pos = 0;
          std::string token;
          while ((pos = s.find(delimiter)) != std::string::npos)
          {
            token = s.substr(0, pos);
            work_vec.push_back(token);
            s.erase(0, pos + delimiter.length());
          }
          work_vec.push_back(s); // work vec has the x,y,z coords

          sElm drop;
          double x_pos = stod(work_vec[1]);
          double y_pos = stod(work_vec[2]);
          double z_pos = stod(work_vec[3]);
          int surf_id = stoi(work_vec[4]);
          double plan_ind = stod(work_vec[5]);

          //iterate elements to see if there's a match
          double eps = 1e-7; //tolerance for checking if same point
          apf::MeshEntity* e;
          apf::Vector3 pt; 
          apf::MeshIterator* itr = m->begin(3);
          while((e = m->iterate(itr))){
            pt = apf::getLinearCentroid(m,e);
            bool match = abs(pt[0] - x_pos) < eps && abs(pt[1] - y_pos) < eps 
                      && abs(pt[2] - z_pos) < eps;
            if(match){
              apf::setScalar(Shock_IDs,e,0,surf_id);
              apf::setScalar(plan_val,e,0,plan_ind);

              //undo shock detection of elements that are just noise
              if(surf_id == 0){
                // apf::setScalar(Shock_Param,e,0,0);
              }
              break;
            }
          }
        }

        cout << "["<< PCU_Comm_Self() << "]\tDone reading and adding to field data" << endl;

        std::cerr << "Exiting shock system ID and segmentation code " << std::endl;
      }
      PCU_Barrier();

      if(in.extendShocks){
        //routine to extend the shock detection result
        extendShocks(m,in);

        /* 
        recalculate planarity after extension
        */
        
        double spacing = 0.1;
        calcPlanarity(m,spacing,in);

        /*
        create second PG field for special handling of 2 short-1 long dimension adaptation
        where interactions are found (high planarity)
        */
        interactionHandling(m,in);
      }


      /* Create vertex level shock indicator for aniso adaptation */
      // Mark vertices that have an adjacent shock containing element 
      PCU_Comm_Begin();
      apf::Field* shock_vert = apf::createFieldOn(m,"shock_vert",apf::SCALAR);
      v_itr = m->begin(0);
      while((v_tmp = m->iterate(v_itr))){
        double v_shock = 0;
        apf::Adjacent elements;
        m->getAdjacent(v_tmp,3,elements);
        for(int i = 0; i < elements.size(); ++i){
          if(apf::getScalar(Shock_Param,elements[i],0)){
            v_shock = 1;
            break;
          }
        }
        apf::setScalar(shock_vert,v_tmp,0,v_shock);

        if(!m->isOwned(v_tmp)){
          apf::Copies remotes;
          m->getRemotes(v_tmp,remotes);
          int owningPart = m->getOwner(v_tmp);

          PCU_COMM_PACK(owningPart,remotes[owningPart]); //send entity
          PCU_COMM_PACK(owningPart,v_shock); //send int 
        }
      }
      m->end(v_itr);

      PCU_Comm_Send();

      double v_shock_rec;
      while(PCU_Comm_Receive()){
        PCU_COMM_UNPACK(ent);
        PCU_COMM_UNPACK(v_shock_rec);

        if(!m->isOwned(ent)){
          std::cerr << "ERROR: Data sent to non-owner entity" << std::endl;
          std::exit(1);
        }

        double v_shock_curr = apf::getScalar(shock_vert,ent,0);
        v_shock_curr = v_shock_curr + v_shock_rec;

        apf::setScalar(shock_vert,ent,0,v_shock_curr);
      }

      apf::synchronize(shock_vert);

      /* Vertex level planarity indicator */

      PCU_Comm_Begin();
      /* not implemented correctly need to fix
      need to average only neighbor shock elements not all elements */
      apf::Field* plan_vert = apf::createFieldOn(m,"plan_vert",apf::SCALAR);
      v_itr = m->begin(0);
      while((v_tmp = m->iterate(v_itr))){
        double v_plan = 0;
        apf::Adjacent elements;
        m->getAdjacent(v_tmp,3,elements);
        int count = 0;
        for(int i = 0; i < elements.size(); ++i){
          if(apf::getScalar(plan_val,elements[i],0) >= 0){
            count++;
            v_plan += apf::getScalar(plan_val,elements[i],0);
          }
        }
        v_plan /= count;
        apf::setScalar(plan_vert,v_tmp,0,v_plan);

        if(!m->isOwned(v_tmp)){
          apf::Copies remotes;
          m->getRemotes(v_tmp,remotes);
          int owningPart = m->getOwner(v_tmp);

          PCU_COMM_PACK(owningPart,remotes[owningPart]); //send entity
          PCU_COMM_PACK(owningPart,v_plan); //send int 
        }
      }
      m->end(v_itr);

      PCU_Comm_Send();

      double v_plan_rec;
      while(PCU_Comm_Receive()){
        PCU_COMM_UNPACK(ent);
        PCU_COMM_UNPACK(v_plan_rec);

        if(!m->isOwned(ent)){
          std::cerr << "ERROR: Data sent to non-owner entity" << std::endl;
          std::exit(1);
        }

        double v_plan_curr = apf::getScalar(plan_vert,ent,0);
        v_plan_curr = std::min(v_plan_curr,v_plan_rec);

        apf::setScalar(plan_vert,ent,0,v_plan_curr);
      }

      apf::synchronize(plan_vert);
      /* End vertex level shock indicator */


      /*
      // Shock System Identification process  //
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

      /*
      End Shock Detection Code
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
      if(step > 1)
        MSA_adapt(adapter, progress);
      // MSA_adapt(adapter, progress); //temporary to test solution migration
      MSA_delete(adapter);


      if(step > 1) {
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
      }


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
