#ifndef PC_SHOCK_H
#define PC_SHOCK_H

#include "pcWriteFiles.h"
#include <SimField.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfSIM.h>
#include <apfMDS.h>
#include <chef.h>
#include <phasta.h>
#include <MeshSimAdapt.h>

// for shock surface ID
#include <armadillo>
#include <unordered_set>
#include <vector>


typedef unordered_set<int>::iterator nbd_iter;
typedef unordered_set<int>& set_ref;

namespace pc {
  struct sElm{
    int ElmID;
    int SurfID = 0;

    double x_pos; 
    double y_pos; 
    double z_pos; //xyz coords
    double svd_val = 0.0; // planarity indicator

    std::unordered_set<int> nbs; //neighbours

    apf::MeshEntity* mesh_ent = nullptr;
  };


  void extendShocks(apf::Mesh2*& m,ph::Input& in);

  void interactionHandling(apf::Mesh2*& m, ph::Input& in);

  int identifyExitFace(apf::Mesh2*& m, apf::MeshEntity* element, apf::Vector3 xt, apf::Vector3 t);

  void ReadData(std::vector<sElm> &AllElms, int proc, int cum_elements);

  void GetNbs(std::vector<sElm>& Elms, double spacing);

  void get_n_layers(int n, set_ref cur_set, set_ref work_set,set_ref ret_set, std::vector<sElm>& Elms);

  void CalcPlanarity(std::vector<sElm>& Elms);

  void calcPlanarity(apf::Mesh2*& m,double spacing,ph::Input& in);

  std::unordered_set<int> outliner(std::unordered_set<int> &start, std::vector<sElm> &All_Elms, int cur_ID);

  void SurfTracer(std::unordered_set<int> &cur_set, std::vector<sElm> &All_Elms, int cur_ID);

  std::vector<int> SurfaceSorter(std::vector<sElm> &All_Elms);

  void WriteData(std::vector<sElm> &Elms);
}


#endif
