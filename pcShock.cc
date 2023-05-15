#include "pcShock.h"
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

// for shock surface ID
#include <armadillo>
#include <unordered_set>
#include <vector>
#include <queue>
#include <map>

namespace pc{

	void extendShocks(apf::Mesh2*& m,ph::Input& in){
    std::cout << "Enter shock extension code" << std::endl;
		int nsd = m->getDimension();

		//find fields that are need from mesh
		apf::Field* p_filt = m->findField("P_Filt");
		apf::Field* shock_param = m->findField("Shock_Param");
		apf::Field* planarity = m->findField("planarity");
		apf::Field* shock_id = m->findField("Shock_ID");
    apf::Field* pg_avg = m->findField("PG_avg");
    
    apf::Field* ray_dist = apf::createField(m,"ray_dist",apf::SCALAR,apf::getConstant(nsd));


    //initialize ray_dist field and find max planarity in mesh
    double plan_max = 0;
    double plan_curr = 0;
    apf::MeshEntity* elm;
		apf::MeshIterator* elm_itr = m->begin(nsd);
		while((elm = m->iterate(elm_itr))){
      apf::setScalar(ray_dist,elm,0,0);
      plan_curr = apf::getScalar(planarity,elm,0);
      if(plan_curr > plan_max){
        plan_max = plan_curr;
      }
    }
    m->end(elm_itr);

    PCU_Max_Double(plan_max);
    
    std::cout << "Ray distance field created" << std::endl;

		//useful variables
		double elm_plan;
		double elm_id;
		apf::Vector3 pg;
		apf::Vector3 xt;
		apf::Vector3 xtdt;
		apf::Vector3 n;
		apf::Vector3 t1;
		apf::Vector3 t2;
		apf::Vector3 e;
		apf::Vector3 t_hat;
		apf::Vector3 t;

		//parameters for shooting rays
		double plan_thres = 0.8 * plan_max; 
    std::cout << plan_thres << std::endl;
		int num_rays = 180;
		double angle_increment = 360./(double)num_rays;
		double ray_length = 8; //[m], should be larger than larger diagonal of bounding box
    double sin_theta,cos_theta;

    int count_elements = 0;

    /*
    Iterate elements 
    */

    elm_itr = m->begin(nsd);
    while((elm = m->iterate(elm_itr))){
      //shock param == 3 --> red zone
      //shock param == 2 --> blue zone
      //shock param == 1 --> green zone

      elm_plan = apf::getScalar(planarity,elm,0);

      if(elm_plan > plan_thres){
        apf::setScalar(shock_param,elm,0,3); //red zone
      }
      else{
        apf::Adjacent neighbors;
        apf::getBridgeAdjacent(m,elm,0,3,neighbors);
        for(int j = 0; j < neighbors.size(); j++){
          if(apf::getScalar(shock_param,elm,0) && 
             apf::getScalar(planarity,neighbors[j],0) > plan_thres){
            apf::setScalar(shock_param,elm,0,2); //blue zone
          }
        }
      }
    }
    m->end(elm_itr);

    
    std::cout << "Begin element iteration" << std::endl;

		//iterate over mesh elements 
    elm_itr = m->begin(nsd);
		while((elm = m->iterate(elm_itr))){
			elm_plan = apf::getScalar(planarity,elm,0);
			xt = apf::getLinearCentroid(m,elm); //start of ray

			//begin extension routine if exceeding threshold
      //added geometric filter as well inlet BC giving issues with PG dir
			// if(elm_plan > plan_thres && xt[1] < 0.5){
      if(apf::getScalar(shock_param,elm,0) == 2 && xt[1] < 0.5){
        ++count_elements;

        //perform averaging of pg to find a smoothed p_filt
        apf::Downward pg_vertices;
        m->getDownward(elm,0,pg_vertices);
        apf::Vector3 pg_averaged(0,0,0);
        apf::Vector3 pg_tmp;
        for(int i = 0; i < 4; i++){
          apf::getComponents(pg_avg,pg_vertices[i],0,&pg_tmp[0]);
          pg_averaged = pg_averaged + pg_tmp;
        }
        pg = pg_averaged / 4.;

				// apf::getComponents(p_filt,elm,0,&pg[0]);
				n = pg.normalize();
				e = apf::Vector3(1,0,0);
				if(n[1] < n[0]){
          if(n[2] < n[1]){
            e[0] = 0;
            e[2] = 1;
          }
          else{
            e[0] = 0;
            e[1] = 1;
          } 
        }
        else if (n[2] < n[0]){
          e[0] = 0;
          e[2] = 1;
        }
				t1 = apf::cross(n,e);
				t1 = t1.normalize();
        t2 = apf::cross(n,t1);
        t2 = t2.normalize();

        double shockID = apf::getScalar(shock_id,elm,0);

        //loop all rays
        for(int i_theta = 0; i_theta < num_rays; i_theta++){
          // std::cout << "Element: " << count_elements << "\tray: " << i_theta << "/" << num_rays << std::endl;

          //rays will be shot by moving around unit circle
          sin_theta = sin(angle_increment*i_theta);
          cos_theta = cos(angle_increment*i_theta);

          t_hat = t1*cos_theta + t2*sin_theta; //unit vector for direction of travel
          t = t_hat*ray_length; //trajectory vector
          xtdt = xt + t; //final position
        
          //find which face the ray intersects
          //convention for the ordering of faces and vertices found in apf documentation
          apf::Downward elm_faces;
          m->getDownward(elm,2,elm_faces); //ordered set of faces

          int exit_face_ind = identifyExitFace(m,elm,xt,t);

          if(exit_face_ind == -1){
            continue;
          }

          apf::MeshEntity* exit_face = elm_faces[exit_face_ind];
          apf::MeshEntity* in_face;

          //only enter new elements if they are not classified on model faces
          apf::ModelEntity* me = m->toModel(exit_face);

          int count = 0;
          int max_rounds = 15;

          while(m->getModelType(me) == nsd){ //mesh interior faces classified on model region
            count++;
            //get element we are now entering
            apf::Adjacent face_elms;
            m->getAdjacent(exit_face,3,face_elms);

            apf::MeshEntity* new_elm;
            new_elm = face_elms[0];
            if(new_elm == elm){
              new_elm = face_elms[1];
            }

            in_face = exit_face;

            if(apf::getScalar(shock_param,new_elm,0)==1 && count > 6){
              //break ray tracing if real shock detection result
              break;
            }

            //find distance from start of ray
            apf::Vector3 new_elm_cent = apf::getLinearCentroid(m,new_elm);
            double distance = (new_elm_cent - xt).getLength();
            //
            if(apf::getScalar(shock_id,new_elm,0) <= 0 || apf::getScalar(shock_param,new_elm,0) == 3){
              double cur_dist = apf::getScalar(ray_dist,new_elm,0);
              if(cur_dist == 0 || distance < cur_dist){
                if(apf::getScalar(shock_param,new_elm,0) != 3){
                  apf::setScalar(shock_id,new_elm,0,-1*shockID);
                  apf::setScalar(ray_dist,new_elm,0,distance);
                }
                apf::setVector(p_filt,new_elm,0,pg);

                //after pressure gradient synchronization step so just assigning
                //vertexes to have the element's pressure gradient for mesh adapt
                apf::Downward pg_vtxs;
                m->getDownward(new_elm,0,pg_vtxs);
                for(int v_ind = 0; v_ind < 4; v_ind++){
                  apf::setVector(pg_avg,pg_vtxs[v_ind],0,pg);
                }
              }
            }

            elm = new_elm;
            //determine the next face the ray passes through
            //find which face the ray intersects
            //convention for the ordering of faces and vertices found in apf documentation
            m->getDownward(elm,2,elm_faces); //ordered set of faces

            exit_face_ind = identifyExitFace(m,elm,xt,t);

            if(exit_face_ind == -1){
              break;
            }

            exit_face = elm_faces[exit_face_ind];
            me = m->toModel(exit_face);
          }
        }//end loop over i_theta (rays)
			}
		}
    m->end(elm_itr);

    /*

    code here to backtrack along new shock results and clip to the right 
    distance in the case of interactions
    */


    //iterate over mesh and mark new elements as shocks
    elm_itr = m->begin(nsd);
    while((elm = m->iterate(elm_itr))){
      double read_shock_id = apf::getScalar(shock_id,elm,0);
      if(read_shock_id < 0){
        if(apf::getScalar(shock_param,elm,0) != 3){
          apf::setScalar(shock_param,elm,0,4);
          apf::setScalar(shock_id,elm,0,-1*read_shock_id);
        }
      }
    }
    m->end(elm_itr);

		if(!PCU_Comm_Self()){
			std::cout << "extendShocks complete" << std::endl;
		}
	}

  int identifyExitFace(apf::Mesh2*& m, apf::MeshEntity* elm, apf::Vector3 xt, apf::Vector3 t){
    apf::Vector3 P1;
    apf::Vector3 P2; 
    apf::Vector3 P3;
    apf::Vector3 vec1;
    apf::Vector3 vec2;
    apf::Vector3 cross;
    double dot;
    bool pos;
    bool neg;

    apf::Downward elm_faces;
    apf::Downward elm_vtxs;

    
    if(m->getType(elm) != 4){
      std::cerr << "Current only tets are supported\nEXITING" << std::endl;
      exit(1);
    }

    m->getDownward(elm,2,elm_faces); //ordered set of faces
    m->getDownward(elm,0,elm_vtxs); //ordered set of vtxs

    int exit_face_ind = -1;
    //loop over candidate faces - doing without for loop to make it easier
    //F0
    m->getPoint(elm_vtxs[0],0,P1);
    m->getPoint(elm_vtxs[2],0,P2);
    m->getPoint(elm_vtxs[1],0,P3);
    pos = false;
    neg = false;

    vec1 = P1 - xt;
    vec2 = P2 - xt;
    cross = apf::cross(vec1,vec2);
    dot = cross*t;
    pos = dot >= 0 || pos;
    neg = dot < 0  || neg;

    vec1 = P2 - xt;
    vec2 = P3 - xt;
    cross = apf::cross(vec1,vec2);
    dot = cross*t;
    pos = dot >= 0 || pos;
    neg = dot < 0  || neg;

    vec1 = P3 - xt;
    vec2 = P1 - xt;
    cross = apf::cross(vec1,vec2);
    dot = cross*t;
    pos = dot >= 0 || pos;
    neg = dot < 0  || neg;

    if(pos && !neg){ //all dot products had same sign
      exit_face_ind = 0;
    }

    //F1
    m->getPoint(elm_vtxs[0],0,P1);
    m->getPoint(elm_vtxs[1],0,P2);
    m->getPoint(elm_vtxs[3],0,P3);
    pos = false;
    neg = false;

    vec1 = P1 - xt;
    vec2 = P2 - xt;
    cross = apf::cross(vec1,vec2);
    dot = cross*t;
    pos = dot >= 0 || pos;
    neg = dot < 0  || neg;

    vec1 = P2 - xt;
    vec2 = P3 - xt;
    cross = apf::cross(vec1,vec2);
    dot = cross*t;
    pos = dot >= 0 || pos;
    neg = dot < 0  || neg;

    vec1 = P3 - xt;
    vec2 = P1 - xt;
    cross = apf::cross(vec1,vec2);
    dot = cross*t;
    pos = dot >= 0 || pos;
    neg = dot < 0  || neg;

    if(pos && !neg){
      exit_face_ind = 1;
    }

    //F2
    m->getPoint(elm_vtxs[1],0,P1);
    m->getPoint(elm_vtxs[2],0,P2);
    m->getPoint(elm_vtxs[3],0,P3);
    pos = false;
    neg = false;

    vec1 = P1 - xt;
    vec2 = P2 - xt;
    cross = apf::cross(vec1,vec2);
    dot = cross*t;
    pos = dot >= 0 || pos;
    neg = dot < 0  || neg;

    vec1 = P2 - xt;
    vec2 = P3 - xt;
    cross = apf::cross(vec1,vec2);
    dot = cross*t;
    pos = dot >= 0 || pos;
    neg = dot < 0  || neg;

    vec1 = P3 - xt;
    vec2 = P1 - xt;
    cross = apf::cross(vec1,vec2);
    dot = cross*t;
    pos = dot >= 0 || pos;
    neg = dot < 0  || neg;

    if(pos && !neg){
      exit_face_ind = 2;
    }

    //F3
    m->getPoint(elm_vtxs[0],0,P1);
    m->getPoint(elm_vtxs[3],0,P2);
    m->getPoint(elm_vtxs[2],0,P3);
    pos = false;
    neg = false;

    vec1 = P1 - xt;
    vec2 = P2 - xt;
    cross = apf::cross(vec1,vec2);
    dot = cross*t;
    pos = dot >= 0 || pos;
    neg = dot < 0  || neg;

    vec1 = P2 - xt;
    vec2 = P3 - xt;
    cross = apf::cross(vec1,vec2);
    dot = cross*t;
    pos = dot >= 0 || pos;
    neg = dot < 0  || neg;

    vec1 = P3 - xt;
    vec2 = P1 - xt;
    cross = apf::cross(vec1,vec2);
    dot = cross*t;
    pos = dot >= 0 || pos;
    neg = dot < 0  || neg;

    if(pos && !neg){
      exit_face_ind = 3;
    }

    return exit_face_ind;
  }

  void calcPlanarity(apf::Mesh2*& m, double spacing, ph::Input& in){
    /*
    ONLY WORKS IN SERIAL SO DON'T BE DUMB
    */

    int nsd = m->getDimension();
    apf::MeshIterator* elm_itr = m->begin(nsd);
    apf::MeshEntity* elm;

    apf::Field* shock_param = m->findField("Shock_Param");
		apf::Field* planarity = m->findField("planarity");

    //fill vector of shock elements
    std::vector<sElm> all_shock_elements;
    int elmID = 0;
    while((elm = m->iterate(elm_itr))){
      if(apf::getScalar(shock_param,elm,0) > 0){
        sElm element;
        apf::Vector3 CellCent = apf::getLinearCentroid(m, elm);

        element.x_pos = CellCent[0];
        element.y_pos = CellCent[1];
        element.z_pos = CellCent[2];
        element.ElmID = elmID++;

        element.mesh_ent = elm;

        all_shock_elements.push_back(element);
      }
    }
    m->end(elm_itr);


    //calculate neighbors and planarity
    GetNbs(all_shock_elements,spacing);
    std::vector<int> surf_sizes = SurfaceSorter(all_shock_elements);
    CalcPlanarity(all_shock_elements);

    //loop over vector and write back to the mesh
    for(int i = 0; i < all_shock_elements.size(); i++){
      elm = all_shock_elements[i].mesh_ent;
      double plan = all_shock_elements[i].svd_val;

      apf::setScalar(planarity,elm,0,plan);
    }
  }

  void interactionHandling(apf::Mesh2*& m, ph::Input& in){
    //find fields that are need from mesh
		apf::Field* p_filt = m->findField("P_Filt");
		apf::Field* shock_param = m->findField("Shock_Param");
		apf::Field* planarity = m->findField("planarity");
		apf::Field* shock_id = m->findField("Shock_ID");
    apf::Field* pg_avg = m->findField("PG_avg");
    int nsd = m->getDimension();

    double threshold = 0.4;


    apf::Field* shock_line = apf::createFieldOn(m,"shock_line",apf::VECTOR);
    apf::Field* shock_line_marker = apf::createFieldOn(m,"shock_line_marker",apf::SCALAR);
    apf::MeshIterator* vtx_itr = m->begin(0);
    apf::MeshEntity* vtx;
    while((vtx = m->iterate(vtx_itr))){
      apf::Vector3 temp(0.0,0.0,0.0);
      apf::setVector(shock_line,vtx,0,temp);
      apf::setScalar(shock_line_marker,vtx,0,0);
    }
    m->end(vtx_itr);

    /*
    idea, for each shock elment, search if there is a nearby element with 
    a pressure gradient who's dot product is pretty far off from its own
    - p_filt has already been smoothed and extended for the shocks
    - if condition met, add a vector field to hold vector for line of shock interaction
      and set a marker for this
    - when adapting, check if this marker is true and handle accordingly
      - size unit vectors are PG, this vector, and vector in cross product direction
    */
    apf::MeshIterator* elm_itr = m->begin(nsd);
    apf::MeshEntity* elm;
    while((elm = m->iterate(elm_itr))){
      if(apf::getScalar(shock_param,elm,0) > 0){
        //pull two gets of adjacent elements via vertex bridge
        bool shock_line_bool = false;
        apf::Vector3 elm_pg;
        apf::getVector(p_filt,elm,0,elm_pg);
        apf::Adjacent elements1;
        apf::Vector3 temp_pg;
        apf::Vector3 pg_min;
        pg_min = elm_pg;
        apf::getBridgeAdjacent(m,elm,0,3,elements1);
        double dot = 1;
        double dot_temp;
        for(int i = 0; i < elements1.size() && !shock_line_bool; i++){
          apf::getVector(p_filt,elements1[i],0,temp_pg);
          dot_temp = (elm_pg.normalize()) * (temp_pg.normalize());
          if(dot_temp < dot){
            dot = dot_temp;
            pg_min = temp_pg;
          }
          
          apf::Adjacent elements2;
          apf::getBridgeAdjacent(m,elements1[i],0,3,elements2);
          for(int j = 0; j < elements2.size(); j++){
            apf::getVector(p_filt,elements1[i],0,temp_pg);
            dot_temp = (elm_pg.normalize()) * (temp_pg.normalize());
            
            if(dot_temp < dot){
              dot = dot_temp;
              pg_min = temp_pg;
            }
          }
          
        }
        apf::Adjacent adj_vtx;
        m->getAdjacent(elm,0,adj_vtx);
        for(int k = 0; k < adj_vtx.size(); k++){
          apf::Vector3 cross = apf::cross(elm_pg,pg_min);
          apf::setVector(shock_line,adj_vtx[k],0,cross.normalize());
          apf::setScalar(shock_line_marker,adj_vtx[k],0,dot);
        }
        
      }
    }
    m->end(elm_itr);
  }



	void ReadData(std::vector<sElm> &AllElms, int proc, int cum_elements)
  {
    std::string file = "ShockElms-" + std::to_string(proc) + ".txt";
    std::ifstream in_str(file);
    if (!in_str)
    {
      cerr << "Could not open " << file << " to read\n";
      exit(1);
    }

    string parse;
    cout << "Reading " << proc << endl;
    string delimiter = ",";
    int id_iter = cum_elements;

    while (in_str >> parse)
    {
      vector<string> work_vec;

      string s = parse;
      size_t pos = 0;
      string token;
      while ((pos = s.find(delimiter)) != std::string::npos)
      {
        token = s.substr(0, pos);
        work_vec.push_back(token);
        s.erase(0, pos + delimiter.length());
      }
      work_vec.push_back(s); // work vec has the x,y,z coords

      sElm drop;
      drop.x_pos = stod(work_vec[1]);
      drop.y_pos = stod(work_vec[2]);
      drop.z_pos = stod(work_vec[3]);

      drop.ElmID = id_iter;
      id_iter++;

      AllElms.push_back(drop);
    }

    cout << "Done reading " << proc << endl;
    return;
  }

  void GetNbs(std::vector<sElm>& Elms, double spacing){
    //get neighbors = elements within a set distance between them
    float cur2 = 0.0;
    int comp1 = 0;

    std::cout << "Elements size: " << Elms.size() << std::endl;
    for (int i = 0; i < Elms.size(); i++)
    {
      for (int j = i+1; j < Elms.size(); j++)
      { // account for symmetry?
        if (i == j)
        {
          continue;
        }
        double dist = (Elms[i].x_pos - Elms[j].x_pos) * (Elms[i].x_pos - Elms[j].x_pos) 
                    + (Elms[i].y_pos - Elms[j].y_pos) * (Elms[i].y_pos - Elms[j].y_pos) 
                    + (Elms[i].z_pos - Elms[j].z_pos) * (Elms[i].z_pos - Elms[j].z_pos);
        dist = sqrt(dist); // absolute dist between elms
        if (dist < spacing)
        {
          Elms[i].nbs.insert(Elms[j].ElmID);
          Elms[j].nbs.insert(Elms[i].ElmID);
        }
      }
      // cur2 = 10 * i / Elms.size();
      // cur2 = std::ceil(cur2);
      // if (comp1 < cur2)
      // {
      //   comp1++;
      //   std::cout << 10 * cur2 << "%" << std::endl;
      // }
    }

    std::cout << "Done finding neighbors" << std::endl;
    return;
  }

  void get_n_layers(int n, set_ref cur_set, set_ref work_set,set_ref ret_set, std::vector<sElm>& Elms){
    if (n == 0)
      return; // recursive stopping condition

    for (nbd_iter i = cur_set.begin(); i != cur_set.end(); i++)
    {
      //**********IDK what's going on here ************
      // insert all of cur_set into ret, since we ignore it in work
      if (ret_set.find(*i) == ret_set.end()) ; ret_set.insert(*i);
      for (nbd_iter j = Elms[*i].nbs.begin(); j != Elms[*i].nbs.end(); j++)
      {
        // insert cur's nbs into work and ret
        if (ret_set.find(*j) == ret_set.end())
        {
          work_set.insert(*j);
          ret_set.insert(*j);
        }
      }
    }

    n--; // decrement layer counter
    cur_set = work_set;
    work_set.clear(); // generally set up for the recursion
    get_n_layers(n, cur_set, work_set, ret_set, Elms);

    return;
  }
  
  void CalcPlanarity(std::vector<sElm>& Elms){
    std::cout << "Begin planarity calculation" << std::endl;

    float cur3 = 0.0;
    int percent_done = 0;
    int comp2 = 0;
    for (int i = 0; i < Elms.size(); i++)
    {
      int tenth = Elms.size()/10;
      
      if(i%tenth==0){
        std::cout << "Planarity calc " << i/tenth*10 << "\% completed" << std::endl;
      }


      // cur3=100*i/Elms.size();
      // cur3=ceil(cur3);
      // if (i % 100 == 0)
      // {
      //   // comp2++;
      //   cout << i << " so far" << endl;
      // }

      std::unordered_set<int> tmp1;
      std::unordered_set<int> tmp2;
      std::unordered_set<int> friends;
      tmp1.insert(Elms[i].ElmID);
      get_n_layers(7, tmp1, tmp2, friends, Elms); // friends includes itself

      arma::mat locMeme(3, friends.size());
      int iter = 0; // for the matrix insertion memes

      double x_c = 0.0; // Elms[i].x_pos;
      double y_c = 0.0; // Elms[i].y_pos;
      double z_c = 0.0; // Elms[i].z_pos;

      for (nbd_iter j = friends.begin(); j != friends.end(); j++)
      { // populate matrix for svd
        locMeme(0, iter) = Elms[*j].x_pos - x_c;
        locMeme(1, iter) = Elms[*j].y_pos - y_c;
        locMeme(2, iter) = Elms[*j].z_pos - z_c;
        iter++;
      }

      for (int j = 0; j < friends.size(); j++)
      { // comment loop out to remove centering
        x_c += locMeme(0, j);
        y_c += locMeme(1, j);
        z_c += locMeme(2, j);
      }
      x_c = x_c / friends.size();
      y_c = y_c / friends.size();
      z_c = z_c / friends.size(); // this too
      for (int k = 0; k < friends.size(); k++)
      {
        locMeme(0, k) = locMeme(0, k) - x_c;
        locMeme(1, k) = locMeme(1, k) - y_c;
        locMeme(2, k) = locMeme(2, k) - z_c;
      }

      arma::vec s = arma::svd(locMeme);
      double pln = 0.0; // default to as planar as possible
      if (friends.size() > 3)
      { // need more than 3 points for a guess
        pln = (s(2) * s(2)) / ((s(0) * s(0)) + (s(1) * s(1)) + (s(2) * s(2)));
      }
      Elms[i].svd_val = pln; // svd_val is the normalized 3rd singular value
      // cout << friends.size() << " svd val: " <<pln<<endl;
    }
    std::cout << "Finished planarity calculation" << std::endl;

    return;
  }

  std::unordered_set<int> outliner(std::unordered_set<int> &start, std::vector<sElm> &All_Elms, int cur_ID)
  {
    // int start_size = 0;
    // //count number of neighbors of each element in set and sum
    // for (nbd_iter i = start.begin(); i != start.end(); i++)
    // {
    //   start_size += All_Elms[*i].nbs.size();
    // }
    std::unordered_set<int> ret;
    // ret.rehash(start_size); //start with the max size, save rehash time

    //iterate over neighboring elements
    for (nbd_iter i = start.begin(); i != start.end(); i++)
    { 
      //grab all of current element's neighbors and iterate
      std::unordered_set<int> work = All_Elms[*i].nbs;
      for (nbd_iter j = work.begin(); j != work.end(); j++)
      { 
        //search current set of neighbors and exclude if contained
        if (start.find(*j) == start.end())
        { 
          //check that element is unnumbered
          if (All_Elms[*j].SurfID == 0)
          { 
            //insert element to numbering set
            ret.insert(*j);
          }
        }
      }
    }
    // a set full of the unnumbered neighbours
    return ret;
  }

  void SurfTracer(std::unordered_set<int> &cur_set, std::vector<sElm> &All_Elms, int cur_ID)
  {
    //loop neighbors and set surface ID to current ID
    for (nbd_iter i = cur_set.begin(); i != cur_set.end(); i++)
    {
      All_Elms[*i].SurfID = cur_ID; // set the ID
    }
    //fill set with next group of elements to add to current shock system 
    std::unordered_set<int> next_set = outliner(cur_set, All_Elms, cur_ID);
    if (next_set.size() > 0)
    { 
      // recurse if set is not empty
      SurfTracer(next_set, All_Elms, cur_ID);
    }
    return;
  }

  std::vector<int> SurfaceSorter(std::vector<sElm> &All_Elms)
  {
    int SurfID_iter = 1;
    //loop over all elements
    for (int i = 0; i < All_Elms.size(); i++)
    {
      //if non surface ID has been assigned
      if (All_Elms[i].SurfID == 0)
      {
        //create set for all elements in shock system
        std::unordered_set<int> begin;
        begin.insert(i); // initialize the function w/ a size 1 set

        SurfTracer(begin, All_Elms, SurfID_iter);
        SurfID_iter++;
      }
    }

    //vector of number of elements in each shock system
    std::vector<int> ret;
    ret.resize(SurfID_iter);
    for (int i = 0; i < All_Elms.size(); i++)
    {
      ret[All_Elms[i].SurfID]++;
    }

    return ret;
  }

  void WriteData(std::vector<sElm> &Elms){
    std::ofstream out_str("processed_shock_elements.txt");
    if (!out_str.good())
    {
      cerr << "Can't open " << "processed_shock_elements.txt" << " to write." << endl;
      exit(1);
    }

    out_str << setprecision(std::numeric_limits<double>::digits10 + 1);
    // out_str << "\"vtkOriginalPointIds\",\"Points:0\",\"Points:1\",\"Points:2\",\"Surface_ID\",\"Planarity_Ind\" " << std::endl;

    for (int i = 0; i < Elms.size(); i++)
    {
      out_str << i << ",";
      out_str << Elms[i].x_pos << ",";
      out_str << Elms[i].y_pos << ",";
      out_str << Elms[i].z_pos << ",";
      out_str << Elms[i].SurfID << ",";
      out_str << Elms[i].svd_val << ",";
      out_str << "dummy" << std::endl;
    }

    return;
  }
}
