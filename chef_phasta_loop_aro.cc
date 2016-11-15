#include "chefPhasta.h"
#include "samSz.h"
#include <PCU.h>
#include <chef.h>
#include <spr.h>
#include <phasta.h>
#include "phIO.h"
#include <phstream.h>
#include <sam.h>
#include <apfMDS.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

#include <apfSIM.h>
#include <gmi_sim.h>
#include <SimPartitionedMesh.h>
#include <MeshSimAdapt.h>
#include <SimField.h>
#include <cstring>

#ifndef WRITE_VTK
#define WRITE_VTK
#endif

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }

  static apf::Field* getSprSF(apf::Mesh2* m) {
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

  static FILE* openstream_read(ph::Input& in, const char* path) {
    std::string fname(path);
    std::string restartStr("restart");
    FILE* f = NULL;
    if( fname.find(restartStr) != std::string::npos )
      f = openRStreamRead(in.rs);
    else {
      fprintf(stderr,
        "ERROR %s type of stream %s is unknown... exiting\n",
        __func__, fname.c_str());
      exit(1);
    }
    return f;
  }

  void setupChef(ph::Input& ctrl, int step) {
    //don't split or tetrahedronize
    ctrl.splitFactor = 1;
    ctrl.tetrahedronize = 0;
    ctrl.timeStepNumber = step;
    ctrl.solutionMigration = 1;
    if(step>0) {
      if(!PCU_Comm_Self()) {
        fprintf(stderr, "STATUS error based adapt %d\n", step);
        fprintf(stderr, "STATUS ctrl.attributeFileName %s step %d\n",
            ctrl.attributeFileName.c_str(), step);
      }
      ctrl.adaptStrategy = 1; //error field adapt
      ctrl.adaptFlag = 1;
    }
  }
  
  bool overwriteMeshCoord(apf::Mesh2* m) { 
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

  bool isMeshqGood(apf::Mesh* m, double crtn) { 
    apf::Field* meshq = m->findField("meshQ");
    if (!meshq) {
      fprintf(stderr, "Not find meshQ field.");
      return true;  
    }
    apf::MeshEntity* elm; 
    apf::MeshIterator* itr = m->begin(m->getDimension());
    while( (elm = m->iterate(itr)) ) {
      if (apf::getScalar(meshq, elm, 0) < crtn) {
        apf::destroyField(meshq);
        return false; 
      } 
    }
    m->end(itr);
    apf::destroyField(meshq);
    return true; 
  }
  
  void writeSequence (apf::Mesh2* m, int step, const char* filename) {
    std::ostringstream oss; 
    oss << filename << step;
    const std::string tmp = oss.str();
#ifdef WRITE_VTK
    apf::writeVtkFiles(tmp.c_str(),m);
#endif
  }

  void writePHTfiles (int step, int nstep, int nproc) {
    std::ostringstream oss;
    oss << "solution_" << step << ".pht";
    const std::string tp = oss.str();
    const char* filename = tp.c_str();
    FILE* sFile = fopen (filename, "w");
    fprintf (sFile, "<?xml version=\"1.0\" ?>\n");
    fprintf (sFile, "<PhastaMetaFile number_of_pieces=\"%d\">\n", nproc);
    fprintf (sFile, "  <GeometryFileNamePattern pattern=\"%d/%d-procs_case/geombc.%d.%%d\"\n",step,nproc,step);
    fprintf (sFile, "                           has_piece_entry=\"1\"\n");
    fprintf (sFile, "                           has_time_entry=\"0\"/>\n");
    fprintf (sFile, "  <FieldFileNamePattern pattern=\"%d-procs_case/restart.%%d.%%d\"\n",nproc);
    fprintf (sFile, "                        has_piece_entry=\"1\"\n");
    fprintf (sFile, "                        has_time_entry=\"1\"/>\n");
    fprintf (sFile, "  <TimeSteps number_of_steps=\"%d\"\n", nstep);
    fprintf (sFile, "             auto_generate_indices=\"1\"\n");
    fprintf (sFile, "             start_index=\"%d\"\n", step+1);
    fprintf (sFile, "             increment_index_by=\"1\"\n");
    fprintf (sFile, "             start_value=\"0.0\"\n");
    fprintf (sFile, "             increment_value_by=\"1.0e-6\">\n");
    fprintf (sFile, "  </TimeSteps>\n");
    fprintf (sFile, "  <Fields number_of_fields=\"7\">\n");
    fprintf (sFile, "    <Field paraview_field_tag=\"pressure\"\n");
    fprintf (sFile, "           phasta_field_tag=\"solution\"\n");
    fprintf (sFile, "           start_index_in_phasta_array=\"0\"\n");
    fprintf (sFile, "           number_of_components=\"1\"/>\n");
    fprintf (sFile, "    <Field paraview_field_tag=\"velocity\"\n");
    fprintf (sFile, "           phasta_field_tag=\"solution\"\n");
    fprintf (sFile, "           start_index_in_phasta_array=\"1\"\n");
    fprintf (sFile, "           number_of_components=\"3\"\n");
    fprintf (sFile, "           data_dependency=\"0\"\n");
    fprintf (sFile, "           data_type=\"double\"/>\n");
    fprintf (sFile, "    <Field paraview_field_tag=\"temperature\"\n");
    fprintf (sFile, "           phasta_field_tag=\"solution\"\n");
    fprintf (sFile, "           start_index_in_phasta_array=\"4\"\n");
    fprintf (sFile, "           number_of_components=\"1\"/>\n");
    fprintf (sFile, "    <Field paraview_field_tag=\"motion_coords\"\n");
    fprintf (sFile, "           phasta_field_tag=\"motion_coords\"\n");
    fprintf (sFile, "           start_index_in_phasta_array=\"0\"\n");
    fprintf (sFile, "           number_of_components=\"3\"\n");
    fprintf (sFile, "           data_dependency=\"0\"\n");
    fprintf (sFile, "           data_type=\"double\"/>\n");
    fprintf (sFile, "    <Field paraview_field_tag=\"mesh_vel\"\n");
    fprintf (sFile, "           phasta_field_tag=\"mesh_vel\"\n");
    fprintf (sFile, "           start_index_in_phasta_array=\"0\"\n");
    fprintf (sFile, "           number_of_components=\"3\"\n");
    fprintf (sFile, "           data_dependency=\"0\"\n");
    fprintf (sFile, "           data_type=\"double\"/>\n");
    fprintf (sFile, "    <Field paraview_field_tag=\"material_type\"\n");
    fprintf (sFile, "           phasta_field_tag=\"material_type\"\n");           
    fprintf (sFile, "           start_index_in_phasta_array=\"0\"\n");           
    fprintf (sFile, "           number_of_components=\"1\"\n");           
    fprintf (sFile, "           data_dependency=\"1\"/>\n"); 
    fprintf (sFile, "    <Field paraview_field_tag=\"meshQ\"\n");
    fprintf (sFile, "           phasta_field_tag=\"meshQ\"\n");          
    fprintf (sFile, "           start_index_in_phasta_array=\"0\"\n");          
    fprintf (sFile, "           number_of_components=\"1\"\n");          
    fprintf (sFile, "           data_dependency=\"1\"/>\n");            
    fprintf (sFile, "  </Fields>\n");
    fprintf (sFile, "</PhastaMetaFile>\n");
    fclose (sFile);
  } 

  void runMeshAdapter(ph::Input& in, apf::Mesh2*& m, apf::Field*& szFld) {
    if (m->findField("material_type"))
      apf::destroyField(m->findField("material_type"));
    if (m->findField("meshQ"))
      apf::destroyField(m->findField("meshQ"));

    /* Or obtain size field based on a certain field
       use temperature field for spr error estimation */
//      apf::Field* szFld = getSprSF(m);
 
    if(in.simmetrixMesh == 1) {
      apf::MeshSIM* sim_m = dynamic_cast<apf::MeshSIM*>(m);
      pParMesh sim_pm = sim_m->getMesh();

      /* create the Simmetrix adapter */
      pMSAdapt adapter = MSA_new(sim_pm, 1);

      /* copy the size field from APF to the Simmetrix adapter */
      apf::MeshEntity* v;
      apf::MeshIterator* it = m->begin(0);
      while ((v = m->iterate(it))) {
        double size = apf::getScalar(szFld, v, 0);
        MSA_setVertexSize(adapter, (pVertex) v, size);
      }
      m->end(it);
      apf::destroyField(szFld);

      /* unpacked solution into serveral fields */
      apf::Field* pre_field = chef::extractField(m,"solution","pressure",1,apf::SCALAR,in.simmetrixMesh);
      apf::Field* vel_field = chef::extractField(m,"solution","velocity",2,apf::VECTOR,in.simmetrixMesh);
      apf::Field* tem_field = chef::extractField(m,"solution","temperature",5,apf::SCALAR,in.simmetrixMesh);

      /* put these field explicitly into pPList */
      int num_flds = 3; // for now 
      pField* sim_flds = new pField[num_flds];
      sim_flds[0] = apf::getSIMField(pre_field);
      sim_flds[1] = apf::getSIMField(vel_field);
      sim_flds[2] = apf::getSIMField(tem_field);
      pPList sim_fld_lst = PList_new();
      for (int i = 0; i < num_flds; i++)
        PList_append(sim_fld_lst, sim_flds[i]);
      assert(num_flds == PList_size(sim_fld_lst));

      /* set fields to be mapped */
      MSA_setMapFields(adapter, sim_fld_lst);
      PList_delete(sim_fld_lst);

      /* run the adapter */
      pProgress progress = Progress_new();
      MSA_adapt(adapter, progress);
      Progress_delete(progress);
      MSA_delete(adapter);

      /* transfer data back to apf */
      apf::Field* solution = chef::combineField(m,"solution","pressure","velocity","temperature");
    }
    else {
      assert(szFld);
      apf::synchronize(szFld);
      apf::synchronize(m->getCoordinateField());
      /* do SCOREC mesh adaptation */
      chef::adapt(m,szFld,in);
    }
  }

} //end namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_Protect();
//debugging for simmetrix mesh
  Sim_readLicenseFile(0);
  SimPartitionedMesh_start(0, 0);
  Sim_logOn("loopDriver.log");
  SimUtil_start();
  SimModel_start();
  SimField_start();
  gmi_sim_start();
  gmi_register_sim();
  MS_init();
//end debugging
  if( argc != 2 ) {
    if(!PCU_Comm_Self())
      fprintf(stderr, "Usage: %s <maxTimeStep>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  int maxStep = atoi(argv[1]);
  chefPhasta::initModelers();
  rstream rs = makeRStream();
  grstream grs = makeGRStream();
  ph::Input ctrl;
  ctrl.load("samAdaptLoop.inp");
  /* load the model and mesh */
  gmi_register_mesh();
  gmi_model* g = 0;
  apf::Mesh2* m = 0;
  chef::cook(g, m, ctrl, rs, grs);
  /* load input file for solver */
  phSolver::Input inp("solver.inp", "input.config");
  int step = 0; int phtStep = 0; int seq  = 0;
  writeSequence(m,seq,"test_"); seq++;
  do {
    m->verify();
    /* take the initial mesh as size field */
    apf::Field* isoSF = samSz::isoSize(m);
    apf::Field* szFld = multipleSF(m, isoSF, 2.0);
//    apf::Field* szFld = multipleSF(m, isoSF, 0.5);
    step = phasta(inp,grs,rs);
    ctrl.rs = rs; 
    clearGRStream(grs);
    if(!PCU_Comm_Self())
      fprintf(stderr, "STATUS ran to step %d\n", step);
    if( step >= maxStep )
      break;
    setupChef(ctrl,step);
//    bool doAdaptation = !isMeshqGood(m, ctrl.meshqCrtn);
// make the adaptaion run anyway
//    doAdaptation = false;
    bool doAdaptation = true;
// delele above when finish debug
    m->verify();
    writePHTfiles(phtStep, step-phtStep, PCU_Comm_Peers()); phtStep = step; 
    if ( doAdaptation ) {
      chef::readAndAttachFields(ctrl,m);
      overwriteMeshCoord(m);
      m->verify();
      writeSequence(m,seq,"test_"); seq++;
      /* do mesh adaptation */ 
      runMeshAdapter(ctrl,m,szFld);
      writeSequence(m,seq,"test_"); seq++;
      m->verify(); 
      chef::balance(ctrl,m);
    }
    chef::preprocess(m,ctrl,grs);
    if ( doAdaptation )
      writeSequence(m,seq,"test_"); seq++;
    clearRStream(rs);
  } while( step < maxStep );
  destroyGRStream(grs);
  destroyRStream(rs);
  freeMesh(m);
  chefPhasta::finalizeModelers();
//debugging for simmetrix mesh
  gmi_sim_stop();
  SimField_stop();
  SimModel_stop();
  SimUtil_stop();
  Sim_logOff();
  SimPartitionedMesh_stop();
  Sim_unregisterAllKeys();
//end debugging
  PCU_Comm_Free();
  MPI_Finalize();
}
