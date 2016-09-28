#include "chefPhasta.h"
#include "samSz.h"
#include <PCU.h>
#include <chef.h>
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

#ifndef WRITE_VTK
#define WRITE_VTK
#endif

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
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

  apf::Field* getField(apf::Mesh* m) {
    /* if the value of the fldIdx'th index from the fldName
     * field is greater than fldLimit then multiply the current
     * isotropic mesh size at the vertex by szFactor */
    const unsigned fldIdx = 1;
    const double fldLimit = 1.0;
    const double szFactor = 0.5;
    const char* fldName = "motion_coords";
    return sam::errorThreshold(m,fldName,fldIdx,fldLimit,szFactor);
  }
  
  static FILE* openfile_read(ph::Input&, const char* path) {
    return fopen(path, "r");
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
    if(step>1) {
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

  void printMeshCoord (apf::Mesh2* m) { 
    apf::MeshEntity* vtx;
    apf::Vector3 points;
    apf::MeshIterator* itr = m->begin(0);
    int countNode = 0; 
    while( (vtx = m->iterate(itr)) ) {
      m->getPoint(vtx, 0, points);
      fprintf(stderr, "Coord of node %d: (%f, %f, %f)\n", 
              countNode, points[0], points[1], points[2]);
      countNode++;
    }
    m->end(itr);
  }

  void printMotionCoordField (apf::Mesh2* m) { 
    apf::Field* f = m->findField("motion_coords");
    assert(f);
    double* vals = new double[apf::countComponents(f)];
    assert(apf::countComponents(f) == 3);
    apf::MeshEntity* vtx;
    apf::Vector3 points; 
    apf::MeshIterator* itr = m->begin(0);
    int countNode = 0; 
    while( (vtx = m->iterate(itr)) ) {
      apf::getComponents(f, vtx, 0, vals);
      fprintf(stderr, "Motion_Coords of %d: (%f, %f, %f)\n", 
              countNode, vals[0], vals[1], vals[2]);
      countNode++;
    }
    m->end(itr); 
    delete [] vals;
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

  void attachSizeField (apf::Field* sf, apf::Mesh2* m) {
    int out_size = 1;
    apf::Field* f = apf::createPackedField(m, "Size Field", out_size);
    double* data = new double[out_size]; 
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(0);
    while ((e = m->iterate(it))) {
      apf::getComponents(sf,e, 0, data);
      apf::setComponents(f, e, 0, data);
    }
    m->end(it);
  }
/*
  bool readAndCheckMeshqField (ph::Input& in, apf::Mesh2* m) { 
    double* data;
    int nodes, vars, step;
    char hname[1024];
    const char* anyfield = "";
    setupInputSubdir(in.restartFileName);
    std::string filename = buildRestartFileName(in.restartFileName, in.timeStepNumber);
    FILE* f = in.openfile_read(in, filename.c_str());
    if (!f) {
      fprintf(stderr,"failed to open \"%s\"!\n", filename.c_str());
      abort();
    }
    int swap = ph_should_swap(f);
//    int ret = ph_read_field(f, anyfield, swap,
//        &data, &nodes, &vars, &step, "MeshQ");
    if(ret==0 || ret==1)
      return true;
    assert(step == in.timeStepNumber);

    fclose(f);
  }
*/
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

}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_Protect();
  if( argc != 2 ) {
    if(!PCU_Comm_Self())
      fprintf(stderr, "Usage: %s <maxTimeStep>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  int maxStep = atoi(argv[1]);
  chefPhasta::initModelers();
  grstream grs = makeGRStream();
  ph::Input ctrl;
  ctrl.load("samAdaptLoop.inp");
  /* setup file reading */
//  ctrl.openfile_read = openfile_read;
  /* load the model and mesh */
  apf::Mesh2* m = apf::loadMdsMesh(
      ctrl.modelFileName.c_str(),ctrl.meshFileName.c_str());
  chef::preprocess(m,ctrl);
  chef::preprocess(m,ctrl,grs);
  rstream rs = makeRStream();
  /* setup stream reading */
  ctrl.openfile_read = openstream_read;
  ctrl.rs = rs;
  phSolver::Input inp("solver.inp", "input.config");
  int step = 0; int phtStep = 0; 
  int loop = 0;
  int seq  = 0;
  writeSequence(m,seq,"test_"); seq++; 
  do {
    m->verify();
    /* take the initial mesh as size field */
    apf::Field* isoSF = samSz::isoSize(m);
    apf::Field* szFld = multipleSF(m, isoSF, 1.1);
    step = phasta(inp,grs,rs);
    ctrl.rs = rs; 
    clearGRStream(grs);
    if(!PCU_Comm_Self())
      fprintf(stderr, "STATUS ran to step %d\n", step);
    setupChef(ctrl,step);
    chef::readAndAttachFields(ctrl,m);
    m->verify();
    overwriteMeshCoord(m);
    m->verify();
    bool doAdaptation = !isMeshqGood(m, ctrl.meshqCrtn);
    m->verify();
// make the adaptaion run anyway
//    doAdaptation = false; 
    doAdaptation = true; 
// delele above when finish debug
    apf::destroyField(m->findField("material_type"));
    m->verify();
    if ( doAdaptation ) {
      writePHTfiles(phtStep, step-phtStep, PCU_Comm_Peers()); phtStep = step; 
      writeSequence(m,seq,"test_"); seq++; 
    }
    /* Or obtain size field based on a certain field*/
//    apf::Field* szFld = getField(m);
    apf::synchronize(szFld);
    apf::synchronize(m->getCoordinateField());
    assert(szFld);
    if ( doAdaptation ) {
      chef::adapt(m,szFld,ctrl);
      writeSequence(m,seq,"test_"); seq++;
      m->verify();
    } 
    apf::destroyField(szFld);
    chef::balanceAndReorder(ctrl,m);
    chef::preprocess(m,ctrl,grs);
    clearRStream(rs);
    loop++; 
  } while( loop < maxStep );
  destroyGRStream(grs);
  destroyRStream(rs);
  freeMesh(m);
  chefPhasta::finalizeModelers();
  PCU_Comm_Free();
  MPI_Finalize();
}
