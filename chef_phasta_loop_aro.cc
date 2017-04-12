#include "chefPhasta.h"
#include "sam.h"
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
      ctrl.writeGeomBCFiles = 0;
    }
  }

  bool overwriteAPFCoord(apf::Mesh2* m) {
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

  bool overwriteSIMCoord(apf::Mesh2* m) {
    apf::MeshSIM* apf_msim = dynamic_cast<apf::MeshSIM*>(m);
    pParMesh ppm = apf_msim->getMesh();
    pMesh pm = PM_mesh(ppm,0);

    M_write(pm, "before_overwrite.sms", 0, NULL);

    VIter vi = M_vertexIter(pm);
    const char* filename = "temp_sim_coord";
    double* loc = new double[3];
    FILE* fp = fopen (filename, "w");
    while (pVertex v = VIter_next(vi)) {
      V_coord(v, loc);
      fprintf(fp, "%.16lg %.16lg %.16lg\n", loc[0], loc[1], loc[2]);
    }
    VIter_delete(vi);
    fclose (fp);

    printf("using MeshMover_run to update the simmetrix mesh coords\n");
    assert(0);
    return true;
  }

  void overwriteMeshCoord(ph::Input& in, apf::Mesh2* m) {
    bool done = false;
    if (in.simmetrixMesh)
      done = overwriteSIMCoord(m);
    else
      done = overwriteAPFCoord(m);
    assert(done);
  }

  int isMeshqGood(apf::Mesh* m, double crtn) {
    apf::Field* meshq = m->findField("meshQ");
    if (!meshq) {
      if (!PCU_Comm_Self())
        fprintf(stderr, "Not find meshQ field.\n");
      assert(meshq);
    }
    int meshGood = 1;
    apf::MeshEntity* elm;
    apf::MeshIterator* itr = m->begin(m->getDimension());
    while( (elm = m->iterate(itr)) ) {
      if (apf::getScalar(meshq, elm, 0) < crtn) {
        meshGood = 0;
        break;
      }
    }
    m->end(itr);
    PCU_Barrier();
    if (PCU_Min_Int(meshGood) && !PCU_Comm_Self())
      printf("Mesh is Good. No need for adaptation!\n");
    return PCU_Min_Int(meshGood);
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

  void runMeshAdapter(ph::Input& in, apf::Mesh2*& m, apf::Field*& orgSF, int step) {
    if (m->findField("material_type"))
      apf::destroyField(m->findField("material_type"));
    if (m->findField("meshQ"))
      apf::destroyField(m->findField("meshQ"));
    in.writeGeomBCFiles = 1; //write GeomBC file for visualization

    /* use the size field of the mesh before mesh motion */
    apf::Field* szFld = orgSF;

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

      PM_write(sim_pm, "before_adapt.sms", NULL);

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
    m->verify();
    chef::balance(in,m);
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
    apf::Field* szFld = samSz::isoSize(m);
    step = phasta(inp,grs,rs);
    ctrl.rs = rs;
    clearGRStream(grs);
    if(!PCU_Comm_Self())
      fprintf(stderr, "STATUS ran to step %d\n", step);
    if( step >= maxStep )
      break;
    setupChef(ctrl,step);
    chef::readAndAttachFields(ctrl,m);
    overwriteMeshCoord(ctrl,m);
//    int doAdaptation = !isMeshqGood(m, ctrl.meshqCrtn);
    int doAdaptation = 1; // make it run anyway

    m->verify();
    if ( doAdaptation ) {
      writePHTfiles(phtStep, step-phtStep, PCU_Comm_Peers()); phtStep = step;
      writeSequence(m,seq,"test_"); seq++;
      /* do mesh adaptation */
      runMeshAdapter(ctrl,m,szFld,step);
      writeSequence(m,seq,"test_"); seq++;
    }
    chef::preprocess(m,ctrl,grs);
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
