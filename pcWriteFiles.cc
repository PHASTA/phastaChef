#include "pcWriteFiles.h"
#include "SimModel.h"
#include "apfSIM.h"
#include "apfMDS.h"
#include <PCU.h>
#include <sstream>
#include <cstdio>
#include <cassert>
#include <algorithm>
#include <math.h>
#include <SimPartitionedMesh.h>
#include "SimModel.h"
#include "SimUtil.h"
#include "SimParasolidKrnl.h"
#include "SimMeshTools.h"

namespace pc {
  void writeSequence (apf::Mesh2* m, int step, const char* filename) {
    std::ostringstream oss;
    oss << filename << step;
    const std::string tmp = oss.str();
    apf::writeVtkFiles(tmp.c_str(),m);
  }

  void writeSIMModel (pGModel model, int step, const char* filename) {
    std::ostringstream oss;
    oss << filename << step << ".smd";
    const std::string tmp = oss.str();
    GM_write(model,tmp.c_str(),0,NULL);
  }

  void writeSIMMesh (pParMesh mesh, int step, const char* filename) {
    std::ostringstream oss;
    oss << filename << step << ".sms";
    const std::string tmp = oss.str();
    PM_write(mesh,tmp.c_str(),NULL);
  }

  void writePHTfiles (int old_step, int step, phSolver::Input& inp) {
    int nfields = 7;
    int ntout = min((int)inp.GetValue("Number of Timesteps between Restarts"),
                    (int)inp.GetValue("Number of Timesteps"));
    double dt = (double)inp.GetValue("Time Step Size");
    if((string)inp.GetValue("Write non-linear residual to restart") == "Yes")
      nfields = nfields + 3;
    try {
      if((string)inp.GetValue("Error Estimation Option") != "False")
        nfields = nfields + 3;
    }
    catch(...){}
    int nproc = PCU_Comm_Peers();
    int nstep = max((step - old_step) / ntout, 1);
    int start_step = ceil((double)old_step / (double)ntout) * ntout;
    std::ostringstream oss;
    oss << "solution_" << old_step << ".pht";
    const std::string tp = oss.str();
    const char* filename = tp.c_str();
    FILE* sFile = fopen (filename, "w");
    fprintf (sFile, "<?xml version=\"1.0\" ?>\n");
    fprintf (sFile, "<PhastaMetaFile number_of_pieces=\"%d\">\n", nproc);
    if (old_step == 0)
      fprintf (sFile, "  <GeometryFileNamePattern pattern=\"%d-procs_case/geombc.dat.%%d\"\n",nproc);
    else
      fprintf (sFile, "  <GeometryFileNamePattern pattern=\"%d/%d-procs_case/geombc.%d.%%d\"\n",old_step,nproc,old_step);
    fprintf (sFile, "                           has_piece_entry=\"1\"\n");
    fprintf (sFile, "                           has_time_entry=\"0\"/>\n");
    fprintf (sFile, "  <FieldFileNamePattern pattern=\"%d-procs_case/restart.%%d.%%d\"\n",nproc);
    fprintf (sFile, "                        has_piece_entry=\"1\"\n");
    fprintf (sFile, "                        has_time_entry=\"1\"/>\n");
    fprintf (sFile, "  <TimeSteps number_of_steps=\"%d\"\n", nstep);
    fprintf (sFile, "             auto_generate_indices=\"1\"\n");
    fprintf (sFile, "             start_index=\"%d\"\n", start_step+ntout);
    fprintf (sFile, "             increment_index_by=\"%d\"\n",ntout);
    fprintf (sFile, "             start_value=\"%12.16e\"\n",(double)((start_step+ntout)*dt));
    fprintf (sFile, "             increment_value_by=\"%12.16e\">\n",dt);
    fprintf (sFile, "  </TimeSteps>\n");
    fprintf (sFile, "  <Fields number_of_fields=\"%d\">\n",nfields);
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
    if((string)inp.GetValue("Write non-linear residual to restart") == "Yes") {
      fprintf (sFile, "    <Field paraview_field_tag=\"residual_mass\"\n");
      fprintf (sFile, "           phasta_field_tag=\"residual\"\n");
      fprintf (sFile, "           start_index_in_phasta_array=\"0\"\n");
      fprintf (sFile, "           number_of_components=\"1\"\n");
      fprintf (sFile, "           data_dependency=\"0\"\n");
      fprintf (sFile, "           data_type=\"double\"/>\n");
      fprintf (sFile, "    <Field paraview_field_tag=\"residual_momentum\"\n");
      fprintf (sFile, "           phasta_field_tag=\"residual\"\n");
      fprintf (sFile, "           start_index_in_phasta_array=\"1\"\n");
      fprintf (sFile, "           number_of_components=\"3\"\n");
      fprintf (sFile, "           data_dependency=\"0\"\n");
      fprintf (sFile, "           data_type=\"double\"/>\n");
      fprintf (sFile, "    <Field paraview_field_tag=\"residual_energy\"\n");
      fprintf (sFile, "           phasta_field_tag=\"residual\"\n");
      fprintf (sFile, "           start_index_in_phasta_array=\"4\"\n");
      fprintf (sFile, "           number_of_components=\"1\"\n");
      fprintf (sFile, "           data_dependency=\"0\"\n");
      fprintf (sFile, "           data_type=\"double\"/>\n");
    }
    try {
     if((string)inp.GetValue("Error Estimation Option") != "False") {
      fprintf (sFile, "    <Field paraview_field_tag=\"error_mass\"\n");
      fprintf (sFile, "           phasta_field_tag=\"VMS_error\"\n");
      fprintf (sFile, "           start_index_in_phasta_array=\"0\"\n");
      fprintf (sFile, "           number_of_components=\"1\"\n");
      fprintf (sFile, "           data_dependency=\"1\"\n");
      fprintf (sFile, "           data_type=\"double\"/>\n");
      fprintf (sFile, "    <Field paraview_field_tag=\"error_momt\"\n");
      fprintf (sFile, "           phasta_field_tag=\"VMS_error\"\n");
      fprintf (sFile, "           start_index_in_phasta_array=\"1\"\n");
      fprintf (sFile, "           number_of_components=\"1\"\n");
      fprintf (sFile, "           data_dependency=\"1\"\n");
      fprintf (sFile, "           data_type=\"double\"/>\n");
      fprintf (sFile, "    <Field paraview_field_tag=\"error_engy\"\n");
      fprintf (sFile, "           phasta_field_tag=\"VMS_error\"\n");
      fprintf (sFile, "           start_index_in_phasta_array=\"2\"\n");
      fprintf (sFile, "           number_of_components=\"1\"\n");
      fprintf (sFile, "           data_dependency=\"1\"\n");
      fprintf (sFile, "           data_type=\"double\"/>\n");
     }
    }
    catch(...){}
    fprintf (sFile, "  </Fields>\n");
    fprintf (sFile, "</PhastaMetaFile>\n");
    fclose (sFile);
  }

}
