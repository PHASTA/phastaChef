#include "pcWriteFiles.h"
#include "SimModel.h"
#include "apfSIM.h"
#include "apfMDS.h"
#include <sstream>
#include <cstdio>
#include <cassert>

namespace pc {
  void writeSequence (apf::Mesh2* m, int step, const char* filename) {
    std::ostringstream oss;
    oss << filename << step;
    const std::string tmp = oss.str();
    apf::writeVtkFiles(tmp.c_str(),m);
  }

  void writePHTfiles (int step, int nstep, int nproc) {
    std::ostringstream oss;
    oss << "solution_" << step << ".pht";
    const std::string tp = oss.str();
    const char* filename = tp.c_str();
    FILE* sFile = fopen (filename, "w");
    fprintf (sFile, "<?xml version=\"1.0\" ?>\n");
    fprintf (sFile, "<PhastaMetaFile number_of_pieces=\"%d\">\n", nproc);
    fprintf (sFile, "  <GeometryFileNamePattern pattern=\"%d-procs_case/geombc.%d.%%d\"\n",nproc,step);
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

  void writeStats(ph::Input& in, gmi_model* g, apf::Mesh2* m, int step) {
    // create anisotropic size field
    apf::Field* sizefld  = apf::createFieldOn(m, "size",  apf::VECTOR);
    apf::Field* framefld = apf::createFieldOn(m, "frame", apf::MATRIX);

    // declaration
    double isoSize;
    double size[1];
    double anisosize[3][3];

    // initial anisosize
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        anisosize[i][j] = 0.0;
      }
    }

    // query the size field and convert to size and frame
    apf::MeshEntity* v;
    apf::MeshIterator* vit = m->begin(0);
    while ((v = m->iterate(vit))) {
      if (in.simmetrixMesh == 1) {
        pVertex meshVertex = reinterpret_cast<pVertex>(v);
        int sztype = V_size(meshVertex, size, anisosize);

        // consider iso size field as anisotropic size field
        if (sztype == 1) {
          anisosize[0][0] = size[0];
          anisosize[1][1] = size[0];
          anisosize[2][2] = size[0];
        }
      }
      else {
        // for now, we only use isotropic size field
        assert(m->findField("isoSize"));
        apf::Field* apfSz = m->findField("isoSize");
        size[0] = apf::getScalar(apfSz,v,0);
        anisosize[0][0] = size[0];
        anisosize[1][1] = size[0];
        anisosize[2][2] = size[0];
      }

      // transfer to apf field sf_mag and sf_dir
      /* note that the frame in Simmetrix is stored by row
         while in PUMI, it is stored by column */
      apf::Vector3 v_mag;
      v_mag[0] = sqrt(anisosize[0][0]*anisosize[0][0]
                    + anisosize[0][1]*anisosize[0][1]
                    + anisosize[0][2]*anisosize[0][2]);
      v_mag[1] = sqrt(anisosize[1][0]*anisosize[1][0]
                    + anisosize[1][1]*anisosize[1][1]
                    + anisosize[1][2]*anisosize[1][2]);
      v_mag[2] = sqrt(anisosize[2][0]*anisosize[2][0]
                    + anisosize[2][1]*anisosize[2][1]
                    + anisosize[2][2]*anisosize[2][2]);
      apf::Matrix3x3 v_dir;
      v_dir = apf::Matrix3x3(
               anisosize[0][0]/v_mag[0], anisosize[1][0]/v_mag[1], anisosize[2][0]/v_mag[2],
               anisosize[0][1]/v_mag[0], anisosize[1][1]/v_mag[1], anisosize[2][1]/v_mag[2],
               anisosize[0][2]/v_mag[0], anisosize[1][2]/v_mag[1], anisosize[2][2]/v_mag[2]);
      apf::setVector(sizefld,  v, 0, v_mag);
      apf::setMatrix(framefld, v, 0, v_dir);
    }
    m->end(vit);

    // build filename for smb file
    std::ostringstream oss;
    oss << "for_stats_" << step << "_.smb";
//    oss << "for_stats_" << step << ".sms";
    const std::string tmp = oss.str();

    apf::Mesh2* smb_mesh = apf::createMdsMesh(g, m);
//    printf("verify apf mesh\n");
//    m->verify();
//    printf("verify MDS mesh\n");
//    smb_mesh->verify();
    smb_mesh->writeNative(tmp.c_str());
    smb_mesh->destroyNative();
    apf::destroyMesh(smb_mesh);

/* we cannot write sms mesh right now.
   When we do convert,
   there is no field on the new mesh */
//    m->writeNative(tmp.c_str());

    // delete size and frame fields
    apf::destroyField(sizefld);
    apf::destroyField(framefld);
  }

}
