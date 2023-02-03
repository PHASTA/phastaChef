#Setting Version Number, Project Name
cmake_minimum_required (VERSION 2.8)
project(PHASTACHEF Fortran C CXX)
enable_language(Fortran)
enable_language(CXX)
enable_language(C)

execute_process(
  COMMAND git rev-parse HEAD
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  OUTPUT_FILE "${CMAKE_BINARY_DIR}/VERSION.phastaChef"
)

# Set options for doxygen documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY
)
add_custom_target(doc
  ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif()

set(PHASTA_CHEF_ENABLED ON)

enable_testing()
include(CTest)
set(MPIRUN "mpirun"
  CACHE string
  "the mpirun or srun executable")
set(MPIRUN_PROCFLAG "-np"
  CACHE string
  "the command line flag to give process count to MPIRUN")
set(PHASTA_SRC_DIR phasta CACHE FILEPATH "path to phasta source code")
add_subdirectory(${PHASTA_SRC_DIR} ${CMAKE_BINARY_DIR}/phasta)

#test to see if simmetrix models are supported
find_package(SCOREC 2.1.0 REQUIRED CONFIG PATHS ${SCOREC_PREFIX} NO_DEFAULT_PATH)
if(TARGET SCOREC::gmi_sim)
  set(GMI_SIM_FOUND 1)
  message(STATUS "Found gmi_sim package - enabling Simmetrix modeler support")
  add_definitions(-DGMI_SIM_FOUND)
endif()

find_package(phastaIC PATHS ${CMAKE_BINARY_DIR})
message(STATUS "PHASTAIC_LIBS ${PHASTAIC_LIBS}")
find_package(phastaC PATHS ${CMAKE_BINARY_DIR})
message(STATUS "PHASTAC_LIBS ${PHASTAC_LIBS}")

if( NOT phastaIC_FOUND AND NOT phastaC_FOUND )
 message(FATAL_ERROR "Neither the incompressible or compressible solver was found")
endif()

find_library(ACUSOLVE_LIB libles)

macro(setup_exe exename srcname IC)
  set(src ${srcname}
    pcWriteFiles.cc
    pcUpdateMesh.cc
    pcAdapter.cc
    pcTimeDepMesh.cc
    pcSmooth.cc
    pcError.cc
  )

  add_executable(${exename} ${src})
  set_target_properties(${exename} PROPERTIES HAS_CXX TRUE)
  set_target_properties(${exename} PROPERTIES HAS_CXX TRUE)
  set_target_properties(${exename} PROPERTIES LINKER_LANGUAGE Fortran)
  set_target_properties(${exename} PROPERTIES LINKER_LANGUAGE Fortran)
  #chef
  target_link_libraries(${exename} PRIVATE SCOREC::core)
  #phasta
  if( ${IC} )
    include_directories(${PHASTAIC_INCLUDE_DIRS})
    target_link_libraries(${exename} PRIVATE ${PHASTAIC_LIBS})
  else()
    include_directories(${PHASTAC_INCLUDE_DIRS})
    target_link_libraries(${exename} PRIVATE ${PHASTAC_LIBS})
  endif()
endmacro(setup_exe)

setup_exe(chefPhasta_posix chef_phasta_posix.cc ${phastaIC_FOUND})
setup_exe(chefPhasta_stream chef_phasta_stream.cc ${phastaIC_FOUND})
setup_exe(chefPhastaLoop_stream_ur chef_phasta_loop_stream_ur.cc ${phastaIC_FOUND})
setup_exe(chefPhastaLoop_stream_adapt chef_phasta_adaptLoop.cc ${phastaIC_FOUND})
setup_exe(chefPhastaLoop_files_adapt chef_phasta_adaptLoop_files.cc ${phastaIC_FOUND})
setup_exe(chefPhastaLoop_sam_stream_adapt chef_phasta_sam_adaptLoop.cc ${phastaIC_FOUND})
setup_exe(loopChefPhasta loopChefPhasta.cc ${phastaIC_FOUND})
setup_exe(transferAndAdapter transferAndAdapter.cc ${phastaIC_FOUND})
setup_exe(solutionProjection solutionProjection.cc ${phastaIC_FOUND})
setup_exe(calcEfficiency calcEfficiency.cc ${phastaIC_FOUND})
setup_exe(meshGrading meshGrading.cc ${phastaIC_FOUND})

setup_exe(_posix chef_phasta_posix.cc)
setup_exe(_stream chef_phasta_stream.cc)
setup_exe(Loop_ur chef_phasta_loop_stream_ur.cc)
setup_exe(Loop_adapt chef_phasta_adaptLoop.cc)
setup_exe(Loop_adapt_posix chef_phasta_adaptLoop_posix.cc)
setup_exe(Loop_adapt_ramdisk chef_phasta_adaptLoop_ramdisk.cc)
setup_exe(Loop_sam_adapt chef_phasta_sam_adaptLoop.cc)
setup_exe(Loop_sam_adapt_sz chef_phasta_sam_adaptLoop_sz.cc)
setup_exe(Loop_sam_adapt_sz_posix chef_phasta_sam_adaptLoop_sz_posix.cc)
setup_exe(Loop_sam_adapt_sz_ramdisk chef_phasta_sam_adaptLoop_sz_ramdisk.cc)
setup_exe(Loop_sam_adapt_ramdisk chef_phasta_sam_adaptLoop_ramdisk.cc)
create_exe(chefPhastaICLoop_aro chef_phasta_loop_aro.cc
  "${PHASTAIC_INCLUDE_DIRS}" "${PHASTAIC_LIBS}")
add_subdirectory(test)
