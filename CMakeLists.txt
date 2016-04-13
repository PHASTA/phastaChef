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

set(PHASTA_CHEF_ENABLED ON)

enable_testing()
include(CTest)
set(MPIRUN "mpirun"
  CACHE string
  "the mpirun or srun executable")
set(MPIRUN_PROCFLAG "-np"
  CACHE string
  "the command line flag to give process count to MPIRUN")
add_subdirectory(phstream)
set(CORE_SRC_DIR core CACHE FILEPATH "path to pumi core source code")
add_subdirectory(${CORE_SRC_DIR} ${CMAKE_BINARY_DIR}/core)
set(PHASTA_SRC_DIR phasta CACHE FILEPATH "path to phasta source code")
add_subdirectory(${PHASTA_SRC_DIR} ${CMAKE_BINARY_DIR}/phasta)

#test to see if simmetrix models are supported
find_package(gmi_sim QUIET PATHS ${CMAKE_BINARY_DIR})
if( ${GMI_SIM_FOUND} )
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

find_package(ph REQUIRED PATHS ${CMAKE_BINARY_DIR}) #chef pkg
find_library(ACUSOLVE_LIB libles)

macro(setup_exe exename srcname IC)
  add_executable(${exename} ${srcname})
  set_target_properties(${exename} PROPERTIES HAS_CXX TRUE)
  set_target_properties(${exename} PROPERTIES HAS_CXX TRUE)
  set_target_properties(${exename} PROPERTIES LINKER_LANGUAGE Fortran)
  set_target_properties(${exename} PROPERTIES LINKER_LANGUAGE Fortran)

  #chef
  include_directories(${PH_INCLUDE_DIRS})
  target_link_libraries(${exename} ${PH_LIBS})

  #phasta
  if( ${IC} )
    include_directories(${PHASTAIC_INCLUDE_DIRS})
    target_link_libraries(${exename} ${PHASTAIC_LIBS})
  else()
    include_directories(${PHASTAC_INCLUDE_DIRS})
    target_link_libraries(${exename} ${PHASTAC_LIBS})
  endif()
endmacro(setup_exe)

setup_exe(chefPhasta_posix chef_phasta_posix.cc ${phastaIC_FOUND})
setup_exe(chefPhasta_stream chef_phasta_stream.cc ${phastaIC_FOUND})
setup_exe(chefPhastaLoop_stream_ur chef_phasta_loop_stream_ur.cc ${phastaIC_FOUND})
setup_exe(chefPhastaLoop_stream_adapt chef_phasta_adaptLoop.cc ${phastaIC_FOUND})
setup_exe(chefPhastaLoop_files_adapt chef_phasta_adaptLoop_files.cc ${phastaIC_FOUND})
setup_exe(chefPhastaLoop_sam_stream_adapt chef_phasta_sam_adaptLoop.cc ${phastaIC_FOUND})

add_subdirectory(test)
