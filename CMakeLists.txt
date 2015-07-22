#Setting Version Number, Project Name
cmake_minimum_required (VERSION 2.8)
project(PHASTACHEF Fortran C CXX)
enable_language(Fortran)
enable_language(CXX)
enable_language(C)

enable_testing()
include(CTest)
set(MPIRUN "mpirun"
  CACHE string
  "the mpirun or srun executable")
set(MPIRUN_PROCFLAG "-np"
  CACHE string
  "the command line flag to give process count to MPIRUN")

add_subdirectory(core)
add_subdirectory(phasta)

find_package(ph REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_library(ACUSOLVE_LIB libles)

macro(setup_exe exename srcname)
  add_executable(${exename} ${srcname})
  target_link_libraries(${exename} common phastaIO)
  set_target_properties(${exename} PROPERTIES HAS_CXX TRUE)
  set_target_properties(${exename} PROPERTIES HAS_CXX TRUE)
  set_target_properties(${exename} PROPERTIES LINKER_LANGUAGE Fortran)
  set_target_properties(${exename} PROPERTIES LINKER_LANGUAGE Fortran)

  #chef
  include_directories(${PH_INCLUDE_DIRS})
  target_link_libraries(${exename} ${PH_LIBS})

  #phasta
  add_dependencies(incompressible common)
  include_directories(${PHASTA_SOURCE_DIR}/phSolver/common)
  if(PHASTA_USE_MPI)
    target_link_libraries(${exename}
      ${MPI_LIBRARIES} ${MPI_Fortran_LIBRARIES}
      ${CMAKE_THREAD_LIBS_INIT})
  endif(PHASTA_USE_MPI)
  target_link_libraries(${exename} ${ACUSOLVE_LIB})
  target_link_libraries(${exename} common incompressible common phastaIO)
  if(PHASTA_BUILD_SHAPEFUNCTION)
    target_link_libraries(${exename} shapeFunction)
  endif(PHASTA_BUILD_SHAPEFUNCTION)
  if(PHASTA_BUILD_PHSHAPE)
    target_link_libraries(${exename} phshape)
  endif(PHASTA_BUILD_PHSHAPE)
  target_link_libraries(${exename} ${ACUSOLVE_LIB} incompressible)
endmacro(setup_exe)

setup_exe(chef_phasta chef_phasta.cc)

add_subdirectory(test)
