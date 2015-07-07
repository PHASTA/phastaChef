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

add_executable(chef_phasta chef_phasta.cc)
set_target_properties(phastaIC.exe PROPERTIES HAS_CXX TRUE)
set_target_properties(incompressible PROPERTIES HAS_CXX TRUE)
set_target_properties(phastaIC.exe PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(incompressible PROPERTIES LINKER_LANGUAGE Fortran)

#chef
find_package(ph REQUIRED PATHS ${CMAKE_BINARY_DIR})
include_directories(${PH_INCLUDE_DIRS})
target_link_libraries(chef_phasta ${PH_LIBS})

#phSolver
find_library(ACUSOLVE_LIB libles)
add_dependencies(incompressible common)
include_directories(${PHASTA_SOURCE_DIR}/phSolver/common)
if(PHASTA_USE_MPI)
  target_link_libraries(chef_phasta 
    ${MPI_LIBRARIES} ${MPI_Fortran_LIBRARIES} 
    ${CMAKE_THREAD_LIBS_INIT}) 
endif(PHASTA_USE_MPI)
target_link_libraries(chef_phasta ${ACUSOLVE_LIB})
target_link_libraries(chef_phasta common incompressible common phastaIO)
if(PHASTA_BUILD_SHAPEFUNCTION)
  target_link_libraries(chef_phasta shapeFunction)
endif(PHASTA_BUILD_SHAPEFUNCTION)
if(PHASTA_BUILD_PHSHAPE)
  target_link_libraries(chef_phasta phshape)
endif(PHASTA_BUILD_PHSHAPE)
target_link_libraries(chef_phasta ${ACUSOLVE_LIB} incompressible)

add_subdirectory(test)
