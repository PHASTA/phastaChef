set(testLabel "chefphasta")
if( ${phastaIC_FOUND} )
  set(CDIR ${CASES}/incompressible)
  set(casename ${testLabel}_posix_incompressible)
  add_test(NAME ${casename}
    COMMAND ${CMAKE_COMMAND}
    -DNAME=${casename}
    -DWORKDIR=${CDIR}
    -DINPFILE=${PHASTA_SOURCE_DIR}/phSolver/common/input.config
    -DMPIRUN=${MPIRUN}
    -DMPIRUN_PROCFLAG=${MPIRUN_PROCFLAG}
    -DEXE=${PHASTACHEF_BINARY_DIR}/chefPhastaIC_posix
    -DNUMPROCS=4
    -DSRCDIRS=4-procs_case
    -P ${CMAKE_CURRENT_SOURCE_DIR}/runphasta.cmake
    )
  
  set(casename ${testLabel}_stream_incompressible)
  add_test(NAME ${casename}
    COMMAND ${CMAKE_COMMAND}
    -DNAME=${casename}
    -DWORKDIR=${CDIR}
    -DINPFILE=${PHASTA_SOURCE_DIR}/phSolver/common/input.config
    -DMPIRUN=${MPIRUN}
    -DMPIRUN_PROCFLAG=${MPIRUN_PROCFLAG}
    -DEXE=${PHASTACHEF_BINARY_DIR}/chefPhastaIC_stream
    -DNUMPROCS=4
    -DSRCDIRS=4-procs_case
    -P ${CMAKE_CURRENT_SOURCE_DIR}/runphasta.cmake
    )
  
  set(casename ${testLabel}_loopStreamUR_incompressible)
  add_test(NAME ${casename}
    COMMAND ${CMAKE_COMMAND}
    -DNAME=${casename}
    -DWORKDIR=${CDIR}
    -DINPFILE=${PHASTA_SOURCE_DIR}/phSolver/common/input.config
    -DMPIRUN=${MPIRUN}
    -DMPIRUN_PROCFLAG=${MPIRUN_PROCFLAG}
    -DEXE=${PHASTACHEF_BINARY_DIR}/chefPhastaICLoop_stream_ur
    -DNUMPROCS=4
    -DMAXSTEPS=12
    -DSRCDIRS=4-procs_case$<SEMICOLON>4
    -P ${CMAKE_CURRENT_SOURCE_DIR}/runphasta.cmake
    )

  
  if( ${GMI_SIM_FOUND} )
    set(casename ${testLabel}_loopStreamAdapt_incompressible)
    add_test(NAME ${casename}
      COMMAND ${CMAKE_COMMAND}
      -DNAME=${casename}
      -DWORKDIR=${CASES}/incompressibleAdapt
      -DINPFILE=${PHASTA_SOURCE_DIR}/phSolver/common/input.config
      -DMPIRUN=${MPIRUN}
      -DMPIRUN_PROCFLAG=${MPIRUN_PROCFLAG}
      -DEXE=${PHASTACHEF_BINARY_DIR}/chefPhastaICLoop_stream_adapt
      -DNUMPROCS=2
      -DMAXSTEPS=8
      -DSRCDIRS=2-procs_case$<SEMICOLON>4$<SEMICOLON>8
      -P ${CMAKE_CURRENT_SOURCE_DIR}/runphasta.cmake
      )

    set(casename ${testLabel}_loopSamStreamAdapt_incompressible)
    add_test(NAME ${casename}
      COMMAND ${CMAKE_COMMAND}
      -DNAME=${casename}
      -DWORKDIR=${CASES}/incompressibleAdapt
      -DINPFILE=${PHASTA_SOURCE_DIR}/phSolver/common/input.config
      -DMPIRUN=${MPIRUN}
      -DMPIRUN_PROCFLAG=${MPIRUN_PROCFLAG}
      -DEXE=${PHASTACHEF_BINARY_DIR}/chefPhastaICLoop_sam_stream_adapt
      -DNUMPROCS=2
      -DMAXSTEPS=8
      -DSRCDIRS=2-procs_case$<SEMICOLON>4
      -P ${CMAKE_CURRENT_SOURCE_DIR}/runphasta.cmake
      )

    set(casename ${testLabel}_loopFilesAdapt_incompressible)
    add_test(NAME ${casename}
      COMMAND ${CMAKE_COMMAND}
      -DNAME=${casename}
      -DWORKDIR=${CASES}/incompressibleAdapt
      -DINPFILE=${PHASTA_SOURCE_DIR}/phSolver/common/input.config
      -DMPIRUN=${MPIRUN}
      -DMPIRUN_PROCFLAG=${MPIRUN_PROCFLAG}
      -DEXE=${PHASTACHEF_BINARY_DIR}/chefPhastaICLoop_files_adapt
      -DNUMPROCS=2
      -DMAXSTEPS=8
      -DSRCDIRS=2-procs_case$<SEMICOLON>4
      -P ${CMAKE_CURRENT_SOURCE_DIR}/runphasta.cmake
      )
  endif()
endif()

if( ${phastaC_FOUND} )
  set(casename ${testLabel}_compressibleShockTube)
  add_test(NAME ${casename}
    COMMAND ${CMAKE_COMMAND}
    -DNAME=${casename}
    -DWORKDIR=${CASES}/compressibleShockTube/SAM_small
    -DINPFILE=${PHASTA_SOURCE_DIR}/phSolver/common/input.config
    -DMPIRUN=${MPIRUN}
    -DMPIRUN_PROCFLAG=${MPIRUN_PROCFLAG}
    -DEXE=${PHASTACHEF_BINARY_DIR}/chefPhastaCLoop_sam_stream_adapt
    -DNUMPROCS=4
    -DMAXSTEPS=100
    -DSRCDIRS=4-procs_case$<SEMICOLON>adapt-4-mdsMesh_bz2$<SEMICOLON>50
    -P ${CMAKE_CURRENT_SOURCE_DIR}/runphasta.cmake
    )
endif()
