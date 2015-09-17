macro(moveDir name work src)
  set(tgtdir ${work}/${src}_${name})
  add_test(
    NAME ${name}_rm${src}
    COMMAND rm -rf ${tgtdir} && true #don't report a failure of rm
    WORKING_DIRECTORY ${work}
  )
  add_test(
    NAME ${name}_mv${src}
    COMMAND mv ${work}/${src} ${tgtdir}
    WORKING_DIRECTORY ${work}
  )
endmacro(moveDir)

set(CDIR ${CASES}/incompressible)

set(chefPhasta_posix ${PHASTACHEF_BINARY_DIR}/chefPhasta_posix)
add_test(chefPhastaStream_copyInpCfg
  cp ${PHASTA_SOURCE_DIR}/phSolver/common/input.config ${CDIR})

add_test(
  NAME chefPhasta_incompressible
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${chefPhasta_posix}
  WORKING_DIRECTORY ${CDIR}
)
set(cmd
  ${PHASTA_BINARY_DIR}/bin/checkphasta
  ${CDIR}/4-procs_case/
  ${CDIR}/4-procs_case-SyncIO-2_ref/
  2 1e-6)
add_test(
  NAME chefPhasta_compareIncompressible
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${cmd}
  WORKING_DIRECTORY ${CDIR}
)
moveDir(chefPhasta ${CDIR} 4-procs_case)

set(chefPhasta_stream ${PHASTACHEF_BINARY_DIR}/chefPhasta_stream)
add_test(
  NAME chefPhastaStream_incompressible
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${chefPhasta_stream}
  WORKING_DIRECTORY ${CDIR}
)
set(cmd
  ${PHASTA_BINARY_DIR}/bin/checkphasta
  ${CDIR}/4-procs_case/
  ${CDIR}/4-procs_case-SyncIO-2_ref/
  2 1e-6)
add_test(
  NAME chefPhastaStream_compareIncompressible
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${cmd}
  WORKING_DIRECTORY ${CDIR}
)
moveDir(chefPhastaStream ${CDIR} 4-procs_case)

set(chefPhastaChef_stream ${PHASTACHEF_BINARY_DIR}/chefPhastaChef_stream)
add_test(
  NAME chefPhastaChefStream_incompressible
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${chefPhastaChef_stream}
  WORKING_DIRECTORY ${CDIR}
)
moveDir(chefPhastaChefStream ${CDIR} 4-procs_case)

set(chefPhastaLoop_stream ${PHASTACHEF_BINARY_DIR}/chefPhastaLoop_stream)
set(maxTimeStep 8)
add_test(
  NAME chefPhastaLoopStream_incompressible
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${chefPhastaLoop_stream} ${maxTimeStep}
  WORKING_DIRECTORY ${CDIR}
)
moveDir(chefPhastaLoopStream ${CDIR} 4-procs_case)

set(chefPhastaLoop_stream_ur ${PHASTACHEF_BINARY_DIR}/chefPhastaLoop_stream_ur)
set(maxTimeStep 12)
add_test(
  NAME chefPhastaLoopStreamUR_incompressible
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${chefPhastaLoop_stream_ur} ${maxTimeStep}
  WORKING_DIRECTORY ${CDIR}
)
moveDir(chefPhastaLoopStreamUR ${CDIR} 4-procs_case)
moveDir(chefPhastaLoopStreamUR ${CDIR} 4)
