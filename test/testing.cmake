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
add_test(
  NAME chefPhasta_mvProcsCaseDir
  COMMAND mv ${CDIR}/4-procs_case ${CDIR}/4-procs_case_chefPhasta
  WORKING_DIRECTORY ${CDIR}
)

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
add_test(
  NAME chefPhastaStream_mvProcsCaseDir
  COMMAND mv ${CDIR}/4-procs_case ${CDIR}/4-procs_case_chefPhastaStream
  WORKING_DIRECTORY ${CDIR}
)

set(chefPhastaChef_stream ${PHASTACHEF_BINARY_DIR}/chefPhastaChef_stream)
add_test(
  NAME chefPhastaChefStream_incompressible
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${chefPhastaChef_stream}
  WORKING_DIRECTORY ${CDIR}
)
add_test(
  NAME chefPhastaChefStream_mvProcsCaseDir
  COMMAND mv ${CDIR}/4-procs_case ${CDIR}/4-procs_case_chefPhastaChefStream
  WORKING_DIRECTORY ${CDIR}
)
