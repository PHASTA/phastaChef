set(CDIR ${CASES}/incompressible)
set(chef_phasta_stream ${PHASTACHEF_BINARY_DIR}/chef_phasta_stream)
add_test(chefPhastaStream_copyInpCfg
  cp ${PHASTA_SOURCE_DIR}/phSolver/common/input.config ${CDIR})
add_test(
  NAME chefPhastaStream_incompressible-sync
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${chef_phasta_stream}
  WORKING_DIRECTORY ${CDIR}
)
set(cmd
  ${PHASTA_BINARY_DIR}/bin/checkphasta
  ${CDIR}/4-procs_case/
  ${CDIR}/4-procs_case-SyncIO-2_ref/
  2 1e-6)
add_test(
  NAME chefPhastaStream_compareIncompressibleSync
  COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4 ${cmd}
  WORKING_DIRECTORY ${CDIR}
)
