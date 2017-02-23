macro(cmd dir exe)
  message("${exe} ${ARGN}")
  execute_process(
    COMMAND ${exe} ${ARGN}
    WORKING_DIRECTORY ${dir}
    OUTPUT_VARIABLE out
    ERROR_VARIABLE out
    RESULT_VARIABLE res
    )
  message("${out}")
  if(res)
    message(FATAL_ERROR "Error running ${exe}")
  else()
    message("Success")
  endif()
endmacro()

cmd(${WORKDIR} cp ${INPFILE} ${WORKDIR})
cmd(${WORKDIR} ${MPIRUN} ${MPIRUN_PROCFLAG} ${NUMPROCS} ${EXE} ${MAXSTEPS})
message("running rm and mv")
foreach(src ${SRCDIRS})
  set(tgtdir ${WORKDIR}/${src}_${NAME})
  cmd(${WORKDIR} rm -rf ${tgtdir})
  cmd(${WORKDIR} mv ${WORKDIR}/${src} ${tgtdir})
endforeach()
