macro(serve TESTNAME PARTS FACTOR WORKDIR)
  set(exe ${PHASTACHEF_BINARY_DIR}/chef_phasta)
  math(EXPR OUTPARTS "${PARTS} * ${FACTOR}")
  add_test(NAME "${TESTNAME}"
    COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} ${OUTPARTS} ${exe}
    WORKING_DIRECTORY "${WORKDIR}")
endmacro()

if (PCU_COMPRESS)
  set(MDIR ${CASES}/crossflow/1-1-Chef-Tet-Part/run)
  serve(chef_phasta0 1 1 ${MDIR})
  set(MDIR ${CASES}/crossflow/1-1-Chef-Tet-Part)
  add_test(NAME chef_phasta1
    COMMAND diff -r -x .svn out_mesh/ good_mesh/
    WORKING_DIRECTORY ${MDIR})
  set(MDIR ${CASES}/crossflow/2-1-Chef-Tet-Part/run)
  serve(chef_phasta2 1 2 ${MDIR})
  set(MDIR ${CASES}/crossflow/2-1-Chef-Tet-Part/4-2-Chef-Part/run)
  serve(chef_phasta3 2 2 ${MDIR})
  set(MDIR ${CASES}/crossflow/4-1-Chef-Tet-Part/run)
  serve(chef_phasta4 1 4 ${MDIR})
  set(MDIR ${CASES}/crossflow/4-1-Chef-Tet-Part/4-4-Chef-Part-ts20/run)
  serve(chef_phasta5 4 1 ${MDIR})
  set(MDIR ${CASES}/crossflow/4-1-Chef-Tet-Part/4-4-Chef-Part-ts20)
  add_test(NAME chef_phasta6
    COMMAND diff -r -x .svn out_mesh/ good_mesh/
    WORKING_DIRECTORY ${MDIR})
else()
  message(WARNING "Testing disabled. Build with PCU_COMPRESS=ON.")
endif()
