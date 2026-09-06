# set the name of the test
set( TESTNAME sparse )

# list the test sources
set( SOURCES
        test_SpMatrix.cpp
        test_Solver.cpp
        )

# add the test
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )

# Tier 2: the parallel nested-dissection wiring ( ParMETIS / PT-Scotch through
# DistMatrix and the PETSc path ). A separate binary because it needs the Tier 2
# main ( launcher sentinel + verdict fold ) and must never be registered
# serially - see the Tier 1 / Tier 2 note in Add_Test.cmake. Gated on the
# libraries at CMake level: Add_Test.cmake registers on USE_MPI and TESTRANKS
# alone, and an #ifdef inside the sources would leave a backend-less binary
# green. No `fast` label: this must not leak into check-fast.
if( USE_METIS OR USE_SCOTCH )

    set( TESTNAME sparsempi )

    set( SOURCES
            test_SparseMPI.cpp
            )

    set( TESTRANKS 2 4 )

    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )

endif()
