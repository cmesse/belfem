# Tier 1: always-run tests (no MPI required)
set( TESTNAME comm )

set( SOURCES
        test_CommUtils.cpp
        test_Communicator.cpp
        )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )

# Tier 2: the MPI tests. A separate binary because it needs a different main()
# ( rank sentinel + Allreduce over the result ) and because it must never be
# registered serially - see the Tier 1 / Tier 2 note in Add_Test.cmake.
#
# TESTRANKS 2 4: two ranks is the floor almost every test here guards on, and
# four is required by RingPassVector4Ranks. Dropping 4 would leave that test
# permanently skipped, which is the failure this tier exists to close. On a
# host with fewer than four cores the launcher oversubscribes; these are
# correctness runs, not performance runs.
set( TESTNAME commmpi )

set( SOURCES
        test_CommMPI.cpp
        )

set( TESTRANKS 2 4 )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )
