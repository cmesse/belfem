# set the name of the test
set( TESTNAME containers )
set( TESTLABELS fast )

# list the test sources
set( SOURCES
        test_Cell.cpp
        test_DynamicBitset.cpp
        test_Bitset.cpp
        test_ShiftRegister.cpp
        test_Map.cpp
        test_Set.cpp
        test_Queue.cpp
        test_StringList.cpp
        test_Genome.cpp
        )

# add the test
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )
