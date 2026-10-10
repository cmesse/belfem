# set the name of the test
set( TESTNAME homology )
set( TESTLABELS fast )

include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/homology )

# list the test sources
set( SOURCES
        test_Cohomology.cpp
        test_CohomologyPeriodic.cpp
        test_CutData.cpp
        test_CutProcessor.cpp
        test_CutFactory_closed_rim.cpp
        test_CutFactory_unwind.cpp
        test_CutFactory_terminal_loop.cpp
        test_InterfaceProcessor.cpp
        )

# Add_Test reverses this list, and dependents must precede their
# dependencies on the link line: homology written last lands first
set( LIBLIST graph bezier interpolation integration mesh homology )

# add the test
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )
