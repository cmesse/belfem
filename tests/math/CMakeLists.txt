# set the name of the test
set( TESTNAME math )
set( TESTLABELS fast )

include_directories( ${BELFEM_SOURCE_DIR}/math/quaternion )
include_directories( ${BELFEM_SOURCE_DIR}/math/tensor )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )

# list the test sources
set( SOURCES
        test_Quaternion.cpp
        test_Tensor.cpp
        test_TensorKernels.cpp
        test_MathTools.cpp
        test_GraphVertex.cpp
        test_GraphAlgorithms.cpp
        test_GraphTools.cpp
        test_GraphSpfa.cpp
        test_Spline.cpp
        test_Polynomials.cpp
        )

set( LIBLIST spline )

# add the test
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )
