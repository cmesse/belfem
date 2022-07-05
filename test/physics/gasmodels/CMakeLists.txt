# List source files
set( TESTNAME gasmodels )

set( SOURCES
        cl_GM_Gas_Airprop.cpp
        cl_GM_Gas_Gibbs.cpp
        cl_GM_EoS_AlphaFunction.cpp
        cl_GM_EoS_Cubic_State.cpp
        cl_GM_EoS_Cubic_Departure.cpp
        cl_GM_EoS_Cubic_AlphaBetaKappa.cpp
        cl_GM_Eos_Cubic_Entropy.cpp
        cl_GM_Gas_Highpressure_Transport.cpp
        cl_GM_Gas_PrandtlMeyer.cpp
        cl_GM_EoS_Hydrogen_Vapor.cpp
        cl_GM_EoS_Hydrogen_Caloric.cpp
        cl_GM_EoS_Oxygen_Vapor.cpp
        cl_GM_EoS_Oxygen_Caloric.cpp
        cl_GM_EoS_Methane_Vapor.cpp
        cl_GM_EoS_Methane_Caloric.cpp
      )

include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )
include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/gastables )
include_directories( ${BELFEM_SOURCE_DIR}/physics/gasmodels )
set ( LIBLIST
        gastables
        gasmodels )

# add the test
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Test.cmake )