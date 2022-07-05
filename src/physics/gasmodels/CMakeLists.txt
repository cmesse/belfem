set( LIBNAME gasmodels )

set( SOURCES
        cl_Gas.cpp
        cl_GM_EoS.cpp
        cl_GM_EoS_AlphaFunction.cpp
        cl_GM_EoS_AlphaFunctionFactory.cpp
        cl_GM_EoS_Idgas.cpp
        cl_GM_EoS_Cubic.cpp
        cl_GM_Helmholtz.cpp
        cl_GM_EoS_Hydrogen.cpp
        cl_GM_EoS_Oxygen.cpp
        cl_GM_EoS_Methane.cpp
        cl_GM_HelmholtzTransport.cpp
        cl_GM_HelmholtzTransport_Methane.cpp
        fn_GM_Helmholtz_DerivTest.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )
include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/gastables )
include_directories( ${BELFEM_SOURCE_DIR}/physics/gasmodels )
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

set ( LIBLIST
        sparse
        gastables
        gasmodels )

if( USE_EXAMPLES )
    set( EXECNAME gasmodel )
    set( MAIN     main.cpp )
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

    set( EXECNAME formation )
    set( MAIN     formation.cpp )
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

    set( EXECNAME helmholtz )
    set( MAIN     helmholtz.cpp )
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

endif()