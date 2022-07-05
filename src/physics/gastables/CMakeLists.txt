set( LIBNAME gastables )

set( SOURCES
        fn_GT_data_path.cpp
        fn_GT_fix_label.cpp
        fn_GT_fix_capitals.cpp
        cl_GT_HeatPoly.cpp
        cl_GT_HeatPolyCustom.cpp
        cl_GT_HeatPolyGlue.cpp
        cl_GT_TransportPoly.cpp
        cl_GT_TransportPolyCustom.cpp
        cl_GT_TransportPolyGlue.cpp
        cl_GT_TransportPolyEmpty.cpp
        cl_GT_GasData.cpp
        cl_GT_InputThermo.cpp
        cl_GT_InputTransport.cpp
        cl_GT_InputData.cpp
        cl_GT_InputAlpha.cpp
        fn_GT_idgas_mu.cpp
        fn_GT_idgas_lambda.cpp
        cl_GT_RefGas.cpp
        cl_GT_RefGasFactory.cpp
        cl_GT_Arguments.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )
include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/gastables )
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

set ( LIBLIST
        sparse
        gastables
         )

if( USE_EXAMPLES )
    set( EXECNAME gastable )
    set( MAIN     main.cpp )
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
endif()