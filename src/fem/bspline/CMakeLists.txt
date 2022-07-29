set( LIBNAME bspline )

set( SOURCES
        cl_BS_TMatrix.cpp
        cl_BS_Basis.cpp
        cl_BS_Element.cpp
        cl_BS_Mapper.cpp
        cl_BS_LookupTable.cpp

        )

if ( USE_GASMODELS )
    set(SOURCES
            ${SOURCES}
            cl_GM_EoS_TableGas.cpp
            cl_TableGas.cpp)
endif()

include_directories( ${BELFEM_SOURCE_DIR}/math/graph )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include_directories( ${BELFEM_SOURCE_DIR}/fem/interpolation )
include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/materials )
if ( USE_GASMODELS )
    include_directories( ${BELFEM_SOURCE_DIR}/physics/gastables )
    include_directories( ${BELFEM_SOURCE_DIR}/physics/gasmodels )
    include_directories( ${BELFEM_SOURCE_DIR}/physics/atmosphere )
    include_directories( ${BELFEM_SOURCE_DIR}/physics/boundarylayer )
endif()

set ( LIBLIST
        graph
        integration
        interpolation
        mesh
        bspline
        materials
        )

if ( USE_GASMODELS )
    set(LIBLIST
            ${LIBLIST}
            gastables
            gasmodels
            atmosphere
            boundarylayer )
endif()

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

if( USE_EXAMPLES AND USE_GASMODELS )
    set( EXECNAME makehotair )
    set( MAIN makehotair.cpp)
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

    set( EXECNAME hotairtable )
    set( MAIN hotairtable.cpp)
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

    set( EXECNAME stagnation )
    set( MAIN stagnation.cpp)
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
endif()

set( EXECNAME makematdb )
set( MAIN makematdb.cpp)
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

set( EXECNAME makematsp )
set( MAIN makematsp.cpp)
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
