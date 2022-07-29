set( LIBNAME materials )

set( SOURCES
        en_Materials.cpp
        cl_Material.cpp
        cl_IsotropicMaterial.cpp
        cl_OrthotropicMaterial.cpp
        cl_Material_AltraMat80.cpp
        cl_Material_Aluminum.cpp
        cl_Material_Copper.cpp
        cl_Material_CuCrZr.cpp
        cl_Material_Ti6Al4V.cpp
        cl_Material_CCSIC.cpp
        cl_Material_Inconel718.cpp
        cl_Material_Rohacell51.cpp
        cl_Material_Zirconia.cpp
        cl_Material_Inconel750X.cpp
        cl_Material_HastelloyC276.cpp
        cl_Material_Silver.cpp
        cl_Material_YBCO.cpp
        cl_Material_Pb40Sn60.cpp
        cl_Material_SAE301.cpp
        cl_Material_SAE314.cpp
	    cl_MaterialFactory.cpp
	    nist_functions.cpp
        )

if ( USE_GASMODELS )
    set(SOURCES
            ${SOURCES}
            cl_Material_Air.cpp)
endif()

include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/physics )


if ( USE_GASMODELS )
    include_directories( ${BELFEM_SOURCE_DIR}/physics/gastables )
    include_directories( ${BELFEM_SOURCE_DIR}/physics/gasmodels )
    set ( LIBLIST
            gastables
            gasmodels
            materials
            )
else()
    set ( LIBLIST
            materials
            )
endif()

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )



if( USE_EXAMPLES )
    set( EXECNAME material )
    set( MAIN     main.cpp )
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
endif()
