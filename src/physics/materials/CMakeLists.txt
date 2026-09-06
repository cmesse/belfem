set( LIBNAME materials )

set( SOURCES
        cl_Material_Copper.cpp
        cl_Material_Aluminum.cpp
        cl_Material_Chromium.cpp
        cl_Material_Silver.cpp
        cl_Material_Indium.cpp
        cl_Material_YBCO.cpp
        cl_Material_HastelloyC276.cpp
        cl_Material_Magnesia.cpp
        cl_Material_Nickel.cpp
        cl_Material_Iron.cpp
        cl_Material_Alloy.cpp
        cl_Material_Ferromagnetic.cpp
        cl_Material_Metal.cpp
        cl_Material_UserDefined.cpp
        cl_Material_SplineLookupTable.cpp
        cl_Material_BhSplineCurve.cpp
        cl_Material.cpp
        cl_BhCurve.cpp
        cl_JcFunction.cpp
        cl_Material_Abundance.cpp
        cl_MaterialFactory.cpp
        cl_Material_WhiteTin.cpp
        cl_Material_Lead.cpp
        fn_material_data_path.cpp
        fn_embed_python_guide.cpp
        debye.f90
        )

include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/bezier )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include_directories( ${CMAKE_BINARY_DIR}/src/numerics/integration )
include_directories( ${BELFEM_SOURCE_DIR}/math/tensor )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/database )
include_directories( ${BELFEM_SOURCE_DIR}/mesh )
set ( LIBLIST
        mesh
        integration
        interpolation
        bezier
        database
        materials
)

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

# Ensure integration library is built first (needed for Fortran module intpoints.mod)
add_dependencies( ${BELFEM_MODULE_TARGET} ${LIBPREFIX}_integration )

set( EXECNAME material )
set( MAIN     main.cpp )
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
# only needed for development purposes
set( EXECNAME mattest )
set( MAIN mattest.cpp)
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
