set( LIBNAME thermal )

set( SOURCES
        cl_ThermalFactory.cpp
        cl_ThermalBoundaryConditionFactory.cpp
        cl_IWG_MaxwellThermal.cpp
        matrices/mt_thermal_phi.cpp
        matrices/mt_thermal_h.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/fem/interpolation )
include_directories( ${BELFEM_SOURCE_DIR}/fem/iwg )
include_directories( ${BELFEM_SOURCE_DIR}/fem/thermal )
include_directories( ${BELFEM_SOURCE_DIR}/fem/kernel )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/math/graph )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/sources )
include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/materials )
include_directories( ${BELFEM_SOURCE_DIR}/physics/database )

set ( LIBLIST
        sparse
        mesh
        integration
        interpolation
        materials
        iwg
        kernel
        )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )