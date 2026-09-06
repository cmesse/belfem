set( LIBNAME maxwell )

set( SOURCES
        cl_MaxwellFactory.cpp
        cl_IWG_Maxwell.cpp
        cl_IWG_MaxwellPostproc.cpp
        cl_MaxwellPostprocessor.cpp
        cl_MaxwellBoundaryConditionFactory.cpp
        cl_Maxwell_FieldList.cpp
        cl_Maxwell_TMatrix.cpp
        matrices/mt_maxwell_phi.cpp
        matrices/mt_maxwell_h.cpp
        matrices/mt_maxwell_symmetry.cpp
        matrices/mt_maxwell_background.cpp
        matrices/mt_maxwell_l2_phi.cpp
        matrices/mt_maxwell_l2_h.cpp
        matrices/mt_maxwell_l2_b.cpp
)

include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/homology )
include_directories( ${BELFEM_SOURCE_DIR}/fem/interpolation )
include_directories( ${BELFEM_SOURCE_DIR}/fem/iwg )
include_directories( ${BELFEM_SOURCE_DIR}/fem/maxwell )
include_directories( ${BELFEM_SOURCE_DIR}/fem/kernel )
include_directories( ${BELFEM_SOURCE_DIR}/fem/postproc )
include_directories( ${BELFEM_SOURCE_DIR}/homology )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/math/graph )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/sources )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/ode )
include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/materials )
include_directories( ${BELFEM_SOURCE_DIR}/physics/database )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )


    set ( LIBLIST
            sparse
            bezier
            mesh
            sources
            integration
            ode
            database
            materials
            interpolation
            integration
            iwg
            homology
            kernel
            postproc
            maxwell
            thermal
    )

