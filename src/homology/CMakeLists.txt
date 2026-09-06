set( LIBNAME homology )

set( SOURCES
        cl_CutFactory.cpp
        cl_Chain.cpp
        cl_Cochain.cpp
        cl_SimplicialComplex.cpp
        cl_Homology.cpp
        cl_Cohomology.cpp
        cl_BeltedTree.cpp
        cl_CutProcessorManual.cpp
        cl_CutProcessor.cpp
        cl_InterfaceProcessor.cpp
        cl_CutData.cpp
        cl_CutSet.cpp
        cl_Topology.cpp
        fn_Smith.cpp
)

include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/materials )
include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/spline )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/sources )
include_directories( ${BELFEM_SOURCE_DIR}/fem/interpolation )
include_directories( ${BELFEM_SOURCE_DIR}/fem/iwg )
include_directories( ${BELFEM_SOURCE_DIR}/fem/kernel )
include_directories( ${BELFEM_SOURCE_DIR}/fem/postproc )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/math/graph )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )