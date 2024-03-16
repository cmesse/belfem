set( LIBNAME kernel )

set( SOURCES
        en_FEM_DomainType.cpp
        cl_FEM_Dof.cpp
        cl_FEM_Element.cpp
        cl_FEM_Calculator.cpp
        cl_FEM_Group.cpp
        cl_FEM_Block.cpp
        cl_FEM_SideSet.cpp
        cl_FEM_Shell.cpp
        cl_FEM_Cut.cpp
        cl_FEM_Bearing.cpp
        cl_FEM_DofManagerBase.cpp
        cl_FEM_DofMgr_Parameters.cpp
        cl_FEM_DofMgr_DofData.cpp
        cl_FEM_DofMgr_BlockData.cpp
        cl_FEM_DofMgr_SideSetData.cpp
        cl_FEM_DofMgr_BearingData.cpp
        cl_FEM_DofMgr_FieldData.cpp
        cl_FEM_DofMgr_SolverData.cpp
        cl_FEM_DofManager.cpp
        cl_FEM_BoundaryCondition.cpp
        cl_FEM_KernelParameters.cpp
        cl_FEM_Kernel.cpp
        cl_FEM_DomainGroup.cpp
        FEM_geometry.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/materials )


include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include_directories( ${BELFEM_SOURCE_DIR}/fem/interpolation )
include_directories( ${BELFEM_SOURCE_DIR}/fem/edgeinterpolation )
include_directories( ${BELFEM_SOURCE_DIR}/fem/iwg )
include_directories( ${BELFEM_SOURCE_DIR}/fem/postproc )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/math/graph )


if( USE_GASMODELS )
    include_directories( ${BELFEM_SOURCE_DIR}/physics/gastables )
    include_directories( ${BELFEM_SOURCE_DIR}/physics/gasmodels )

    set ( LIBLIST
            sparse
            mesh
            gastables
            gasmodels
            integration
            interpolation
            materials
            postproc
            kernel
            iwg
            kernel
            )
else()
    set ( LIBLIST
            sparse
            mesh
            integration
            interpolation
            materials
            postproc
            kernel
            iwg
            kernel
            )
endif()

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

if( USE_EXAMPLES )
    set( EXECNAME tensiletest )
    set( MAIN tensiletest.cpp)
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )


    set( EXECNAME petsctest )
    set( MAIN petsctest.cpp)
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

    set( EXECNAME staticheat )
    set( MAIN staticheat.cpp)
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

    set( EXECNAME transientheat )
    set( MAIN transientheat.cpp)
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
endif()