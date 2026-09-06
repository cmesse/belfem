set( LIBNAME kernel )

set( SOURCES
        en_FEM_BoundaryConditionType.cpp
        cl_FEM_Dof.cpp
        cl_FEM_Element.cpp
        cl_FEM_Tmatrix.cpp
        cl_FEM_Calculator.cpp
        cl_FEM_Group.cpp
        cl_FEM_Block.cpp
        cl_FEM_SideSet.cpp
        cl_FEM_Bearing.cpp
        cl_FEM_Controller.cpp
        cl_FEM_DofManagerBase.cpp
        cl_FEM_DofMgr_Parameters.cpp
        cl_FEM_DofMgr_DofData.cpp
        cl_FEM_DofMgr_BlockData.cpp
        cl_FEM_DofMgr_SideSetData.cpp
        cl_FEM_DofMgr_BearingData.cpp
        cl_FEM_DofMgr_FieldData.cpp
        cl_FEM_DofMgr_SolverData.cpp
        cl_FEM_DofMgr_EigenValues.cpp
        cl_FEM_DofManager.cpp
        cl_FEM_PhysicalBoundaryCondition.cpp
        cl_FEM_KernelParameters.cpp
        cl_FEM_Kernel.cpp
        cl_FEM_Domain.cpp
        cl_FEM_Postprocessor.cpp
        cl_MeshChecker.cpp
        cl_Pipette.cpp
        cl_ThinShellFactory.cpp
        cl_EdgeCutter.cpp
)

include_directories( ${BELFEM_SOURCE_DIR}/physics )
include_directories( ${BELFEM_SOURCE_DIR}/physics/materials )


include_directories( ${BELFEM_SOURCE_DIR}/mesh )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/math/graph )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/sources )

include_directories( ${BELFEM_SOURCE_DIR}/fem/interpolation )
include_directories( ${BELFEM_SOURCE_DIR}/fem/iwg )
include_directories( ${BELFEM_SOURCE_DIR}/fem/kernel )
include_directories( ${BELFEM_SOURCE_DIR}/fem/postproc )



# LIBLIST feeds Add_Executable link lines only ( Add_Library ignores it );
# circuit removed 2026-08-15: no kernel source references circuit symbols —
# the real executables link circuit in src/executables/CMakeLists.txt
set ( LIBLIST
        sparse
        mesh
        integration
        sources
        interpolation
        materials
        postproc
        kernel
        iwg
        kernel
        )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

if( USE_EXAMPLES )

    set( EXECNAME petsctest )
    set( MAIN petsctest.cpp)
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

    set( EXECNAME poisson )
    set( MAIN poisson.cpp)
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

    #set( EXECNAME staticheat )
    #set( MAIN staticheat.cpp)
    #include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

    set( EXECNAME hangingnodes )
    set( MAIN hangingnodes.cpp)
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
endif()
