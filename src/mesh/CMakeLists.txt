set( LIBNAME mesh )

set( SOURCES
        cl_ReferenceElement.cpp
        meshtools.cpp
        vtktools.cpp
        geometrytools.cpp
        fn_Mesh_ratio.cpp
        fn_entity_type.cpp
        fn_read_node_tmatrix_from_hdf5.cpp
        fn_to_master_orientation.cpp
        fn_to_master_orientation.cpp
        fn_Mesh_symrcm_nodes.cpp
        cl_Mesh_Basis.cpp
        cl_Vertex.cpp
        cl_Node.cpp
        cl_Edge.cpp
        cl_Face.cpp
        cl_Facet.cpp
        cl_Element.cpp
        cl_ControlPoint.cpp
        cl_Element_Factory.cpp
        cl_Block.cpp
        cl_SideSet.cpp
        cl_ThinShell.cpp
        cl_Mesh_GlobalVariable.cpp
        cl_Mesh_Field.cpp
        cl_Mesh.cpp
        cl_Mesh_GmshReader.cpp
        cl_Mesh_ExodusWriter.cpp
        cl_Mesh_SourceExpander.cpp
        cl_ProtoMesh.cpp
        cl_Mesh_BfmFile.cpp
        cl_Mesh_VtkWriter.cpp
        cl_Mesh_OrderConverter.cpp
        cl_Mesh_Distributor.cpp
        cl_Mesh_Partitioner.cpp
        cl_Mesh_CurvedElementChecker.cpp
        cl_Mesh_ConnectivityCalculator.cpp
        cl_Mesh_Periodicity.cpp
        cl_Mesh_PeriodicityFactory.cpp
        cl_TensorMeshConfig.cpp
        cl_TensorMeshFactory.cpp
        cl_EdgeFactory.cpp
        cl_FaceFactory.cpp
        cl_Segment.cpp
        cl_Curve.cpp
        cl_CurveFactory.cpp
        en_DomainType.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/core )
include_directories( ${BELFEM_SOURCE_DIR}/comm )
include_directories( ${BELFEM_SOURCE_DIR}/io )
include_directories( ${BELFEM_SOURCE_DIR}/containers )
include_directories( ${BELFEM_SOURCE_DIR}/math/tools )
include_directories( ${BELFEM_SOURCE_DIR}/numerics/integration )
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

set ( LIBLIST
        integration
        mesh
         )

# msh2exo writes exodus. Without that backend it would still build and link --
# ExodusWriter compiles to stubs -- and then produce nothing, so the target is
# gated instead. A name in BELFEM_INSTALL_EXECUTABLES with no target is ignored
# ( config/globals.cmake ), so the install entry needs no guard
if( USE_EXODUS )
    set( EXECNAME msh2exo )
    set( MAIN     msh2exo.cpp )
    include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
endif()

#set( EXECNAME meshtest )
#set( MAIN     main.cpp )
#include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )

