set( LIBNAME mesh )

set( SOURCES
        meshtools.cpp
        vtktools.cpp
        geometrytools.cpp
        fn_Mesh_ratio.cpp
        fn_entity_type.cpp
        cl_Vertex.cpp
        cl_Node.cpp
        cl_Edge.cpp
        cl_Face.cpp
        cl_Facet.cpp
        cl_Element.cpp
        cl_Element_Factory.cpp
        cl_Block.cpp
        cl_SideSet.cpp
        cl_Mesh_GlobalVariable.cpp
        cl_Mesh_Field.cpp
        cl_Mesh.cpp
        cl_Mesh_GmshReader.cpp
        cl_Mesh_ExodusWriter.cpp
        cl_Mesh_HDF5Writer.cpp
        cl_Mesh_HDF5Reader.cpp
        cl_Mesh_VtkWriter.cpp
        cl_Mesh_OrderConverter.cpp
        cl_Mesh_Partitioner.cpp
        cl_Mesh_Scissors.cpp
        cl_Mesh_TapeRoller.cpp
        cl_Mesh_CurvedElementChecker.cpp
        cl_TensorMeshFactory.cpp
        cl_EdgeFactory.cpp
        cl_FaceFactory.cpp
        )

include_directories( ${BELFEM_SOURCE_DIR}/core )
include_directories( ${BELFEM_SOURCE_DIR}/comm )
include_directories( ${BELFEM_SOURCE_DIR}/containers )

include( ${BELFEM_CONFIG_DIR}/scripts/Add_Library.cmake )

set ( LIBLIST
        mesh )

set( EXECNAME meshtest )
set( MAIN     main.cpp )
include( ${BELFEM_CONFIG_DIR}/scripts/Add_Executable.cmake )
