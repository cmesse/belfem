//
// Created by Christian Messe on 2019-07-25.
//

#ifndef BELFEM_CL_MESH_GMSHREADER_HPP
#define BELFEM_CL_MESH_GMSHREADER_HPP

#include "typedefs.hpp"
#include "cl_Ascii.hpp"
#include "cl_Cell.hpp"
#include "cl_Node.hpp"
#include "cl_Element.hpp"
#include "cl_Facet.hpp"
#include "cl_Map.hpp"
#include "cl_Block.hpp"
#include "cl_SideSet.hpp"
#include "cl_Mesh.hpp"

namespace belfem
{
    namespace mesh
    {
        /*struct GmshPhysicalGroup
        {
            string mLabel = "untitled";
            uint mID = 0;
            uint mDimension = 0;
        };*/

//------------------------------------------------------------------------------

        class GmshReader : public Ascii
        {
            string mFilename;

            // scale factor to convert spatial unit of mesh into m
            const real mMeshScale = 1.0 ;

            real   mVersion = 0.0;

            uint   mNumberOfNodes = 0;
            uint   mNumberOfElements = 0;
            uint   mNumberOfPhysicalGroups = 0;

            // tags in the file
            size_t mNodesTag = 0;
            size_t mElementsTag = 0;
            size_t mPhysicalTag = 0;

            // Pointer to mesh object
            Mesh * mMesh;

            // dimensions
            uint   mNumberOfDimensions = 0;

            // cell containing Nodes
            Cell< Node * > mNodes;

            // cell containing Elements
            Cell< Element * > mElements;

            // map managing Nodes by index
            Map< id_t, Node* > mNodeMap;

            // number of elements per dimension
            Vector< index_t > mNumberOfElementsPerDimension;

            Cell< Element * > mVertices;
            Cell< Element * > mEdges;
            Cell< Element * > mFaces;
            Cell< Element * > mVolumes;

            // groups ids
            Vector< uint > mNodeGroupIDs;
            Vector< uint > mEdgeGroupIDs;
            Vector< uint > mFaceGroupIDs;
            Vector< uint > mVolumeGroupIDs;

            /**
             * flag telling if mesh is deleted on destruction
             * once the mesh is called using this->mesh()
             * or passed in the constructor, it is the
             *   responsibility of the user to delete it
             */
            bool mOwnMesh = false;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            GmshReader( const string & aPath, Mesh * aMesh = nullptr, const real aMeshScale = 1.0 );

//------------------------------------------------------------------------------

            ~GmshReader();

//------------------------------------------------------------------------------

            Mesh *
            mesh()
            {

                mOwnMesh = false;

                return mMesh;
            }

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            tidy_up_buffer();

//------------------------------------------------------------------------------

            void
            check_tag_existence();

//------------------------------------------------------------------------------

            void
            read_version();

//------------------------------------------------------------------------------

            void
            read_tags_v22();

//------------------------------------------------------------------------------

            void
            read_nodes_tag_v22( size_t & aCount );

//------------------------------------------------------------------------------

            void
            read_elements_tag_v22( size_t & aCount );

//------------------------------------------------------------------------------

            void
            read_physical_tag( size_t & aCount );

//------------------------------------------------------------------------------

            void
            read_nodes_v22();

//------------------------------------------------------------------------------

            void
            read_elements_v22();

//------------------------------------------------------------------------------

            void
            read_mesh_v41();

//------------------------------------------------------------------------------

            void
            read_entity_tag_v41( size_t & aCount );

//------------------------------------------------------------------------------

            void
            read_nodes_v41( size_t & aCount );

//------------------------------------------------------------------------------

            void
            read_elements_v41( size_t & aCount );

//------------------------------------------------------------------------------

            void
            convert_orders_for_quadratic_volume_elements_to_exo();

//------------------------------------------------------------------------------

            void
            create_node_map();

//------------------------------------------------------------------------------

            /**
             * determine if this is a 2D or a 3D mesh
             */
            void
            read_mesh_dimension();

//------------------------------------------------------------------------------

            void
            create_element_containers_per_dimension();

//------------------------------------------------------------------------------

            void
            create_group_ids();

//------------------------------------------------------------------------------

            void
            create_blocks();

//------------------------------------------------------------------------------

            void
            check_element_orientation();

//------------------------------------------------------------------------------

            void
            create_sidesets();

//------------------------------------------------------------------------------

            void
            delete_unused_nodes_and_elements();

//------------------------------------------------------------------------------

            void
            create_mesh();

//------------------------------------------------------------------------------

            void
            create_vertices();

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_MESH_GMSHREADER_HPP
