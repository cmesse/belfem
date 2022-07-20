//
// Created by christian on 3/30/22.
//

#ifndef BELFEM_CL_MESH_TAPEROLLER_HPP
#define BELFEM_CL_MESH_TAPEROLLER_HPP

#include "cl_Mesh.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Map.hpp"
#include "st_Mesh_Layer.hpp"

namespace belfem
{
    namespace mesh
    {
        // todo: fix ownership of new elements after creation!
        class TapeRoller
        {
            const proc_t mMyRank ;
            Mesh * mMesh ;

            const uint mNumberOfBlocks ; // physical layers
            const uint mElementOrder ;
            const uint mNumberOfGhostLayers ; // dof layers

            Cell< id_t > mSelectedSideSets ;
            Cell< id_t > mMasterBlocks ;

            id_t mMaxNodeId ;
            id_t mMaxEdgeId ;
            id_t mMaxFaceId ;
            id_t mMaxElementID ;

            id_t mMaxSideSetId ;

            // ids with elements that are to be cloned
            Vector< index_t > mElementIDs ;

            Cell< Layer * > mGhostLayers ;

            // map that links old node ids to position in layer
            Map< id_t, index_t > mNodeMap ;
            Map< id_t, index_t > mEdgeMap ;
            Map< id_t, index_t > mElementMap ;

            Vector< id_t > mGhostSideSetIDs ;

            // backup container for orientation flipping
            Cell< mesh::Node * > mOrientationBackup ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            TapeRoller( Mesh * aMesh, const uint aNumberOfLayers, const uint aElementOrder );

//------------------------------------------------------------------------------

            ~TapeRoller();

//------------------------------------------------------------------------------

            void
            add_sidesets( const Vector< id_t > & aSideSetIDs );

//------------------------------------------------------------------------------

            void
            get_sidesets( Vector< id_t > & aSideSetIDs );

//------------------------------------------------------------------------------

            void
            add_master_blocks( const Vector< id_t > & aMasterBlockIDs );

//------------------------------------------------------------------------------

            /**
             * returns the max block id before tape creation
             */
            id_t
            run();

//------------------------------------------------------------------------------

            /**
             * returns the ids of the ghost sidesets (call after this->run() )
             */
             const Vector< id_t > &
             ghost_sideset_ids() const;

//------------------------------------------------------------------------------

            void
            flip_element_orientation();

//------------------------------------------------------------------------------

            void
            revert_element_orientation();

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
            void
            get_max_ids();

//------------------------------------------------------------------------------

            void
            create_layers();

//------------------------------------------------------------------------------

            void
            clone_nodes();

//------------------------------------------------------------------------------

            void
            clone_elements();

//------------------------------------------------------------------------------

            void
            clone_edges();

//------------------------------------------------------------------------------

            void
            clone_faces();

//------------------------------------------------------------------------------

            void
            fix_tape_to_node_connectivities( const id_t aMaxBlockID );

//------------------------------------------------------------------------------

            void
            shift_nodes();

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        inline const Vector< id_t > &
        TapeRoller::ghost_sideset_ids() const
        {
            return mGhostSideSetIDs ;
        }

//------------------------------------------------------------------------------

    }
}
#endif //BELFEM_CL_MESH_TAPEROLLER_HPP
