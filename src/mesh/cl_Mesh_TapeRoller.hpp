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

            Cell< id_t > mSelectedSideSets ;

            const uint mNumberOfLayers ;

            id_t mMaxNodeId ;
            id_t mMaxEdgeId ;
            id_t mMaxFaceId ;
            id_t mMaxElementID ;

            id_t mMaxSideSetId ;

            // ids with elements that are to be cloned
            Vector< index_t > mElementIDs ;

            Cell< Layer * > mLayers ;

            // map that links old node ids to position in layer
            Map< id_t, index_t > mNodeMap ;
            Map< id_t, index_t > mEdgeMap ;
            Map< id_t, index_t > mElementMap ;

            Vector< id_t > mGhostSideSetIDs ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            TapeRoller( Mesh * aMesh, const uint aNumberOfLayers );

//------------------------------------------------------------------------------

            ~TapeRoller();

//------------------------------------------------------------------------------

            void
            add_sideset( const id_t aSideSetID );

//------------------------------------------------------------------------------

            void
            add_sidesets( const Vector< id_t > & aSideSetID );

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
        private:
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

//-------------------------------------------------------------------------

            void
            clone_edges();

//-------------------------------------------------------------------------

            void
            clone_faces();

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        inline const Vector< id_t > &
        TapeRoller::ghost_sideset_ids() const
        {
            return mGhostSideSetIDs ;
        }

//---------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_MESH_TAPEROLLER_HPP
