//
// Created by christian on 9/15/21.
//

#ifndef BELFEM_CL_MESH_SCISSORS_HPP
#define BELFEM_CL_MESH_SCISSORS_HPP

#include "cl_Mesh.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Map.hpp"

namespace belfem
{
    namespace mesh
    {
        namespace scissors
        {
            class CutData
            {
                const id_t mID ;

                // vector contains sideset id, plus side and minus side
                Cell< Vector< id_t > > mData ;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                CutData( const id_t aCutID,
                         const id_t aSideSetID,
                         const id_t aPlusBlock,
                         const id_t aMinusBlock );

                CutData( const id_t aCutID,
                         const Vector< id_t > & aSideSetIDs,
                         const id_t aPlusBlock,
                         const id_t aMinusBlock );

                CutData( const id_t aCutID,
                         const Vector< id_t > & aSideSetIDs,
                         const Vector< id_t > & aPlusBlocks,
                         const Vector< id_t > & aMinusBlocks );

//------------------------------------------------------------------------------

                ~CutData() = default ;

//------------------------------------------------------------------------------

                id_t
                id() const;

//------------------------------------------------------------------------------

                Vector< id_t > &
                data( const index_t aSideSetCount );

//------------------------------------------------------------------------------

                index_t
                num_sidesets() const ;

//------------------------------------------------------------------------------

            };

//------------------------------------------------------------------------------

            inline id_t
            CutData::id() const
            {
                return mID ;
            }

//------------------------------------------------------------------------------

            inline index_t
            CutData::num_sidesets() const
            {
                return mData.size() ;
            }

//------------------------------------------------------------------------------

            inline Vector< id_t > &
            CutData::data( const index_t aSideSetCount )
            {
                return mData( aSideSetCount );
            }

//------------------------------------------------------------------------------
        }

        class Scissors
        {
            Mesh * mMesh ;
            const proc_t mMyRank ;

            const index_t mNumberOfOriginalNodes ;

            Cell< Node * >    & mNodes ;
            Cell< Facet * >   & mFacets ;
            Cell< Element * > & mElements ;

            Cell< Node * >      mOriginalNodes ;
            Cell< Node * >      mDuplicateNodes ;

            const Vector< id_t > mAirBlockIDs ;

            index_t mMaxNodeID = 0 ;
            index_t mMaxElementID = 0 ;
            index_t mMaxSideSetID = 0 ;

            index_t mNumberOfCuts = 0 ;
            index_t mNumberOfTapes = 0 ;

            Cell< scissors::CutData * > mCuts ;
            Cell< scissors::CutData * > mTapes ;

            // flattened dataset
            Matrix< id_t > mTapeData ;
            Matrix< id_t > mCutData ;

            index_t mNumberOfNodes = 0 ;
            index_t mNumberOfFacets = 0 ;

            // map that links minus node ids with plus nodes
            //Map< id_t, mesh::Node * > mNodeMapTapes ;
            //Map< id_t, mesh::Node * > mNodeMapCuts ;

            // map that links minus facet ids with plus facets
            //Map< id_t, mesh::Facet * > mFacetMap ;

            // list of connector sidesets
            Vector< id_t > mConnectorSetIDs ;
            Map< id_t, id_t > mSideSetToCutIDs ;


            // check if a block is minus
            Map< id_t, index_t > mMinus ;

            // check if a block is plus
            Map< id_t, index_t > mPlus ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Scissors( Mesh * aMesh, const Vector< id_t > & aAirBlockIDs );

            ~Scissors();

//------------------------------------------------------------------------------

            void
            cut( const id_t aSidesetID,
                 const id_t aPlusBlock,
                 const id_t aMinusBlock,
                 const bool aIsSheet=false );

//------------------------------------------------------------------------------

            void
            cut( const Vector< id_t > & aSidesetIDs,
                 const id_t aPlusBlock,
                 const id_t aMinusBlock,
                 const bool aIsSheet=false );

//------------------------------------------------------------------------------

            void
            cut( const Vector< id_t > & aSidesetIDs,
                 const Vector< id_t > & aMinusBlocks,
                 const Vector< id_t > & aPusBlocks,
                 const bool aIsSheet=false );

//------------------------------------------------------------------------------

            void
            finalize();

//------------------------------------------------------------------------------

            void
            get_cut_ids( const Vector< id_t > & aSideSetIDs, Vector< id_t > & aCutIDs );

//------------------------------------------------------------------------------

            const Matrix< id_t > &
            cut_data() const ;

//------------------------------------------------------------------------------

            void
            collect_cut_data( const bool aTapeMode = false );

//------------------------------------------------------------------------------
        private:

//------------------------------------------------------------------------------

            void
            compute_max_sideset_id();

//------------------------------------------------------------------------------

            void
            compute_max_node_id();

//------------------------------------------------------------------------------

            void
            compute_max_element_id();

//------------------------------------------------------------------------------

            void
            count_nodes( const bool aTapeMode  );

//------------------------------------------------------------------------------

            void
            count_facets( const bool aTapeMode  );

//------------------------------------------------------------------------------

            void
            duplicate_nodes( const bool aTapeMode );

//------------------------------------------------------------------------------

            void
            create_cut_table( const bool aTapeMode  );

//------------------------------------------------------------------------------

            void
            relink_facets();

//------------------------------------------------------------------------------

            void
            relink_elements( const bool aTapeMode );

//------------------------------------------------------------------------------

            void
            create_connectors();

//------------------------------------------------------------------------------

            // help function to check input sanity
            void
            check_compatibility( const bool aTapeMode );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline const Matrix< id_t > &
        Scissors::cut_data() const
        {
            return mCutData ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_MESH_SCISSORS_HPP
