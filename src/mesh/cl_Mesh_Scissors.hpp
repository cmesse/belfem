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
                const bool mIsTape ;

                // vector contains sideset id, plus side and minus side
                Cell< Vector< id_t > > mData ;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                CutData( const id_t aCutID,
                         const id_t aSideSetID,
                         const id_t aPlusBlock,
                         const id_t aMinusBlock,
                         const bool aIsTape );

                CutData( const id_t aCutID,
                         const Vector< id_t > & aSideSetIDs,
                         const id_t aPlusBlock,
                         const id_t aMinusBlock,
                         const bool aIsTape );

                CutData( const id_t aCutID,
                         const Vector< id_t > & aSideSetIDs,
                         const Vector< id_t > & aPlusBlocks,
                         const Vector< id_t > & aMinusBlocks,
                         const bool aIsTape );

//------------------------------------------------------------------------------

                ~CutData() = default ;

//------------------------------------------------------------------------------

                id_t
                id() const;

//------------------------------------------------------------------------------

                bool
                is_tape() const;

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

            inline bool
            CutData::is_tape() const
            {
                return mIsTape ;
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

            Cell< scissors::CutData * > mCutsAndTapes ;

            // flattened dataset
            Matrix< id_t > mTapeData ;
            Matrix< id_t > mCutData ;

            index_t mNumberOfNodes = 0 ;
            index_t mNumberOfFacets = 0 ;


            // list of connector sidesets
            Vector< id_t > mConnectorSetIDs ;
            Map< id_t, id_t > mSideSetToCutIDs ;


            // check if a block is minus
            Map< id_t, index_t > mMinus ;

            // check if a block is plus
            Map< id_t, index_t > mPlus ;

            // new data from here
            //const Matrix< uint > mCase =  { { 0, 1, 2 }, { 3, 4, 5 }, { 6, 7, 8 } } ;
            //const Vector< uint > mLookup = { 0, 1, 1, 1, 1, 1, 1, 2, 3 };

            Cell< SideSet * > mCuts ;

            // - - - - - -  - - - - - - - - - - - - - - - - - - - -
            // new from here
            // - - - - - -  - - - - - - - - - - - - - - - - - - - -

            //! all used sidesets on the original mesh
            Cell< SideSet * > mAllSideSets ;

            //! sidesets on the original mesh that become tapes
            Cell< SideSet * > mTapeSideSets ;

            //! sidesets on the original mesh that become cuts
            Cell< SideSet * > mCutSideSets ;

            //! lookup table for number of nodes to be crated
            //! rows: num tapes, cols: num cuts
            const Matrix< uint > mLookup =  { { 0, 1, 1 }, { 1, 1, 1 }, { 1, 2, 3 } } ;

            //! number of tape sidesets per node
            Vector< uint > mNumTapesPerNode ;

            //! number of cut sidesets per node
            Vector< uint > mNumCutsPerNode ;

            //! multi purpose counting container
            Vector< uint > mNodeCounter ;

            //! list of plus blocks per node
            Cell< Vector< id_t > > mPlusBlocksTable ;
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Scissors( Mesh * aMesh, const Vector< id_t > & aAirBlockIDs );

//----------------------------------------------------------------------------

            ~Scissors();

//------------------------------------------------------------------------------

            void
            cut( const id_t aSidesetID,
                 const id_t aPlusBlock,
                 const id_t aMinusBlock,
                 const bool aIsTape=false );

//------------------------------------------------------------------------------

            void
            cut( const Vector< id_t > & aSidesetIDs,
                 const id_t aPlusBlock,
                 const id_t aMinusBlock,
                 const bool aIsTape=false );

//------------------------------------------------------------------------------

            void
            cut( const Vector< id_t > & aSidesetIDs,
                 const Vector< id_t > & aMinusBlocks,
                 const Vector< id_t > & aPlusBlocks,
                 const bool aIsTape=false );

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

            // help function to check input sanity
            void
            check_compatibility( const bool aTapeMode );

//------------------------------------------------------------------------------
            // NEW FROM HERE
//------------------------------------------------------------------------------

            void
            collect_sidesets();

//------------------------------------------------------------------------------

            void
            collect_original_nodes();

//------------------------------------------------------------------------------

            void
            count_duplicate_nodes();

//------------------------------------------------------------------------------

            void
            create_node_duplicates();

//------------------------------------------------------------------------------

            void
            relink_elements();

//------------------------------------------------------------------------------

            void
            create_connectors();

//------------------------------------------------------------------------------

            void
            unflag_nodes();

//------------------------------------------------------------------------------

            void
            flag_nodes();

//------------------------------------------------------------------------------

            void
            add_nodes_to_mesh();

//------------------------------------------------------------------------------

            void
            rearrange_nodes( Cell< Node * > & aMaster,
                             Cell< Node * > & aSlaveIn,
                             Cell< Node * > & aSlaveOut,
                             Vector< uint > & aFlags );


//------------------------------------------------------------------------------

            void
            relink_facets();

//------------------------------------------------------------------------------

            uint
            inside_node( const ElementType aElementType );

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
