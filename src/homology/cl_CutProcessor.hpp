/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef CL_COHOMOLOGYPROCESSOR_HPP
#define CL_COHOMOLOGYPROCESSOR_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Map.hpp"
#include "cl_OrderedMap.hpp"
#include "cl_DynamicBitset.hpp"
#include "cl_Mesh.hpp"
#include "cl_Cohomology.hpp"
#include "cl_CutData.hpp"
#include "cl_CutSet.hpp"

namespace belfem
{
    namespace mesh
    {

        /**
         * @brief Converts thick cuts into the thin cuts the FEM assembly consumes.
         *
         * @ingroup grp_homology
         * @see @ref homology_thick_thin_cuts_and_conjugate_edges
         */
        class CutProcessor
        {
            Mesh       * mMesh ;
            bool         mIs2D ;

            const uint   mNumberOfDimensions ;
            const uint   mNumberOfCuts ;
            const ElementType mElementType ;

            Vector< id_t > mPhiBlocks ;
            Vector< id_t > mNonPhiBlocks ;

            //! list of all elements that contribute to the cohomology
            Cell< Element * > mElements ;

            //! list of all edges that contribute to the cohomology
            Cell< Edge * > mCohomologyEdges ;

            Cell< DynamicBitset * > mCohomologyEdgesPlus ;
            Cell< DynamicBitset * > mCohomologyEdgesMinus ;

            Vector< id_t > mPhiBoundaries ;
            Vector< id_t > mPhiBoundariesAndPeriodic ;

            //! list of all original nodes that contribute to the cohomology
            Cell< Node * > mNodes ;

            //! container that translates cohomology data
            Cell< CutData * > mCutData ;

            //! keyed by cut-pattern hex; OrderedMap ( std::map ) so iteration is
            //! deterministic by key — duplicate node IDs are then reproducible
            //! across runs ( aligns with the Step 3 deterministic-pairing design )
            OrderedMap< string, CutSet * > mCutSets ;

            //! string with hex patterns
            Cell< string > mPatterns ;

            //! maximum ID for elements
            id_t mMaxElementID ;

            //! maximum ID for nodes
            id_t mMaxNodeID ;

            //! maximum sideset ID
            id_t mMaxSidesetID ;

            //! abstract nodes for current dofs
            Cell< Node * > mAbstractNodes ;

            // for more efficient bitset identification
            Map< size_t, Cell< DynamicBitset* > > mHashToBitsets;

//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            CutProcessor(
                    Mesh * aMesh,
                    Cohomology * aCohomology,
                    const Vector< id_t > & aPhiBlocks,
                    const Vector< id_t > & aNonPhiBlocks,
                    const Vector< id_t > & aPhiInterfaces,
                    const Vector< id_t > & aPhiBoundaries,
                    const Vector< id_t > & aPhiPeriodic );

            ~CutProcessor();


            Cell< Node * > &
            abstract_nodes();

            id_t
            max_node_id();

            id_t
            max_sideset_id();

            void
            save_debug_meshes();

//-----------------------------------------------------------------------------
        private:
//-----------------------------------------------------------------------------

            void
            create_thin_cut_sidesets();

//-----------------------------------------------------------------------------

            void
            collect_elements();

//-----------------------------------------------------------------------------
            
            void
            collect_edges();

//-----------------------------------------------------------------------------

            void
            collect_facets();

//-----------------------------------------------------------------------------

            void
            compute_edge_bitsets();

//-----------------------------------------------------------------------------

            void
            compute_node_bitsets();

//-----------------------------------------------------------------------------

            void
            collect_nodes();

//-----------------------------------------------------------------------------

            void
            determine_cut_sets();

//-----------------------------------------------------------------------------

            index_t
            pattern_index( const DynamicBitset * aBitset );

//-----------------------------------------------------------------------------

            void
            create_abstract_nodes();

//-----------------------------------------------------------------------------

            void
            create_cut_sets();

//-----------------------------------------------------------------------------

            void
            flip_node_bitsets_tri3(
                    Element * aElement,
                    Cell< DynamicBitset * > & aBitsets );

            void
            check_edges_tri3(
                    Element * aElement,
                    Cell< DynamicBitset * > & aBitsets );

//-----------------------------------------------------------------------------

            void
            flip_node_bitsets_tri6(
                    Element * aElement,
                    Cell< DynamicBitset * > & aBitsets );

            void
            check_edges_tri6(
                    Element * aElement,
                    Cell< DynamicBitset * > & aBitsets );

//-----------------------------------------------------------------------------

            void
            flip_node_bitsets_tet4(
                    Element * aElement,
                    Cell< DynamicBitset * > & aBitsets );

            void
            check_edges_tet4(
                    Element * aElement,
                    Cell< DynamicBitset * > & aBitsets );

//-----------------------------------------------------------------------------

            void
            flip_node_bitsets_tet10(
                    Element * aElement,
                    Cell< DynamicBitset * > & aBitsets );

            void
            check_edges_tet10(
                    Element * aElement,
                    Cell< DynamicBitset * > & aBitsets );

            void
            check_node_bitsets( Cell< DynamicBitset * > & aBitsets );

//-----------------------------------------------------------------------------

            void
            check_edge(
                    Element * aElement,
                    const uint aEdgeIndex,
                    const uint aMinusNodeIndex,
                    const uint aPlusNodeIndex,
                    Cell< DynamicBitset * > & aBitsets );

//-----------------------------------------------------------------------------

            void
            check_midside(
                    const uint aMinusNodeIndex,
                    const uint aPlusNodeIndex,
                    const uint aMidsideNodeIndex,
                    Cell< DynamicBitset * > & aBitsets );

//-----------------------------------------------------------------------------

            void
            duplicate_nodes();

//-----------------------------------------------------------------------------

            //void
            //flip_node_bitsets();

//-----------------------------------------------------------------------------

            void
            relink_elements();

//-----------------------------------------------------------------------------

            void
            relink_element( Element * aElement,
                            Cell< DynamicBitset * > & aBitsets );

//-----------------------------------------------------------------------------

            void
            collect_duplicates();

//-----------------------------------------------------------------------------
        };

        inline Cell< Node * > &
        CutProcessor::abstract_nodes()
        {
            return mAbstractNodes;
        }

//-----------------------------------------------------------------------------

        inline id_t
        CutProcessor::max_node_id()
        {
            return mMaxNodeID;
        }

//-----------------------------------------------------------------------------

        inline id_t
        CutProcessor::max_sideset_id()
        {
            return mMaxSidesetID;
        }

//-----------------------------------------------------------------------------


        inline void
        CutProcessor::check_edge(
                Element   * aElement,
                const uint  aEdgeIndex,
                const uint  aMinusNodeIndex,
                const uint  aPlusNodeIndex,
                Cell< DynamicBitset * > & aBitsets )
        {


            if ( ! aElement->edge( aEdgeIndex )->is_flagged() ) return ;

            for ( uint c=0; c<mNumberOfCuts; ++c )
            {
                int tValue = mCutData( c )->weight( aElement->edge( aEdgeIndex ) ) ;

                if ( aElement->edge_direction( aEdgeIndex ) == mIs2D ) tValue *= -1 ;

                if ( tValue == 1 )
                {
                    if ( aElement->node( aPlusNodeIndex )->is_flagged() ) aBitsets( aPlusNodeIndex )->set( c ) ;
                }
                else if ( tValue == -1 )
                {
                    if ( aElement->node( aMinusNodeIndex )->is_flagged() ) aBitsets( aMinusNodeIndex )->set( c ) ;
                }
            }
        }

//-----------------------------------------------------------------------------

        inline void
        CutProcessor::check_midside(
                const uint aMinusNodeIndex,
                const uint aPlusNodeIndex,
                const uint aMidsideNodeIndex,
                Cell< DynamicBitset * > & aBitsets )
        {
            if (  aBitsets( aPlusNodeIndex )->count() == 0 ) return ;
            if (  aBitsets( aMinusNodeIndex )->count() == 0 ) return ;

            for ( uint c=0; c<mNumberOfCuts; ++c )
            {
                if ( aBitsets( aPlusNodeIndex )->test( c ) and aBitsets( aMinusNodeIndex )->test( c ) )
                {
                    aBitsets( aMidsideNodeIndex )->set( c );
                }
            }
        }

//-----------------------------------------------------------------------------


    }
}
#endif //CL_COHOMOLOGYPROCESSOR_HPP
