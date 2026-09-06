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

#include "assert.hpp"
#include "stringtools.hpp"
#include "cl_CutData.hpp"
#include "cl_Element_Factory.hpp"
#include "en_SolverEnums.hpp"


namespace belfem
{
    namespace mesh
    {
//-----------------------------------------------------------------------------

        CutData::CutData(
                  Mesh           * aMesh,
                  Cohomology     * aCohomology,
            const index_t          aCohomologyIndex ) :
                mMesh( aMesh ),
                mCochain( aCohomology->get_Generators()( 1 )( aCohomologyIndex ) ),
                mIndex( aCohomologyIndex ),
                mID( aCohomologyIndex + 1 ),
                mNumberOfCuts( aCohomology->get_Generators()( 1 ).size() )
        {

        }

//-----------------------------------------------------------------------------

        CutData::~CutData()
        {
            if ( mCohomologyPlus != nullptr )
            {
                delete mCohomologyPlus ;
            }
            if ( mCohomologyMinus != nullptr )
            {
                delete mCohomologyMinus ;
            }
        }

//-----------------------------------------------------------------------------

        void
        CutData::add_thin_cut_sidesets_to_mesh( id_t & aMaxSideSetID, id_t & aMaxElementID )
        {
            ElementFactory tFactory ;


            string tFormat = "cut_" + format_with_leading_zeros( mNumberOfCuts );

            SideSet * tSideSet =  new SideSet( ++aMaxSideSetID, 0 ) ;

            tSideSet->label() = sprint( tFormat.c_str(), mID );


            Cell< Facet * > & tFacets = tSideSet->facets() ;

            if( mMesh->number_of_dimensions() == 2 )
            {
                tFacets.set_size( mThinCutEdges.size(), nullptr );

                mMesh->unflag_all_elements() ;
                for ( Element * tElement : mCutElements )
                {
                    tElement->flag();
                }

                index_t tSide = gNoIndex ;


                DynamicBitset tBitset(3);

                index_t tCount = 0 ;

                for( Edge * tEdge : mThinCutEdges )
                {

                    BELFEM_ASSERT( tEdge->number_of_elements() == 2, "Edge needs two elements" );
                    BELFEM_ASSERT( tEdge->element( 0 )->is_flagged() xor tEdge->element( 1 )->is_flagged(), "Invalid cut pattern");

                    Element * tMaster = tEdge->element( 0 )->is_flagged() ? tEdge->element( 0 ) : tEdge->element( 1 ) ;

                    // identify side, note: something is off with the original pointers

                    tBitset.reset();

                    for ( uint i=0; i<2; ++i )
                    {
                        for ( uint k=0; k<tMaster->number_of_corner_nodes(); ++k )
                        {

                            real dx = tMaster->node( k )->x() - tEdge->node( i )->x() ;
                            real dy= tMaster->node( k )->y() - tEdge->node( i )->y() ;
                            if ( std::sqrt( dx*dx + dy*dy ) < BELFEM_EPSILON )
                            {
                                tBitset.set( k );
                                break;
                            }
                        }
                    }

                    switch ( tBitset.to_int() )
                    {
                        case( 3 ) :
                        {
                            tSide = 0 ;
                            break ;
                        }
                        case( 5 ) :
                        {
                            tSide = 2 ;
                            break ;
                        }
                        case( 6 ) :
                        {
                            tSide = 1 ;
                            break;
                        }
                        default:
                        {
                            BELFEM_ERROR( false, "Invalid side : %s", tBitset.to_string().c_str() );
                        }
                    }

                    ElementType tType = element_type_from_numnodes( 1, tEdge->number_of_nodes() );

                    Facet * tFacet = new Facet( tFactory.create_element( tType, ++aMaxElementID ) );
                    tFacet->set_master( tMaster, tSide );
                    tFacets( tCount++ ) = tFacet ;
                }
            }
            else
            {
                ElementType tType =  element_type_from_numnodes( mMesh->number_of_dimensions() - 1, mThinCutFaces( 0 )->number_of_nodes() );

                // allocate memory
                index_t tCount = 0 ;
                tFacets.set_size( mThinCutFaces.size(), nullptr );
                for ( Face * tFace : mThinCutFaces )
                {
                    Facet * tFacet = new Facet( tFactory.create_element( tType, ++aMaxElementID ) );

                    if ( tFace->master() != nullptr )
                    {
                        tFacet->set_master( tFace->master(), tFace->index_on_master(), true );
                        tFacet->set_slave( tFace->slave(), tFace->index_on_slave(), tFace->orientation_on_slave() );
                    }
                    else
                    {
                        // Step 5a: a slave-only seam face ( target periodic plane;
                        // fix_face_slaves cleared its master ) would crash the
                        // unconditional master() deref above. Emit a one-sided
                        // facet from the slave element instead, so the thin-cut
                        // sideset still covers it ( one facet per mThinCutFaces
                        // entry, keeping the preallocated tFacets count exact ).
                        BELFEM_ERROR( tFace->slave() != nullptr,
                            "Thin-cut face %lu has neither master nor slave element",
                            ( long unsigned int ) tFace->id() );

                        tFacet->set_master( tFace->slave(), tFace->index_on_slave(), true );
                    }

                    tFacets( tCount++ ) = tFacet ;
                }


            }

            mMesh->add_sideset( tSideSet ) ;
        }

//-----------------------------------------------------------------------------

        void
        CutData::collect_edges()
        {
            // grab the edge map
            OrderedMap< index_t, int > & tMap   = mCochain->getSimplicesMap() ;
            Cell< Edge * >   & tEdges = mMesh->edges() ;

            // count edges
            // note: in theory, all zero coefficients should be gone here,
            //       we do this test just as a precaution
            index_t tCount = 0 ;
            for ( auto & tPair : tMap )
            {
                if ( tPair.second != 0 )
                {
                    tEdges( tPair.first )->flag() ;
                    ++tCount ;
                    if ( tEdges( tPair.first )->is_periodic() )
                    {
                        tEdges( tPair.first )->periodic()->flag() ;
                        ++tCount ;
                    }
                }
            }

            mCohomologyEdges.set_size( tCount, nullptr ) ;
            tCount = 0 ;
            for ( auto & tPair : tMap )
            {
                if ( tPair.second != 0 )
                {
                    mCohomologyEdges( tCount++ ) = tEdges( tPair.first ) ;
                    if ( tEdges( tPair.first )->is_periodic() )
                    {
                        mCohomologyEdges( tCount++ ) = tEdges( tPair.first )->periodic() ;
                    }
                }
            }

            for ( Edge * tEdge : mCohomologyEdges )
            {
                tEdge->unflag() ;
            }
        }

        void
        CutData::collect_coefficients( const index_t aNumEdges )
        {
            mCohomologyPlus  = new DynamicBitset( aNumEdges );
            mCohomologyMinus = new DynamicBitset( aNumEdges );

            Cell< Edge * >   & tEdges = mMesh->edges() ;

            OrderedMap< index_t, int > & tMap   = mCochain->getSimplicesMap() ;
            for ( auto & tPair : tMap )
            {
                if ( tPair.second == 1 )
                {
                    mCohomologyPlus->set( tEdges( tPair.first )->index() );
                }
                else if ( tPair.second == -1 )
                {
                    mCohomologyMinus->set( tEdges( tPair.first )->index() );
                }
                if ( tEdges( tPair.first )->is_periodic() )
                {
                    if ( tPair.second == 1 )
                    {
                        mCohomologyPlus->set( tEdges( tPair.first )->periodic()->index() );
                    }
                    else if ( tPair.second == -1 )
                    {
                        mCohomologyMinus->set( tEdges( tPair.first )->periodic()->index() );
                    }
                }
            }
        }

        void
        CutData::collect_elements( const Vector< id_t > & aNonPhiDomains )
        {
            BELFEM_ASSERT( mCohomologyEdges.size() > 0, "No edges found for cohomology %u .", ( unsigned int ) mID );

            for ( Edge * tEdge : mCohomologyEdges )
            {
                tEdge->flag() ;
            }

            // flag candidates
            for ( Edge * tEdge : mCohomologyEdges )
            {
                for ( uint e=0; e<tEdge->number_of_elements(); ++e )
                {
                    tEdge->element( e )->flag() ;
                }
            }

            uint tDim = mMesh->number_of_dimensions() ;

            // count candidates
            for ( Element * tElement : mMesh->elements() )
            {
                if ( tElement->is_flagged() )
                {
                    // check number of flagged edges
                    uint tEdgeCount = 0 ;

                    tElement->unflag() ;

                    for ( uint e=0; e<tElement->number_of_edges(); ++e )
                    {
                        if ( tElement->edge( e )->is_flagged() )
                        {
                            ++tEdgeCount ;
                        }
                        if ( tEdgeCount >= tDim )
                        {
                            tElement->flag() ;
                            break ;
                        }
                    }
                }
            }

            // unflag elements on forbidden domains
            for ( id_t tID : aNonPhiDomains )
            {
                mMesh->block( tID )->unflag_elements() ;
            }

            // count elements
            index_t tCount = 0 ;
            for ( Element * tElement : mMesh->elements() )
            {
                if ( tElement->is_flagged() )
                {
                    tElement->set_index( tCount++ ) ;
                }
            }

            // sanity check
            BELFEM_ERROR( tCount > 0, "No elements found for cohomology %u .",
             ( unsigned int ) mID );

            // collect elements
            mCutElements.set_size( tCount, nullptr ) ;
            mCutCases.set_size( tCount, 0 ) ;
            tCount = 0 ;
            if ( mMesh->number_of_dimensions() == 2 )
            {
                for ( Element * tElement : mMesh->elements() )
                {
                    if ( tElement->is_flagged() )
                    {
                        tElement->unflag() ;
                        mCutElements( tCount ) = tElement ;
                        mCutCases( tCount++ )  = this->determine_cut_case_2d( tElement );
                    }
                }
            }
            else
            {
                for ( Element * tElement : mMesh->elements() )
                {
                    if ( tElement->is_flagged() )
                    {
                        tElement->unflag() ;
                        mCutElements( tCount ) = tElement ;
                        mCutCases( tCount++ )  = this->determine_cut_case_3d( tElement );
                    }
                }
            }

            for ( Edge * tEdge : mCohomologyEdges )
            {
                tEdge->unflag() ;
            }
        }

        int
        CutData::determine_cut_case_2d( Element * aElement )
        {

            uint tPattern   = 0 ;

            if( aElement->edge( 0 )->is_flagged() )
            {
                tPattern += 1 ;
            }
            if( aElement->edge( 1 )->is_flagged() )
            {
                tPattern += 2 ;
            }
            if( aElement->edge( 2 )->is_flagged() )
            {
                tPattern += 4 ;
            }

            uint tEdgeIn  = BELFEM_UINT_MAX ;
            uint tEdgeOut = BELFEM_UINT_MAX ;
            int aCase = 0 ;

            switch( tPattern )
            {
                case( 3 ) :
                {
                    // second corner / third edge
                    tEdgeIn  = 0 ;
                    tEdgeOut = 1 ;
                    aCase    = 3 ;
                    break ;
                }
                case( 5 ) :
                {
                    // first corner / second edge
                    tEdgeIn  = 2 ;
                    tEdgeOut = 0 ;
                    aCase    = 2 ;
                    break ;
                }
                case( 6 ) :
                {
                    // third corner / first edge
                    tEdgeIn  = 1 ;
                    tEdgeOut = 2 ;
                    aCase    = 1 ;
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid cut orientation in Element (%u)",
                                   ( long unsigned int ) aElement->id() );
                }
            } // end switch

            int tCoeffIn  = this->weight( aElement->edge( tEdgeIn ) );
            if( ! aElement->edge_direction( tEdgeIn ) )
            {
                tCoeffIn *= -1 ;
            }
            int tCoeffOut = this->weight( aElement->edge( tEdgeOut ) );
            if( ! aElement->edge_direction( tEdgeOut ) )
            {
                tCoeffOut *= -1 ;
            }

            BELFEM_ASSERT( tCoeffIn + tCoeffOut == 0, "Invalid cut coefficients in Element (%u)",
                           ( long unsigned int ) aElement->id() );

            if( tCoeffIn == -1 && tCoeffOut == 1 )
            {
                aCase *= -1 ;
            }

            return aCase ;
        }

        int
        CutData::determine_cut_case_3d( Element * aElement )
        {
            uint tPattern = 0 ;

            // determine value based on flagged pattern
            // this segment is deliberately unrolled
            if( aElement->edge( 0 )->is_flagged() )
            {
                tPattern += 1 ;
            }
            if( aElement->edge( 1 )->is_flagged() )
            {
                tPattern += 2 ;
            }
            if( aElement->edge( 2 )->is_flagged() )
            {
                tPattern += 4 ;
            }
            if( aElement->edge( 3 )->is_flagged() )
            {
                tPattern += 8 ;
            }
            if( aElement->edge( 4 )->is_flagged() )
            {
                tPattern += 16 ;
            }
            if( aElement->edge( 5 )->is_flagged() )
            {
                tPattern += 32 ;
            }

            uint tNumEdges = 0 ;
            int aCase = 0 ;
            uint tIndex[4];
            int  tValue[4];

            // with the flag pattern, we can determine the cut case
            switch( tPattern )
            {
                case( 13 ) :
                {
                    // second surface
                    aCase  = 2 ;
                    tNumEdges = 3 ;

                    tIndex[ 0 ] =  0 ;
                    tIndex[ 1 ] =  2 ;
                    tIndex[ 2 ] =  3 ;

                    tValue[ 0 ] =  1 ;
                    tValue[ 1 ] = -1 ;
                    tValue[ 2 ] =  1 ;

                    break ;
                }
                case( 19 ) :
                {
                    // third surface
                    aCase  = 3 ;
                    tNumEdges = 3 ;

                    tIndex[ 0 ] =  0 ;
                    tIndex[ 1 ] =  4 ;
                    tIndex[ 2 ] =  1 ;

                    tValue[ 0 ] = -1 ;
                    tValue[ 1 ] =  1 ;
                    tValue[ 2 ] =  1 ;

                    break ;
                }
                case( 38 ) :
                {
                    // first surface
                    aCase = 1 ;
                    tNumEdges = 3 ;

                    tIndex[ 0 ] =  2 ;
                    tIndex[ 1 ] =  1 ;
                    tIndex[ 2 ] =  5 ;

                    tValue[ 0 ] =  1 ;
                    tValue[ 1 ] = -1 ;
                    tValue[ 2 ] =  1 ;

                    break ;
                }
                case( 56 ) :
                {
                    // fourth surface
                    aCase  = 4 ;
                    tNumEdges = 3 ;

                    tIndex[ 0 ] =  3 ;
                    tIndex[ 1 ] =  5 ;
                    tIndex[ 2 ] =  4 ;

                    tValue[ 0 ] =  -1 ;
                    tValue[ 1 ] =  -1 ;
                    tValue[ 2 ] =  -1 ;

                    break ;
                }
                case( 30 ) :
                {
                    // second diagonal: edge 1 -> edge 3
                    aCase = 6;
                    tNumEdges = 4;

                    tIndex[ 0 ] =  1;
                    tIndex[ 1 ] =  2;
                    tIndex[ 2 ] =  3;
                    tIndex[ 3 ] =  4;

                    tValue[ 0 ] =  1;
                    tValue[ 1 ] = -1;
                    tValue[ 2 ] =  1;
                    tValue[ 3 ] =  1;

                    break;
                }
                case( 53 ) :
                {
                    // third diagonal: edge 2 -> edge 4
                    aCase = 7 ;
                    tNumEdges = 4 ;

                    tIndex[ 0 ] =  2;
                    tIndex[ 1 ] =  0;
                    tIndex[ 2 ] =  4;
                    tIndex[ 3 ] =  5;

                    tValue[ 0 ] =  1;
                    tValue[ 1 ] = -1;
                    tValue[ 2 ] =  1;
                    tValue[ 3 ] =  1;

                    break ;
                }
                case( 43 ) :
                {
                    // first diagonal: edge 0 -> edge 5
                    aCase = 5 ;
                    tNumEdges = 4 ;

                    tIndex[ 0 ] =  0;
                    tIndex[ 1 ] =  1;
                    tIndex[ 2 ] =  5;
                    tIndex[ 3 ] =  3;

                    tValue[ 0 ] =  1;
                    tValue[ 1 ] = -1;
                    tValue[ 2 ] =  1;
                    tValue[ 3 ] =  1;

                    break ;
                }
                default :
                {
                    // an inadmissible flagged-edge pattern points at a generator or
                    // closedness defect upstream; a too-coarse mesh is caught earlier,
                    // in clean_spfa(), with an edge-ID certificate
                    BELFEM_ERROR( false,
                        "Inadmissible cut edge pattern ( code %u ) on Element %lu: the cleaned generator should not produce this. This indicates a cohomology generator or closedness defect, not mesh coarseness.",
                           ( unsigned int ) tPattern,
                           ( long unsigned int ) aElement->id() );
                }
            }

            // apply coefficient and adjust orientation
            for( uint e=0; e<tNumEdges; ++e )
            {
                if( aElement->edge_direction( tIndex[ e ] ) )
                {
                    tValue[ e ] *=  this->weight( aElement->edge( tIndex[ e ] ) );
                }
                else
                {
                    tValue[ e ] *= -this->weight( aElement->edge( tIndex[ e ] ) );
                }
            }

            if( aCase <= 4 )
            {
                if ( tValue[ 0 ] == 1 && tValue[ 1 ] == 1 && tValue[ 2 ] == 1 )
                {
                    return aCase ;
                }
                else if ( tValue[ 0 ] == -1 && tValue[ 1 ] == -1 && tValue[ 2 ] == -1 )
                {
                    return - aCase ;
                }
                else
                {
                    // unit but not sign-coherent ( or zero: edge flagged without a
                    // coefficient ) -> generator orientation defect, not coarseness
                    BELFEM_ERROR( false,
                        "Cut coefficients on Element %lu are not sign-coherent ( %d %d %d ). This indicates a generator orientation defect, not mesh coarseness.",
                        ( long unsigned int ) aElement->id(),
                        tValue[ 0 ], tValue[ 1 ], tValue[ 2 ] );
                }
            }
            else
            {
                if ( tValue[ 0 ] == 1 && tValue[ 1 ] == 1 && tValue[ 2 ] == 1 && tValue[ 3 ] == 1 )
                {
                    return aCase ;
                }
                else if ( tValue[ 0 ] == -1 && tValue[ 1 ] == -1 && tValue[ 2 ] == -1  && tValue[ 3 ] == -1 )
                {
                    return - aCase ;

                }
                else
                {
                    BELFEM_ERROR( false,
                        "Cut coefficients on Element %lu are not sign-coherent ( %d %d %d %d ). This indicates a generator orientation defect, not mesh coarseness.",
                        ( long unsigned int ) aElement->id(),
                        tValue[ 0 ], tValue[ 1 ], tValue[ 2 ], tValue[ 3 ] );
                }
            }
            return 0 ;
        }

//-----------------------------------------------------------------------------

        void
        CutData::flag_nodes()
        {
            if ( mMesh->number_of_dimensions() == 2 )
            {
                if ( mMesh->has_periodicity() )
                {
                    for ( Edge * tEdge : mThinCutEdges )
                    {
                        for ( uint e=0; e<tEdge->number_of_nodes(); ++e )
                        {
                            Node * tNode = tEdge->node( e ) ;
                            tNode->flag() ;
                            if ( tNode->is_periodic() )
                            {
                                tNode->periodic()->flag() ;
                            }
                        }
                    }
                }
                else
                {
                    for ( Edge * tEdge : mThinCutEdges )
                    {
                        tEdge->flag_nodes() ;
                    }
                }
            }
            else
            {
                if ( mMesh->has_periodicity() )
                {
                    for  ( Face * tFace : mThinCutFaces )
                    {
                        for ( uint e=0; e<tFace->number_of_nodes(); ++e )
                        {
                            Node * tNode = tFace->node( e ) ;
                            tNode->flag() ;
                            if ( tNode->is_periodic() )
                            {
                                tNode->periodic()->flag() ;
                            }
                        }
                    }
                }
                else
                {
                    for  ( Face * tFace : mThinCutFaces )
                    {
                        tFace->flag_nodes() ;
                    }
                }
            }
        }

//-----------------------------------------------------------------------------

        Mesh *
        CutData::create_debug_mesh()
        {
            mMesh->unflag_all_nodes() ;
            for ( Edge * tEdge : mCohomologyEdges )
            {
                tEdge->flag_nodes() ;
            }

            if ( mMesh->number_of_dimensions() == 2 )
            {
                for ( Edge * tEdge : mThinCutEdges )
                {
                    tEdge->flag_nodes() ;
                }
            }
            else
            {
                for ( Face * tFace : mThinCutFaces )
                {
                    tFace->flag_nodes() ;
                }
            }

            // count nodes
            index_t tCount = 0 ;
            for ( Node * tNode : mMesh->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    tNode->set_index( tCount++ ) ;
                }
            }

            Mesh * aMesh = new Mesh( mMesh->number_of_dimensions() ) ;

            Cell< Node * > & tNodes = aMesh->nodes() ;
            tNodes.set_size( tCount, nullptr ) ;
            tCount = 0 ;

            for ( Node * tNode : mMesh->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    tNodes( tCount++ ) = new Node( tNode->id(), tNode->x(), tNode->y(), tNode->z() ) ;
                }
            }

            ElementFactory tFactory ;

            Block * tHomology = new Block( mID, mCohomologyEdges.size() ) ;

            for ( Edge * tEdge : mCohomologyEdges )
            {
                Element * tElement = tFactory.create_element( ElementType::LINE2, tEdge->id() );
                tElement->insert_node( tNodes( tEdge->node( 0 )->index() ) , 0 );
                tElement->insert_node( tNodes( tEdge->node( 1 )->index() ) , 1 );
                tHomology->insert_element( tElement );
            }

            aMesh->add_block( tHomology );


            if ( mMesh->number_of_dimensions() == 2 )
            {
                if ( mThinCutEdges.size() > 0 )
                {
                    Block * tCut = new Block( 2* mID, mThinCutEdges.size() ) ;


                    for ( Edge * tEdge : mThinCutEdges )
                    {
                        Element * tElement = tFactory.create_element( ElementType::LINE2, tEdge->id() );
                        tElement->insert_node( tNodes( tEdge->node( 0 )->index() ) , 0 );
                        tElement->insert_node( tNodes( tEdge->node( 1 )->index() ) , 1 );

                        tCut->insert_element( tElement );
                    }

                    aMesh->add_block( tCut ) ;
                }
            }
            else
            {
                if ( mThinCutFaces.size() > 0 )
                {
                    Block * tCut = new Block( 2* mID, mThinCutFaces.size() ) ;

                    ElementType tType =  element_type_from_numnodes( mMesh->number_of_dimensions() - 1, mThinCutFaces( 0 )->number_of_nodes() );

                    for ( Face * tFace : mThinCutFaces )
                    {

                        Element * tElement = tFactory.create_element( tType, tFace->id() );
                        for ( uint k=0; k<tFace->number_of_nodes(); ++k )
                        {
                            tElement->insert_node( tNodes( tFace->node( k )->index() ) , k );
                        }
                        tCut->insert_element( tElement );
                    }

                    aMesh->add_block( tCut ) ;
                }
            }

            aMesh->finalize();
            return aMesh ;
        }

//-----------------------------------------------------------------------------

        void
        CutData::unflag_all_elements()
        {
            for ( Element * tElement : mCutElements )
            {
                tElement->unflag() ;
            }
        }

//-----------------------------------------------------------------------------

        void
        CutData::flag_all_elements()
        {
            for ( Element * tElement : mCutElements )
            {
                tElement->flag() ;
            }
        }

    }
}
