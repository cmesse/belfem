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
#include "cl_CutProcessor.hpp"

namespace belfem
{
    namespace mesh
    {
//-----------------------------------------------------------------------------

        CutProcessor::CutProcessor(
                Mesh * aMesh,
                Cohomology * aCohomology,
                const Vector< id_t > & aPhiBlocks,
                const Vector< id_t > & aNonPhiBlocks,
                const Vector< id_t > & aPhiInterfaces,
                    const Vector< id_t > & aPhiBoundaries,
                    const Vector< id_t > & aPhiPeriodic ) :
                mMesh( aMesh ),
                mIs2D( aMesh->number_of_dimensions() == 2 ),
                mNumberOfDimensions( aMesh->number_of_dimensions() ),
                mNumberOfCuts( aCohomology->get_Generators()( 1 ).size() ),
                mElementType(  aMesh->max_element_order() == 1 ?
                               ( aMesh->number_of_dimensions() == 2 ? ElementType::TRI3 : ElementType::TET4 ) :
                               ( aMesh->number_of_dimensions() == 2 ? ElementType::TRI6 : ElementType::TET10 )  ),
                mPhiBlocks( aPhiBlocks ),
                mNonPhiBlocks( aNonPhiBlocks )
        {
            // get maximum IDs from mesh
            mMaxNodeID    = mMesh->max_node_id() ;


            mCutData.set_size( mNumberOfCuts, nullptr );

            for ( uint c=0; c<mCutData.size(); ++c )
            {
                mCutData( c ) = new CutData( mMesh, aCohomology, c );
            }
            this->collect_edges();

            this->collect_elements();

            mPhiBoundaries.set_size( aPhiInterfaces.length() + aPhiBoundaries.length() );
            mPhiBoundariesAndPeriodic.set_size( mPhiBoundaries.length() + aPhiPeriodic.length() );
            index_t tCount = 0 ;
            for ( id_t tID : aPhiInterfaces )
            {
                mPhiBoundaries( tCount ) = tID ;
                mPhiBoundariesAndPeriodic( tCount++ ) = tID ;
            }
            for ( id_t tID : aPhiBoundaries )
            {
                mPhiBoundaries( tCount ) = tID ;
                mPhiBoundariesAndPeriodic( tCount++ ) = tID ;
            }
            for ( id_t tID : aPhiPeriodic )
            {
                mPhiBoundariesAndPeriodic( tCount++ ) = tID ;
            }

            this->collect_facets();

            this->collect_nodes();

            this->determine_cut_sets();
            this->create_cut_sets();

            this->compute_edge_bitsets();

            this->compute_node_bitsets();

            this->create_abstract_nodes();

            this->duplicate_nodes();

            this->relink_elements();

            this->collect_duplicates() ;

            this->create_thin_cut_sidesets();

            mMesh->finalize();
        }

//-----------------------------------------------------------------------------

        CutProcessor::~CutProcessor()
        {
            for ( DynamicBitset * tBitset : mCohomologyEdgesPlus )
            {
                delete tBitset ;
            }
            for ( DynamicBitset * tBitset : mCohomologyEdgesMinus )
            {
                delete tBitset ;
            }

            for ( CutData * tCutData : mCutData )
            {
                delete tCutData ;
            }

            for ( auto tPair : mCutSets )
            {
                delete tPair.second ;
            }

            for ( auto tPair : mHashToBitsets )
            {
                Cell< DynamicBitset * > & tBitsets = tPair.second ;

                for ( DynamicBitset * tBitset : tBitsets )
                {
                    delete tBitset ;
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        CutProcessor::save_debug_meshes()
        {
            for ( CutData * tCutData : mCutData )
            {
                Mesh * tMesh = tCutData->create_debug_mesh();
                tMesh->save( "cut_" + std::to_string( tCutData->index() ) + ".vtk" );
                delete tMesh;
            }
        }

//-----------------------------------------------------------------------------

        void
        CutProcessor::collect_edges()
        {
            // select edges
            mMesh->unflag_all_edges();

            for ( CutData * tCutData : mCutData )
            {
                tCutData->collect_edges() ;
            }

            for ( CutData * tCutData : mCutData )
            {
                for ( Edge * tEdge : tCutData->cohomology_edges() )
                {
                    tEdge->flag() ;
                }
            }

            // count edges
            index_t tCount = 0 ;
            for ( Edge * tEdge : mMesh->edges() )
            {
                if ( tEdge->is_flagged() )
                {
                    tEdge->set_index( tCount++ ) ;
                }
            }

            // collect edges
            mCohomologyEdges.set_size( tCount, nullptr );
            tCount = 0 ;
            for ( Edge * tEdge : mMesh->edges() )
            {
                if ( tEdge->is_flagged() )
                {
                    mCohomologyEdges( tCount++ ) = tEdge ;
                }
            }

            // once we know the total size of the edges, we can collect the coefficients
            for ( CutData * tCutData : mCutData )
            {
                tCutData->collect_coefficients( tCount );
            }

            for ( Edge * tEdge : mCohomologyEdges )
            {
                tEdge->unflag();
            }
        }

//-----------------------------------------------------------------------------

        void
        CutProcessor::compute_edge_bitsets()
        {
            index_t tNumEdges = mCohomologyEdges.size() ;

            mCohomologyEdgesPlus.set_size( tNumEdges, nullptr );
            for ( index_t e=0; e<tNumEdges; ++e )
            {
                mCohomologyEdgesPlus( e ) = new DynamicBitset( mNumberOfCuts ) ;
            }

            mCohomologyEdgesMinus.set_size( tNumEdges, nullptr );
            for ( index_t e=0; e<tNumEdges; ++e )
            {
                mCohomologyEdgesMinus( e ) = new DynamicBitset( mNumberOfCuts ) ;
            }

            uint c = 0 ;
            for ( CutData * tCutData : mCutData )
            {
                for ( Edge * tEdge : tCutData->cohomology_edges() )
                {
                    int tVal = tCutData->weight( tEdge );

                    if ( tVal == 1 )
                    {
                        mCohomologyEdgesPlus( tEdge->index() )->set( c ) ;
                    }
                    else if ( tVal == -1 )
                    {
                        mCohomologyEdgesMinus( tEdge->index() )->set( c ) ;
                    }
                }
                ++c ;
            }
        }

//-----------------------------------------------------------------------------

        void
        CutProcessor::collect_elements()
        {
            // comparison object for sorting
            struct
            {
                inline bool
                operator()( const Element * aA, const Element * aB )
                {
                    return aA->id() < aB->id();
                }
            } opElementId;

            // make sure that all flags are reset
            mMesh->unflag_all_elements() ;

            // collect elements per cut data
            for ( CutData * tCutData : mCutData )
            {
                tCutData->collect_elements( mNonPhiBlocks );
            }

            // count total number of elements
            index_t tCount = 0 ;
            for ( CutData * tCutData : mCutData )
            {
                tCount += tCutData->elements().size() ;
            }

            mElements.set_size( tCount, nullptr ) ;

            tCount = 0 ;

            for ( CutData * tCutData : mCutData )
            {
                for ( Element * tElement : tCutData->elements() )
                {
                    mElements( tCount++ ) = tElement ;
                }
            }

            // make sure that elements are unique
            unique( mElements ) ;

            // order elements after their ID
            sort( mElements , opElementId );

            // set the element indices
            tCount = 0 ;
            for ( Element * tElement : mElements )
            {
                tElement->set_index( tCount++ ) ;
                tElement->unflag();
            }
        }

//-----------------------------------------------------------------------------

        void
        CutProcessor::collect_facets()
        {
            for ( CutData * tCutData : mCutData )
            {
                for ( Element * tElement : tCutData->elements() )
                {
                    tElement->unflag();
                }
            }

            if ( mMesh->number_of_dimensions() == 2 )
            {
                mMesh->unflag_all_edges();

                for ( CutData * tCutData : mCutData )
                {

                    for ( Element * tElement : tCutData->elements() )
                    {
                        tElement->flag();
                    }

                    // flag edges
                    index_t tCount = 0 ;
                    for ( Element * tElement : tCutData->elements() )
                    {
                        int tCase = tCutData->cut_case( tCount++ );

                        if ( tCase > 0 )
                        {
                            tElement->edge( tCase-1 )->flag() ;
                        }
                    }

                    // count edges
                    tCount = 0 ;
                    for ( Edge * tEdge : mMesh->edges() )
                    {
                        if ( tEdge->is_flagged() )
                        {
                            ++tCount ;
                        }
                    }

                    // collect edges in temporary container
                    Cell< Edge * > tEdges( tCount, nullptr ) ;
                    tCount = 0 ;

                    for ( Edge * tEdge : mMesh->edges() )
                    {
                        if ( tEdge->is_flagged() )
                        {
                            tEdges( tCount++ ) = tEdge ;
                        }
                    }

                    // remove self canceling edges
                    for ( Edge * tEdge : tEdges )
                    {
                        if ( tEdge->number_of_elements() > 1 )
                        {
                            if ( tEdge->element( 0 )->is_flagged() && tEdge->element( 1 )->is_flagged() )
                            {
                                tEdge->unflag() ;
                            }
                        }
                    }

                    // unflag edges on boundaries
                    for ( id_t tID : mPhiBoundaries )
                    {
                        Cell< Facet * > & tFacets = mMesh->sideset( tID )->facets() ;
                        for ( Facet * tFacet : tFacets )
                        {
                            tFacet->master()->edge( tFacet->index_on_master() )->unflag() ;
                        }
                    }

                    // count edges
                    tCount = 0 ;
                    for ( Edge * tEdge : tEdges )
                    {
                        if ( tEdge->is_flagged() )
                        {
                            ++tCount ;
                        }
                    }

                    BELFEM_ERROR( tCount > 0,
                        "Thin cut %u is empty. The solver supports only thin cuts; thick cuts are not implemented yet. No conjugated edge survived the self-cancellation and boundary passes. The thick-cut cohomology generator is not a tight representative: locally, it is not the boundary of a one-sided region (see src/homology/doc/thick_thin_cuts_and_conjugate_edges.md). No user-side remedy exists yet: the representative the cohomology engine chose cannot be pushed onto element edges on this mesh.",
                        ( unsigned int ) tCutData->index() );

                    // collect edges
                    Cell< Edge * > & tThinCut = tCutData->thin_cut_edges() ;
                    tThinCut.set_size( tCount, nullptr ) ;

                    tCount = 0 ;
                    for ( Edge * tEdge : tEdges )
                    {
                        if ( tEdge->is_flagged() )
                        {
                            tThinCut( tCount++ ) = tEdge ;
                            tEdge->unflag() ;
                        }
                    }

                    for ( Element * tElement : tCutData->elements() )
                    {
                        tElement->unflag();
                    }
                }
            }
            else
            {
                mMesh->unflag_all_faces();
                mMesh->unflag_all_edges();

                for ( id_t tID : mPhiBoundariesAndPeriodic )
                {
                    mMesh->sideset( tID )->flag_edges() ;
                }

                for ( CutData * tCutData : mCutData )
                {
                    // flag facets
                    index_t tCount = 0 ;
                    for ( Element * tElement : tCutData->elements() )
                    {
                        int tCase = tCutData->cut_case( tCount++ );

                        if ( 0 < tCase && tCase < 5 )
                        {
                            tElement->face( tCase-1 )->flag() ;
                            tElement->flag();
                        }
                    }

                    // quotient self-cancel ( todo/closed/periodic_cap_cut_emission.md ):
                    // the interior self-cancel below works on shared Face
                    // objects and cannot see that two periodically identified
                    // cap faces are the same quotient face. If both copies
                    // claim the cut, the claims cancel exactly like an
                    // interior double-flag; a single claim stays where it
                    // fired ( one-sided emission + straight ties imposes the
                    // jump exactly once in the quotient )
                    for ( Face * tFace : mMesh->faces() )
                    {
                        if ( tFace->is_flagged() && tFace->is_periodic() )
                        {
                            Face * tPartner = tFace->periodic() ;
                            if ( tPartner != tFace && tPartner->is_flagged() )
                            {
                                tFace->unflag() ;
                                tPartner->unflag() ;
                            }
                        }
                    }

                    // count facets
                    tCount = 0 ;
                    for ( Face * tFace : mMesh->faces() )
                    {
                        if ( tFace->is_flagged() )
                        {
                            ++tCount ;
                        }
                    }

                    // collect facets in temporary container
                    Cell< Face * > tFaces( tCount, nullptr ) ;
                    tCount = 0 ;
                    for ( Face * tFace : mMesh->faces() )
                    {
                        if ( tFace->is_flagged() )
                        {
                            tFaces( tCount++ ) = tFace ;
                        }
                    }

                    // remove self canceling faces
                    for ( Face * tFace : tFaces )
                    {
                        if ( tFace->slave() != nullptr && tFace->master() != nullptr )
                        {
                            if ( tFace->master()->is_flagged() && tFace->slave()->is_flagged() )
                            {
                                tFace->unflag() ;
                            }
                        }
                    }
                    // unflag faces on boundaries
                    for ( id_t tID : mPhiBoundaries )
                    {
                        Cell< Facet * > & tFacets = mMesh->sideset( tID )->facets() ;
                        for ( Facet * tFacet : tFacets )
                        {
                            tFacet->master()->face( tFacet->index_on_master() )->unflag() ;
                        }
                    }

                    // collect faces
                    tCount = 0 ;
                    for ( Face * tFace : tFaces )
                    {   if ( tFace->is_flagged() )
                        {
                            ++tCount ;
                        }
                    }

                    Cell< Face * > tTempCut( tCount, nullptr );


                    tTempCut.set_size( tCount, nullptr ) ;

                    tCount = 0 ;
                    for ( Face * tFace : tFaces )
                    {
                        if ( tFace->is_flagged() )
                        {
                            tTempCut( tCount++ ) = tFace ;
                        }
                    }

                    // faces emitted before the dangling-face pruning below;
                    // a healthy one-sheet cut loses only a few of them
                    const index_t tNumEmitted = tTempCut.size() ;
                    index_t tNumIterations = 0 ;

                    index_t tCountOld = gNoIndex ;
                    tCount = 0 ;

                    for ( index_t k=0; tCountOld!=tCount; ++k )
                    {
                        tCountOld = tCount ;
                        tNumIterations = k ;

                        tCount = 0 ;
                        for ( Face * tFace : tTempCut )
                        {
                            if ( ! tFace->is_flagged() ) continue;

                            for ( uint e=0; e<tFace->number_of_edges(); ++e )
                            {
                                uint tNumFaces = 0 ;

                                // grab edge
                                Edge * tEdge = tFace->edge( e ) ;

                                // check if edge sits on boundary
                                if ( tEdge->is_flagged() ) continue ;

                                for ( uint f=0; f<tEdge->number_of_faces(); ++f )
                                {
                                    if ( tEdge->face( f )->is_flagged() )
                                    {
                                        ++tNumFaces ;
                                    }
                                }

                                if ( tNumFaces < 2 )
                                {
                                    tFace->unflag() ;
                                    break ;
                                }
                            }

                            if ( tFace->is_flagged() )
                            {
                                ++tCount ;
                            }


                        }

                        // note, we want to prevent an infinite loop here
                        // this number might have to be adjusted
                        BELFEM_ERROR( k<1000, "Too many iterations" );
                    }

                    // the pruning is a no-op when the emitted faces close into
                    // one sheet ( a tight generator ); if the generator is not
                    // tight the set has holes and the loop eats the surface
                    // inward from every hole, leaving fragments that impose no
                    // jump. Abort rather than hand the solver a cut that is not
                    // there ( observed on corc: 7k faces -> 242 in 70 passes )
                    const index_t tNumPruned = tNumEmitted - tCount ;

                    BELFEM_ERROR( tCount > 0,
                        "Thin cut %u is empty. The solver supports only thin cuts; thick cuts are not implemented yet. No conjugated face survived emission, self-cancellation, boundary removal, or dangling-face pruning (%lu faces entered pruning; %lu passes). The thick-cut cohomology generator is not a tight one-sheet representative: locally, it is not the boundary of a one-sided region (see src/homology/doc/thick_thin_cuts_and_conjugate_edges.md). No user-side remedy exists yet: the representative the cohomology engine chose cannot be pushed onto element faces on this mesh.",
                        ( unsigned int ) tCutData->index(),
                        ( long unsigned int ) tNumEmitted,
                        ( long unsigned int ) tNumIterations );

                    BELFEM_ERROR( tNumPruned <= tNumEmitted / 4,
                        "Thin cut %u collapsed. The solver supports only thin cuts; thick cuts are not implemented yet. Dangling-face pruning removed %lu of %lu emitted faces in %lu passes (a healthy cut loses a few faces in one or two passes). The thick-cut cohomology generator is not a tight one-sheet representative: locally, it is not the boundary of a one-sided region, so the thin cut would not impose the jump (see src/homology/doc/thick_thin_cuts_and_conjugate_edges.md). No user-side remedy exists yet: the representative the cohomology engine chose cannot be pushed onto element faces on this mesh.",
                        ( unsigned int ) tCutData->index(),
                        ( long unsigned int ) tNumPruned,
                        ( long unsigned int ) tNumEmitted,
                        ( long unsigned int ) tNumIterations );

                    Cell< Face * > & tThinCut = tCutData->thin_cut_faces() ;
                    tThinCut.set_size( tCount, nullptr ) ;
                    tCount = 0 ;
                    for ( Face * tFace : tTempCut )
                    {
                        if ( tFace->is_flagged() )
                        {
                            tThinCut( tCount++ ) = tFace ;
                        }
                    }

                    for ( Face * tFace : tTempCut )
                    {
                        tFace->unflag() ;
                    }
                    for ( Element * tElement : tCutData->elements() )
                    {
                        tElement->unflag();
                    }

                }

                for ( id_t tID : mPhiBoundariesAndPeriodic )
                {
                    mMesh->sideset( tID )->unflag_edges() ;
                }

            }
        }

//-----------------------------------------------------------------------------


        void
        CutProcessor::compute_node_bitsets()
        {

            uint tNumNodesPerElement = number_of_nodes( mElementType );
            mMesh->unflag_all_edges() ;
            for ( Edge * tEdge : mCohomologyEdges )
            {
                tEdge->flag() ;
            }

            // create the temporary bitsets
            Cell< DynamicBitset * > tBitsets( tNumNodesPerElement, nullptr );
            for ( uint k=0; k<tNumNodesPerElement; ++k )
            {
                tBitsets( k ) = new DynamicBitset( mNumberOfCuts ) ;
            }

            switch ( mElementType )
            {
                case ElementType::TRI3 :
                {
                    for ( Element * tElement : mElements )
                    {
                        this->flip_node_bitsets_tri3( tElement, tBitsets ) ;
                    }
                    break ;
                }
                case ElementType::TRI6 :
                {
                    for ( Element * tElement : mElements )
                    {
                        this->flip_node_bitsets_tri6( tElement, tBitsets ) ;
                    }
                    break ;
                }
                case ElementType::TET4 :
                {
                    for ( Element * tElement : mElements )
                    {
                        this->flip_node_bitsets_tet4( tElement, tBitsets ) ;
                    }
                    break ;
                }
                case ElementType::TET10 :
                {
                    for ( Element * tElement : mElements )
                    {
                        this->flip_node_bitsets_tet10( tElement, tBitsets ) ;
                    }
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "This should not happen");
                }
            }

            // delete the bitset
            for ( DynamicBitset * tBitset : tBitsets )
            {
                delete tBitset ;
            }
        }

//-----------------------------------------------------------------------------

        void
        CutProcessor::collect_nodes()
        {
            mMesh->unflag_all_nodes() ;

            for ( CutData * tCutData : mCutData )
            {
                tCutData->flag_nodes();
            }

            // count nodes
            index_t tCount = 0 ;
            for ( Node * tNode : mMesh->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    ++tCount ;
                }
                else
                {
                    tNode->set_index( gNoIndex );
                }
            }

            // collect nodes
            mNodes.set_size( tCount, nullptr );
            tCount = 0 ;
            for ( Node * tNode : mMesh->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    mNodes( tCount++ ) = tNode ;
                }
            }

            // the sorting is not neccessary, but it's fast and easier to debug
            sort( mNodes.begin(), mNodes.end(), []( Node * a, Node * b ) { return a->id() < b->id(); } ) ;

            tCount = 0 ;
            for ( Node * tNode : mNodes )
            {
                tNode->set_index( tCount++ ) ;
            }

        }

//-----------------------------------------------------------------------------

        void
        CutProcessor::determine_cut_sets()
        {
            mMesh->unflag_all_edges() ;
            for ( Edge * tEdge : mCohomologyEdges )
            {
                tEdge->flag() ;
            }

            uint tNumNodesPerElement = number_of_nodes( mElementType );

            // create the temporary bitsets
            Cell< DynamicBitset * > tBitsets( tNumNodesPerElement, nullptr );
            for ( uint k=0; k<tNumNodesPerElement; ++k )
            {
                tBitsets( k ) = new DynamicBitset( mNumberOfCuts ) ;
            }

            switch ( mElementType )
            {
                case ElementType::TRI3 :
                {
                    for ( Element * tElement : mElements )
                    {
                        this->check_edges_tri3( tElement, tBitsets ) ;
                        this->check_node_bitsets( tBitsets ) ;
                    }
                    break ;
                }
                case ElementType::TRI6 :
                {
                    for ( Element * tElement : mElements )
                    {
                        this->check_edges_tri6( tElement, tBitsets ) ;
                        this->check_node_bitsets( tBitsets ) ;
                    }
                    break ;
                }
                case ElementType::TET4 :
                {
                    for ( Element * tElement : mElements )
                    {
                        this->check_edges_tet4( tElement, tBitsets ) ;
                        this->check_node_bitsets( tBitsets ) ;
                    }
                    break ;
                }
                case ElementType::TET10 :
                {
                    for ( Element * tElement : mElements )
                    {
                        this->check_edges_tet10( tElement, tBitsets ) ;
                        this->check_node_bitsets( tBitsets ) ;
                    }
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "This should not happen");
                }
            }

            // delete the bitset
            for ( DynamicBitset * tBitset : tBitsets )
            {
                delete tBitset ;
            }

            index_t tCount = 0 ;
            mPatterns.clear();
            for ( auto tPair : mHashToBitsets )
            {
                Cell< DynamicBitset * > & tBitsets = tPair.second ;

                for ( DynamicBitset * tBitset : tBitsets )
                {
                    tBitset->set_index( tCount++ ) ;
                    mPatterns.push( tBitset->to_hex() );
                }
            }
        }

        index_t
        CutProcessor::pattern_index( const DynamicBitset * aBitset )
        {
            // end delete
            Cell< DynamicBitset * > & tBitsetWithSameHash
                    = mHashToBitsets( aBitset->hash() );

            for ( DynamicBitset * tOtherBitset : tBitsetWithSameHash )
            {
                if ( *tOtherBitset == *aBitset )
                {
                    return tOtherBitset->index() ;
                }
            }
            return gNoIndex ;
        }

//-----------------------------------------------------------------------------

        void
        CutProcessor::create_abstract_nodes()
        {
            // create abstract nodes
            mAbstractNodes.set_size( mNumberOfCuts, nullptr ) ;
            for ( uint c=0; c<mNumberOfCuts; ++c )
            {
                mAbstractNodes( c ) = new Node( ++mMaxNodeID ) ;
            }
        }

//-----------------------------------------------------------------------------

        void
        CutProcessor::create_cut_sets()
        {
            index_t tNumSets = mPatterns.size() ;


            for ( index_t s=0; s<tNumSets; ++s )
            {
                mCutSets[ mPatterns( s ) ] = new CutSet( mMesh, mNodes, mPatterns( s ), mNumberOfCuts ) ;
            }

        }
//-----------------------------------------------------------------------------

        void
        CutProcessor::create_thin_cut_sidesets()
        {
            mMaxElementID = mMesh->max_element_id();
            mMaxSidesetID = mMesh->max_block_and_sideset_id() ;

            if ( mMesh->number_of_dimensions() == 2 )
            {
                mMesh->unflag_all_edges() ;
            }
            else
            {
                mMesh->unflag_all_faces() ;
            }

            for ( CutData * tCutData : mCutData )
            {
                tCutData->add_thin_cut_sidesets_to_mesh(
                        mMaxSidesetID, mMaxElementID );
            }
        }


//-----------------------------------------------------------------------------

        void
        CutProcessor::check_edges_tri3(
                Element * aElement,
                Cell< DynamicBitset * > & aBitsets )
        {
            for ( DynamicBitset * tBitset : aBitsets )
            {
                tBitset->reset() ;
            }

            this->check_edge( aElement, 0, 0, 1, aBitsets ) ;
            this->check_edge( aElement, 1, 1, 2, aBitsets ) ;
            this->check_edge( aElement, 2, 2, 0, aBitsets ) ;
        }

        void
        CutProcessor::flip_node_bitsets_tri3(
                Element * aElement,
                Cell< DynamicBitset * > & aBitsets )
        {

            this->check_edges_tri3( aElement, aBitsets );

            for ( uint k=0; k<3; ++k )
            {
                // skip empty bitsets
                if ( aBitsets( k )->count() == 0 ) continue;

                // lock the bitset
                aBitsets( k )->lock();

                // get the hex value
                const string & tHex = mPatterns( this->pattern_index( aBitsets( k ) ) );

                // set the flag of the node bitset
                mCutSets( tHex )->node_bitset()->set( aElement->node( k )->index() ) ;

                // unlock the bitset
                aBitsets( k )->unlock();
            }
        }

//-----------------------------------------------------------------------------

        void
        CutProcessor::check_edges_tet4(
                Element * aElement,
                Cell< DynamicBitset * > & aBitsets )
        {
            for ( DynamicBitset * tBitset : aBitsets )
            {
                tBitset->reset() ;
            }

            this->check_edge( aElement, 0, 0, 1, aBitsets ) ;
            this->check_edge( aElement, 1, 1, 2, aBitsets ) ;
            this->check_edge( aElement, 2, 2, 0, aBitsets ) ;
            this->check_edge( aElement, 3, 0, 3, aBitsets ) ;
            this->check_edge( aElement, 4, 1, 3, aBitsets ) ;
            this->check_edge( aElement, 5, 2, 3, aBitsets ) ;
        }

        void
        CutProcessor::flip_node_bitsets_tet4(
                Element * aElement,
                Cell< DynamicBitset * > & aBitsets )
        {
            this->check_edges_tet4( aElement, aBitsets );

            for ( uint k=0; k<4; ++k )
            {
                // skip empty bitsets
                if ( aBitsets( k )->count() == 0 ) continue;

                // lock the bitset
                aBitsets( k )->lock();

                // get the hex value
                const string & tHex = mPatterns( this->pattern_index( aBitsets( k ) ) );

                // set the flag of the node bitset
                mCutSets( tHex )->node_bitset()->set( aElement->node( k )->index() ) ;

                // unlock the bitset
                aBitsets( k )->unlock();
            }

        }

//-----------------------------------------------------------------------------

        void
        CutProcessor::check_edges_tri6(
                Element * aElement,
                Cell< DynamicBitset * > & aBitsets )
        {
            for ( DynamicBitset * tBitset : aBitsets )
            {
                tBitset->reset() ;
            }

            this->check_edge( aElement, 0, 0, 1, aBitsets ) ;
            this->check_edge( aElement, 1, 1, 2, aBitsets ) ;
            this->check_edge( aElement, 2, 2, 0, aBitsets ) ;
            this->check_midside( 0, 1, 3, aBitsets );
            this->check_midside( 1, 2, 4, aBitsets );
            this->check_midside( 2, 0, 5, aBitsets );

        }

        void
        CutProcessor::flip_node_bitsets_tri6(
                Element * aElement,
                Cell< DynamicBitset * > & aBitsets )
        {
            this->check_edges_tri6( aElement, aBitsets );

            for ( uint k=0; k<6; ++k )
            {
                // skip empty bitsets
                if ( aBitsets( k )->count() == 0 ) continue;

                // lock the bitset
                aBitsets( k )->lock();

                // get the hex value
                const string & tHex = mPatterns( this->pattern_index( aBitsets( k ) ) );

                // set the flag of the node bitset
                mCutSets( tHex )->node_bitset()->set( aElement->node( k )->index() ) ;

                // unlock the bitset
                aBitsets( k )->unlock();
            }
        }

//-----------------------------------------------------------------------------

        void
        CutProcessor::check_edges_tet10(
                Element * aElement,
                Cell< DynamicBitset * > & aBitsets )
        {
            for ( DynamicBitset * tBitset : aBitsets )
            {
                tBitset->reset() ;
            }

            this->check_edge( aElement, 0, 0, 1, aBitsets ) ;
            this->check_edge( aElement, 1, 1, 2, aBitsets ) ;
            this->check_edge( aElement, 2, 2, 0, aBitsets ) ;
            this->check_edge( aElement, 3, 0, 3, aBitsets ) ;
            this->check_edge( aElement, 4, 1, 3, aBitsets ) ;
            this->check_edge( aElement, 5, 2, 3, aBitsets ) ;
            this->check_midside( 0, 1, 4, aBitsets );
            this->check_midside( 1, 2, 5, aBitsets );
            this->check_midside( 2, 0, 6, aBitsets );
            this->check_midside( 0, 3, 7, aBitsets );
            this->check_midside( 1, 3, 8, aBitsets );
            this->check_midside( 2, 3, 9, aBitsets );
        }

        void
        CutProcessor::flip_node_bitsets_tet10(
                Element * aElement,
                Cell< DynamicBitset * > & aBitsets )
        {
            this->check_edges_tet10( aElement, aBitsets );

            for ( uint k=0; k<10; ++k )
            {
                // skip empty bitsets
                if ( aBitsets( k )->count() == 0 ) continue;

                // lock the bitset
                aBitsets( k )->lock();

                // get the hex value
                const string & tHex = mPatterns( this->pattern_index( aBitsets( k ) ) );

                // set the flag of the node bitset
                mCutSets( tHex )->node_bitset()->set( aElement->node( k )->index() ) ;

                // unlock the bitset
                aBitsets( k )->unlock();
            }
        }

        void
        CutProcessor::check_node_bitsets(
                Cell< DynamicBitset * > & aBitsets )
        {
            for ( DynamicBitset * tBitset : aBitsets )
            {
                tBitset->lock() ;
                if ( tBitset->count() == 0 ) continue;

                // get the hash of the bitset
                size_t tHash =  tBitset->hash() ;

                // test if has exists in map
                if ( ! mHashToBitsets.key_exists( tHash ) )
                {
                    // create a new Cell
                    Cell< DynamicBitset * > tBitsetsWithSameHash( 1, nullptr );

                    // create a new bitset
                    DynamicBitset * tNewBitset = new DynamicBitset( tBitset->size() );

                    // set bitsets
                    *tNewBitset = *tBitset ;

                    tNewBitset->lock() ;

                    // add bitset to container
                    tBitsetsWithSameHash( 0 ) = tNewBitset ;

                    // add cell to map
                    mHashToBitsets[ tHash ] = tBitsetsWithSameHash ;
                }
                else
                {
                    // grab the bitsets with the same hash
                    Cell< DynamicBitset * > & tBitsetWithSameHash = mHashToBitsets( tHash );

                    bool tIsNew = true ;

                    // if the hash exists, we have to check if the bitsets already exists
                    for ( DynamicBitset * tOtherBitset : tBitsetWithSameHash )
                    {
                        if ( *tOtherBitset == *tBitset )
                        {
                            tIsNew = false ;
                            break ;
                        }
                    }

                    if ( tIsNew )
                    {

                        // create a new bitset
                        DynamicBitset * tNewBitset = new DynamicBitset( tBitset->size() );

                        // set bitsets
                        *tNewBitset = *tBitset ;

                        tNewBitset->lock() ;

                        tBitsetWithSameHash.push( tNewBitset );
                    }
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        CutProcessor::duplicate_nodes()
        {
            for ( auto tPair : mCutSets )
            {
                tPair.second->create_duplicates( mMaxNodeID, mAbstractNodes ) ;
            }
        }

//-----------------------------------------------------------------------------

        void
        CutProcessor::relink_elements()
        {
            ElementType tType = mMesh->max_element_order() == 1 ?
                                ( mMesh->number_of_dimensions() == 2 ? ElementType::TRI3 : ElementType::TET4 ) :
                                ( mMesh->number_of_dimensions() == 2 ? ElementType::TRI6 : ElementType::TET10 ) ;

            uint tNumNodesPerElement = number_of_nodes( tType );

            // create the temporary bitsets
            Cell< DynamicBitset * > tBitsets( tNumNodesPerElement, nullptr );
            for ( uint k=0; k<tNumNodesPerElement; ++k )
            {
                tBitsets( k ) = new DynamicBitset( mNumberOfCuts ) ;
            }

            switch ( tType )
            {
                case ElementType::TRI3 :
                {
                    for ( Element * tElement : mElements )
                    {
                        this->flip_node_bitsets_tri3( tElement, tBitsets ) ;
                        this->relink_element( tElement, tBitsets ) ;
                    }
                    break ;
                }
                case ElementType::TRI6 :
                {
                    for ( Element * tElement : mElements )
                    {
                        this->flip_node_bitsets_tri6( tElement, tBitsets ) ;
                        this->relink_element( tElement, tBitsets ) ;
                    }
                    break ;
                }
                case ElementType::TET4 :
                {
                    for ( Element * tElement : mElements )
                    {
                        this->flip_node_bitsets_tet4( tElement, tBitsets ) ;
                        this->relink_element( tElement, tBitsets ) ;
                    }
                    break ;
                }
                case ElementType::TET10 :
                {
                    for ( Element * tElement : mElements )
                    {
                        this->flip_node_bitsets_tet10( tElement, tBitsets ) ;
                        this->relink_element( tElement, tBitsets ) ;
                    }
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "This should not happen");
                }
            }

            // delete the bitset
            for ( DynamicBitset * tBitset : tBitsets )
            {
                delete tBitset ;
            }
        }

//-----------------------------------------------------------------------------

        void
        CutProcessor::relink_element(
                Element * aElement,
                Cell< DynamicBitset * > & aBitsets )
        {
            for ( uint k=0; k<aElement->number_of_nodes(); ++k )
            {
                if ( aBitsets( k )->count() > 0 )
                {
                    // get the duplicate
                    Node * tDup = mCutSets( aBitsets( k )->to_hex() )->duplicate( aElement->node( k ) ) ;
                    BELFEM_ASSERT( tDup != nullptr , "Could not find duplicate") ;

                    // relink the element
                    aElement->insert_node( tDup, k ) ;
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        CutProcessor::collect_duplicates()
        {
            mMesh->unflag_all_nodes();

            for ( id_t tID : mPhiBlocks )
            {
                mMesh->block( tID )->flag_nodes() ;
            }

            // count nodes
            index_t tCount = mAbstractNodes.size() ;

            for ( auto & tPair : mCutSets )
            {
                for ( auto & tPair2 : tPair.second->duplicate_map() )
                {
                    if ( tPair2.second->is_flagged() )
                    {
                        ++tCount ;
                    }
                }
            }

            Cell< Node * > tNodes( tCount, nullptr ) ;

            tCount = 0 ;
            for ( Node * tAbstractNode : mAbstractNodes )
            {
                tNodes( tCount++ ) = tAbstractNode ;
            }

            for ( auto & tPair : mCutSets )
            {
                for ( auto & tPair2 : tPair.second->duplicate_map() )
                {
                    if ( tPair2.second->is_flagged() )
                    {
                        tNodes( tCount++ ) = tPair2.second ;
                    }
                    else
                    {
                        delete tPair2.second ;
                    }
                }
                tPair.second->duplicate_map().clear() ;
            }

            // add nodes to mesh
            append( mMesh->nodes(), tNodes );
        }

//-----------------------------------------------------------------------------
    }
}
