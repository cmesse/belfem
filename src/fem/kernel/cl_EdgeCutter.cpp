/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include "cl_EdgeCutter.hpp"
#include "cl_Mesh_ConnectivityCalculator.hpp"
#include "fn_norm.hpp"
#include "fn_cross.hpp"
#include "fn_dot.hpp"
#include "cl_Element_Factory.hpp"

namespace belfem
{
    namespace mesh
    {
        EdgeCutter::EdgeCutter( Mesh * aMesh )
            : mMesh( aMesh ), mElementMapper( new fem::ElementMapper( ) )
        {
            mMesh->unflag_all_facets();

            if ( ! mMesh->test_connectivity( Connectivity::FacetToFacet ) )
            {
                ConnectivityCalculator tCalc( aMesh );
                tCalc.connect_facets_to_facets();
            }
            mP.set_size( 3 );
            mQ.set_size( 3 );
            mU.set_size( 3 );
            mV.set_size( 3 );
            mW.set_size( 3 );
        }

        EdgeCutter::~EdgeCutter()
        {
            delete mFacetBitset;
            delete mElementMapper;
        }

        void
        EdgeCutter::select_sidesets( const Vector< id_t > & aSideSets )
        {
            mSideSets.reserve( aSideSets.length() );
            mFacets.clear() ;

            for ( id_t tID : aSideSets )
            {
                SideSet * tSideSet = mMesh->sideset( tID );
                append(mFacets, tSideSet->facets() );
                mSideSets.push( tSideSet );
            }

            delete mFacetBitset ;
            mFacetBitset = new DynamicBitset( mFacets.size() );

        }


        void
        EdgeCutter::process_curve( Curve * aCurve )
        {
            BELFEM_ASSERT( mSideSets.size() > 0, "No sidesets selected" );

            this->facet_bfs( aCurve );

            Matrix< real > tNormals, tTangents, tBinomials;
            Vector< real > tDistances;
            this->compute_node_vectors( aCurve, tNormals, tTangents, tBinomials, tDistances );

            Cell< Node * > tNodes;
            this->compute_temporary_nodes( aCurve, tBinomials, tNodes );
            this->project_temporary_nodes( mFacets, tNodes );

            // for debugging
            this->save_testmesh( tNodes );

            for ( Node * tNode : tNodes )
            {
                delete tNode ;
            }
        }

        void
        EdgeCutter::facet_bfs( Curve * aCurve )
        {
            BELFEM_ASSERT( mMesh->test_connectivity( Connectivity::FacetToFacet ), "Facet connectivity is required for facet BFS" );

            index_t tCount = 0 ;
            for ( Facet * tFacet : mFacets )
            {
                tFacet->set_level( gNoIndex );
                tFacet->flag( 1 );
                tFacet->set_index( tCount++ );
                for ( uint k=0; k<tFacet->number_of_nodes(); ++k )
                {
                    tFacet->node( k )->unflag( 1 );
                }
            }

            Cell< Node * > & tNodes = aCurve->nodes();
            for ( Node * tNode : tNodes )
            {
                tNode->flag( 1 );
            }

            mFacetBitset->reset();
            for ( Facet * tFacet : mFacets )
            {
                for ( uint f=0; f<tFacet->number_of_facets(); ++f )
                {
                    Facet * tOther = tFacet->facet( f );
                    if ( tOther->is_flagged( 1 ) )
                    {
                        for ( uint k=0; k<tOther->number_of_corner_nodes(); ++k )
                        {
                            if ( tOther->node( k )->original()->is_flagged( 1 ) )
                            {
                                mFacetBitset->set( tFacet->index() );
                                break;
                            }
                        }
                    }
                }
            }
            Cell< index_t > tIndices ;
            index_t tLevel = 0 ;


            tIndices.set_size( 1, 0 );
            while ( tIndices.size() != 0 )
            {
                mFacetBitset->where( tIndices );
                mFacetBitset->reset();

                for ( index_t tIndex : tIndices )
                {
                    Facet * tFacet = mFacets( tIndex );
                    tFacet->set_level( tLevel );
                    ++tCount ;
                    for ( uint f=0; f<tFacet->number_of_facets(); ++f )
                    {
                        Facet * tOther = tFacet->facet( f );
                        if ( ! tOther->is_flagged( 1 ) ) continue;
                        if ( tOther->level() != gNoIndex ) continue;
                        mFacetBitset->set( tOther->index() );
                    }
                }
                ++tLevel ;
            }
            std::cout << "Found " << tCount << " facets on " << tLevel << " levels" << std::endl;
            std::sort( mFacets.begin(), mFacets.end(), []( Facet * a, Facet * b ) { return a->level() < b->level(); } );
            mFacetsPerLayerBegin.set_size( tLevel, 0 );
            mFacetsPerLayerEnd.set_size( tLevel, mFacets.size() );
            tCount = 0 ;
            tLevel = 0 ;

            for ( Facet * tFacet : mFacets )
            {
                if ( tFacet->level() != tLevel )
                {
                    mFacetsPerLayerEnd( tLevel ) = tCount ;
                    tFacet->set_index( tCount++ );
                    mFacetsPerLayerBegin( ++tLevel ) = tCount ;
                }
                else
                {
                    tFacet->set_index( tCount++ );
                }
            }
            ++tLevel ;
            std::cout << "Found " << tCount << " facets on " << tLevel << " levels" << std::endl;

            for ( index_t l=0; l<tLevel; ++l )
            {
                std::cout << "   " << mFacetsPerLayerBegin( l ) << " - " << mFacetsPerLayerEnd( l ) << std::endl;
            }
        }

        void
        EdgeCutter::compute_node_vectors(
            Curve * aCurve,
            Matrix< real > & aNormals,
            Matrix< real > & aTangents,
            Matrix< real > & aBinomials,
            Vector< real > & aDistances )
        {
            Cell< Node * > & tNodes = aCurve->nodes();

            index_t tCount = tNodes.size() ;
            aNormals.set_size( 3, tCount );
            aTangents.set_size( 3, tCount );
            aBinomials.set_size( 3, tCount );
            aDistances.set_size( tCount );
            tCount = 0 ;
            for ( Node * tNode : tNodes )
            {
                if ( tNode == tNodes.first() )
                {
                    tNodes( 0 )->get_coords( mU );
                    tNodes( 1 )->get_coords( mV );
                    tNodes( 2 )->get_coords( mW );
                }
                else if ( tNode == tNodes.last() )
                {
                    tNodes( tNodes.size() - 3 )->get_coords( mU );
                    tNodes( tNodes.size() - 2 )->get_coords( mV );
                    tNodes( tNodes.size() - 1 )->get_coords( mW );
                }
                else
                {
                    tNodes( tCount - 1 )->get_coords( mU );
                    tNodes( tCount )->get_coords( mV );
                    tNodes( tCount + 1 )->get_coords( mW );
                }

                mW -= mV ;
                mV -= mU ;
                mU = mV + mW ;
                mU /= norm( mU );
                aTangents.set_col( tCount, mU );

                for ( uint f=0; f<tNode->number_of_facets(); ++f )
                {
                    Facet * tFacet = tNode->facet( f );
                    BELFEM_ASSERT( tFacet->element()->type() == ElementType::TRI3 ,"need TRI3 facet" );
                    uint c = 0 ;
                    for ( uint k=0; k<tFacet->number_of_corner_nodes(); ++k )
                    {
                        if ( tFacet->node( k )->original()->is_flagged( 1 ) ) ++c ;
                    }
                    if ( c == 2 )
                    {
                        tFacet->node( 0 )->original()->get_coords( mW );
                        tFacet->node( 1 )->original()->get_coords( mU );
                        tFacet->node( 2 )->original()->get_coords( mV );

                        mU -= mW ;
                        mV -= mW ;

                        mW = cross( mU, mV );
                        mW /= norm( mW );

                        aNormals.set_col( tCount, mW );

                        mP.fill( 0. );
                        for ( uint k=0; k<tFacet->number_of_corner_nodes(); ++k )
                        {
                            tFacet->node( k )->get_coords( mV );
                            mP += mV ;
                        }
                        mP /= tFacet->number_of_corner_nodes();
                        break;
                    }
                }

                // binomial = tangent × normal
                mU = aTangents.col( tCount );
                mV = aNormals.col( tCount );
                mW = cross( mU, mV );
                mW /= norm( mW );
                aBinomials.set_col( tCount, mW );

                aDistances( tCount ) = dot( mV, mP );
                ++tCount ;
            }
        }

        void
        EdgeCutter::compute_temporary_nodes( Curve * aCurve, Matrix< real > & aBinomials, Cell< Node * > & aNodes )
        {
            real tScale = this->determine_sign( aCurve, aBinomials ) * mConnectorWitdh ;

            aNodes.set_size(  aCurve->nodes().size(), nullptr );
            index_t tCount = 0 ;
            for ( Node * tOrg : aCurve->nodes() )
            {
                tOrg->get_coords( mU );
                mV = aBinomials.col( tCount );
                mW = mU + tScale * mV ;
                Node * tDup = new Node( tOrg->id() );
                tDup->set_coords( mW );
                aNodes( tCount++ ) = tDup ;
            }
        }


        real
        EdgeCutter::determine_sign( Curve * aCurve, Matrix< real > & aBinomials )
        {
            Node * tNode = aCurve->nodes().first();

            Facet * tFacet = nullptr ;
            Node * tOther = nullptr ;

            uint c = 0;
            for ( uint f=0; f<tNode->number_of_facets(); ++f )
            {
                c=0;
                tFacet = tNode->facet( f );
                if ( ! tFacet->is_flagged( 1 ) ) continue;
                for ( uint k=0; k<tFacet->number_of_corner_nodes(); ++k )
                {
                    if ( tFacet->node( k )->original()->is_flagged( 1 ) ) ++c ;
                }
                if ( c == 2 ) break;

            }

            for ( uint k=0; k<tFacet->number_of_corner_nodes(); ++k )
            {
                if ( ! tFacet->node( k )->original()->is_flagged( 1 ) )
                {
                    tOther = tFacet->node( k )->original();
                    break;
                }
            }
            BELFEM_ASSERT( tOther != nullptr, "No facet found" );

            tNode->get_coords( mU );
            mV = aBinomials.col( 0 );
            tOther->get_coords( mW );

            mU -= mW ;

            // find the nb (flip sign deliberately)
            return  dot( mU, mV ) > 0. ? -1 : 1. ;
        }

        void
        EdgeCutter::save_testmesh( Cell< Node * > & aNodes )
        {
            Mesh * tMesh = new Mesh( 3 );

            Cell< Node * > & tNodes = tMesh->nodes();
            tNodes.set_size( aNodes.size(), nullptr );
            index_t tCount = 0 ;
            for ( Node * tNode : aNodes )
            {
                tNodes( tCount++ ) = new Node( tNode->id(), tNode->x(), tNode->y(), tNode->z() );
            }

            ElementFactory tFactory ;
            index_t n = aNodes.size() - 1;

            Block * tBlock = new Block( 1, n );

            for ( index_t k=0; k<n; ++k )
            {
                Element * tElement = tFactory.create_element( ElementType::LINE2, k+1 );

                tElement->insert_node( tNodes( k ), 0 );
                tElement->insert_node( tNodes( k+1 ), 1 );
                tBlock->insert_element( tElement );
            }

            tMesh->add_block( tBlock );
            tMesh->finalize();
            tMesh->save( "testmesh.vtk" );
            delete tMesh ;
        }


        void
        EdgeCutter::project_temporary_nodes(
            Cell< Facet * > & aFacets,
            Cell< Node * > & aNodes )
        {
            Matrix< real > tXmax, tXmin ;
            tXmax.set_size( 3, aFacets.size());
            tXmin.set_size( 3, aFacets.size() );

            // create the bounding boxes
            index_t tCount = 0 ;
            for ( Facet * tFacet : aFacets )
            {
                mU.fill( -BELFEM_REAL_MAX );
                mV.fill( BELFEM_REAL_MAX );

                Element * tM = tFacet->master() ;
                if ( tM != nullptr )
                {
                    for ( uint k=0; k<tM->number_of_nodes(); ++k )
                    {
                        tM->node( k )->get_coords( mW );
                        for ( uint d=0; d<3; ++d )
                        {
                            if ( mW( d ) > mU( d ) ) mU( d ) = mW( d );
                            if ( mW( d ) < mV( d ) ) mV( d ) = mW( d );
                        }
                    }
                }

                Element * tS = tFacet->slave() ;
                if ( tS != nullptr )
                {
                    for ( uint k=0; k<tS->number_of_nodes(); ++k )
                    {
                        tS->node( k )->get_coords( mW );
                        for ( uint d=0; d<3; ++d )
                        {
                            if ( mW( d ) > mU( d ) ) mU( d ) = mW( d );
                            if ( mW( d ) < mV( d ) ) mV( d ) = mW( d );
                        }
                    }
                }

                mU += BELFEM_EPSILON ;
                mV -= BELFEM_EPSILON ;

                tXmax.set_col( tCount, mU );
                tXmin.set_col( tCount, mV );
                ++tCount ;
            }

            tCount = 0 ;
            for ( Node * tNode : aNodes )
            {
                tNode->get_coords( mP );
                std::cout << "projecting node " << tCount << std::endl ;
                bool tProjected = false ;
                for ( Facet * tFacet : aFacets )
                {
                    if ( ! this->inside_bounding_box( mP, tXmin.col( tFacet->index() ), tXmax.col( tFacet->index() ) ) )
                    {
                        continue ;
                    }

                    // compute the facet plane: normal mW and offset d = n̂·centroid
                    tFacet->node( 1 )->get_coords( mU );
                    tFacet->node( 2 )->get_coords( mV );
                    tFacet->node( 0 )->get_coords( mW );
                    mQ = mU + mV + mW ;
                    mQ /= 3. ;
                    mU -= mW ;
                    mV -= mW ;
                    mW = cross( mU, mV );
                    mW /= norm( mW );
                    real d = dot( mW, mQ ) ;

                    // project mP onto the facet plane n̂·x = d
                    mU = mP - ( dot( mW, mP ) - d ) * mW ;

                    // verify the projection lands inside the facet's parent
                    // element on either side
                    Element * tM = tFacet->master() ;
                    if ( tM != nullptr )
                    {
                        mElementMapper->link( tM );
                        if ( mElementMapper->evaluate( mU, mV ) )
                        {
                            tNode->set_coords( mU );
                            std::cout << "#shift node " << tCount << std::endl ;
                            tNode->set_index( tFacet->index() );
                            tProjected = true ;
                            break ;
                        }
                    }

                    Element * tS = tFacet->slave() ;
                    if ( tS != nullptr )
                    {
                        mElementMapper->link( tS );
                        if ( mElementMapper->evaluate( mU, mV ) )
                        {
                            tNode->set_coords( mU );
                            std::cout << "#shift node " << tCount << std::endl ;
                            tNode->set_index( tFacet->index() );
                            tProjected = true ;
                            break ;
                        }
                    }
                }
                if ( ! tProjected )
                {
                    std::cout << "  no facet matched for node " << tCount << std::endl ;
                }
                ++tCount ;
            }
        }

        bool
        EdgeCutter::inside_bounding_box( const Vector< real > & aPoint, const Vector< real > & aXmin, const Vector< real > & aXmax ) const
        {
            for ( uint d=0; d<3; ++d )
            {
                if ( aPoint( d ) < aXmin( d ) || aPoint( d ) > aXmax( d ) ) return false;
            }
            return true;
        }

    }
}