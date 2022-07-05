//
// Created by christian on 11/1/21.
//
#include "typedefs.hpp"
#include "commtools.hpp"

#include "fn_FEM_compute_normb.hpp"

namespace belfem
{
    namespace fem
    {
//----------------------------------------------------------------------------

        void
        compute_normb( DofManager * aField, const bool aBiotSavartFlag )
        {
            BELFEM_ERROR( ! ( aBiotSavartFlag && comm_rank() > 0 ),
                         "If the Biot-Savart flag is set, compute_normb() must only be called by the master proc" );

            // get the mesh
            Mesh * tMesh = aField->mesh() ;

            // get the nodes
            Cell< mesh::Node * > & tNodes = tMesh->nodes();

            string tLabel = aBiotSavartFlag ? "normBiotSavart" : "normb" ;
            Vector< real > & tB = tMesh->field_exists( tLabel ) ?
                                  tMesh->field_data(   tLabel ) :
                                  tMesh->create_field( tLabel ) ;

            // get my rank
            proc_t tMyRank = comm_rank();

            if ( aField->mesh()->number_of_dimensions() == 2 )
            {
                Vector< real > & tBx = aBiotSavartFlag ? tMesh->field_data( "BiotSavartx" ) :
                                       tMesh->field_data( "bx" );

                Vector< real > & tBy = aBiotSavartFlag ? tMesh->field_data( "BiotSavarty" ) :
                                       tMesh->field_data( "by" );

                index_t k ;

                if( aBiotSavartFlag || comm_size() == 1 )
                {
                    for ( mesh::Node * tNode : tNodes )
                    {
                        k = tNode->index() ;
                        tB( k ) = std::sqrt( tBx( k ) * tBx( k ) + tBy( k ) * tBy( k ) );
                    }
                }
                else
                {
                    for ( mesh::Node * tNode: tNodes )
                    {
                        if ( tNode->owner() == tMyRank )
                        {
                            k = tNode->index();

                            tB( k ) = std::sqrt( tBx( k ) * tBx( k ) + tBy( k ) * tBy( k ) );
                        }
                    }
                }

            }
            else if ( aField->mesh()->number_of_dimensions() == 3 )
            {
                Vector< real > & tBx = aBiotSavartFlag ? tMesh->field_data( "BiotSavartx" ) :
                                       tMesh->field_data( "bx" );

                Vector< real > & tBy = aBiotSavartFlag ? tMesh->field_data( "BiotSavarty" ) :
                                       tMesh->field_data( "by" );

                Vector< real > & tBz = aBiotSavartFlag ? tMesh->field_data( "BiotSavartz" ) :
                                       tMesh->field_data( "bz" );

                index_t k ;

                if( aBiotSavartFlag || comm_size() == 1 )
                {
                    for ( mesh::Node * tNode: tNodes )
                    {
                        k = tNode->index();

                        tB( k ) = std::sqrt(
                                tBx( k ) * tBx( k )
                                + tBy( k ) * tBy( k )
                                + tBz( k ) * tBz( k ));
                    }
                }
                else
                {
                    for ( mesh::Node * tNode: tNodes )
                    {
                        if ( tNode->owner() == tMyRank )
                        {
                            k = tNode->index();

                            tB( k ) = std::sqrt(
                                    tBx( k ) * tBx( k )
                                + tBy( k ) * tBy( k )
                                    + tBz( k ) * tBz( k ));
                        }
                    }
                }
            }

            if( ! aBiotSavartFlag )
            {
                comm_barrier();
                Cell< string > tLabels( 1, tLabel );
                aField->collect_fields( tLabels );

                // get cuts
                Cell< mesh::SideSet * > & tCuts =  tMesh->cuts() ;

                // fix cuts
                if( tMyRank == 0 && tCuts.size() > 0 )
                {
                    for( mesh::SideSet * tCut : tCuts )
                    {
                        tCut->unflag_all_nodes() ;

                        // loop over all elements
                        for( mesh::Facet * tFacet : tCut->facets() )
                        {
                            // first node
                            mesh::Node * tP = tFacet->element()->node( 0 ) ;

                            if( ! tP->is_flagged() )
                            {
                                mesh::Node * tQ = tFacet->element()->node( 1 );

                                tB( tQ->index() ) = tB( tP->index() );
                                tP->flag() ;
                            }
                        }
                    }

                    // tidy up
                    tMesh->unflag_all_nodes() ;
                }
            }
        }

//----------------------------------------------------------------------------
    }
}