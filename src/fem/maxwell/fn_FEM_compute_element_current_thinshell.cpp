//
// Created by christian on 8/1/22.
//
#include "assert.hpp"
#include "fn_FEM_compute_element_current_thinshell.hpp"
#include "cl_IWG_Maxwell.hpp"

namespace belfem
{
    namespace fem
    {
        void
        compute_element_current_thinshell_tri3( DofManager * aMagfield  )
        {

            // get equation

            BELFEM_ERROR(  aMagfield->iwg()->type() == IwgType::MAXWELL_HPHI_TRI3 ||
                           aMagfield->iwg()->type() == IwgType::MAXWELL_HA_TRI3,
                           "this function requires a tri3 base type iwg" );

            // get equation
            // IWG_Maxwell * tIWG = reinterpret_cast< IWG_Maxwell * >( aMagfield->iwg() );

            // get mesh
            Mesh * tMesh = aMagfield->mesh() ;

            // number of ghost blocks
            uint tNumGhostSideSets = tMesh->ghost_sideset_ids().length();

            if( tNumGhostSideSets < 2 )
            {
                return ;
            }
            if( comm_rank() > 0 )
            {
                comm_barrier() ;
                return ;
            }

            // collect blocks with ghost blocks
            uint tNumGhostBlocks = tNumGhostSideSets - 1 ;
            uint tNumBlocks = tMesh->number_of_blocks() ;
            Cell< mesh::Block * > tBlocks( tNumGhostBlocks, nullptr );
            uint tOff = tNumBlocks - tNumGhostBlocks ;
            for( uint b=0; b<tNumGhostBlocks; ++b )
            {
                tBlocks( b ) = tMesh->blocks()( tOff++ );
            }

            // collect ghost sidesets
            Cell< mesh::SideSet * > tSideSets( tNumGhostSideSets, nullptr );
            for( uint s=0; s<tNumGhostSideSets; ++s )
            {
                tSideSets( s ) = tMesh->sideset( tMesh->ghost_sideset_ids()( s ) );
            }

            // create the element field
            Vector< real >& tJ = tMesh->field_exists( "elementJ") ?
                   tMesh->field_data( "elementJ") : tMesh->create_field( "elementJ", EntityType::ELEMENT );

            Vector< real > & tH = tMesh->field_data( "edge_h");

            // get number of elements
            index_t tNumElems = tBlocks( 0 )->number_of_elements() ;


            // loop over all layers
            for( uint l=0; l<tNumGhostBlocks; ++l )
            {
                // get block
                Cell< mesh::Element * > & tElements = tBlocks( l )->elements() ;

                // get lower sideset
                mesh::SideSet * tBottom = tSideSets( l );

                // get upper sideset
                mesh::SideSet * tTop = tSideSets( l + 1 );

                // compute thickness
                real tDX = tElements( 0 )->node( 0 )->x()
                          -tElements( 0 )->node( 3 )->x() ;

                real tDY =  tElements( 0 )->node( 0 )->y()
                           -tElements( 0 )->node( 3 )->y() ;


                real tT = std::sqrt( tDX *tDX + tDY * tDY ) ;

                for( index_t e=0; e<tNumElems; ++e )
                {
                    real tH0 = tH( tBottom->facet_by_index( e )->edge(0)->index()  );
                    real tH1 = tH( tTop->facet_by_index( e )->edge(0)->index()  );

                    tJ( tElements( e )->index() ) = ( tH1 - tH0 ) / tT ;
                }

            }

            comm_barrier() ;
        }
    }
}
