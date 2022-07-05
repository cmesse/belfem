//
// Created by christian on 10/28/21.
//
#include "commtools.hpp"
#include "assert.hpp"
#include "cl_Mesh.hpp"
#include "cl_FEM_Kernel.hpp"
#include "en_FEM_DomainType.hpp"
#include "fn_FEM_compute_biot_savart.hpp"
#include "cl_IWG_Maxwell.hpp"
#include "fn_dot.hpp"
#include "fn_det.hpp"
#include "fn_norm.hpp"
#include "cl_Timer.hpp"
#include "fn_cross.hpp"


namespace belfem
{
    namespace fem
    {

//------------------------------------------------------------------------------

        namespace biotsavart
        {
            void
            identify_current_blocks(  DofManager * aField,
                                      Vector< id_t > & aAirBlocks,
                                      Vector< id_t > & aCoilBlocks,
                                      Vector< id_t > & aScBlocks )
            {
                Vector< id_t > tAllCoilBlocks ;
                Vector< id_t > tAllScBlocks ;

                if( aField->parent()->is_master() )
                {
                    Cell< id_t > tAirBlocks ;
                    Cell< id_t > tCoilBlocks ;
                    Cell< id_t > tScBlocks ;

                    // loop over all selected blocks
                    for( id_t tID : aField->iwg()->selected_blocks() )
                    {
                        // get type of element
                        if( aField->mesh()->block( tID )->number_of_elements() > 0 )
                        {
                            DomainType tType = static_cast< DomainType > (
                                aField->mesh()->block( tID )->element( 0 )->physical_tag() ) ;

                            if( tType == DomainType::Coil )
                            {
                                tCoilBlocks.push( tID );
                            }
                            else if ( tType == DomainType::SuperConductor )
                            {
                                tScBlocks.push( tID );
                            }
                            else if ( tType == DomainType::Air )
                            {
                                tAirBlocks.push( tID );
                            }
                        }
                    }

                    // convert cells to vectors
                    aAirBlocks.set_size( tAirBlocks.size() );
                    for( uint k=0; k<tAirBlocks.size(); ++k )
                    {
                        aAirBlocks( k ) = tAirBlocks( k );
                    }
                    tAllCoilBlocks.set_size( tCoilBlocks.size() );
                    for( uint k=0; k<tCoilBlocks.size(); ++k )
                    {
                        tAllCoilBlocks( k ) = tCoilBlocks( k );
                    }

                    tAllScBlocks.set_size( tScBlocks.size() );
                    for( uint k=0; k<tScBlocks.size(); ++k )
                    {
                        tAllScBlocks( k ) = tScBlocks( k );
                    }

                    comm_barrier() ;
                    send_same( aField->parent()->comm_table(), aAirBlocks );
                    send_same( aField->parent()->comm_table(), tAllCoilBlocks );
                    send_same( aField->parent()->comm_table(), tAllScBlocks );

                }
                else
                {
                    comm_barrier();
                    receive( aField->parent()->master(), aAirBlocks );
                    receive( aField->parent()->master(), tAllCoilBlocks );
                    receive( aField->parent()->master(), tAllScBlocks );
                }

                // check which blocks exist on this proc
                uint tCount = 0 ;
                for( id_t tID : tAllCoilBlocks )
                {
                    if( aField->block_exists( tID ) )
                    {
                        ++tCount ;
                    }
                }

                // allocate memory
                aCoilBlocks.set_size( tCount );

                // populate vector
                tCount = 0 ;
                for( id_t tID : tAllCoilBlocks )
                {
                    if( aField->block_exists( tID ) )
                    {
                        aCoilBlocks( tCount++ ) = tID ;
                    }
                }

                tCount = 0 ;
                for( id_t tID : tAllScBlocks )
                {
                    if( aField->block_exists( tID ) )
                    {
                        ++tCount ;
                    }
                }

                // allocate memory
                aScBlocks.set_size( tCount );

                // populate vector
                tCount = 0 ;
                for( id_t tID : tAllScBlocks )
                {
                    if( aField->block_exists( tID ) )
                    {
                        aScBlocks( tCount++ ) = tID ;
                    }
                }

            }

//------------------------------------------------------------------------------

            void
            get_node_coords(
                    DofManager * aField,
                    Vector< index_t >    & aNodeIndices,
                    Matrix< real >       & aNodeCoords )
            {
                if( aField->parent()->is_master() )
                {
                    Mesh * tMesh = aField->mesh();


                    // get nodes
                    Cell< mesh::Node * > & tNodes = tMesh->nodes();

                    // initialize counter
                    index_t tCount = tNodes.size() ;

                    if( aField->mesh()->number_of_dimensions() == 2 )
                    {
                        // allocate memory
                        aNodeIndices.set_size( tCount );
                        aNodeCoords.set_size( 2, tCount );

                        // reset counter
                        tCount = 0;
                        for ( mesh::Node * tNode: tMesh->nodes())
                        {
                            aNodeCoords( 0, tCount ) = tNode->x();
                            aNodeCoords( 1, tCount ) = tNode->y();
                            aNodeIndices( tCount++ ) = tNode->index();
                        }
                    }
                    else
                    {
                        // allocate memory
                        aNodeIndices.set_size( tCount );
                        aNodeCoords.set_size( 3, tCount );

                        // reset counter
                        tCount = 0;
                        for ( mesh::Node * tNode: tMesh->nodes() )
                        {
                            aNodeCoords( 0, tCount ) = tNode->x();
                            aNodeCoords( 1, tCount ) = tNode->y();
                            aNodeCoords( 2, tCount ) = tNode->z();
                            aNodeIndices( tCount++ ) = tNode->index();
                        }
                    }
                    comm_barrier();
                    send_same( aField->parent()->comm_table(), aNodeCoords );
                }
                else
                {
                    comm_barrier();
                    receive( aField->parent()->master(), aNodeCoords );
                }
            }

//------------------------------------------------------------------------------

            void
            get_sideset_node_coords(
                    DofManager * aField,
                    const Vector< id_t > & aSideSetIDs,
                    Vector< index_t >    & aNodeIndices,
                    Matrix< real >       & aNodeCoords )
            {
                if( aField->parent()->is_master() )
                {
                    Mesh * tMesh = aField->mesh();

                    // unflag all nodes on the mesh
                    tMesh->unflag_all_nodes();

                    // flag all air blocks
                    for ( id_t tID: aSideSetIDs )
                    {
                        tMesh->sideset( tID )->flag_all_nodes() ;
                    }

                    // initialize counter
                    index_t tCount = 0;

                    // get nodes
                    Cell< mesh::Node * > & tNodes = tMesh->nodes();

                    // count nodes
                    for ( mesh::Node * tNode: tNodes )
                    {
                        if ( tNode->is_flagged())
                        {
                            ++tCount;
                        }
                    }

                    if( aField->mesh()->number_of_dimensions() == 2 )
                    {
                        // allocate memory
                        aNodeIndices.set_size( tCount );
                        aNodeCoords.set_size( 2, tCount );

                        // reset counter
                        tCount = 0;
                        for ( mesh::Node * tNode: tMesh->nodes())
                        {
                            if ( tNode->is_flagged())
                            {
                                aNodeCoords( 0, tCount ) = tNode->x();
                                aNodeCoords( 1, tCount ) = tNode->y();
                                aNodeIndices( tCount++ ) = tNode->index();
                            }
                        }
                    }
                    else if( aField->mesh()->number_of_dimensions() == 3 )
                    {
                        // allocate memory
                        aNodeIndices.set_size( tCount );
                        aNodeCoords.set_size( 3, tCount );

                        // reset counter
                        tCount = 0;
                        for ( mesh::Node * tNode: tMesh->nodes())
                        {
                            if ( tNode->is_flagged())
                            {
                                aNodeCoords( 0, tCount ) = tNode->x();
                                aNodeCoords( 1, tCount ) = tNode->y();
                                aNodeCoords( 2, tCount ) = tNode->z();
                                aNodeIndices( tCount++ ) = tNode->index();
                            }
                        }
                    }
                    comm_barrier();
                    send_same( aField->parent()->comm_table(), aNodeCoords );
                }
                else
                {
                    comm_barrier();
                    receive( aField->parent()->master(), aNodeCoords );
                }
            }

//------------------------------------------------------------------------------

            void
            compute_biot_savart_2d(
                    DofManager * aField,
                    const Cell< string >       & aBFields,
                    const Vector< id_t >       & aCoilBlocks,
                    const Vector< id_t >       & aScBlocks,
                    const Vector< index_t >    & aNodeIndices,
                    const Matrix< real >       & aNodeCoords )
            {
                // get the fields
                index_t tNumNodes = aNodeCoords.n_cols() ;
                Matrix< real > tB( 2, tNumNodes, 0.0 );

                Vector< real > tP( 2 );
                Vector< real > tR( 2 );

                IWG_Maxwell * tIWG = reinterpret_cast< IWG_Maxwell * > ( aField->iwg() );

                // loop over all nodes
                for( index_t k=0; k<tNumNodes; ++k )
                {
                    // get the x-coordinate of the node
                    tP( 0 ) = aNodeCoords( 0, k );
                    tP( 1 ) = aNodeCoords( 1, k );

                    // loop over current blocks
                    for( id_t tID : aCoilBlocks )
                    {
                        // get block
                        Block * tBlock  = aField->block( tID ) ;

                        uint tNumNodesPerElement = mesh::number_of_nodes( tBlock->element_type() );
                        Vector< real > tJz( tNumNodesPerElement );

                        Vector< real > tX( tNumNodesPerElement );
                        Vector< real > tY( tNumNodesPerElement );

                        // get integration weights
                        const Vector< real > & tW = tBlock->integration_weights();

                        // get number of integration points
                        uint tNumIntPoints = tW.length() ;

                        real tJ ;
                        real tScale ;

                        // loop over all elements of this block
                        for( Element * tElement : tBlock->elements() )
                        {
                            // get current densities
                            tIWG->collect_node_data( tElement, "jz", tJz );

                            // collect node coords
                            for( uint i=0; i<tNumNodesPerElement; ++i )
                            {
                                tX( i ) = tElement->element()->node( i )->x();
                                tY( i ) = tElement->element()->node( i )->y();
                            }

                            // loop over all integration points
                            for( uint j=0; j<tNumIntPoints; ++j )
                            {
                                // r-vector of point
                                tR( 0 ) = tP( 0 ) - dot( tBlock->n( j ), tX ) ;
                                tR( 1 ) = tP( 1 ) - dot( tBlock->n( j ), tY );

                                // interpolate j
                                tJ = dot( tBlock->n( j ), tJz );

                                // determinant
                                tScale = tW( j ) * std::abs( det( tElement->J( j ) ) ) /
                                         std::pow( norm( tR ), 2 );

                                // assemble contribution
                                tB( 0, k ) -= tJ * tR( 1 ) * tScale ;
                                tB( 1, k ) += tJ * tR( 0 ) * tScale ;
                            }
                        }
                    }

                    // loop over all SC blocks
                    for( id_t tID : aScBlocks )
                    {
                        // get block
                        Block * tBlock  = aField->block( tID ) ;

                        uint tNumNodesPerElement = mesh::number_of_nodes( tBlock->element_type() );
                        Vector< real > tX( tNumNodesPerElement );
                        Vector< real > tY( tNumNodesPerElement );

                        // get integration weights
                        //const Vector< real > & tW = tBlock->integration_weights();

                       // uint tNumIntPoints = tW.length() ;

                        // h-field
                        Vector< real > & tH = tBlock->work_phi() ;

                        Vector< real > tJ( 1 ) ;

                        //real tScale ;

                        // loop over all elements of this block
                        for( Element * tElement : tBlock->elements() )
                        {
                            // get edge field
                            tIWG->collect_edge_data( tElement, "edge_h", tH );

                            // collect node coords
                            for( uint i=0; i<tNumNodesPerElement; ++i )
                            {
                                tX( i ) = tElement->element()->node( i )->x();
                                tY( i ) = tElement->element()->node( i )->y();
                            }

                            BELFEM_ERROR( false, "fix this!");
                            // loop over all integration points
                            /*for( uint j=0; j<tNumIntPoints; ++j )
                            {
                                tIWG->update_nabla( tElement, j );

                                // r-vector of point
                                tR( 0 ) = tP( 0 ) - dot( tBlock->n( j ), tX ) ;
                                tR( 1 ) = tP( 1 ) - dot( tBlock->n( j ), tY );

                                // interpolate j
                                tJ = tIWG->curl_h( tElement, j ) * tH ;

                                // determinant
                                tScale = tW( j ) * tIWG->domain_increment() /
                                         std::pow( norm( tR ), 2 );

                                // assemble contribution
                                tB( 0, k ) -= tJ( 0 ) * tR( 1 ) * tScale ;
                                tB( 1, k ) += tJ( 0 ) * tR( 0 ) * tScale ;
                            }*/

                        }
                    }
                }

                comm_barrier() ;

                // collocate data
                if( aField->parent()->is_master() )
                {
                    Cell< Matrix< real > > tData( comm_size(),Matrix< real >() );

                    receive( aField->parent()->comm_table(), tData );

                    // get the b-fields
                    Vector< real > & tBx =
                            aField->mesh()->field_exists( aBFields( 0 ) ) ?
                            aField->mesh()->field_data( aBFields( 0 ) ) :
                            aField->mesh()->create_field( aBFields( 0 ) );


                    Vector< real > & tBy =
                            aField->mesh()->field_exists( aBFields( 1 ) ) ?
                            aField->mesh()->field_data( aBFields( 1 ) ) :
                            aField->mesh()->create_field( aBFields( 1 )) ;

                    proc_t tCommSize = aField->parent()->comm_table().length() ;

                    // reset data
                    tBx.fill( 0.0 );
                    tBy.fill( 0.0 );

                    for( proc_t p=0; p<tCommSize; ++p )
                    {
                        Matrix< real > & tValues = p == aField->parent()->master() ? tB : tData( p ) ;

                        index_t tCount = 0 ;

                        // loop over all nodes
                        for( index_t k : aNodeIndices )
                        {
                            tBx( k ) += tValues( 0, tCount );
                            tBy( k ) += tValues( 1, tCount++ );
                        }
                    }

                    tBx *= 0.5 * constant::mu0 / constant::pi ;
                    tBy *= 0.5 * constant::mu0 / constant::pi ;
                }
                else
                {
                    send( aField->parent()->master(), tB );
                }

                // wait
                comm_barrier() ;
            }
//------------------------------------------------------------------------------

            void
            compute_biot_savart_3d(
                    DofManager * aField,
                    const Cell< string >       & aBFields,
                    const Vector< id_t >       & aCoilBlocks,
                    const Vector< id_t >       & aScBlocks,
                    const Vector< index_t >    & aNodeIndices,
                    const Matrix< real >       & aNodeCoords )
            {
                // get the fields
                index_t tNumNodes = aNodeCoords.n_cols() ;
                Matrix< real > tB( 3, tNumNodes, 0.0 );

                Vector< real > tP( 3 );
                Vector< real > tR( 3 );
                Vector< real > tJ( 3 );
                Vector< real > tJxR( 3 );

                IWG_Maxwell * tIWG = reinterpret_cast< IWG_Maxwell * > ( aField->iwg() );

                // loop over all nodes
                for( index_t k=0; k<tNumNodes; ++k )
                {
                    // get the x-coordinate of the node
                    tP( 0 ) = aNodeCoords( 0, k );
                    tP( 1 ) = aNodeCoords( 1, k );
                    tP( 2 ) = aNodeCoords( 2, k );

                    // loop over current blocks
                    for( id_t tID : aCoilBlocks )
                    {
                        // get block
                        Block * tBlock  = aField->block( tID ) ;

                        uint tNumNodesPerElement = mesh::number_of_nodes( tBlock->element_type() );

                        Vector< real > tX( tNumNodesPerElement );
                        Vector< real > tY( tNumNodesPerElement );
                        Vector< real > tZ( tNumNodesPerElement );

                        Vector< real > tJx( tNumNodesPerElement );
                        Vector< real > tJy( tNumNodesPerElement );
                        Vector< real > tJz( tNumNodesPerElement );

                        // get integration weights
                        const Vector< real > & tW = tBlock->integration_weights();

                        // get number of integration points
                        uint tNumIntPoints = tW.length() ;

                        real tScale ;

                        // loop over all elements of this block
                        for( Element * tElement : tBlock->elements() )
                        {
                            // get current densities
                            tIWG->collect_node_data( tElement, "jx", tJx );
                            tIWG->collect_node_data( tElement, "jy", tJy );
                            tIWG->collect_node_data( tElement, "jz", tJz );

                            // collect node coords
                            for( uint i=0; i<tNumNodesPerElement; ++i )
                            {
                                tX( i ) = tElement->element()->node( i )->x();
                                tY( i ) = tElement->element()->node( i )->y();
                                tZ( i ) = tElement->element()->node( i )->z();
                            }

                            // loop over all integration points
                            for( uint j=0; j<tNumIntPoints; ++j )
                            {
                                // r-vector of point
                                tR( 0 ) = tP( 0 ) - dot( tBlock->n( j ), tX ) ;
                                tR( 1 ) = tP( 1 ) - dot( tBlock->n( j ), tY );
                                tR( 2 ) = tP( 2 ) - dot( tBlock->n( j ), tZ );

                                // interpolate j
                                tJ( 0 ) = dot( tBlock->n( j ), tJx );
                                tJ( 1 ) = dot( tBlock->n( j ), tJy );
                                tJ( 2 ) = dot( tBlock->n( j ), tJz );

                                // determinant
                                tScale = tW( j ) * std::abs( det( tElement->J( j ) ) ) /
                                         std::pow( norm( tR ), 3 );

                                tJxR = cross( tJ, tR );

                                // assemble contribution
                                tB( 0, k ) += tJxR( 0 ) * tScale ;
                                tB( 1, k ) += tJxR( 1 ) * tScale ;
                                tB( 2, k ) += tJxR( 2 ) * tScale ;
                            }
                        }
                    }

                    // loop over all SC blocks
                    for( id_t tID : aScBlocks )
                    {
                        // get block
                        Block * tBlock  = aField->block( tID ) ;

                        uint tNumNodesPerElement = mesh::number_of_nodes( tBlock->element_type() );
                        Vector< real > tX( tNumNodesPerElement );
                        Vector< real > tY( tNumNodesPerElement );
                        Vector< real > tZ( tNumNodesPerElement );

                        // get integration weights
                        //const Vector< real > & tW = tBlock->integration_weights();

                        //uint tNumIntPoints = tW.length() ;

                        // h-field
                        Vector< real > & tH = tBlock->work_phi() ;

                        //real tScale ;

                        // loop over all elements of this block
                        for( Element * tElement : tBlock->elements() )
                        {
                            // get edge field
                            tIWG->collect_edge_data( tElement, "edge_h", tH );

                            // collect node coords
                            for( uint i=0; i<tNumNodesPerElement; ++i )
                            {
                                tX( i ) = tElement->element()->node( i )->x();
                                tY( i ) = tElement->element()->node( i )->y();
                                tZ( i ) = tElement->element()->node( i )->z();
                            }

                            BELFEM_ERROR( false, "fix this!");

                            // loop over all integration points
                            /*for( uint j=0; j<tNumIntPoints; ++j )
                            {
                                tIWG->update_nabla( tElement, j );

                                // r-vector of point
                                tR( 0 ) = tP( 0 ) - dot( tBlock->n( j ), tX );
                                tR( 1 ) = tP( 1 ) - dot( tBlock->n( j ), tY );
                                tR( 2 ) = tP( 2 ) - dot( tBlock->n( j ), tZ );

                                // interpolate j
                                tJ = tIWG->curl_h( tElement, j ) * tH ;
                                // determinant etc
                                tScale = tW( j ) * tIWG->domain_increment()/
                                         std::pow( norm( tR ), 3 );

                                tJxR = cross( tJ, tR );

                                // assemble contribution
                                tB( 0, k ) += tJxR( 0 ) * tScale ;
                                tB( 1, k ) += tJxR( 1 ) * tScale ;
                                tB( 2, k ) += tJxR( 2 ) * tScale ;
                            } */
                        }
                    }

                }

                comm_barrier() ;

                // collocate data
                if( aField->parent()->is_master() )
                {
                    Cell< Matrix< real > > tData( comm_size(),Matrix< real >() );

                    receive( aField->parent()->comm_table(), tData );

                    // get the b-fields
                    Vector< real > & tBx =
                            aField->mesh()->field_exists( aBFields( 0 ) ) ?
                            aField->mesh()->field_data( aBFields( 0 ) ) :
                            aField->mesh()->create_field( aBFields( 0 ) );


                    Vector< real > & tBy =
                            aField->mesh()->field_exists( aBFields( 1 ) ) ?
                            aField->mesh()->field_data( aBFields( 1 ) ) :
                            aField->mesh()->create_field( aBFields( 1 )) ;

                    Vector< real > & tBz =
                            aField->mesh()->field_exists( aBFields( 2 ) ) ?
                            aField->mesh()->field_data( aBFields( 2 ) ) :
                            aField->mesh()->create_field( aBFields( 2 )) ;

                    proc_t tCommSize = aField->parent()->comm_table().length() ;

                    // reset data
                    tBx.fill( 0.0 );
                    tBy.fill( 0.0 );
                    tBz.fill( 0.0 );

                    for( proc_t p=0; p<tCommSize; ++p )
                    {
                        Matrix< real > & tValues = p == aField->parent()->master() ? tB : tData( p ) ;

                        index_t tCount = 0 ;

                        // loop over all nodes
                        for( index_t k : aNodeIndices )
                        {
                            tBx( k ) += tValues( 0, tCount );
                            tBy( k ) += tValues( 1, tCount );
                            tBz( k ) += tValues( 2, tCount );
                        }
                    }

                    tBx *= 0.25 * constant::mu0 / constant::pi ;
                    tBy *= 0.25 * constant::mu0 / constant::pi ;
                    tBz *= 0.25 * constant::mu0 / constant::pi ;
                }
                else
                {
                    send( aField->parent()->master(), tB );
                }

                // wait
                comm_barrier() ;
            }
        }

//------------------------------------------------------------------------------

        void
        compute_biot_savart( DofManager * aField )
        {
            Timer tTimer ;

            switch( aField->iwg()->type() )
            {
                case( IwgType::MAXWELL_HA_TRI3 ) :
                case( IwgType::MAXWELL_HA_TRI6 ) :
                case( IwgType::MAXWELL_HPHI_TRI3 ) :
                case( IwgType::MAXWELL_HPHI_TRI6 ) :
                {
                    Vector< id_t > tAirBlocks ;
                    Vector< id_t > tCoilBlocks ;
                    Vector< id_t > tScBlocks ;
                    biotsavart::identify_current_blocks(  aField, tAirBlocks, tCoilBlocks, tScBlocks );


                    Vector< index_t > tNodeIndices ;
                    Matrix< real > tNodeCoords ;
                    Cell< string > tFields = { "BiotSavartx", "BiotSavarty" };

                    biotsavart::get_node_coords( aField, tNodeIndices, tNodeCoords );
                    biotsavart::compute_biot_savart_2d( aField, tFields, tCoilBlocks, tScBlocks , tNodeIndices, tNodeCoords );
                    break ;
                }
                case( IwgType::MAXWELL_HA_TET4 ) :
                case( IwgType::MAXWELL_HA_TET10 ) :
                case( IwgType::MAXWELL_HPHI_TET4 ) :
                case( IwgType::MAXWELL_HPHI_TET10 ) :
                {
                    Vector< id_t > tAirBlocks ;
                    Vector< id_t > tCoilBlocks ;
                    Vector< id_t > tScBlocks ;

                    biotsavart::identify_current_blocks(  aField, tAirBlocks, tCoilBlocks, tScBlocks );


                    Vector< index_t > tNodeIndices ;
                    Matrix< real > tNodeCoords ;
                    Cell< string > tFields = { "BiotSavartx", "BiotSavarty", "BiotSavartz" };

                    biotsavart::get_node_coords( aField, tNodeIndices, tNodeCoords );
                    biotsavart::compute_biot_savart_3d( aField, tFields, tCoilBlocks, tScBlocks ,tNodeIndices, tNodeCoords );
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "compute_biot_savart() not implemented fot this IWG type");
                }
            }
            if( comm_rank() == 0 )
            {
                std::cout << "    time for computing Biot-savart: " << tTimer.stop() * 0.001 << " s" << std::endl << std::endl;
            }

        }

//------------------------------------------------------------------------------

    }
}