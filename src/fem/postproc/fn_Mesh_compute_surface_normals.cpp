//
// Created by Christian Messe on 14.06.20.
//
#include "assert.hpp"
#include "fn_Mesh_compute_surface_normals.hpp"
#include "cl_Cell.hpp"
#include "cl_Element.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"
#include "meshtools.hpp"
#include "commtools.hpp"
#include "fn_trans.hpp"
#include "fn_cross.hpp"
#include "cl_Timer.hpp"
#include "cl_Logger.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        void
        compute_surface_normals(
                Mesh * aMesh,
                const Vector< id_t > & aGroupIDs,
                const GroupType aGroupType
        )
        {
            if( comm_rank() == 0 )
            {
                Timer tTimer ;

                switch( aMesh->number_of_dimensions() )
                {
                    case( 2 ) :
                    {
                        normals::compute_surface_normals_2d( aMesh, aGroupIDs, aGroupType );
                        break ;
                    }
                    case( 3 ) :
                    {
                        normals::compute_surface_normals_3d( aMesh, aGroupIDs, aGroupType );
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Mesh must be 2D or 3D in order to compute normals" );
                    }
                }

                message( 2, "    Time for computing surface normals:    %8.1f ms", tTimer.stop() );

                proc_t tCommSize = comm_size() ;

                // the next block is only necessary in parallel mode
                if( tCommSize > 1 )
                {
                    comm_barrier() ;

                    // reset timer
                    tTimer.reset();

                    if( aMesh->number_of_dimensions() == 2 )
                    {
                        normals::send_surface_normals_2d( aMesh );
                    }
                    else if( aMesh->number_of_dimensions() == 3 )
                    {
                        normals::send_surface_normals_3d( aMesh );
                    }

                    comm_barrier() ;
                    message( 2, "    Time for distributing surface normals: %8.1f ms\n", tTimer.stop() );

                }
            }
            else
            {
                if( aMesh->number_of_dimensions() == 2 )
                {
                    normals::receive_surface_normals_2d( aMesh, aGroupIDs, aGroupType );
                }
                else if( aMesh->number_of_dimensions() == 3 )
                {
                    normals::receive_surface_normals_3d( aMesh, aGroupIDs, aGroupType );
                }

                comm_barrier() ;
            }

        }

//------------------------------------------------------------------------------

        namespace normals
        {
//------------------------------------------------------------------------------

            void
            compute_surface_normals_2d(
                    Mesh                 * aMesh,
                    const Vector< id_t > & aGroupIDs,
                    const GroupType        aGroupType
            )
            {
                // create the fields if they have not existed already
                normals::create_normal_fields( aMesh );

                Vector< real > & tNormalX = aMesh->field_data( "SurfaceNormalsx" );
                Vector< real > & tNormalY = aMesh->field_data( "SurfaceNormalsy" );

                index_t tNumberOfNodes = aMesh->number_of_nodes() ;

                aMesh->unflag_all_nodes() ;

                // create the interpolation function factory
                fem::InterpolationFunctionFactory tFactory;

                for ( id_t tID : aGroupIDs )
                {
                    // grab elements form entity
                    Cell< Element * > tElements;
                    normals::collect_elements( aMesh, tID, aGroupType, tElements );

                    // get element type of this group
                    ElementType tType = tElements( 0 )->type();

                    // jump to next if no elements exist here
                    if ( tElements.size() == 0 )
                    {
                        continue;
                    }

                    // get string for error message
                    string tGroupString = aGroupType == GroupType::SIDESET ? "sideset" : "block" ;

                    // make sure that this is a 2D type
                    BELFEM_ERROR(
                            geometry_type( tType ) == GeometryType::LINE,
                            "Error while trying to compute surface normals in %s with ID %lu :\n Element type must be LINE",
                            tGroupString.c_str(), ( long unsigned int ) tID );

                    // get number of nodes per element
                    uint tNumNodesPerElement = number_of_nodes( tType );

                    // create the shape function
                    fem::InterpolationFunction * tShape
                            = tFactory.create_lagrange_function( tType );

                    // get parameter coordinates
                    Matrix< real > tXi;
                    tShape->param_coords( tXi );

                    // Cell with evaluated shape
                    Cell< Matrix< real > > tdNdXi( tNumNodesPerElement, Matrix< real >( 1, tNumNodesPerElement ));

                    // evaluate shape function
                    for( uint k=0; k<tNumNodesPerElement; ++k )
                    {
                        tShape->dNdXi( tXi.col( k ), tdNdXi( k ) );

                        // transpose shape
                        tdNdXi( k ) = trans( tdNdXi( k ) );
                    }

                    // matrix with node coordinates
                    Matrix< real > tNodeCoords( 2, tNumNodesPerElement );

                    // Normal vector
                    Vector< real > tNorm( 3, 0.0 );


                    Matrix< real > tDeriv( 2, 1 );

                    // Loop over all elements
                    for( Element * tElement : tElements )
                    {
                        // collect node coords
                        for( uint k=0; k<tNumNodesPerElement; ++k )
                        {
                            tNodeCoords( 0, k ) = tElement->node( k )->x() ;
                            tNodeCoords( 1, k ) = tElement->node( k )->y() ;
                        }

                        // loop over all nodes of element
                        for( uint k=0; k<tNumNodesPerElement; ++k )
                        {
                            // get index of node
                            index_t tIndex = tElement->node( k )->index() ;

                            // flag this node
                            tElement->node( k )->flag() ;

                            // compute derivative for node
                            //       ( 3 x n ) * ( n x 1 )
                            tDeriv = tNodeCoords * tdNdXi( k ) ;

                            // compute normal
                            tNorm( 0 ) =  tDeriv( 1, 0 ) ;
                            tNorm( 1 ) = -tDeriv( 0, 0 ) ;

                            // add values to fields
                            tNormalX( tIndex ) += tNorm( 0 );
                            tNormalY( tIndex ) += tNorm( 1 );
                        }
                    }

                    delete tShape;
                }

                // get nodes on mesh
                Cell< mesh::Node * > & tNodes = aMesh->nodes() ;

                real tNorm ;

                // loop over all nodes on mesh
                for( index_t k=0; k<tNumberOfNodes; ++k )
                {
                    if( tNodes( k )->is_flagged() )
                    {
                        // compute normal
                        tNorm = std::sqrt(
                                  tNormalX( k ) * tNormalX( k )
                                + tNormalY( k ) * tNormalY( k ) );

                        tNormalX( k ) /= tNorm ;
                        tNormalY( k ) /= tNorm ;
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            compute_surface_normals_3d(
                    Mesh                 * aMesh,
                    const Vector< id_t > & aGroupIDs,
                    const GroupType        aGroupType
            )
            {
                // create the fields if they have not existed already
                normals::create_normal_fields( aMesh );

                Vector< real > & tNormalX = aMesh->field_data( "SurfaceNormalsx" );
                Vector< real > & tNormalY = aMesh->field_data( "SurfaceNormalsy" );
                Vector< real > & tNormalZ = aMesh->field_data( "SurfaceNormalsz" );

                index_t tNumberOfNodes = aMesh->number_of_nodes() ;

                // create the interpolation function factory
                fem::InterpolationFunctionFactory tFactory;

                // loop over all entities
                for ( id_t tID : aGroupIDs )
                {
                    // grab elements form entity
                    Cell< Element * > tElements;
                    normals::collect_elements( aMesh, tID, aGroupType, tElements );

                    // jump to next if no elements exist here
                    if ( tElements.size() == 0 )
                    {
                        continue;
                    }

                    // get element type of this group
                    ElementType tType = tElements( 0 )->type();

                    // get string for error message
                    string tGroupString = aGroupType == GroupType::SIDESET ? "sideset" : "block" ;

                    // make sure that this is a 2D type
                    BELFEM_ERROR(
                            geometry_type( tType ) == GeometryType::TRI ||
                            geometry_type( tType ) == GeometryType::QUAD,
                            "Error while trying to compute surface normals in %s with ID %lu :\n Element type must be TRI or QUAD",
                            tGroupString.c_str(), ( long unsigned int ) tID );

                    // get number of nodes per element
                    uint tNumNodesPerElement = number_of_nodes( tType );

                    // create the shape function
                    fem::InterpolationFunction * tShape
                            = tFactory.create_lagrange_function( tType );

                    // get parameter coordinates
                    Matrix< real > tXi;
                    tShape->param_coords( tXi );

                    // Cell with evaluated shape
                    Cell< Matrix< real > > tdNdXi( tNumNodesPerElement, Matrix< real >( 2, tNumNodesPerElement ));

                    // evaluate shape function
                    for( uint k=0; k<tNumNodesPerElement; ++k )
                    {
                        tShape->dNdXi( tXi.col( k ), tdNdXi( k ) );

                        // transpose shape
                        tdNdXi( k ) = trans( tdNdXi( k ) );
                    }

                    // matrix with node coordinates
                    Matrix< real > tNodeCoords( 3, tNumNodesPerElement );

                    // Normal vector
                    Vector< real > tNorm( 3 );

                    Matrix< real > tDeriv( 3, 2 );

                    // Loop over all elements
                    for( Element * tElement : tElements )
                    {
                        // collect node coords
                        for( uint k=0; k<tNumNodesPerElement; ++k )
                        {
                            tNodeCoords( 0, k ) = tElement->node( k )->x() ;
                            tNodeCoords( 1, k ) = tElement->node( k )->y() ;
                            tNodeCoords( 2, k ) = tElement->node( k )->z() ;
                        }

                        // loop over all nodes of element
                        for( uint k=0; k<tNumNodesPerElement; ++k )
                        {
                            // get index of node
                            index_t tIndex = tElement->node( k )->index() ;

                            // flag this node
                            tElement->node( k )->flag() ;

                            // compute derivative for node
                            //       ( 3 x n ) * ( n x 2 )
                            tDeriv = tNodeCoords * tdNdXi( k ) ;

                            // compute normal
                            tNorm = cross( tDeriv.col( 0 ), tDeriv.col( 1 ) );

                            // add values to fields
                            tNormalX( tIndex ) += tNorm( 0 );
                            tNormalY( tIndex ) += tNorm( 1 );
                            tNormalZ( tIndex ) += tNorm( 2 );
                        }
                    }
                    delete tShape;
                }

                // get nodes on mesh
                Cell< mesh::Node * > & tNodes = aMesh->nodes() ;

                real tNorm ;

                // loop over all nodes on mesh
                for( index_t k=0; k<tNumberOfNodes; ++k )
                {
                    if( tNodes( k )->is_flagged() )
                    {
                        // compute normal
                        tNorm = std::sqrt(
                                  tNormalX( k )*tNormalX( k )
                                + tNormalY( k )*tNormalY( k )
                                + tNormalZ( k )*tNormalZ( k ) );

                        tNormalX( k ) /= tNorm ;
                        tNormalY( k ) /= tNorm ;
                        tNormalZ( k ) /= tNorm ;
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            create_normal_fields( Mesh * aMesh )
            {
                // create fields
                if( ! aMesh->field_exists( "SurfaceNormalsx" ) )
                {
                    aMesh->create_field( "SurfaceNormalsx") ;
                }
                if( ! aMesh->field_exists( "SurfaceNormalsy" ) )
                {
                    aMesh->create_field( "SurfaceNormalsy") ;
                }
                if( ! aMesh->field_exists( "SurfaceNormalsz") && aMesh->number_of_dimensions()==3 )
                {
                    aMesh->create_field( "SurfaceNormalsz") ;
                }
            }

//------------------------------------------------------------------------------

            void
            collect_elements(
                    Mesh * aMesh,
                    const id_t aGroupID,
                    const GroupType aGroupType,
                    Cell< Element * > & aElements )
            {
                switch( aGroupType )
                {
                    case( GroupType::BLOCK ) :
                    {
                        // make sure that block exists
                        if( ! aMesh->block_exists( aGroupID ) )
                        {
                            // reset container
                            aElements.clear() ;

                            // exit this function
                            return ;
                        }

                        // grab Block from Mesh
                        Block * tBlock = aMesh->block( aGroupID );

                        // get number of elements
                        index_t tNumElements = tBlock->number_of_elements() ;

                        // allocate memory
                        aElements.set_size( tNumElements, nullptr );

                        // copy elements into out cell
                        for( index_t e=0; e<tNumElements; ++e )
                        {
                            aElements( e ) = tBlock->element( e );
                        }

                        break;
                    }
                    case( GroupType::SIDESET ):
                    {
                        // make sure that block exists
                        if( ! aMesh->sideset_exists( aGroupID ) )
                        {
                            // reset container
                            aElements.clear() ;

                            // exit this function
                            return ;
                        }

                        // grab sideset
                        SideSet * tSideSet = aMesh->sideset( aGroupID );

                        // get number of elements
                        index_t tNumElements = tSideSet->number_of_facets() ;

                        // allocate memory
                        aElements.set_size( tNumElements, nullptr );

                        // copy elements into out cell
                        for( index_t e=0; e<tNumElements; ++e )
                        {
                            aElements( e ) = tSideSet->facet_by_index( e )->element() ;
                        }

                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Invalid group type");
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            send_surface_normals_2d( Mesh * aMesh )
            {
                BELFEM_ASSERT( comm_rank() == 0, "send_surface_normals_2d() must be called by master proc");

                // create the communication list
                Vector< proc_t > tCommList ;
                create_commlist( tCommList );
                uint tNumProcs = tCommList.length() ;

                // initialize container with node IDs
                Cell< Vector< id_t > > tNodeIDs( tNumProcs, Vector< id_t >() );

                // get data
                receive( tCommList, tNodeIDs );

                Cell< Vector< real > > tData( tNumProcs, Vector< real >() );

                // allocate containers
                for( uint p=0; p<tNumProcs; ++p )
                {
                    tData( p ).set_size( tNodeIDs( p ).length() );
                }

                // link data to mesh
                Vector< real > & tNormalX = aMesh->field_data( "SurfaceNormalsx" );
                Vector< real > & tNormalY = aMesh->field_data( "SurfaceNormalsy" );

                // X-Values
                for( uint p=0; p<tNumProcs; ++p )
                {
                    // reset counter
                    index_t tCount = 0 ;


                    Vector< real > & tNx = tData( p );
                    for( id_t tID : tNodeIDs( p ) )
                    {
                        tNx( tCount++ ) = tNormalX( aMesh->node( tID )->index() );
                    }
                }
                send( tCommList, tData );

                // Y-Values
                for( uint p=0; p<tNumProcs; ++p )
                {
                    // reset counter
                    index_t tCount = 0 ;

                    Vector< real > & tNy = tData( p );
                    for( id_t tID : tNodeIDs( p ) )
                    {
                        tNy( tCount++ ) = tNormalY( aMesh->node( tID )->index() );
                    }
                }
                send( tCommList, tData );
            }

//------------------------------------------------------------------------------

            void
            send_surface_normals_3d( Mesh * aMesh )
            {
                BELFEM_ASSERT( comm_rank() == 0, "send_surface_normals_3d() must be called by master proc");

                // create the communication list
                Vector< proc_t > tCommList ;
                create_commlist( tCommList );
                uint tNumProcs = tCommList.length() ;

                // initialize container with node IDs
                Cell< Vector< id_t > > tNodeIDs( tNumProcs, Vector< id_t >() );

                // get data
                receive( tCommList, tNodeIDs );

                Cell< Vector< real > > tData( tNumProcs, Vector< real >() );

                // allocate containers
                for( uint p=0; p<tNumProcs; ++p )
                {
                    tData( p ).set_size( tNodeIDs( p ).length() );
                }

                // link data to mesh
                Vector< real > & tNormalX = aMesh->field_data( "SurfaceNormalsx" );
                Vector< real > & tNormalY = aMesh->field_data( "SurfaceNormalsy" );
                Vector< real > & tNormalZ = aMesh->field_data( "SurfaceNormalsz" );


                // X-Values
                for( uint p=0; p<tNumProcs; ++p )
                {
                    // reset counter
                    index_t tCount = 0 ;

                    Vector< real > & tNx = tData( p );

                    for( id_t tID : tNodeIDs( p ) )
                    {
                        tNx( tCount++ ) = tNormalX( aMesh->node( tID )->index() );
                    }
                }
                send( tCommList, tData );

                // Y-Values
                for( uint p=0; p<tNumProcs; ++p )
                {
                    // reset counter
                    index_t tCount = 0 ;

                    Vector< real > & tNy = tData( p );
                    for( id_t tID : tNodeIDs( p ) )
                    {
                        tNy( tCount++ ) = tNormalY( aMesh->node( tID )->index() );
                    }
                }
                send( tCommList, tData );

                // Z-Values
                for( uint p=0; p<tNumProcs; ++p )
                {
                    // reset counter
                    index_t tCount = 0 ;

                    Vector< real > & tNz = tData( p );

                    for( id_t tID : tNodeIDs( p ) )
                    {
                        tNz( tCount++ ) = tNormalZ( aMesh->node( tID )->index() );
                    }
                }
                send( tCommList, tData );
            }

//------------------------------------------------------------------------------

            void
            receive_surface_normals_2d( Mesh * aMesh,
                                        const Vector< id_t > & aGroupIDs,
                                        const GroupType aGroupType )
            {
                BELFEM_ASSERT( comm_rank() != 0, "receive_surface_normals_2d() must not be called by master proc");

                // create fields
                normals::create_normal_fields( aMesh );

                // link data to mesh
                Vector< real > & tNormalX = aMesh->field_data( "SurfaceNormalsx" );
                Vector< real > & tNormalY = aMesh->field_data( "SurfaceNormalsy" );

                // select nodes the information is needed for
                Vector< id_t > tNodeIDs ;
                get_node_ids( aMesh,  aGroupIDs, aGroupType, tNodeIDs);

                // wait for master to be ready
                comm_barrier() ;

                // send request list to master
                send( 0, tNodeIDs );

                // vector with X-coordinates
                Vector< real > tX ;
                receive( 0, tX );
                Vector< real > tY ;
                receive( 0, tY );

                // reset field
                tNormalX.fill( 0.0 );
                tNormalY.fill( 0.0 );

                // write data on mesh
                index_t tCount = 0 ;
                for( Node * tNode : aMesh->nodes() )
                {
                    if( tNode->is_flagged() )
                    {
                        index_t tIndex = tNode->index() ;
                        tNormalX( tIndex ) = tX( tCount );
                        tNormalY( tIndex ) = tY( tCount );

                        ++tCount ;
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            receive_surface_normals_3d( Mesh * aMesh,
                                        const Vector< id_t > & aGroupIDs,
                                        const GroupType aGroupType )
            {
                BELFEM_ASSERT( comm_rank() != 0, "receive_surface_normals_3d() must not be called by master proc");

                // create fields
                normals::create_normal_fields( aMesh );

                // link data to mesh
                Vector< real > & tNormalX = aMesh->field_data( "SurfaceNormalsx" );
                Vector< real > & tNormalY = aMesh->field_data( "SurfaceNormalsy" );
                Vector< real > & tNormalZ = aMesh->field_data( "SurfaceNormalsz" );

                // select nodes the information is needed for
                Vector< id_t > tNodeIDs ;
                get_node_ids( aMesh,  aGroupIDs, aGroupType, tNodeIDs);

                // wait for master to be ready
                comm_barrier() ;

                // send request list to master
                send( 0, tNodeIDs );

                // vector with X-coordinates
                Vector< real > tX ;
                receive( 0, tX );
                Vector< real > tY ;
                receive( 0, tY );
                Vector< real > tZ ;
                receive( 0, tZ );

                // reset field
                tNormalX.fill( 0.0 );
                tNormalY.fill( 0.0 );
                tNormalZ.fill( 0.0 );

                // write data on mesh
                index_t tCount = 0 ;
                for( Node * tNode : aMesh->nodes() )
                {
                    if( tNode->is_flagged() )
                    {
                        index_t tIndex = tNode->index() ;
                        tNormalX( tIndex ) = tX( tCount );
                        tNormalY( tIndex ) = tY( tCount );
                        tNormalZ( tIndex ) = tZ( tCount );

                        ++tCount ;
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            create_commlist( Vector< proc_t > & aCommList )
            {
                proc_t tCommSize = comm_size() - 1 ;

                aCommList.set_size( tCommSize );

                // create communication list
                for( proc_t p=0; p<tCommSize; ++p )
                {
                    aCommList( p ) = p + 1 ;
                }
            }

//------------------------------------------------------------------------------

            void
            get_node_ids( Mesh * aMesh,
                          const Vector< id_t > & aGroupIDs,
                          const GroupType aGroupType,
                          Vector< id_t > & aNodeIDs )
            {
                // fix node indices
                aMesh->update_node_indices() ;

                // unflag all nodes on mesh
                aMesh->unflag_all_nodes() ;

                // flag all nodes
                for ( id_t tID : aGroupIDs )
                {
                    // grab elements form entity
                    Cell< Element * > tElements;
                    normals::collect_elements( aMesh, tID, aGroupType, tElements );

                    for( Element * tElement : tElements )
                    {
                        tElement->flag_nodes() ;
                    }
                }

                // count flagged nodes
                index_t tCount = 0 ;

                for( Node * tNode : aMesh->nodes() )
                {
                    if( tNode->is_flagged() )
                    {
                        // increment node counter
                        ++tCount ;
                    }
                }

                // collect node IDs
                aNodeIDs.set_size( tCount );
                tCount = 0 ;
                for( Node * tNode : aMesh->nodes() )
                {
                    if( tNode->is_flagged() )
                    {
                        aNodeIDs( tCount++ ) = tNode->id() ;
                    }
                }
            }

//------------------------------------------------------------------------------
        }
    }
}