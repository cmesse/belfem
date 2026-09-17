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
#include "fn_sum.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        void
        compute_surface_normals(
                Mesh * aMesh,
                const Vector< id_t > & aGroupIDs,
                const GroupType aGroupType,
                const bool aComputeInParallel
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

                message( InfoLevel::Default, "    Time for computing surface normals:    %u ms", tTimer.stop() );
                proc_t tCommSize = comm_size() ;

                // the next block is only necessary in parallel mode
                if( tCommSize > 1 && aComputeInParallel )
                {
                    comm_barrier() ;

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
                    message( InfoLevel::Default, "    Time for distributing surface normals: %u ms\n", tTimer.stop() );

                }
            }
            else if ( aComputeInParallel )
            {
                comm_barrier() ;

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
                normals::create_normal_fields( aMesh );

                Vector< real > & tNormalX = aMesh->field_data( "SurfaceNormalsx" );
                Vector< real > & tNormalY = aMesh->field_data( "SurfaceNormalsy" );

                aMesh->unflag_all_nodes() ;

                fem::InterpolationFunctionFactory tFactory;

                for ( id_t tID : aGroupIDs )
                {
                    Cell< Element * > tElements;
                    normals::collect_elements( aMesh, tID, aGroupType, tElements );

                    ElementType tType = tElements( 0 )->type();

                    if ( tElements.size() == 0 )
                    {
                        continue;
                    }

                    bool tIsSideSet = aGroupType == GroupType::SIDESET ;
                    string tGroupString = tIsSideSet ? "sideset" : "block" ;

                    BELFEM_ERROR(
                            geometry_type( tType ) == GeometryType::LINE,
                            "Error while trying to compute surface normals in %s with ID %lu :\n Element type must be LINE",
                            tGroupString.c_str(), ( long unsigned int ) tID );

                    uint tNumNodesPerElement = number_of_nodes( tType );

                    fem::InterpolationFunction * tShape
                            = tFactory.create_lagrange_function( tType );

                    Matrix< real > tXi;
                    tShape->param_coords( tXi );

                    Cell< Matrix< real > > tdNdXi( tNumNodesPerElement, Matrix< real >( 1, tNumNodesPerElement ));

                    real tXm = BELFEM_QUIET_NAN ;
                    real tYm = BELFEM_QUIET_NAN ;

                    Vector< real > tXmaster ;
                    Vector< real > tYmaster ;

                    uint tNumNodesOnMaster = 0 ;
                    if( tIsSideSet )
                    {
                        if( aMesh->sideset( tID )->number_of_facets() > 0 )
                        {
                            tNumNodesOnMaster = aMesh->sideset( tID )->facets()( 0 )->master()->number_of_nodes() ;
                            tXmaster.set_size( tNumNodesOnMaster );
                            tYmaster.set_size( tNumNodesOnMaster );
                        }
                    }

                    for( uint k=0; k<tNumNodesPerElement; ++k )
                    {
                        tShape->dNdXi( tXi.col( k ), tdNdXi( k ) );

                        tdNdXi( k ) = trans( tdNdXi( k ) );
                    }

                    Matrix< real > tNodeCoords( 2, tNumNodesPerElement );

                    Vector< real > tNorm( 3, 0.0 );

                    Matrix< real > tDeriv( 2, 1 );

                    if( tIsSideSet )
                    {
                        for ( Element * tElement: tElements )
                        {
                            for ( uint k = 0; k < tNumNodesPerElement; ++k )
                            {
                                tNodeCoords( 0, k ) = tElement->node( k )->x();
                                tNodeCoords( 1, k ) = tElement->node( k )->y();
                            }

                            Element * tMaster = aMesh->facet( tElement->id())->master();

                            for ( uint k = 0; k < tNumNodesOnMaster; ++k )
                            {
                                tXmaster( k ) = tMaster->node( k )->x();
                                tYmaster( k ) = tMaster->node( k )->y();
                            }

                            tXm = sum( tXmaster ) / tNumNodesOnMaster;
                            tYm = sum( tYmaster ) / tNumNodesOnMaster;

                            for ( uint k = 0; k < tNumNodesPerElement; ++k )
                            {
                                index_t tIndex = tElement->node( k )->index();

                                tElement->node( k )->flag();

                                //       ( 2 x n ) * ( n x 1 )
                                tDeriv = tNodeCoords * tdNdXi( k );

                                tNorm( 0 ) = tDeriv( 1, 0 );
                                tNorm( 1 ) = -tDeriv( 0, 0 );

                                // check sign.
                                // let g = X + a * tDeriv
                                //     h = Xm + b * n

                                if( (  tNorm( 0 ) * ( tNodeCoords( 0, k ) - tXm )
                                     + tNorm( 1 ) * ( tNodeCoords( 1, k ) - tYm ) ) > 0 )
                                {
                                    tNormalX( tIndex ) += tNorm( 0 );
                                    tNormalY( tIndex ) += tNorm( 1 );
                                }
                                else
                                {
                                    // add values to fields, but flip sign
                                    tNormalX( tIndex ) -= tNorm( 0 );
                                    tNormalY( tIndex ) -= tNorm( 1 );
                                }
                            }

                        }
                    }
                    else
                    {
                        for ( Element * tElement: tElements )
                        {
                            for ( uint k = 0; k < tNumNodesPerElement; ++k )
                            {
                                tNodeCoords( 0, k ) = tElement->node( k )->x();
                                tNodeCoords( 1, k ) = tElement->node( k )->y();
                            }

                            for ( uint k = 0; k < tNumNodesPerElement; ++k )
                            {
                                index_t tIndex = tElement->node( k )->index();

                                tElement->node( k )->flag();

                                //       ( 2 x n ) * ( n x 1 )
                                tDeriv = tNodeCoords * tdNdXi( k );

                                tNorm( 0 ) = tDeriv( 1, 0 );
                                tNorm( 1 ) = -tDeriv( 0, 0 );

                                tNormalX( tIndex ) += tNorm( 0 );
                                tNormalY( tIndex ) += tNorm( 1 );
                            }
                        }
                    }
                    delete tShape;
                }

                Cell< mesh::Node * > & tNodes = aMesh->nodes() ;

                real tNorm ;

                index_t tIndex = 0 ;

                for(  mesh::Node * tNode : tNodes )
                {
                    if( tNode->is_flagged() )
                    {
                        tNorm = std::sqrt(
                                tNormalX( tIndex ) * tNormalX( tIndex )
                                + tNormalY( tIndex ) * tNormalY( tIndex ) );

                        tNormalX( tIndex ) /= tNorm ;
                        tNormalY( tIndex ) /= tNorm ;

                        if( tNode->number_of_duplicates() > 0 )
                        {
                            index_t tOrgIndex = tNode->original()->index() ;
                            tNormalX( tOrgIndex ) = tNormalX( tIndex );
                            tNormalY( tOrgIndex ) = tNormalY( tIndex );
                        }
                    }
                    ++tIndex ;
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
                normals::create_normal_fields( aMesh );

                Vector< real > & tNormalX = aMesh->field_data( "SurfaceNormalsx" );
                Vector< real > & tNormalY = aMesh->field_data( "SurfaceNormalsy" );
                Vector< real > & tNormalZ = aMesh->field_data( "SurfaceNormalsz" );

                fem::InterpolationFunctionFactory tFactory;

                for ( id_t tID : aGroupIDs )
                {
                    Cell< Element * > tElements;
                    normals::collect_elements( aMesh, tID, aGroupType, tElements );

                    if ( tElements.size() == 0 )
                    {
                        continue;
                    }

                    ElementType tType = tElements( 0 )->type();

                    string tGroupString = aGroupType == GroupType::SIDESET ? "sideset" : "block" ;

                    BELFEM_ERROR(
                            geometry_type( tType ) == GeometryType::TRI ||
                            geometry_type( tType ) == GeometryType::QUAD,
                            "Error while trying to compute surface normals in %s with ID %lu :\n Element type must be TRI or QUAD",
                            tGroupString.c_str(), ( long unsigned int ) tID );

                    uint tNumNodesPerElement = number_of_nodes( tType );

                    fem::InterpolationFunction * tShape
                            = tFactory.create_lagrange_function( tType );

                    Matrix< real > tXi;
                    tShape->param_coords( tXi );

                    Cell< Matrix< real > > tdNdXi( tNumNodesPerElement, Matrix< real >( 2, tNumNodesPerElement ));

                    for( uint k=0; k<tNumNodesPerElement; ++k )
                    {
                        tShape->dNdXi( tXi.col( k ), tdNdXi( k ) );

                        tdNdXi( k ) = trans( tdNdXi( k ) );
                    }

                    Matrix< real > tNodeCoords( 3, tNumNodesPerElement );

                    Vector< real > tNorm( 3 );

                    Matrix< real > tDeriv( 3, 2 );

                    for( Element * tElement : tElements )
                    {
                        for( uint k=0; k<tNumNodesPerElement; ++k )
                        {
                            tNodeCoords( 0, k ) = tElement->node( k )->x() ;
                            tNodeCoords( 1, k ) = tElement->node( k )->y() ;
                            tNodeCoords( 2, k ) = tElement->node( k )->z() ;
                        }

                        for( uint k=0; k<tNumNodesPerElement; ++k )
                        {
                            index_t tIndex = tElement->node( k )->index() ;

                            tElement->node( k )->flag() ;

                            //       ( 3 x n ) * ( n x 2 )
                            tDeriv = tNodeCoords * tdNdXi( k ) ;

                            tNorm = cross( tDeriv.col( 0 ), tDeriv.col( 1 ) );

                            tNormalX( tIndex ) += tNorm( 0 );
                            tNormalY( tIndex ) += tNorm( 1 );
                            tNormalZ( tIndex ) += tNorm( 2 );
                        }
                    }
                    delete tShape;
                }

                Cell< mesh::Node * > & tNodes = aMesh->nodes() ;

                real tNorm ;

                index_t tIndex = 0 ;

                for(  mesh::Node * tNode : tNodes )
                {
                    if( tNodes( tIndex )->is_flagged() )
                    {
                        tNorm = std::sqrt(
                                  tNormalX( tIndex )*tNormalX( tIndex )
                                + tNormalY( tIndex )*tNormalY( tIndex )
                                + tNormalZ( tIndex )*tNormalZ( tIndex ) );

                        tNormalX( tIndex ) /= tNorm ;
                        tNormalY( tIndex ) /= tNorm ;
                        tNormalZ( tIndex ) /= tNorm ;

                        if( tNode->number_of_duplicates() > 0 )
                        {
                            index_t tOrgIndex = tNode->original()->index() ;
                            tNormalX( tOrgIndex ) = tNormalX( tIndex );
                            tNormalY( tOrgIndex ) = tNormalY( tIndex );
                            tNormalZ( tOrgIndex ) = tNormalZ( tIndex );
                        }
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            create_normal_fields( Mesh * aMesh )
            {
                if( ! aMesh->field_exists( "SurfaceNormalsx" ) )
                {
                    aMesh->create_field( "SurfaceNormalsx") ;
                }
                if( ! aMesh->field_exists( "SurfaceNormalsy" ) )
                {
                    aMesh->create_field( "SurfaceNormalsy") ;
                }
                if( ! aMesh->field_exists( "SurfaceNormalsz") )
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
                        if( ! aMesh->block_exists( aGroupID ) )
                        {
                            aElements.clear() ;

                            return ;
                        }

                        Block * tBlock = aMesh->block( aGroupID );

                        index_t tNumElements = tBlock->number_of_elements() ;

                        aElements.set_size( tNumElements, nullptr );

                        for( index_t e=0; e<tNumElements; ++e )
                        {
                            aElements( e ) = tBlock->element( e );
                        }

                        break;
                    }
                    case( GroupType::SIDESET ):
                    {
                        if( ! aMesh->sideset_exists( aGroupID ) )
                        {
                            aElements.clear() ;

                            return ;
                        }

                        SideSet * tSideSet = aMesh->sideset( aGroupID );

                        index_t tNumElements = tSideSet->number_of_facets() ;

                        aElements.set_size( tNumElements, nullptr );

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

                uint tNumProcs = comm_size() ;

                Cell< Vector< id_t > > tNodeIDs( tNumProcs, Vector< id_t >() );

                collect( tNodeIDs );

                Cell< Vector< real > > tData( tNumProcs, Vector< real >() );

                for( uint p=0; p<tNumProcs; ++p )
                {
                    tData( p ).set_size( tNodeIDs( p ).length() );
                }

                Vector< real > & tNormalX = aMesh->field_data( "SurfaceNormalsx" );
                Vector< real > & tNormalY = aMesh->field_data( "SurfaceNormalsy" );

                for( uint p=0; p<tNumProcs; ++p )
                {
                    index_t tCount = 0 ;

                    Vector< real > & tNx = tData( p );
                    for( id_t tID : tNodeIDs( p ) )
                    {
                        tNx( tCount++ ) = tNormalX( aMesh->node( tID )->index() );
                    }
                }
                distribute( tData );

                for( uint p=0; p<tNumProcs; ++p )
                {
                    index_t tCount = 0 ;

                    Vector< real > & tNy = tData( p );
                    for( id_t tID : tNodeIDs( p ) )
                    {
                        tNy( tCount++ ) = tNormalY( aMesh->node( tID )->index() );
                    }
                }
                distribute( tData );
            }

//------------------------------------------------------------------------------

            void
            send_surface_normals_3d( Mesh * aMesh )
            {
                BELFEM_ASSERT( comm_rank() == 0, "send_surface_normals_3d() must be called by master proc");

                uint tNumProcs = comm_size();

                Cell< Vector< id_t > > tNodeIDs( tNumProcs, Vector< id_t >() );

                collect( tNodeIDs );

                Cell< Vector< real > > tData( tNumProcs, Vector< real >() );

                for( uint p=0; p<tNumProcs; ++p )
                {
                    tData( p ).set_size( tNodeIDs( p ).length() );
                }

                Vector< real > & tNormalX = aMesh->field_data( "SurfaceNormalsx" );
                Vector< real > & tNormalY = aMesh->field_data( "SurfaceNormalsy" );
                Vector< real > & tNormalZ = aMesh->field_data( "SurfaceNormalsz" );

                for( uint p=0; p<tNumProcs; ++p )
                {
                    index_t tCount = 0 ;

                    Vector< real > & tNx = tData( p );

                    for( id_t tID : tNodeIDs( p ) )
                    {
                        tNx( tCount++ ) = tNormalX( aMesh->node( tID )->index() );
                    }
                }
                distribute( tData );

                for( uint p=0; p<tNumProcs; ++p )
                {
                    index_t tCount = 0 ;

                    Vector< real > & tNy = tData( p );
                    for( id_t tID : tNodeIDs( p ) )
                    {
                        tNy( tCount++ ) = tNormalY( aMesh->node( tID )->index() );
                    }
                }
                distribute( tData );

                for( uint p=0; p<tNumProcs; ++p )
                {
                    index_t tCount = 0 ;

                    Vector< real > & tNz = tData( p );

                    for( id_t tID : tNodeIDs( p ) )
                    {
                        tNz( tCount++ ) = tNormalZ( aMesh->node( tID )->index() );
                    }
                }
                distribute( tData );
            }

//------------------------------------------------------------------------------

            void
            receive_surface_normals_2d( Mesh * aMesh,
                                        const Vector< id_t > & aGroupIDs,
                                        const GroupType aGroupType )
            {
                BELFEM_ASSERT( comm_rank() != 0, "receive_surface_normals_2d() must not be called by master proc");

                normals::create_normal_fields( aMesh );

                Vector< real > & tNormalX = aMesh->field_data( "SurfaceNormalsx" );
                Vector< real > & tNormalY = aMesh->field_data( "SurfaceNormalsy" );

                Vector< id_t > tNodeIDs ;
                get_node_ids( aMesh,  aGroupIDs, aGroupType, tNodeIDs);

                comm_barrier() ;

                send( tNodeIDs );

                Vector< real > tX ;
                receive( tX );
                Vector< real > tY ;
                receive( tY );

                tNormalX.fill( 0.0 );
                tNormalY.fill( 0.0 );

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

                normals::create_normal_fields( aMesh );

                Vector< real > & tNormalX = aMesh->field_data( "SurfaceNormalsx" );
                Vector< real > & tNormalY = aMesh->field_data( "SurfaceNormalsy" );
                Vector< real > & tNormalZ = aMesh->field_data( "SurfaceNormalsz" );

                Vector< id_t > tNodeIDs ;
                get_node_ids( aMesh,  aGroupIDs, aGroupType, tNodeIDs);

                comm_barrier() ;

                send( tNodeIDs );

                Vector< real > tX ;
                receive( tX );

                Vector< real > tY ;
                receive( tY );

                Vector< real > tZ ;
                receive( tZ );

                tNormalX.fill( 0.0 );
                tNormalY.fill( 0.0 );
                tNormalZ.fill( 0.0 );

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
            get_node_ids( Mesh * aMesh,
                          const Vector< id_t > & aGroupIDs,
                          const GroupType aGroupType,
                          Vector< id_t > & aNodeIDs )
            {
                aMesh->update_node_indices() ;

                aMesh->unflag_all_nodes() ;

                for ( id_t tID : aGroupIDs )
                {
                    Cell< Element * > tElements;
                    normals::collect_elements( aMesh, tID, aGroupType, tElements );

                    for( Element * tElement : tElements )
                    {
                        tElement->flag_nodes() ;
                    }
                }

                index_t tCount = 0 ;

                for( Node * tNode : aMesh->nodes() )
                {
                    if( tNode->is_flagged() )
                    {
                        ++tCount ;
                    }
                }

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
