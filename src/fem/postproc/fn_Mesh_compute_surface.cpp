//
// Created by Christian Messe on 27.06.20.
//

#include "fn_Mesh_compute_volume.hpp"
#include "cl_Cell.hpp"
#include "cl_Element.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"
#include "fn_intpoints.hpp"
#include "meshtools.hpp"
#include "fn_intpoints_auto_integration_order.hpp"
#include "fn_intpoints.hpp"
#include "commtools.hpp"
#include "fn_det.hpp"
#include "fn_sum.hpp"
#include "fn_inv.hpp"
#include "fn_dot.hpp"
#include "geometrytools.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        real
        compute_surface( Mesh * aMesh, const Vector< id_t > & aSideSetIDs )
        {
            BELFEM_ASSERT( aMesh->is_kernel_mesh() || comm_size() == 1,
                          "in parralell mode, input argument of compute_volume must be to Kernel->mesh() " );

            // get my rank
            const proc_t tMyRank = comm_rank() ;

            // create the factory
            fem::InterpolationFunctionFactory tFactory ;

            // total surface
            real aSurface = 0.0 ;

            // nuber of dimensions of mesh
            uint tNumDim = aMesh->number_of_dimensions() ;

            // geometry jacobian
            /*Matrix< real > tJ( tNumDim, tNumDim );

            // Container with surface normals
            Cell< Field * > tNormalFields( tNumDim, nullptr );

            BELFEM_ERROR( aMesh->field_exists( "SurfaceNormalsx"),
                "You need to compute the surface normals before calculating the surface" );

            tNormalFields( 0 ) = aMesh->field( "SurfaceNormalsx" );
            tNormalFields( 1 ) = aMesh->field( "SurfaceNormalsy" );
            if( tNumDim == 3 )
            {
                tNormalFields( 2 ) = aMesh->field( "SurfaceNormalsz" );
            } */

            // loop over all blocks
            for( id_t tID : aSideSetIDs )
            {
                // check if block exists
                if( aMesh->sideset_exists( tID ) )
                {
                    // get pointer to facets
                    Cell< mesh::Facet * > & tFacets = aMesh->sideset( tID )->facets() ;

                    // skip this sideset if it is empty
                    if( tFacets.size() == 0 )
                    {
                        continue;
                    }

                    // get the type of the surface element
                    ElementType tElemType = tFacets( 0 )->element()->type() ;

                    // create the shape function
                    fem::InterpolationFunction * tShape = tFactory.create_lagrange_function( tElemType ) ;

                    // number of nodes for the surface element
                    uint tNumNodes = number_of_nodes( tElemType );

                    // integration points for facet
                    Matrix< real > tXi ;

                    // integreation weights for facet
                    Vector< real > tW ;

                    uint tOrder = auto_integration_order( tElemType ) ;

                    // populate data
                    intpoints(
                            IntegrationScheme::GAUSS,
                            geometry_type( tElemType ),
                            tOrder,
                            tW,
                            tXi );

                    // number of integration points
                    uint tNumPoints = tW.length() ;

                    // number of surface dimensions
                    uint tNumSurfaceDim = dimension( tElemType ) ;

                    // container with derivatives
                    Cell< Matrix< real > > tdNdXi( tNumPoints, Matrix< real >( tNumSurfaceDim, tNumNodes ) );

                    // populate data
                    for( uint k=0; k<tNumPoints; ++k )
                    {
                        tShape->dNdXi( tXi.col( k ), tdNdXi( k ) );
                    }

                    // delete the shape function
                    delete tShape ;

                    // container for node coordinates
                    Matrix< real > tNodeCoords( tNumNodes,  tNumDim );

                    // loop over all elements
                    for( mesh::Facet * tFacet : tFacets )
                    {
                        // check of this element is owned
                        if( tFacet->master()->owner() == tMyRank )
                        {
                            real tS = 0.0 ;

                            // collect node coordinates
                            collect_node_coords( tFacet->element(), tNodeCoords, tNumDim );

                            // loop over all integration points
                            for( uint k=0; k<tNumPoints; ++k )
                            {
                                // compute surface increment
                                real tDetJ = mesh::compute_surface_increment(
                                        tdNdXi( k ), tNodeCoords );

                                // add contribution to surface
                                tS += tW( k ) * tDetJ ;
                            }

                            aSurface += tS ;
                        }// end if own element
                    } // end element loop
                } // end sideset exists
            } // end loop over all selected sideset

            proc_t tCommSize = comm_size() ;

            // check if we are in parallel mode
            if( tCommSize > 1 )
            {
                // check if this is the master
                if ( tMyRank == 0 )
                {
                    Vector< proc_t > tCommList ;

                    // create commlist
                    create_master_commlist( tCommList );

                    // data container
                    Vector< real > tData( tCommSize - 1 );

                    // collect data from other procs
                    receive( tCommList, tData );

                    // add data to volume
                    aSurface += sum( tData );

                    // send data back to other procs
                    tData.fill( aSurface );
                    send( tCommList, tData );

                }
                else
                {
                    // send my contribution to master
                    send( 0, aSurface );

                    // receive sum from master
                    receive( 0, aSurface );
                }
            }

            return aSurface ;
        }

//------------------------------------------------------------------------------
    }
}