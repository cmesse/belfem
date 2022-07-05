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
#include "commtools.hpp"
#include "fn_det.hpp"
#include "fn_sum.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        real
        compute_volume( Mesh * aMesh, const id_t  aBlockID )
        {
            return compute_volume( aMesh ,
                                   Vector< id_t >( 1, aBlockID ) );
        }

//------------------------------------------------------------------------------

    real
    compute_volume( Mesh * aMesh, const Vector< id_t > & aBlockIDs )
    {
        BELFEM_ASSERT( aMesh->is_kernel_mesh() || comm_size() == 1,
                      "in parralell mode, input argument of compute_volume must be to Kernel->mesh() " );

        // get my rank
        const proc_t tMyRank = comm_rank() ;

        // create the factory
        fem::InterpolationFunctionFactory tFactory ;

        // total volume
        real aVolume = 0.0 ;

        // nuber of dimensions of mesh
        uint tNumDim = aMesh->number_of_dimensions() ;

        // geometry jacobian
        Matrix< real > tJ( tNumDim, tNumDim );

        // loop over all blocks
        for( id_t tID : aBlockIDs )
        {
            // check if block exists
            if( aMesh->block_exists( tID ) )
            {
                // get pointer to elements
                Cell< Element * > & tElements = aMesh->block( tID )->elements() ;

                // skip this block if it is empty
                if( tElements.size() == 0 )
                {
                    continue;
                }

                // get the element type
                ElementType tElemType = aMesh->block( tID )->element_type() ;

                // create the shape function
                fem::InterpolationFunction * tShape = tFactory.create_lagrange_function( tElemType ) ;

                // determine interpolation order
                uint tOrder = auto_integration_order( tElemType );

                // populate shape function

                // integration points
                Matrix< real > tXi ;

                // integration weights
                Vector< real > tW ;


                // number of nodes per element
                uint tNumNodes = number_of_nodes( tElemType );

                // element coordinates
                Matrix< real > tX( tNumNodes, tNumDim );

                // populate vectors
                intpoints( IntegrationScheme::GAUSS, geometry_type( tElemType ), tOrder, tW, tXi );

                // get number of integration points
                uint tNumIntpoints = tW.length() ;

                // populate shape derivative
                Cell< Matrix< real > > tdNdXi( tNumIntpoints, Matrix< real >( tNumDim, tNumNodes ) );

                for( uint k=0; k<tNumIntpoints; ++k )
                {
                    tShape->dNdXi( tXi.col( k ), tdNdXi( k ) );
                }

                // tidy up
                delete tShape ;

                // loop over all elements on this block
                for( Element * tElement : tElements )
                {

                    // check if this element is owned
                    if( tElement->owner() == tMyRank )
                    {
                        real tElVolume = 0.0 ;

                        // populate element coordinates
                        for( uint i=0; i<tNumDim; ++i )
                        {
                            for( uint k=0; k<tNumNodes; ++k )
                            {
                                tX( k, i ) = tElement->node( k )->x( i );
                            }
                        }

                        // loop over all integration points
                        for( uint k=0; k<tNumIntpoints; ++k )
                        {
                            // compute geometry jacobian
                            tJ = tdNdXi( k ) * tX ;

                            // add value to total volume
                            tElVolume += tW( k ) * std::abs( det( tJ ) );

                        }

                        aVolume += tElVolume ;
                    } // end if own element
                } // end element loop
            } // end block exists
        } // end loop over all selected blocks

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
                aVolume += sum( tData );

                // send data back to other procs
                tData.fill( aVolume );
                send( tCommList, tData );

            }
            else
            {
                // send my contribution to master
                send( 0, aVolume );

                // receive sum from master
                receive( 0, aVolume );
            }
        }

        return aVolume ;
    }

//------------------------------------------------------------------------------
    }
}