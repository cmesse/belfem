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
                      "in parallel mode, input argument of compute_volume must be Kernel->mesh()" );

        const proc_t tMyRank = comm_rank() ;

        fem::InterpolationFunctionFactory tFactory ;

        // total volume
        real aVolume = 0.0 ;

        // nuber of dimensions of mesh
        uint tNumDim = aMesh->number_of_dimensions() ;

        // geometry jacobian
        Matrix< real > tJ( tNumDim, tNumDim );

        for( id_t tID : aBlockIDs )
        {
            if( aMesh->block_exists( tID ) )
            {
                Cell< Element * > & tElements = aMesh->block( tID )->elements() ;

                // skip this block if it is empty
                if( tElements.size() == 0 )
                {
                    continue;
                }

                ElementType tElemType = aMesh->block( tID )->element_type() ;

                fem::InterpolationFunction * tShape = tFactory.create_lagrange_function( tElemType ) ;

                uint tOrder = auto_integration_order( tElemType );


                // integration points
                Matrix< real > tXi ;

                // integration weights
                Vector< real > tW ;


                // number of nodes per element
                uint tNumNodes = number_of_nodes( tElemType );

                // element coordinates
                Matrix< real > tX( tNumNodes, tNumDim );

                intpoints( IntegrationScheme::GAUSS, geometry_type( tElemType ), tOrder, tW, tXi );

                uint tNumIntpoints = tW.length() ;

                Cell< Matrix< real > > tdNdXi( tNumIntpoints, Matrix< real >( tNumDim, tNumNodes ) );

                for( uint k=0; k<tNumIntpoints; ++k )
                {
                    tShape->dNdXi( tXi.col( k ), tdNdXi( k ) );
                }

                delete tShape ;

                for( Element * tElement : tElements )
                {

                    if( tElement->owner() == tMyRank )
                    {
                        real tElVolume = 0.0 ;

                        for( uint i=0; i<tNumDim; ++i )
                        {
                            for( uint k=0; k<tNumNodes; ++k )
                            {
                                tX( k, i ) = tElement->node( k )->x( i );
                            }
                        }

                        for( uint k=0; k<tNumIntpoints; ++k )
                        {
                            tJ = tdNdXi( k ) * tX ;

                            tElVolume += tW( k ) * std::abs( det( tJ ) );

                        }

                        aVolume += tElVolume ;
                    } // end if own element
                } // end element loop
            } // end block exists
        } // end loop over all selected blocks

        proc_t tCommSize = comm_size() ;

        if( tCommSize > 1 )
        {
            if ( tMyRank == 0 )
            {
                Vector< real > tData( tCommSize, 0 );

                collect( tData );

                aVolume += sum( tData );

            }
            else
            {
                send( aVolume );
            }

            broadcast( aVolume );
        }

        return aVolume ;
    }

//------------------------------------------------------------------------------
    }
}
