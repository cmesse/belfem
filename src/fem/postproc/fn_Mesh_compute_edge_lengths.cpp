//
// Created by Christian Messe on 14.06.20.
//

#include "fn_Mesh_compute_edge_lengths.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"
#include "meshtools.hpp"
#include "fn_norm.hpp"
#include "fn_intpoints.hpp"
#include "fn_trans.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        void
        compute_edge_lengths(
                const uint                aDimension,
                Cell< mesh::Element * > & aElements,
                Vector< real >          & aEdgeLength )
        {
            if( aElements.size() > 0 )
            {
                // get type of first element
                ElementType tType = aElements( 0 )->type();

#if !defined( NDEBUG ) || defined( DEBUG )
                // make sure that all elements are from the same type
                for( mesh::Element * tElement : aElements )
                {
                    BELFEM_ERROR( tElement->type() == tType,
                            "All elements in input container must have the same type." );
                }
#endif
                // make sure that the passed elements are of the type line
                BELFEM_ASSERT( geometry_type( tType ) == GeometryType::LINE,
                        "All elements in input container must be of the type LINE."
                        ) ;

                // allocate memory for output vector
                aEdgeLength.set_size( aElements.size() );

                // initialize counter
                index_t tCount = 0 ;

                if( tType == ElementType::LINE2 )
                {
                    Vector< real > tA( aDimension );
                    Vector< real > tB( aDimension );

                    for( Element * tElement : aElements )
                    {
                        // populate coordinate vectors
                        for( uint i=0; i<aDimension; ++i )
                        {
                            tA( i ) = tElement->node( 0 )->x( i ) ;
                        }
                        for( uint i=0; i<aDimension; ++i )
                        {
                            tB( i ) = tElement->node( 1 )->x( i ) ;
                        }

                        aEdgeLength( tCount++ ) = norm( tB - tA );
                    }
                }
                else
                {
                    // create the factory
                    fem::InterpolationFunctionFactory tFactory;

                    // create the shape function
                    fem::InterpolationFunction * tShape = tFactory.create_lagrange_function( tType );

                    // get the nuber of nodes per element
                    uint tNumNodesPerElement = number_of_nodes( tType );

                    // integration weights
                    Vector< real > tW;

                    // integration points
                    Matrix< real > tXi;

                    // compute the integration order
                    uint tOrder = tType == ElementType::LINE3 ? 7 : 9;

                    // get the integration points
                    intpoints(
                            IntegrationScheme::GAUSS,
                            GeometryType::LINE,
                            tOrder,
                            tW,
                            tXi );

                    // get the number of integration points
                    uint tNumIntpoints = tW.length();

                    // populate shape array
                    Cell< Matrix< real > > tdNdXi( tNumIntpoints, Matrix< real >( 1, tNumNodesPerElement ));
                    for ( uint k = 0; k < tNumIntpoints; ++k )
                    {
                        tShape->dNdXi( tXi.col( k ), tdNdXi( k ));
                        tdNdXi( k ) = trans( tdNdXi( k ));
                    }

                    // matrix with node coordinates
                    Matrix< real > tNodeCoords( aDimension, tNumNodesPerElement );

                    // loop over all elements
                    for( Element * tElement : aElements )
                    {
                        // populate node coordinates
                        for( uint k=0; k<tNumNodesPerElement; ++k )
                        {
                            Node * tNode = tElement->node( k );

                            for( uint i=0; i<aDimension; ++i )
                            {
                                tNodeCoords( i, k ) = tNode->x( i );
                            }
                        }

                        // init value for length
                        real tLength = 0 ;


                        // perform integration
                        for( uint k=0; k<tNumIntpoints; ++k )
                        {
                            tLength += tW( k ) * norm( tNodeCoords.matrix_data() * tdNdXi( k ).matrix_data() );
                        }

                        // add length to container
                        aEdgeLength( tCount++ ) = tLength ;
                    }

                    // delete the shape
                    delete tShape;
                }

            }
        }

//------------------------------------------------------------------------------

    }
}