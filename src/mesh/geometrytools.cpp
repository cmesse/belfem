//
// Created by Christian Messe on 28.06.20.
//

#include "geometrytools.hpp"
#include "fn_dot.hpp"
#include "cl_Node.hpp"

namespace belfem
{
    namespace mesh
    {

//------------------------------------------------------------------------------

        void
        collect_node_coords(
                Element *         aElement,
                Matrix< real > & aNodeCoords,
                const uint aNumberOfDimensions,
                const bool aCornerNodesOnly )
        {
            // check dimensions
            uint tNumDim = aNumberOfDimensions == 0 ?
                           dimension( aElement->type() ) : aNumberOfDimensions ;

            // get number of nodes
            uint tNumNodes = aCornerNodesOnly ?
                             aElement->number_of_corner_nodes() :
                             aElement->number_of_nodes() ;

            // check size of matrix
            BELFEM_ASSERT( aNodeCoords.n_cols() == tNumDim ,
                          "Number of colums for input matrix does not match ( is %u, expect %u )",
                          ( unsigned int ) aNodeCoords.n_cols(),
                          ( unsigned int ) tNumDim );

            BELFEM_ASSERT( aNodeCoords.n_rows() == tNumNodes ,
                          "Number of rows for input matrix does not match ( is %u, expect %u )",
                          ( unsigned int ) aNodeCoords.n_rows(),
                          ( unsigned int ) tNumNodes );

            // populate data
            for( uint i=0; i<tNumDim; ++i )
            {
                for( uint k=0; k<tNumNodes; ++k )
                {
                    aNodeCoords( k, i ) = aElement->node( k )->x( i );
                }
            }

        }

//------------------------------------------------------------------------------

        void
        collect_node_data(
                Element              * aElement,
                const Vector< real > & aFieldOnMesh ,
                Vector< real >       & aNodalField,
                const bool             aCornerNodesOnly )
        {

            // get number of nodes
            uint tNumNodes = aCornerNodesOnly ?
                             aElement->number_of_corner_nodes() :
                             aElement->number_of_nodes() ;

            // check size of vector
            BELFEM_ASSERT( tNumNodes == aNodalField.length() ,
                          "Number of nodes on element %lu and output vector do not match ( is %u, expect %u )",
                          ( long unsigned int ) aElement->id(),
                          ( unsigned int ) tNumNodes,
                          ( unsigned int ) aNodalField.length() );

            // collect data
            for( uint k=0; k<tNumNodes; ++k )
            {
                aNodalField( k ) = aFieldOnMesh( aElement->node( k )->index() );
            }
        }

//------------------------------------------------------------------------------

        // Bronstein, Eq. 8.152c
        real
        compute_surface_increment(
                const Matrix< real > & adNdXi,
                const Matrix< real > & aNodeCoords,
                const uint             aNumSpatialDimensions  )
        {
            // get number of dimensions
            uint tNumDim = aNumSpatialDimensions==0 ? adNdXi.n_rows() + 1 : aNumSpatialDimensions ;

            // check size of matrix
            BELFEM_ASSERT( aNodeCoords.n_cols() == tNumDim ,
                          "Number of colums for input matrix does not match ( is %u, expect %u )",
                          ( unsigned int ) aNodeCoords.n_cols(),
                          ( unsigned int ) tNumDim );

            BELFEM_ASSERT( aNodeCoords.n_rows() ==  adNdXi.n_cols() ,
                          "Number of rows for input matrix does not match ( is %u, expect %u )",
                          ( unsigned int ) aNodeCoords.n_rows(),
                          ( unsigned int )  adNdXi.n_cols() );

            // See Bronstein Eq 3.490
            real tE = 0.0 ;
            real tG = 0.0 ;
            real tF = 0.0 ;

            if( tNumDim == 2 )
            {
                for( uint i=0; i<tNumDim; ++i )
                {
                    real tdXdXi = dot( adNdXi.row( 0 ), aNodeCoords.col( i ) );
                    tE += tdXdXi * tdXdXi ;
                }

                return std::sqrt( tE );
            }
            else if( tNumDim == 3 )
            {
                for( uint i=0; i<tNumDim; ++i )
                {
                    real tdXdXi  = dot( adNdXi.row( 0 ), aNodeCoords.col( i ) );
                    real tdXdEta = dot( adNdXi.row( 1 ), aNodeCoords.col( i ));

                    tE += tdXdXi  * tdXdXi ;
                    tG += tdXdEta * tdXdEta ;
                    tF += tdXdXi  * tdXdEta ;
                }

                return std::sqrt( std::abs( tE * tG - tF * tF ) );
            }
            else
            {
                BELFEM_ERROR( false, "invalid dimension: %u", (unsigned  int ) tNumDim );
                return BELFEM_SIGNALING_NAN ;
            }
        }

//------------------------------------------------------------------------------
    }
}