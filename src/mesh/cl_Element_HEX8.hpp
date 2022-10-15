//
// Created by Christian Messe on 2019-07-30.
//

#ifndef BELFEM_CL_ELEMENT_HEX8_HPP
#define BELFEM_CL_ELEMENT_HEX8_HPP

#include "typedefs.hpp"
#include "Mesh_Enums.hpp"
#include "cl_LagrangeElement.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        template <>
        ElementType
        LagrangeElement< 8, 8, 12, 6, 6 >::type() const
        {
            return ElementType::HEX8;
        }

//------------------------------------------------------------------------------

        template <>
        void
        LagrangeElement< 8, 8, 12, 6, 6 >::get_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 4, nullptr );

            switch( aFacetIndex )
            {
                case( 0 ):
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 1 ];
                    aNodes( 2 ) = mNodes[ 5 ];
                    aNodes( 3 ) = mNodes[ 4 ];
                    break;
                }
                case( 1 ):
                {
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    aNodes( 2 ) = mNodes[ 6 ];
                    aNodes( 3 ) = mNodes[ 5 ];
                    break;
                }
                case( 2 ):
                {
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 3 ];
                    aNodes( 2 ) = mNodes[ 7 ];
                    aNodes( 3 ) = mNodes[ 6 ];
                    break;
                }
                case( 3 ):
                {
                    aNodes( 0 ) = mNodes[ 4 ];
                    aNodes( 1 ) = mNodes[ 7 ];
                    aNodes( 2 ) = mNodes[ 3 ];
                    aNodes( 3 ) = mNodes[ 0 ];
                    break;
                }
                case( 4 ):
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 3 ];
                    aNodes( 2 ) = mNodes[ 2 ];
                    aNodes( 3 ) = mNodes[ 1 ];
                    break;
                }
                case( 5 ):
                {
                    aNodes( 0 ) = mNodes[ 4 ];
                    aNodes( 1 ) = mNodes[ 5 ];
                    aNodes( 2 ) = mNodes[ 6 ];
                    aNodes( 3 ) = mNodes[ 7 ];
                    break;
                }
                default:
                {
                    this->throw_facet_error( aFacetIndex );
                }
            }

        }

//------------------------------------------------------------------------------

        template <>
        void
        LagrangeElement< 8, 8, 12, 6, 6 >::get_corner_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
        {
            this->get_nodes_of_facet( aFacetIndex, aNodes );
        }

//------------------------------------------------------------------------------

        template <>
        void
        LagrangeElement< 8, 8, 12, 6, 6 >::get_nodes_of_edge( const uint aEdgeIndex, Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 2, nullptr );

            switch ( aEdgeIndex )
            {
                case ( 0 ):
                {
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 0 ];
                    break;
                }
                case ( 1 ):
                {
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 1 ];
                    break;
                }
                case ( 2 ):
                {
                    aNodes( 0 ) = mNodes[ 3 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    break;
                }
                case ( 3 ):
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 3 ];
                    break;
                }
                case ( 4 ):
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 4 ];
                    break;
                }
                case ( 5 ):
                {
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 5 ];
                    break;
                }
                case ( 6 ):
                {
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 6 ];
                    break;
                }
                case ( 7 ):
                {
                    aNodes( 0 ) = mNodes[ 3 ];
                    aNodes( 1 ) = mNodes[ 7 ];
                    break;
                }
                case ( 8 ):
                {
                    aNodes( 0 ) = mNodes[ 4 ];
                    aNodes( 1 ) = mNodes[ 5 ];
                    break;
                }
                case ( 9 ):
                {
                    aNodes( 0 ) = mNodes[ 5 ];
                    aNodes( 1 ) = mNodes[ 6 ];
                    break;
                }
                case ( 10 ):
                {
                    aNodes( 0 ) = mNodes[ 6 ];
                    aNodes( 1 ) = mNodes[ 7 ];
                    break;
                }
                case ( 11 ):
                {
                    aNodes( 0 ) = mNodes[ 7 ];
                    aNodes( 1 ) = mNodes[ 4 ];
                    break;
                }
                default:
                {
                    this->throw_edge_error( aEdgeIndex );
                }
            }
        }

//------------------------------------------------------------------------------

        template <>
        void
        LagrangeElement< 8, 8, 12, 6, 6 >::get_edges_of_facet(
                const uint aFacetIndex, Cell< Edge * > & aEdges )
        {
            // allocate the node container
            aEdges.set_size( 4, nullptr );

            switch( aFacetIndex )
            {
                case( 0 ):
                {
                    aEdges( 0 ) = mEdges[  0 ];
                    aEdges( 1 ) = mEdges[  5 ];
                    aEdges( 2 ) = mEdges[  8 ];
                    aEdges( 3 ) = mEdges[  4 ];
                    break;
                }
                case( 1 ):
                {
                    aEdges( 0 ) = mEdges[  1 ];
                    aEdges( 1 ) = mEdges[  6 ];
                    aEdges( 2 ) = mEdges[  9 ];
                    aEdges( 3 ) = mEdges[  5 ];
                    break;
                }
                case( 2 ):
                {
                    aEdges( 0 ) = mEdges[  2 ];
                    aEdges( 1 ) = mEdges[  7 ];
                    aEdges( 2 ) = mEdges[ 10 ];
                    aEdges( 3 ) = mEdges[  6 ];
                    break;
                }
                case( 3 ):
                {
                    aEdges( 0 ) = mEdges[ 11 ];
                    aEdges( 1 ) = mEdges[  7 ];
                    aEdges( 2 ) = mEdges[  3 ];
                    aEdges( 3 ) = mEdges[  4 ];
                    break;
                }
                case( 4 ):
                {
                    aEdges( 0 ) = mEdges[  3 ];
                    aEdges( 1 ) = mEdges[  2 ];
                    aEdges( 2 ) = mEdges[  1 ];
                    aEdges( 3 ) = mEdges[  0 ];
                    break;
                }
                case( 5 ):
                {
                    aEdges( 0 ) = mEdges[  8 ];
                    aEdges( 1 ) = mEdges[  9 ];
                    aEdges( 2 ) = mEdges[ 10 ];
                    aEdges( 3 ) = mEdges[ 11 ];
                    break;
                }
                default:
                {
                    this->throw_facet_error( aFacetIndex );
                }
            }
        }

//------------------------------------------------------------------------------
    } /* namespace mesh */
} /* namespace belfem */

#endif //BELFEM_CL_ELEMENT_HEX8_HPP
