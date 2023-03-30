//
// Created by Christian Messe on 30.04.20.
//

#ifndef BELFEM_CL_ELEMENT_HEX64_HPP
#define BELFEM_CL_ELEMENT_HEX64_HPP


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
        LagrangeElement< 64, 8, 12, 6, 6 >::type() const
        {
            return ElementType::HEX64;
        }

//------------------------------------------------------------------------------

        template <>
        void
        LagrangeElement< 64, 8, 12, 6, 6 >::get_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 16, nullptr );

            switch( aFacetIndex )
            {
                case( 0 ):
                {
                    aNodes(  0 ) = mNodes[  0 ];
                    aNodes(  1 ) = mNodes[  1 ];
                    aNodes(  2 ) = mNodes[  5 ];
                    aNodes(  3 ) = mNodes[  4 ];
                    aNodes(  4 ) = mNodes[  8 ];
                    aNodes(  5 ) = mNodes[  9 ];
                    aNodes(  6 ) = mNodes[ 16 ];
                    aNodes(  7 ) = mNodes[ 17 ];
                    aNodes(  8 ) = mNodes[ 25 ];
                    aNodes(  9 ) = mNodes[ 24 ];
                    aNodes( 10 ) = mNodes[ 13 ];
                    aNodes( 11 ) = mNodes[ 12 ];
                    aNodes( 12 ) = mNodes[ 36 ];
                    aNodes( 13 ) = mNodes[ 37 ];
                    aNodes( 14 ) = mNodes[ 38 ];
                    aNodes( 15 ) = mNodes[ 39 ];
                    break;
                }
                case( 1 ):
                {
                    aNodes(  0 ) = mNodes[  1 ];
                    aNodes(  1 ) = mNodes[  2 ];
                    aNodes(  2 ) = mNodes[  6 ];
                    aNodes(  3 ) = mNodes[  5 ];
                    aNodes(  4 ) = mNodes[ 14 ];
                    aNodes(  5 ) = mNodes[ 15 ];
                    aNodes(  6 ) = mNodes[ 20 ];
                    aNodes(  7 ) = mNodes[ 21 ];
                    aNodes(  8 ) = mNodes[ 29 ];
                    aNodes(  9 ) = mNodes[ 28 ];
                    aNodes( 10 ) = mNodes[ 17 ];
                    aNodes( 11 ) = mNodes[ 16 ];
                    aNodes( 12 ) = mNodes[ 44 ];
                    aNodes( 13 ) = mNodes[ 45 ];
                    aNodes( 14 ) = mNodes[ 46 ];
                    aNodes( 15 ) = mNodes[ 47 ];
                    break;
                }
                case( 2 ):
                {
                    aNodes(  0 ) = mNodes[  2 ];
                    aNodes(  1 ) = mNodes[  3 ];
                    aNodes(  2 ) = mNodes[  7 ];
                    aNodes(  3 ) = mNodes[  6 ];
                    aNodes(  4 ) = mNodes[ 18 ];
                    aNodes(  5 ) = mNodes[ 19 ];
                    aNodes(  6 ) = mNodes[ 22 ];
                    aNodes(  7 ) = mNodes[ 23 ];
                    aNodes(  8 ) = mNodes[ 31 ];
                    aNodes(  9 ) = mNodes[ 30 ];
                    aNodes( 10 ) = mNodes[ 21 ];
                    aNodes( 11 ) = mNodes[ 20 ];
                    aNodes( 12 ) = mNodes[ 48 ];
                    aNodes( 13 ) = mNodes[ 49 ];
                    aNodes( 14 ) = mNodes[ 50 ];
                    aNodes( 15 ) = mNodes[ 51 ];
                    break;
                }
                case( 3 ):
                {
                    aNodes(  0 ) = mNodes[  0 ];
                    aNodes(  1 ) = mNodes[  4 ];
                    aNodes(  2 ) = mNodes[  7 ];
                    aNodes(  3 ) = mNodes[  3 ];
                    aNodes(  4 ) = mNodes[ 12 ];
                    aNodes(  5 ) = mNodes[ 13 ];
                    aNodes(  6 ) = mNodes[ 26 ];
                    aNodes(  7 ) = mNodes[ 27 ];
                    aNodes(  8 ) = mNodes[ 23 ];
                    aNodes(  9 ) = mNodes[ 22 ];
                    aNodes( 10 ) = mNodes[ 11 ];
                    aNodes( 11 ) = mNodes[ 10 ];
                    aNodes( 12 ) = mNodes[ 40 ];
                    aNodes( 13 ) = mNodes[ 41 ];
                    aNodes( 14 ) = mNodes[ 42 ];
                    aNodes( 15 ) = mNodes[ 43 ];
                    break;
                }
                case( 4 ):
                {
                    aNodes(  0 ) = mNodes[  0 ];
                    aNodes(  1 ) = mNodes[  3 ];
                    aNodes(  2 ) = mNodes[  2 ];
                    aNodes(  3 ) = mNodes[  1 ];
                    aNodes(  4 ) = mNodes[ 10 ];
                    aNodes(  5 ) = mNodes[ 11 ];
                    aNodes(  6 ) = mNodes[ 19 ];
                    aNodes(  7 ) = mNodes[ 18 ];
                    aNodes(  8 ) = mNodes[ 15 ];
                    aNodes(  9 ) = mNodes[ 14 ];
                    aNodes( 10 ) = mNodes[  9 ];
                    aNodes( 11 ) = mNodes[  8 ];
                    aNodes( 12 ) = mNodes[ 32 ];
                    aNodes( 13 ) = mNodes[ 33 ];
                    aNodes( 14 ) = mNodes[ 34 ];
                    aNodes( 15 ) = mNodes[ 35 ];
                    break;
                }
                case( 5 ):
                {
                    aNodes(  0 ) = mNodes[  4 ];
                    aNodes(  1 ) = mNodes[  5 ];
                    aNodes(  2 ) = mNodes[  6 ];
                    aNodes(  3 ) = mNodes[  7 ];
                    aNodes(  4 ) = mNodes[ 24 ];
                    aNodes(  5 ) = mNodes[ 25 ];
                    aNodes(  6 ) = mNodes[ 28 ];
                    aNodes(  7 ) = mNodes[ 29 ];
                    aNodes(  8 ) = mNodes[ 30 ];
                    aNodes(  9 ) = mNodes[ 31 ];
                    aNodes( 10 ) = mNodes[ 27 ];
                    aNodes( 11 ) = mNodes[ 26 ];
                    aNodes( 12 ) = mNodes[ 52 ];
                    aNodes( 13 ) = mNodes[ 53 ];
                    aNodes( 14 ) = mNodes[ 54 ];
                    aNodes( 15 ) = mNodes[ 55 ];
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
        LagrangeElement< 64, 8, 12, 6, 6 >::get_corner_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
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
                    aNodes( 0 ) = mNodes[  0 ];
                    aNodes( 1 ) = mNodes[  4 ];
                    aNodes( 2 ) = mNodes[  7 ];
                    aNodes( 3 ) = mNodes[  3 ];
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
        LagrangeElement< 64, 8, 12, 6, 6 >::get_nodes_of_edge( const uint aEdgeIndex, Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 4, nullptr );

            switch ( aEdgeIndex )
            {
                case ( 0 ):
                {
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 0 ];
                    aNodes( 2 ) = mNodes[ 9 ];
                    aNodes( 3 ) = mNodes[ 8 ];
                    break;
                }
                case ( 1 ):
                {
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 1 ];
                    aNodes( 2 ) = mNodes[ 15 ];
                    aNodes( 3 ) = mNodes[ 14 ];
                    break;
                }
                case ( 2 ):
                {
                    aNodes( 0 ) = mNodes[ 3 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    aNodes( 2 ) = mNodes[ 19 ];
                    aNodes( 3 ) = mNodes[ 18 ];
                    break;
                }
                case ( 3 ):
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 3 ];
                    aNodes( 2 ) = mNodes[ 10 ];
                    aNodes( 3 ) = mNodes[ 11 ];
                    break;
                }
                case ( 4 ):
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 4 ];
                    aNodes( 2 ) = mNodes[ 12 ];
                    aNodes( 3 ) = mNodes[ 13 ];
                    break;
                }
                case ( 5 ):
                {
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 5 ];
                    aNodes( 2 ) = mNodes[ 16 ];
                    aNodes( 3 ) = mNodes[ 17 ];
                    break;
                }
                case ( 6 ):
                {
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 6 ];
                    aNodes( 2 ) = mNodes[ 20 ];
                    aNodes( 3 ) = mNodes[ 21 ];
                    break;
                }
                case ( 7 ):
                {
                    aNodes( 0 ) = mNodes[ 3 ];
                    aNodes( 1 ) = mNodes[ 7 ];
                    aNodes( 2 ) = mNodes[ 22 ];
                    aNodes( 3 ) = mNodes[ 23 ];
                    break;
                }
                case ( 8 ):
                {
                    aNodes( 0 ) = mNodes[ 4 ];
                    aNodes( 1 ) = mNodes[ 5 ];
                    aNodes( 2 ) = mNodes[ 24 ];
                    aNodes( 3 ) = mNodes[ 25 ];
                    break;
                }
                case ( 9 ):
                {
                    aNodes( 0 ) = mNodes[ 5 ];
                    aNodes( 1 ) = mNodes[ 6 ];
                    aNodes( 2 ) = mNodes[ 28 ];
                    aNodes( 3 ) = mNodes[ 29 ];
                    break;
                }
                case ( 10 ):
                {
                    aNodes( 0 ) = mNodes[ 6 ];
                    aNodes( 1 ) = mNodes[ 7 ];
                    aNodes( 2 ) = mNodes[ 30 ];
                    aNodes( 3 ) = mNodes[ 31 ];
                    break;
                }
                case ( 11 ):
                {
                    aNodes( 0 ) = mNodes[ 7 ];
                    aNodes( 1 ) = mNodes[ 4 ];
                    aNodes( 2 ) = mNodes[ 27 ];
                    aNodes( 3 ) = mNodes[ 26 ];
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
        LagrangeElement< 64, 8, 12, 6, 6 >::get_edges_of_facet(
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
                    aEdges( 0 ) = mEdges[  4 ];
                    aEdges( 1 ) = mEdges[ 11 ];
                    aEdges( 2 ) = mEdges[  7 ];
                    aEdges( 3 ) = mEdges[  3 ];
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

#endif //BELFEM_CL_ELEMENT_HEX64_HPP
