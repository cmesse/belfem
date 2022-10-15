//
// Created by Christian Messe on 2019-08-04.
//

#ifndef BELFEM_CL_ELEMENT_PENTA15_HPP
#define BELFEM_CL_ELEMENT_PENTA15_HPP

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
        LagrangeElement< 15, 6, 9, 5, 5 >::type() const
        {
            return ElementType::PENTA15;
        }

//------------------------------------------------------------------------------

        template <>
        void
        LagrangeElement< 15, 6, 9, 5, 5 >::get_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
        {
            switch( aFacetIndex )
            {
                case( 0 ):
                    aNodes.set_size( 8, nullptr );
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 1 ];
                    aNodes( 2 ) = mNodes[ 4 ];
                    aNodes( 3 ) = mNodes[ 3 ];
                    aNodes( 4 ) = mNodes[ 6 ];
                    aNodes( 5 ) = mNodes[ 10 ];
                    aNodes( 6 ) = mNodes[ 12 ];
                    aNodes( 7 ) = mNodes[ 9 ];
                    break;
                case( 1 ):
                    aNodes.set_size( 8, nullptr );
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    aNodes( 2 ) = mNodes[ 5 ];
                    aNodes( 3 ) = mNodes[ 4 ];
                    aNodes( 4 ) = mNodes[ 7 ];
                    aNodes( 5 ) = mNodes[ 11 ];
                    aNodes( 6 ) = mNodes[ 13 ];
                    aNodes( 7 ) = mNodes[ 10 ];
                    break;
                case( 2 ):
                    aNodes.set_size( 8, nullptr );
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 0 ];
                    aNodes( 2 ) = mNodes[ 3 ];
                    aNodes( 3 ) = mNodes[ 5 ];
                    aNodes( 4 ) = mNodes[ 8 ];
                    aNodes( 5 ) = mNodes[ 9 ];
                    aNodes( 6 ) = mNodes[ 14 ];
                    aNodes( 7 ) = mNodes[ 11 ];
                    break;
                case( 3 ):
                    aNodes.set_size( 6, nullptr );
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    aNodes( 2 ) = mNodes[ 1 ];
                    aNodes( 3 ) = mNodes[ 8 ];
                    aNodes( 4 ) = mNodes[ 7 ];
                    aNodes( 5 ) = mNodes[ 6 ];
                    break;
                case( 4 ):
                    aNodes.set_size( 6, nullptr );
                    aNodes( 0 ) = mNodes[ 3 ];
                    aNodes( 1 ) = mNodes[ 4 ];
                    aNodes( 2 ) = mNodes[ 5 ];
                    aNodes( 3 ) = mNodes[ 12 ];
                    aNodes( 4 ) = mNodes[ 13 ];
                    aNodes( 5 ) = mNodes[ 14 ];
                    break;
                default:
                {
                    this->throw_facet_error( aFacetIndex );
                }
            }
        }

//------------------------------------------------------------------------------

        template <>
        void
        LagrangeElement< 15, 6, 9, 5, 5 >::get_corner_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
        {
            switch( aFacetIndex )
            {
                case( 0 ):
                {
                    aNodes.set_size( 4, nullptr );
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 1 ];
                    aNodes( 2 ) = mNodes[ 4 ];
                    aNodes( 3 ) = mNodes[ 3 ];
                    break;
                }
                case( 1 ):
                {
                    aNodes.set_size( 4, nullptr );
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    aNodes( 2 ) = mNodes[ 5 ];
                    aNodes( 3 ) = mNodes[ 4 ];
                    break;
                }
                case( 2 ):
                {
                    aNodes.set_size( 4, nullptr );
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 0 ];
                    aNodes( 2 ) = mNodes[ 3 ];
                    aNodes( 3 ) = mNodes[ 5 ];
                    break;
                }
                case( 3 ):
                {
                    aNodes.set_size( 3, nullptr );
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    aNodes( 2 ) = mNodes[ 1 ];
                    break;
                }
                case( 4 ):
                {
                    aNodes.set_size( 3, nullptr );
                    aNodes( 0 ) = mNodes[ 3 ];
                    aNodes( 1 ) = mNodes[ 4 ];
                    aNodes( 2 ) = mNodes[ 5 ];
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
        LagrangeElement< 15, 6, 9, 5, 5 >::get_nodes_of_edge( const uint aEdgeIndex, Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 3, nullptr );

            switch( aEdgeIndex )
            {
                case( 0 ):
                {
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 0 ];
                    aNodes( 2 ) = mNodes[ 6 ];
                    break;
                }
                case( 1 ):
                {
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 1 ];
                    aNodes( 2 ) = mNodes[ 7 ];
                    break;
                }
                case( 2 ):
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    aNodes( 2 ) = mNodes[ 8 ];
                    break;
                }
                case( 3 ):
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 3 ];
                    aNodes( 2 ) = mNodes[ 9 ];
                    break;
                }
                case( 4 ):
                {
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 4 ];
                    aNodes( 2 ) = mNodes[ 10 ];
                    break;
                }
                case( 5 ):
                {
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 5 ];
                    aNodes( 2 ) = mNodes[ 11 ];
                    break;
                }
                case( 6 ):
                {
                    aNodes( 0 ) = mNodes[ 3 ];
                    aNodes( 1 ) = mNodes[ 4 ];
                    aNodes( 2 ) = mNodes[ 12 ];
                    break;
                }
                case( 7 ):
                {
                    aNodes( 0 ) = mNodes[ 4 ];
                    aNodes( 1 ) = mNodes[ 5 ];
                    aNodes( 2 ) = mNodes[ 13 ];
                    break;
                }
                case( 8 ):
                {
                    aNodes( 0 ) = mNodes[ 5 ];
                    aNodes( 1 ) = mNodes[ 3 ];
                    aNodes( 2 ) = mNodes[ 14 ];
                    break;
                }
                default:
                {
                    this->throw_facet_error( aEdgeIndex );
                }
            }
        }

//------------------------------------------------------------------------------

        template <>
        void
        LagrangeElement< 15, 6, 9, 5, 5 >::get_edges_of_facet(
                const uint aFacetIndex, Cell< Edge * > & aEdges )
        {
            switch( aFacetIndex )
            {
                case ( 0 ):
                {
                    aEdges.set_size( 4, nullptr );
                    aEdges( 0 ) = mEdges[ 0 ];
                    aEdges( 1 ) = mEdges[ 4 ];
                    aEdges( 2 ) = mEdges[ 6 ];
                    aEdges( 3 ) = mEdges[ 3 ];
                    break;
                }
                case ( 1 ):
                {
                    aEdges.set_size( 4, nullptr );
                    aEdges( 0 ) = mEdges[ 1 ];
                    aEdges( 1 ) = mEdges[ 5 ];
                    aEdges( 2 ) = mEdges[ 7 ];
                    aEdges( 3 ) = mEdges[ 4 ];
                    break;
                }
                case ( 2 ):
                {
                    aEdges.set_size( 4, nullptr );
                    aEdges( 0 ) = mEdges[ 2 ];
                    aEdges( 1 ) = mEdges[ 3 ];
                    aEdges( 2 ) = mEdges[ 8 ];
                    aEdges( 3 ) = mEdges[ 5 ];
                    break;
                }
                case ( 3 ):
                {
                    aEdges.set_size( 3, nullptr );
                    aEdges( 0 ) = mEdges[ 2 ];
                    aEdges( 1 ) = mEdges[ 1 ];
                    aEdges( 2 ) = mEdges[ 0 ];
                    break;
                }
                case ( 4 ):
                {
                    aEdges.set_size( 3, nullptr );
                    aEdges( 0 ) = mEdges[ 6 ];
                    aEdges( 1 ) = mEdges[ 7 ];
                    aEdges( 2 ) = mEdges[ 8 ];
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

#endif //BELFEM_CL_ELEMENT_PENTA15_HPP
