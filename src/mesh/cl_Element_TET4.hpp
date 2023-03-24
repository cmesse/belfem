//
// Created by Christian Messe on 2019-07-28.
//

#ifndef BELFEM_CL_ELEMENT_TET4_HPP
#define BELFEM_CL_ELEMENT_TET4_HPP

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
        LagrangeElement< 4, 4, 6, 4, 4 >::type() const
        {
            return ElementType::TET4;
        }

//------------------------------------------------------------------------------

        template <>
        void
        LagrangeElement< 4, 4, 6, 4, 4 >::get_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 3, nullptr );

            switch( aFacetIndex )
            {
                case( 0 ):
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 1 ];
                    aNodes( 2 ) = mNodes[ 3 ];
                    break;
                }
                case( 1 ):
                {
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    aNodes( 2 ) = mNodes[ 3 ];
                    break;
                }
                case( 2 ):
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 3 ];
                    aNodes( 2 ) = mNodes[ 2 ];
                    break;
                }
                case( 3 ):
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    aNodes( 2 ) = mNodes[ 1 ];
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
        LagrangeElement< 4, 4, 6, 4, 4 >::get_corner_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
        {
            this->get_nodes_of_facet( aFacetIndex, aNodes );
        }

//------------------------------------------------------------------------------

        template <>
        void
        LagrangeElement< 4, 4, 6, 4, 4 >::get_nodes_of_edge( const uint aEdgeIndex, Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 2, nullptr );

            switch ( aEdgeIndex )
            {
                case ( 0 ):
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 1 ];
                    break;
                }
                case ( 1 ):
                {
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    break;
                }
                case ( 2 ):
                {
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 0 ];
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
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 3 ];
                    break;
                }
                case ( 5 ):
                {
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 3 ];
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
        LagrangeElement< 4, 4, 6, 4, 4 >::get_edges_of_facet(
                const uint aFacetIndex, Cell< Edge * > & aEdges )
        {
            // allocate the node container
            aEdges.set_size( 3, nullptr );

            switch( aFacetIndex )
            {
                case( 0 ):
                {
                    aEdges( 0 ) = mEdges[ 0 ];
                    aEdges( 1 ) = mEdges[ 4 ];
                    aEdges( 2 ) = mEdges[ 3 ];
                    break;
                }
                case( 1 ):
                {
                    aEdges( 0 ) = mEdges[ 1 ];
                    aEdges( 1 ) = mEdges[ 5 ];
                    aEdges( 2 ) = mEdges[ 4 ];
                    break;
                }
                case( 2 ):
                {
                    aEdges( 0 ) = mEdges[ 3 ];
                    aEdges( 1 ) = mEdges[ 5 ];
                    aEdges( 2 ) = mEdges[ 2 ];
                    break;
                }
                case( 3 ):
                {
                    aEdges( 0 ) = mEdges[ 2 ];
                    aEdges( 1 ) = mEdges[ 1 ];
                    aEdges( 2 ) = mEdges[ 0 ];
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
#endif //BELFEM_CL_ELEMENT_TET4_HPP
