//
// Created by Christian Messe on 2019-08-04.
//

#ifndef BELFEM_CL_ELEMENT_TET10_HPP
#define BELFEM_CL_ELEMENT_TET10_HPP

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
        LagrangeElement< 10, 4, 6, 4, 4 >::type() const
        {
            return ElementType::TET10;
        }

//------------------------------------------------------------------------------

        template <>
        void
        LagrangeElement< 10, 4, 6, 4, 4 >::get_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 6, nullptr );

            switch( aFacetIndex )
            {
                case( 0 ):
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 1 ];
                    aNodes( 2 ) = mNodes[ 3 ];
                    aNodes( 3 ) = mNodes[ 4 ];
                    aNodes( 4 ) = mNodes[ 8 ];
                    aNodes( 5 ) = mNodes[ 7 ];
                    break;
                }
                case( 1 ):
                {
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    aNodes( 2 ) = mNodes[ 3 ];
                    aNodes( 3 ) = mNodes[ 5 ];
                    aNodes( 4 ) = mNodes[ 9 ];
                    aNodes( 5 ) = mNodes[ 8 ];
                    break;
                }
                case( 2 ):
                {
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 0 ];
                    aNodes( 2 ) = mNodes[ 3 ];
                    aNodes( 3 ) = mNodes[ 6 ];
                    aNodes( 4 ) = mNodes[ 7 ];
                    aNodes( 5 ) = mNodes[ 9 ];
                    break;
                }
                case( 3 ):
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    aNodes( 2 ) = mNodes[ 1 ];
                    aNodes( 3 ) = mNodes[ 6 ];
                    aNodes( 4 ) = mNodes[ 5 ];
                    aNodes( 5 ) = mNodes[ 4 ];
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
        LagrangeElement< 10, 4, 6, 4, 4 >::get_corner_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
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
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 0 ];
                    aNodes( 2 ) = mNodes[ 3 ];
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
        LagrangeElement< 10, 4, 6, 4, 4 >::get_nodes_of_edge( const uint aEdgeIndex, Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 3, nullptr );

            switch ( aEdgeIndex )
            {
                case ( 0 ):
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 1 ];
                    aNodes( 2 ) = mNodes[ 4 ];
                    break;
                }
                case ( 1 ):
                {
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    aNodes( 2 ) = mNodes[ 5 ];
                    break;
                }
                case ( 2 ):
                {
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 0 ];
                    aNodes( 2 ) = mNodes[ 6 ];
                    break;
                }
                case ( 3 ):
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 3 ];
                    aNodes( 2 ) = mNodes[ 7 ];
                    break;
                }
                case ( 4 ):
                {
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 3 ];
                    aNodes( 2 ) = mNodes[ 8 ];
                    break;
                }
                case ( 5 ):
                {
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 3 ];
                    aNodes( 2 ) = mNodes[ 9 ];
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
        LagrangeElement< 10, 4, 6, 4, 4 >::get_edges_of_facet(
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
                    aEdges( 0 ) = mEdges[ 2 ];
                    aEdges( 1 ) = mEdges[ 3 ];
                    aEdges( 2 ) = mEdges[ 5 ];
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

#endif //BELFEM_CL_ELEMENT_TET10_HPP
