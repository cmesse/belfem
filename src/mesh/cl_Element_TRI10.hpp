//
// Created by Christian Messe on 25.07.20.
//

#ifndef BELFEM_CL_ELEMENT_TRI10_HPP
#define BELFEM_CL_ELEMENT_TRI10_HPP


#include "typedefs.hpp"
#include "Mesh_Enums.hpp"
#include "cl_ElementTemplate.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        template <>
        ElementType
        ElementTemplate< 10, 3, 3, 3, 1 >::type() const
        {
            return ElementType::TRI10;
        }

//------------------------------------------------------------------------------

        template <>
        void
        ElementTemplate< 10, 3, 3, 3, 1 >::get_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 4, nullptr );

            switch( aFacetIndex )
            {
                case( 0 ):
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 1 ];
                    aNodes( 2 ) = mNodes[ 3 ];
                    aNodes( 3 ) = mNodes[ 4 ];
                    break;
                }
                case( 1 ):
                {
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    aNodes( 2 ) = mNodes[ 5 ];
                    aNodes( 3 ) = mNodes[ 6 ];
                    break;
                }
                case( 2 ):
                {
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 0 ];
                    aNodes( 2 ) = mNodes[ 7 ];
                    aNodes( 3 ) = mNodes[ 8 ];
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
        ElementTemplate< 10, 3, 3, 3, 1 >::get_corner_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 2, nullptr );

            switch( aFacetIndex )
            {
                case( 0 ):
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 1 ];
                    break;
                }
                case( 1 ):
                {
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    break;
                }
                case( 2 ):
                {
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 0 ];
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
        ElementTemplate< 10, 3, 3, 3, 1 >::get_edges_of_facet(
                const uint aFacetIndex, Cell< Edge * > & aEdges )
        {
            // allocate the node container
            aEdges.set_size( 1, nullptr );

            switch( aFacetIndex )
            {
                case( 0 ):
                {
                    aEdges( 0 ) = mEdges[ 0 ];
                    break ;
                }
                case( 1 ):
                {
                    aEdges( 0 ) = mEdges[ 1 ];
                    break ;
                }
                case( 2 ):
                {
                    aEdges( 0 ) = mEdges[ 2 ];
                    break ;
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

#endif //BELFEM_CL_ELEMENT_TRI10_HPP
