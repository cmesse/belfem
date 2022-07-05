//
// Created by Christian Messe on 25.07.20.
//

#ifndef BELFEM_CL_ELEMENT_TRI15_HPP
#define BELFEM_CL_ELEMENT_TRI15_HPP



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
        LagrangeElement< 15, 3, 3, 3, 1 >::type() const
        {
            return ElementType::TRI15;
        }

//------------------------------------------------------------------------------

        template <>
        void
        LagrangeElement< 15, 3, 3, 3, 1 >::get_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 5, nullptr );

            switch( aFacetIndex )
            {
                case( 0 ):
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 1 ];
                    aNodes( 2 ) = mNodes[ 3 ];
                    aNodes( 3 ) = mNodes[ 4 ];
                    aNodes( 4 ) = mNodes[ 5 ];
                    break;
                }
                case( 1 ):
                {
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    aNodes( 2 ) = mNodes[ 6 ];
                    aNodes( 3 ) = mNodes[ 7 ];
                    aNodes( 4 ) = mNodes[ 8 ];
                    break;
                }
                case( 2 ):
                {
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 0 ];
                    aNodes( 2 ) = mNodes[ 9 ];
                    aNodes( 3 ) = mNodes[ 10 ];
                    aNodes( 4 ) = mNodes[ 11 ];
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
        LagrangeElement< 15, 3, 3, 3, 1 >::get_edges_of_facet(
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


#endif //BELFEM_CL_ELEMENT_TRI15_HPP
