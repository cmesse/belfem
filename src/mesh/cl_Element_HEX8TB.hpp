/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_CL_ELEMENT_HEX8TB_HPP
#define BELFEM_CL_ELEMENT_HEX8TB_HPP


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
        ElementTemplate< 8, 8, 4, 6, 1 >::type() const
        {
            return ElementType::HEX8TB;
        }

//------------------------------------------------------------------------------

        template <>
        uint
        ElementTemplate< 8, 8, 4, 6, 1 >::dimension() const
        {
            return 3 ;
        }

//------------------------------------------------------------------------------

        template <>
        void
        ElementTemplate< 8, 8, 4, 6, 1 >::get_nodes_of_edge( const uint aEdgeIndex, Cell< Node * > & aNodes )
        {
            aNodes.set_size( 2, nullptr );
            switch ( aEdgeIndex )
            {
                case 0 :
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 1 ];
                    break;
                }
                case 1 :
                {
                    aNodes( 0 ) = mNodes[ 3 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    break;
                }
                case 2 :
                {
                    aNodes( 0 ) = mNodes[ 4 ];
                    aNodes( 1 ) = mNodes[ 5 ];
                    break;
                }
                case 3 :
                {
                    aNodes( 0 ) = mNodes[ 7 ];
                    aNodes( 1 ) = mNodes[ 6 ];
                    break;
                }
                default:
                {
                    this->throw_edge_error( aEdgeIndex );
                }
            }
        }

        template <>
        void
        ElementTemplate< 8, 8, 4, 6, 1 >::get_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 4, nullptr );

            switch( aFacetIndex )
            {
                case 0 :
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 1 ];
                    aNodes( 2 ) = mNodes[ 5 ];
                    aNodes( 3 ) = mNodes[ 4 ];
                    break;
                }
                case 1 :
                {
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    aNodes( 2 ) = mNodes[ 6 ];
                    aNodes( 3 ) = mNodes[ 5 ];
                    break;
                }
                case 2 :
                {
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 3 ];
                    aNodes( 2 ) = mNodes[ 7 ];
                    aNodes( 3 ) = mNodes[ 6 ];
                    break;
                }
                case 3 :
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 4 ];
                    aNodes( 2 ) = mNodes[ 7 ];
                    aNodes( 3 ) = mNodes[ 3 ];
                    break;
                }
                case 4 :
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 3 ];
                    aNodes( 2 ) = mNodes[ 2 ];
                    aNodes( 3 ) = mNodes[ 1 ];
                    break;
                }
                case 5 :
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

        template <>
        void
        ElementTemplate< 8, 8, 4, 6, 1 >::get_corner_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
        {
            this->get_nodes_of_facet( aFacetIndex, aNodes );
        }

        template <>
        void
        ElementTemplate< 8, 8, 4, 6, 1 >::get_edges_of_facet(
                const uint aFacetIndex, Cell< Edge * > & aEdges )
        {
            // allocate the node container
            aEdges.set_size( 2, nullptr );

            switch( aFacetIndex )
            {
                case 0 :
                {
                    // lateral face (0,1,5,4): bottom edge first, then top,
                    // matching the recovery facet's edge slot convention
                    aEdges( 0 ) = mEdges[  0 ];
                    aEdges( 1 ) = mEdges[  2 ];
                    break;
                }
                case 2 :
                {
                    // lateral face (2,3,7,6): bottom edge first, then top
                    aEdges( 0 ) = mEdges[  1 ];
                    aEdges( 1 ) = mEdges[  3 ];
                    break;
                }
                case 4 :
                {
                    aEdges( 0 ) = mEdges[  0 ];
                    aEdges( 1 ) = mEdges[  1 ];
                    break;
                }
                case 5 :
                {
                    aEdges( 0 ) = mEdges[  2 ];
                    aEdges( 1 ) = mEdges[  3 ];
                    break;
                }
                default:
                {
                    // faces 1 and 3 carry no longitudinal dof edges
                    this->throw_facet_error( aFacetIndex );
                }
            }
        }

//-----------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_ELEMENT_HEX8TB_HPP
