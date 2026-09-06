/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef CL_ELEMENT_PENTA6TS_HPP
#define CL_ELEMENT_PENTA6TS_HPP

//
// Created by Christian Messe on 2019-08-04.
//

#ifndef BELFEM_CL_ELEMENT_PENTA6TS_HPP
#define BELFEM_CL_ELEMENT_PENTA6TS_HPP

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
        ElementTemplate< 6, 6, 6, 5, 1 >::type() const
        {
            return ElementType::PENTA6TS;
        }

//------------------------------------------------------------------------------

        template <>
        uint
        ElementTemplate< 6, 6, 6, 5, 1 >::dimension() const
        {
            return 3 ;
        }

//------------------------------------------------------------------------------

        template <>
        bool
        ElementTemplate< 6, 6, 6, 5, 1 >::is_thinshell() const
        {
            return true ;
        }

//------------------------------------------------------------------------------

        template <>
        void
        ElementTemplate< 6, 6, 6, 5, 1 >::get_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
        {
            switch( aFacetIndex )
            {
                case 0 :
                {
                    aNodes.set_size( 4, nullptr );
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 1 ];
                    aNodes( 2 ) = mNodes[ 4 ];
                    aNodes( 3 ) = mNodes[ 3 ];
                    break;
                }
                case 1 :
                {
                    aNodes.set_size( 4, nullptr );
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    aNodes( 2 ) = mNodes[ 5 ];
                    aNodes( 3 ) = mNodes[ 4 ];
                    break;
                }
                case 2 :
                {
                    aNodes.set_size( 4, nullptr );
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 3 ];
                    aNodes( 2 ) = mNodes[ 5 ];
                    aNodes( 3 ) = mNodes[ 2 ];
                    break;
                }
                case 3 :
                {
                    aNodes.set_size( 3, nullptr );
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    aNodes( 2 ) = mNodes[ 1 ];
                    break;
                }
                case 4 :
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
        ElementTemplate< 6, 6, 6, 5, 1 >::get_corner_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
        {
            this->get_nodes_of_facet( aFacetIndex, aNodes );
        }

//------------------------------------------------------------------------------

        template <>
        void
        ElementTemplate< 6, 6, 6, 5, 1 >::get_nodes_of_edge( const uint aEdgeIndex, Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 2, nullptr );

            switch( aEdgeIndex )
            {
                case 0 :
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 1 ];
                    break;
                }
                case 1 :
                {
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    break;
                }
                case 2 :
                {
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 0 ];
                    break;
                }
                case 3 :
                {
                    aNodes( 0 ) = mNodes[ 3 ];
                    aNodes( 1 ) = mNodes[ 4 ];
                    break;
                }
                case 4 :
                {
                    aNodes( 0 ) = mNodes[ 4 ];
                    aNodes( 1 ) = mNodes[ 5 ];
                    break;
                }
                case 5 :
                {
                    aNodes( 0 ) = mNodes[ 5 ];
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
        ElementTemplate< 6, 6, 6, 5, 1 >::get_edges_of_facet(
                const uint aFacetIndex, Cell< Edge * > & aEdges )
        {
            switch( aFacetIndex )
            {
                case 3 :
                {
                    aEdges.set_size( 3, nullptr );
                    aEdges( 0 ) = mEdges[ 2 ];
                    aEdges( 1 ) = mEdges[ 1 ];
                    aEdges( 2 ) = mEdges[ 0 ];
                    break;
                }
                case 4 :
                {
                    aEdges.set_size( 3, nullptr );
                    aEdges( 0 ) = mEdges[ 3 ];
                    aEdges( 1 ) = mEdges[ 4 ];
                    aEdges( 2 ) = mEdges[ 5 ];
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

#endif //BELFEM_CL_ELEMENT_PENTA6TS_HPP


#endif //CL_ELEMENT_PENTA6TS_HPP
