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

#ifndef CL_ELEMENT_QUAD4TS_HPP
#define CL_ELEMENT_QUAD4TS_HPP


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
        ElementTemplate< 4, 4, 2, 4, 1 >::type() const
        {
            return ElementType::QUAD4TS;
        }

//------------------------------------------------------------------------------

        template <>
        uint
        ElementTemplate< 4, 4, 2, 4, 1 >::dimension() const
        {
            return 2 ;
        }

//------------------------------------------------------------------------------

        template <>
        bool
        ElementTemplate< 4, 4, 2, 4, 1 >::is_thinshell() const
        {
            return true ;
        }

//------------------------------------------------------------------------------

        template <>
        void
        ElementTemplate< 4, 4, 2, 4, 1 >::get_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 2, nullptr );

            switch( aFacetIndex )
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
                    aNodes( 1 ) = mNodes[ 3 ];
                    break;
                }
                case 3 :
                {
                    aNodes( 0 ) = mNodes[ 3 ];
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
        ElementTemplate< 4, 4, 2, 4, 1 >::get_corner_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
        {
            this->get_nodes_of_facet( aFacetIndex, aNodes );
        }

//------------------------------------------------------------------------------

        template <>
        void
        ElementTemplate< 4, 4, 2, 4, 1 >::get_nodes_of_edge( const uint aEdgeIndex, Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 2, nullptr );

            // QUAD4TS has only 2 edges corresponding to the bottom and top curve
            // (the two Nédélec-carrying edges of the thin shell)
            switch( aEdgeIndex )
            {
                case 0 :
                {
                    // bottom curve edge: node 0 -> node 1
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 1 ];
                    break;
                }
                case 1 :
                {
                    // top curve edge: node 3 -> node 2 (same tangent direction
                    // as bottom, so parallel edges can share Nédélec DOFs)
                    aNodes( 0 ) = mNodes[ 3 ];
                    aNodes( 1 ) = mNodes[ 2 ];
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
        ElementTemplate< 4, 4, 2, 4, 1 >::get_edges_of_facet(
                const uint aFacetIndex, Cell< Edge * > & aEdges )
        {
            switch( aFacetIndex )
            {
                case 0 :
                {
                    // bottom facet uses the bottom curve edge (mEdges[0])
                    aEdges.set_size( 1, nullptr );
                    aEdges( 0 ) = mEdges[ 0 ];
                    break ;
                }
                case 2 :
                {
                    // top facet uses the top curve edge (mEdges[1]).
                    // mEdges[1] is stored as {3,2} but the facet traversal
                    // here is {2,3}; fem::Element::compute_edge_directions resolves
                    // the local sign from the node-ID comparison.
                    aEdges.set_size( 1, nullptr );
                    aEdges( 0 ) = mEdges[ 1 ];
                    break ;
                }
                // Side facets 1 and 3 have no Nédélec edge in the thin-shell
                // reduction (analogous to PENTA6TS side-quad facets 0/1/2).
                default:
                {
                    this->throw_facet_error( aFacetIndex );
                }
            }
        }

//------------------------------------------------------------------------------
    } /* namespace mesh */
} /* namespace belfem */

#endif //CL_ELEMENT_QUAD4TS_HPP
