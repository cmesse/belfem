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

#ifndef CL_ELEMENT_QUAD9TS_HPP
#define CL_ELEMENT_QUAD9TS_HPP

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
        ElementTemplate< 9, 4, 3, 4, 1 >::type() const
        {
            return ElementType::QUAD9TS;
        }

//------------------------------------------------------------------------------

        template <>
        uint
        ElementTemplate< 9, 4, 3, 4, 1 >::dimension() const
        {
            return 2 ;
        }

//------------------------------------------------------------------------------

        template <>
        bool
        ElementTemplate< 9, 4, 3, 4, 1 >::is_thinshell() const
        {
            return true ;
        }


//------------------------------------------------------------------------------

        template <>
        void
        ElementTemplate< 9, 4, 3, 4, 1 >::get_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
        {
            // allocate the node container
            aNodes.set_size( 3, nullptr );

            switch( aFacetIndex )
            {
                case 0 :
                {
                    aNodes( 0 ) = mNodes[ 0 ];
                    aNodes( 1 ) = mNodes[ 1 ];
                    aNodes( 2 ) = mNodes[ 4 ];
                    break;
                }
                case 1 :
                {
                    aNodes( 0 ) = mNodes[ 1 ];
                    aNodes( 1 ) = mNodes[ 2 ];
                    aNodes( 2 ) = mNodes[ 5 ];
                    break;
                }
                case 2 :
                {
                    aNodes( 0 ) = mNodes[ 2 ];
                    aNodes( 1 ) = mNodes[ 3 ];
                    aNodes( 2 ) = mNodes[ 6 ];
                    break;
                }
                case 3 :
                {
                    aNodes( 0 ) = mNodes[ 3 ];
                    aNodes( 1 ) = mNodes[ 0 ];
                    aNodes( 2 ) = mNodes[ 7 ];
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
        ElementTemplate< 9, 4, 3, 4, 1 >::get_corner_nodes_of_facet( const uint aFacetIndex, Cell< Node * > & aNodes )
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
        ElementTemplate< 9, 4, 3, 4, 1 >::get_edges_of_facet(
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
                    // mEdges[1] is stored with tangent parallel to mEdges[0];
                    // fem::Element::compute_edge_directions resolves the local sign.
                    aEdges.set_size( 1, nullptr );
                    aEdges( 0 ) = mEdges[ 1 ];
                    break ;
                }
                // Side facets 1 and 3 have no Nédélec edge in the thin-shell
                // reduction (mEdges[2] is the center/mid curve, not a side).
                default:
                {
                    this->throw_facet_error( aFacetIndex );
                }
            }
        }

//------------------------------------------------------------------------------
    } /* namespace mesh */
} /* namespace belfem */

#endif //CL_ELEMENT_QUAD9TS_HPP
