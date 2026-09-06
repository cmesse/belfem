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

#ifndef BELFEM_FN_CHECK_FACET_ORIENTATION_HPP
#define BELFEM_FN_CHECK_FACET_ORIENTATION_HPP
#include "assert.hpp"
#include "cl_Facet.hpp"
#include "cl_Element.hpp"
namespace belfem
{
    namespace mesh
    {
        /**
         * Checks whether two adjacent facets have consistent orientation.
         *
         * In 3D, facets are polygons and adjacency means a shared edge:
         * consistent outward normals traverse that edge in opposite
         * directions:
         *
         *     Facet A: ... → P → Q → ...
         *     Facet B: ... → Q → P → ...   (opposite → consistent)
         *
         * If both traverse the shared edge in the same direction,
         * the normals point in opposite directions (inverted).
         *
         * For two-node facets ( 2D line segments ), two distinct facets
         * share a single node, never an edge — consistency there means one
         * traversal continues the other through the shared node:
         *
         *     Facet A: P → Q
         *     Facet B: Q → R          (chain continues → consistent)
         *
         * If both segments start, or both end, at the shared node, the
         * normals are inverted. Coincident twins reduce to the 3D rule:
         * reversed twin → consistent, identical twin → inverted. Nodes
         * are compared via original() because the 2D facet-to-facet
         * connectivity matches duplicates onto their originals.
         *
         * @param[in] aA  first facet
         * @param[in] aB  second facet (must share an edge, or in 2D a node, with aA)
         * @return true if normals are consistent, false if inverted
         */
        inline bool
        check_facet_orientation( const Facet * aA, const Facet * aB )
        {
            const Element * tA = aA->element();
            const Element * tB = aB->element();
            uint n = tA->number_of_corner_nodes();
            uint m = tB->number_of_corner_nodes();

            if ( n == 2 && m == 2 )
            {
                id_t tA0 = tA->node( 0 )->original()->id();
                id_t tA1 = tA->node( 1 )->original()->id();
                id_t tB0 = tB->node( 0 )->original()->id();
                id_t tB1 = tB->node( 1 )->original()->id();

                // B continues A's traversal ( or reversed twin ): consistent
                if ( tA1 == tB0 || tA0 == tB1 )
                {
                    return true;
                }

                // both start or both end at the shared node
                // ( or identical twin ): inverted
                if ( tA0 == tB0 || tA1 == tB1 )
                {
                    return false;
                }
            }
            else
            {
                // loop over all edge pairs (p→q) from A and (r→s) from B
                id_t p ;
                id_t q = tA->node( n-1 )->id();
                id_t r ;
                id_t s ;

                for ( uint i = 0; i < n; ++i )
                {
                    p = q ;
                    q = tA->node( i )->id();
                    s = tB->node( m-1 )->id();
                    for ( uint j = 0; j < m; ++j )
                    {
                        r = s ;
                        s = tB->node( j )->id();

                        // A traverses P→Q, B traverses Q→P: consistent normals
                        if ( p == s && q == r )
                        {
                            return true;
                        }

                        // A traverses P→Q, B traverses P→Q: inverted normals
                        if ( p == r && q == s )
                        {
                            return false;
                        }
                    }
                }
            }

            BELFEM_ERROR( false, "Could not determine orientation of facets." );
            return false;
        }
    }
}
#endif //BELFEM_FN_CHECK_FACET_ORIENTATION_HPP