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

#include "fn_Graph_spfa.hpp"

#include "assert.hpp"
#include "cl_DynamicBitset.hpp"

namespace belfem
{
    namespace graph
    {
//------------------------------------------------------------------------------

        namespace spfa
        {
            /**
             * search for a strictly negative cycle in the predecessor graph
             * ( vertex -> tail of its parent arc ) and collect its arcs
             * in cycle order. Returns true if one was found. A cycle of
             * total weight zero is NOT an infeasibility proof ( the
             * invariant only gives theta(head) >= theta(tail) + w, so a
             * parent cycle sums to <= 0 ) and is skipped.
             */
            static bool
            extract_parent_cycle(
                    const index_t           aStartVertex,
                    const index_t           aNumVertices,
                    const Cell< index_t > & aArcTails,
                    const Cell< int64_t > & aArcWeights,
                    const Cell< index_t > & aParentArcs,
                    Cell< index_t >       & aNegativeCycle )
            {
                DynamicBitset tOnWalk( aNumVertices );
                DynamicBitset tDone( aNumVertices );

                // scratch to unmark the walk afterwards
                Cell< index_t > tWalk ;

                // fast path first: walk back from the triggering vertex.
                // fallback: walk back from every parented vertex.
                for ( index_t s = 0 ; s <= aNumVertices; ++s )
                {
                    const index_t tStart = ( s == 0 ) ? aStartVertex : s - 1 ;

                    if ( tDone.test( tStart ) || aParentArcs( tStart ) == gNoIndex )
                    {
                        continue ;
                    }

                    tWalk.clear() ;
                    index_t tVertex = tStart ;

                    while ( true )
                    {
                        if ( tOnWalk.test( tVertex ) )
                        {
                            // found a cycle: collect arcs starting at tVertex
                            const index_t tCycleStart = tVertex ;
                            int64_t tCycleWeight = 0 ;
                            do
                            {
                                const index_t tArc = aParentArcs( tVertex );
                                aNegativeCycle.push( tArc );
                                tCycleWeight += aArcWeights( tArc );
                                tVertex = aArcTails( tArc );
                            }
                            while ( tVertex != tCycleStart );

                            if ( tCycleWeight >= 0 )
                            {
                                // zero-weight parent cycle: no proof, treat
                                // this walk as a dead end and keep scanning
                                aNegativeCycle.clear() ;
                                break ;
                            }

                            // arcs were collected against arc direction: reverse
                            const index_t tNum = aNegativeCycle.size() ;
                            for ( index_t k = 0; k < tNum / 2; ++k )
                            {
                                const index_t tSwap = aNegativeCycle( k );
                                aNegativeCycle( k ) = aNegativeCycle( tNum - 1 - k );
                                aNegativeCycle( tNum - 1 - k ) = tSwap ;
                            }
                            return true ;
                        }

                        if ( tDone.test( tVertex ) || aParentArcs( tVertex ) == gNoIndex )
                        {
                            // dead end: no cycle through this walk
                            break ;
                        }

                        tOnWalk.set( tVertex );
                        tWalk.push( tVertex );
                        tVertex = aArcTails( aParentArcs( tVertex ) );
                    }

                    // unmark the walk and never visit these vertices again
                    for ( index_t tV : tWalk )
                    {
                        tOnWalk.reset( tV );
                        tDone.set( tV );
                    }
                }

                return false ;
            }
        }

//------------------------------------------------------------------------------

        bool
        spfa_difference_constraints(
                const index_t           aNumVertices,
                const Cell< index_t > & aArcTails,
                const Cell< index_t > & aArcHeads,
                const Cell< int64_t > & aArcWeights,
                Cell< int64_t >       & aTheta,
                Cell< index_t >       & aNegativeCycle )
        {
            const index_t tNumArcs = aArcTails.size() ;

            BELFEM_ASSERT( aArcHeads.size() == tNumArcs && aArcWeights.size() == tNumArcs,
                    "spfa_difference_constraints: arc array sizes do not match ( %lu heads, %lu weights, %lu tails )",
                    ( long unsigned int ) aArcHeads.size(),
                    ( long unsigned int ) aArcWeights.size(),
                    ( long unsigned int ) tNumArcs );

            aTheta.set_size( aNumVertices, 0 );
            aNegativeCycle.clear() ;

            if ( tNumArcs == 0 || aNumVertices == 0 )
            {
                return true ;
            }

            // a negative self arc is a length-one negative cycle
            for ( index_t a = 0; a < tNumArcs; ++a )
            {
                BELFEM_ASSERT( aArcTails( a ) < aNumVertices && aArcHeads( a ) < aNumVertices,
                        "spfa_difference_constraints: arc %lu references vertex out of range",
                        ( long unsigned int ) a );

                if ( aArcTails( a ) == aArcHeads( a ) && aArcWeights( a ) < 0 )
                {
                    aNegativeCycle.push( a );
                    return false ;
                }
            }

            // bucket the arcs by tail vertex ( CSR layout )
            Cell< index_t > tOffsets( aNumVertices + 1, 0 );
            for ( index_t a = 0; a < tNumArcs; ++a )
            {
                ++tOffsets( aArcTails( a ) + 1 );
            }
            for ( index_t v = 0; v < aNumVertices; ++v )
            {
                tOffsets( v + 1 ) += tOffsets( v );
            }
            Cell< index_t > tArcsOfTail( tNumArcs, 0 );
            {
                Cell< index_t > tCursor( tOffsets );
                for ( index_t a = 0; a < tNumArcs; ++a )
                {
                    tArcsOfTail( tCursor( aArcTails( a ) )++ ) = a ;
                }
            }

            // vertices that participate in at least one arc
            DynamicBitset tActive( aNumVertices );
            for ( index_t a = 0; a < tNumArcs; ++a )
            {
                tActive.set( aArcTails( a ) );
                tActive.set( aArcHeads( a ) );
            }

            // solver state
            Cell< index_t > tParentArcs( aNumVertices, gNoIndex );
            Cell< index_t > tUpdateCount( aNumVertices, 0 );
            DynamicBitset   tInQueue( aNumVertices );

            // active vertices; their count is also the ring buffer capacity,
            // since only arc endpoints can ever ( re )enter the queue
            Cell< index_t > tActiveVerts ;
            tActive.where( tActiveVerts, false ); // dense: most vertices carry arcs

            const index_t tNumActive = tActiveVerts.size() ;

            Cell< index_t > tQueue( tNumActive, 0 );
            index_t tQueueHead  = 0 ;
            index_t tQueueTail  = 0 ;
            index_t tQueueCount = 0 ;

            // warm start: instead of the cold virtual super source ( all
            // theta = 0, every arc of the excess landscape violated ), grow a
            // BFS spanning forest over the out-arc adjacency and set
            // theta( head ) = theta( tail ) + w along tree arcs. Both arcs of
            // a difference-constraint pair are then satisfied at init, so
            // only off-tree arcs carry violations and the relax work scales
            // with the actual excess instead of the potential range. Tree
            // arcs enter the predecessor graph with tight invariant
            // theta( head ) = theta( tail ) + w.
            {
                DynamicBitset tVisited( aNumVertices );

                for ( index_t tRoot : tActiveVerts )
                {
                    if ( tVisited.test( tRoot ) )
                    {
                        continue ;
                    }
                    tVisited.set( tRoot );

                    // the ring buffer doubles as the BFS queue ( it is empty
                    // here, and BFS holds each active vertex at most once )
                    tQueue( tQueueTail ) = tRoot ;
                    if ( ++tQueueTail == tNumActive ) tQueueTail = 0 ;
                    ++tQueueCount ;

                    while ( tQueueCount > 0 )
                    {
                        const index_t tU = tQueue( tQueueHead );
                        if ( ++tQueueHead == tNumActive ) tQueueHead = 0 ;
                        --tQueueCount ;

                        for ( index_t k = tOffsets( tU ); k < tOffsets( tU + 1 ); ++k )
                        {
                            const index_t tArc = tArcsOfTail( k );
                            const index_t tV   = aArcHeads( tArc );

                            if ( ! tVisited.test( tV ) )
                            {
                                tVisited.set( tV );
                                aTheta( tV ) = aTheta( tU ) + aArcWeights( tArc );
                                tParentArcs( tV ) = tArc ;

                                tQueue( tQueueTail ) = tV ;
                                if ( ++tQueueTail == tNumActive ) tQueueTail = 0 ;
                                ++tQueueCount ;
                            }
                        }
                    }
                }
            }

            // enqueue every vertex with a violated outgoing arc
            for ( index_t a = 0; a < tNumArcs; ++a )
            {
                const index_t tU = aArcTails( a );

                if ( aTheta( tU ) + aArcWeights( a ) < aTheta( aArcHeads( a ) )
                     && ! tInQueue.test( tU ) )
                {
                    tQueue( tQueueTail ) = tU ;
                    if ( ++tQueueTail == tNumActive ) tQueueTail = 0 ;
                    ++tQueueCount ;
                    tInQueue.set( tU );
                }
            }

            // amortized negative-cycle detection: a cycle in the predecessor
            // graph proves infeasibility the moment it forms, so scan for one
            // every ~4 * tNumActive updates ( O( 1 ) amortized per update )
            // instead of waiting for one vertex to accumulate tNumActive
            // updates. The per-vertex trigger below stays as a backstop.
            const uint64_t tScanInterval =
                    4 * ( uint64_t ) tNumActive + 1024 ;
            uint64_t tUpdatesSinceScan = 0 ;
            index_t  tLastUpdated      = 0 ;

            // hard failsafe sized past SPFA's O( V * E ) worst case; hitting it
            // means an implementation defect, not a hard instance
            const uint64_t tOpCap =
                    ( uint64_t ) ( tNumActive + 1 ) * ( uint64_t ) ( tNumArcs + 1 )
                    + 1000 ;
            uint64_t tNumDequeues = 0 ;

            while ( tQueueCount > 0 )
            {
                const index_t tU = tQueue( tQueueHead );
                if ( ++tQueueHead == tNumActive ) tQueueHead = 0 ;
                --tQueueCount ;
                tInQueue.reset( tU );

                BELFEM_ERROR( ++tNumDequeues <= tOpCap,
                        "spfa_difference_constraints: operation cap exceeded ( implementation defect? )" );

                const int64_t tThetaU = aTheta( tU );

                for ( index_t k = tOffsets( tU ); k < tOffsets( tU + 1 ); ++k )
                {
                    const index_t tArc = tArcsOfTail( k );
                    const index_t tV   = aArcHeads( tArc );
                    const int64_t tNew = tThetaU + aArcWeights( tArc );

                    if ( tNew < aTheta( tV ) )
                    {
                        aTheta( tV ) = tNew ;
                        tParentArcs( tV ) = tArc ;
                        tLastUpdated = tV ;
                        ++tUpdatesSinceScan ;

                        if ( ++tUpdateCount( tV ) > tNumActive )
                        {
                            // backstop trigger: verify-or-continue ( only a
                            // predecessor-graph cycle proves infeasibility )
                            if ( spfa::extract_parent_cycle(
                                    tV, aNumVertices, aArcTails, aArcWeights,
                                    tParentArcs, aNegativeCycle ) )
                            {
                                return false ;
                            }

                            tUpdateCount( tV ) = 0 ;
                        }

                        if ( ! tInQueue.test( tV ) )
                        {
                            tQueue( tQueueTail ) = tV ;
                            if ( ++tQueueTail == tNumActive ) tQueueTail = 0 ;
                            ++tQueueCount ;
                            tInQueue.set( tV );
                        }
                    }
                }

                if ( tUpdatesSinceScan >= tScanInterval )
                {
                    tUpdatesSinceScan = 0 ;

                    if ( spfa::extract_parent_cycle(
                            tLastUpdated, aNumVertices, aArcTails, aArcWeights,
                            tParentArcs, aNegativeCycle ) )
                    {
                        return false ;
                    }
                }
            }

            return true ;
        }

//------------------------------------------------------------------------------
    }
}
