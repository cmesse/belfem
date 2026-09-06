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

#include <algorithm>
#include "typedefs.hpp"
#include "graph_typedefs.hpp"
#include "graphtools.hpp"

#include "commtools.hpp"
#include "assert.hpp"
#include "cl_Timer.hpp"
#include "cl_Logger.hpp"
#include "fn_Graph_METIS.hpp"
#include "fn_Graph_ParMETIS.hpp"

namespace belfem
{

    namespace graph
    {
        void
        parmetis_nd( Graph & aGraph )
        {
#ifdef BELFEM_PARMETIS
            if ( comm_size() < 2 )
            {
                metis_ndp( aGraph, comm_size() );
                return;
            }
            Vector< metis_t > tDistribution ;
            Vector< metis_t > tVertices ;
            Vector< metis_t > tEdges ;

            Timer tTimer;

            metis_t tSeed = 0 ;

            proc_t tCommRank = comm_rank();
            proc_t tCommSize = comm_size();

            // ParMETIS_V3_NodeND refuses a processor that owns no vertex
            // ( libparmetis/weird.c, CheckInputsNodeND ). Root knows the
            // distribution; the verdict is broadcast in BOTH rank branches
            // before any share/receive, so every rank takes the same path
            int tOk = 1 ;

            if ( tCommRank == 0 )
            {
                BELFEM_ERROR( aGraph.size() > 0, "Graph is empty" );


                Cell< Vector< metis_t > > tAllVertices ;
                Cell< Vector< metis_t > > tAllEdges ;

                build_pargraph_adjacency( aGraph, tDistribution, tAllVertices, tAllEdges );

                for ( proc_t p = 0; p < tCommSize; ++p )
                {
                    if ( tDistribution( p + 1 ) - tDistribution( p ) < 1 )
                    {
                        tOk = 0 ;
                    }
                }

                tSeed = static_cast< metis_t >( aGraph.size() );
                comm_barrier();
                broadcast( tOk );
                if ( tOk == 0 )
                {
                    message( InfoLevel::Minimal,
                        "    WARNING: parmetis_nd: %lu vertices leave at least one of %d ranks empty, "
                        "which ParMETIS refuses. Falling back to serial METIS nested dissection.",
                        ( long unsigned int ) aGraph.size(), ( int ) tCommSize );
                    metis_ndp( aGraph, tCommSize );
                    return;
                }
                broadcast( tSeed );
                share( tDistribution );
                distribute(  tAllVertices );
                distribute(  tAllEdges );
                tVertices = std::move( tAllVertices( 0 ) );
                tEdges = std::move( tAllEdges( 0 ) );
            }
            else
            {
                comm_barrier();
                broadcast( tOk );
                if ( tOk == 0 )
                {
                    return;
                }
                broadcast( tSeed );
                receive( tDistribution );
                receive( tVertices );
                receive( tEdges );
            }

            // Set ParMETIS options
            metis_t tOptions[3] = { 1, 0, 0 };

            if ( gLog.info_level() >= static_cast< uint >( InfoLevel::Everything ) )
            {
                tOptions[ 0 ] = 3;
            }
            // random seed
            tOptions[ 2 ] = tSeed ;

            metis_t tNumFlag = 0;  // C-style numbering

            Vector< metis_t > tMyPermutation ( tVertices.length() - 1 );
            // ParMETIS documents sizes as an array of length 2*npes (subdomain
            // sizes, then the separator sizes of every level); a shorter buffer
            // is written past its end
            Vector< metis_t > tMySizes ( 2 * static_cast< size_t >( tCommSize ), 0 );

            auto tWorld = gComm.world();

            // Call ParMETIS nested dissection
            int tStatus = ParMETIS_V3_NodeND(
                tDistribution.data(),       // Vertex distribution
                tVertices.data(),          // Local row pointers
                tEdges.data(),        // Local column indices (global numbering)
                &tNumFlag,             // 0 = C numbering
                tOptions,              // Options array
                tMyPermutation.data(), // Output: local permutation
                tMySizes.data(),       // Output: separator sizes
                &tWorld                // MPI communicator
            );

            BELFEM_ERROR( tStatus == METIS_OK,
             "ParMETIS_NodeND failed with status %d : %s",
             ( unsigned int ) tStatus,
             metis_status( tStatus ).c_str() );

            if ( comm_rank() == 0 )
            {
                message( InfoLevel::Verbose, "    ... time for reordering             : %u ms",
                  ( unsigned int ) tTimer.stop() );

                Cell< Vector< metis_t > > tAllPermutations( comm_size(), {} );

                collect(  tAllPermutations );
                tAllPermutations( 0 ) = std::move( tMyPermutation );

                tMyPermutation.set_size( aGraph.size() );

                Vector< index_t > tCount( gComm.size(), 0 );

                for ( Vertex * tVertex : aGraph )
                {
                    tMyPermutation( tVertex->index() ) = tAllPermutations( tVertex->owner() )( tCount( tVertex->owner() )++);
                }
                apply_graph_permutation( aGraph, tMyPermutation );
            }
            else
            {
                send( tMyPermutation );
            }

#else
            BELFEM_ERROR( false, "We are not linked against ParMETIS");
#endif
        }
    }
}
