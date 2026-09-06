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


#include "fn_compute_permutation.hpp"
#include "fn_create_graph_from_matrix.hpp"
#include "fn_Graph_METIS.hpp"
#include "assert.hpp"
#include "cl_Communicator.hpp"

namespace belfem
{
    namespace sparse
    {
        void
        compute_permutation(
            const SpMatrix & aMatrix ,
            Graph & aGraph ,
            Cell< int_t > & aForwardPermutation ,
            Cell< int_t > & aBackwardPermutation ,
            Cell< int_t > & aIndexPermutation )
        {



            BELFEM_ASSERT( aGraph.size() == 0, "graph has already been allocated" );

            create_graph_from_matrix( aMatrix, aGraph );

            graph::metis_ndp( aGraph, gComm.size() );

            // Validate graph structure after the METIS nested-dissection permutation
            int_t tN = aGraph.size();

            aForwardPermutation.set_size( tN );
            aBackwardPermutation.set_size( tN );
            aIndexPermutation.set_size( aMatrix.number_of_nonzeros() );

            // Build forward and backward permutations for vertices
            for( graph::Vertex * tVertex : aGraph )
            {
                aForwardPermutation( tVertex->index() ) = tVertex->id();
                aBackwardPermutation( tVertex->id() ) = tVertex->index() ;
            }

            // Build index permutation by walking through the permuted graph structure
            // The permuted matrix (created from aGraph) has its nonzeros in a different order
            // We need to map: new_position -> old_position in the values array

            bool tIsCSR = aMatrix.type() == SpMatrixType::CSR;
            const int_t * tOldPointers = aMatrix.pointers();
            const int_t * tOldIndices = aMatrix.indices();

            index_t tNewPos = 0;
            for( graph::Vertex * tVertex : aGraph )
            {
                int_t iOld = tVertex->id();  // Old row/col index (original matrix)

                for ( uint k=0; k<tVertex->number_of_vertices(); ++k )
                {
                    graph::Vertex * tNeighbor = tVertex->vertex( k );
                    int_t jOld = tNeighbor->id();  // Old col/row index

                    // Find position in original matrix
                    // For CSR: search row iOld for column jOld
                    // For CSC: search column iOld for row jOld
                    int_t iSearch = tIsCSR ? iOld : jOld;
                    int_t jSearch = tIsCSR ? jOld : iOld;

                    int_t start = tOldPointers[iSearch];
                    int_t end = tOldPointers[iSearch + 1];

                    int_t tOldPos = 0;
                    bool found = false;
                    for( int_t pos = start; pos < end; ++pos )
                    {
                        if( tOldIndices[pos] == jSearch )
                        {
                            tOldPos = pos;
                            found = true;
                            break;
                        }
                    }

                    BELFEM_ERROR( found, "could not find old position" );
                    aIndexPermutation( tNewPos ) = tOldPos;
                    ++tNewPos;
                }
            }
        }
    }
}
