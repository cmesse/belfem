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
#include "typedefs.hpp"
#include "assert.hpp"

#include "fn_create_graph_from_matrix.hpp"
namespace belfem
{
    namespace sparse
    {
        void
            create_graph_from_matrix(
                const SpMatrix      & aMatrix,
                Graph & aGraph )
        {
            BELFEM_ASSERT(
                aMatrix.type() == SpMatrixType::CSR || aMatrix.type() == SpMatrixType::CSC,
                "Invalid matrix type" );

            BELFEM_ERROR( aMatrix.indexing_base() == 0,
                "set indexing base of this matrix to Cpp (0-based) before calling create_graph_from_matrix()" );
            BELFEM_ASSERT( aGraph.size() == 0, "Graph has already been initialized" );

            id_t tN = aMatrix.type() == SpMatrixType::CSR ?
                         aMatrix.n_rows() : aMatrix.n_cols();

            aGraph.set_size( tN, nullptr );
            const int_t * tPointers = aMatrix.pointers();
            for ( id_t k=0; k<tN; ++k )
            {
                graph::Vertex * tVertex = new graph::Vertex();
                tVertex->set_id( k );
                tVertex->set_index( k );
                int_t s = tPointers[ k+1 ] - tPointers[ k ];
                if ( s > 0 )
                {
                    tVertex->init_vertex_container( s );
                }
                aGraph( k ) = tVertex;
            }

            const int_t * tIndices  = aMatrix.indices();
            for ( id_t k=0; k<tN; ++k )
            {
                for ( int_t j=tPointers[ k ]; j<tPointers[ k+1 ]; ++j )
                {
                    aGraph( k )->insert_vertex( aGraph( tIndices[ j ] ) );
                }
            }

        }
    }
}