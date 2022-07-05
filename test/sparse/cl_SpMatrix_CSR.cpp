//
// Created by Christian Messe on 06.11.19.
//
#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "fn_r2.hpp"
#include "cl_Graph_Vertex.hpp"
#define protected public
#define private   public
#include "cl_SpMatrix.hpp"
#include "cl_Solver.hpp"
#undef protected
#undef private


using namespace belfem;

extern belfem::Cell< belfem::graph::Vertex * > gGraph;

TEST( SPARSE, CSR )
{
    // create the matrix
    SpMatrix tMatrix( gGraph, SpMatrixType::CSR );

    // Vector with pointers
    Vector<int> tPointers( tMatrix.n_rows() + 1 );

    // Vector with indices
    Vector<int> tIndices( tMatrix.number_of_nonzeros());

    // swap indices
    tMatrix.set_indexing_base( SpMatrixIndexingBase::Fortran );
    std::memcpy( tPointers.data(), tMatrix.pointers(), ( tMatrix.n_rows() + 1 ) * sizeof( int ));
    std::memcpy( tIndices.data(), tMatrix.indices(), ( tMatrix.number_of_nonzeros()) * sizeof( int ));
    EXPECT_TRUE( tPointers == Vector<int>( { 1, 4, 6, 9, 12, 14 } ));
    EXPECT_TRUE( tIndices == Vector<int>( { 1, 2, 4, 1, 2, 3, 4, 5, 1, 3, 4, 2, 5 } ));

    tMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );
    std::memcpy( tPointers.data(), tMatrix.pointers(), ( tMatrix.n_rows() + 1 ) * sizeof( int ));
    std::memcpy( tIndices.data(), tMatrix.indices(), ( tMatrix.number_of_nonzeros()) * sizeof( int ));
    EXPECT_TRUE( tPointers == Vector<int>( { 0, 3, 5, 8, 11, 13 } ));
    EXPECT_TRUE( tIndices == Vector<int>( { 0, 1, 3, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4 } ));

    for ( uint k = 0; k < 2; ++k ) // make two runs, one in zero and one in one-based indexing
    {
        // write data
        tMatrix( 0, 0 ) =  1.0;
        tMatrix( 1, 0 ) = -2.0;
        tMatrix( 3, 0 ) = -4.0;
        tMatrix( 0, 1 ) = -1.0;
        tMatrix( 1, 1 ) =  5.0,
        tMatrix( 4, 1 ) =  8.0;
        tMatrix( 2, 2 ) =  4.0;
        tMatrix( 3, 2 ) =  2.0;
        tMatrix( 0, 3 ) = -3.0;
        tMatrix( 2, 3 ) =  6.0;
        tMatrix( 3, 3 ) =  7.0;
        tMatrix( 2, 4 ) =  4.0;
        tMatrix( 4, 4 ) = -5.0;

        // copy values
        Vector<real> tValues = tMatrix.number_of_nonzeros();
        std::memcpy( tValues.data(), tMatrix.data(), ( tMatrix.number_of_nonzeros()) * sizeof( real ));
        EXPECT_TRUE( tValues == Vector<real>( { 1., -1., -3., -2., 5., 4., 6., 4., -4., 2., 7., 8., -5. } ));

        tMatrix.set_indexing_base( SpMatrixIndexingBase::Fortran );

        // Multiplication test
        Vector< real > tX = { 1, 2, 3, 4, 5 };

        Vector< real > tY( 5 );
        tMatrix.multiply( tX, tY );

        EXPECT_TRUE( tY == Vector<real>( { -13., 8., 56., 30., -9., } ));

        tMatrix.set_indexing_base( SpMatrixIndexingBase::Fortran );
    }

    Vector<real> tY = { -13., 8., 56., 30., -9. };
    Vector<real> tX( 5 );
    Vector<real> tExpect = { 1, 2, 3, 4, 5 };

#ifdef BELFEM_SUITESPARSE
    tX.fill( 0 );
    Solver tUMFPACK = Solver( SolverType::UMFPACK ) ;
    tUMFPACK.solve( tMatrix, tX, tY );
    tUMFPACK.free();
    EXPECT_NEAR( r2( tX, tExpect ), 1.0, BELFEM_EPSILON );
#endif

#ifdef BELFEM_MKL
    tX.fill( 0 );
    Solver tPARDISO = Solver( SolverType::PARDISO ) ;
    tPARDISO.solve( tMatrix, tX, tY );
    tPARDISO.free();
    EXPECT_NEAR( r2( tX, tExpect ), 1.0, BELFEM_EPSILON );
#endif

#ifdef BELFEM_MUMPS
tX.fill( 0 );
    Solver tMUMPS = Solver( SolverType::MUMPS ) ;
    tMUMPS.solve( tMatrix, tX, tY );
    tMUMPS.free();
    EXPECT_NEAR( r2( tX, tExpect ), 1.0, BELFEM_EPSILON );
#endif
}
