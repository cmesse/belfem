//
// Created by Christian Messe on 14.07.20.
//

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "commtools.hpp"

#include "cl_Logger.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_SpMatrix.hpp"
#include "cl_Solver.hpp"
#include "fn_r2.hpp"

using namespace belfem;

extern belfem::Cell< belfem::graph::Vertex * > gGraph;

TEST( SOLVER, PETSC )
{
    if( comm_size() <= 4 )
    {
        Matrix <real> tA( 5, 5, 0.0 );

        tA( 0, 0 ) =  1.0;
        tA( 0, 1 ) = -1.0;
        tA( 0, 3 ) = -3.0;
        tA( 1, 0 ) = -2.0;
        tA( 1, 1 ) =  5.0;
        tA( 2, 2 ) =  4.0;
        tA( 2, 3 ) =  6.0;
        tA( 2, 4 ) =  4.0;
        tA( 3, 0 ) = -4.0;
        tA( 3, 2 ) =  2.0;
        tA( 3, 3 ) =  7.0;
        tA( 4, 1 ) =  8.0;
        tA( 4, 4 ) = -5.0;


        SpMatrix tM( tA, SpMatrixType::CSR );

        Vector <real> tY = { -13., 8., 56., 30., -9. };
        Vector <real> tX( 5, 0.1 );

        Vector< real > tExpect = { 1., 2., 3., 4., 5. };

        Solver tSolver( SolverType::PETSC );

        tSolver.solve( tM, tX, tY );

        tSolver.free();

        if( comm_rank() == 0 )
        {
            EXPECT_NEAR( r2( tX, tExpect ), 1.0, 1e-7 );
        }
    }
}
