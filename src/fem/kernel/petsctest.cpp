//
// Created by Christian Messe on 12.07.20.
//


// thermostructural analysis of rocked engines and spacecraft

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "commtools.hpp"

#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_SpMatrix.hpp"
#include "petsctools.hpp"
#include "cl_Solver.hpp"

using namespace belfem;

Communicator gComm;
Logger       gLog( 5 );


int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm = Communicator( &argc, &argv );

    Matrix< real > tB( 5, 5, 0.0 ) ;

    tB( 0, 0 ) =  1.0 ;
    tB( 0, 1 ) = -1.0 ;
    tB( 0, 3 ) = -3.0 ;
    tB( 1, 0 ) = -2.0 ;
    tB( 1, 1 ) =  5.0 ;
    tB( 2, 2 ) =  4.0 ;
    tB( 2, 3 ) =  6.0 ;
    tB( 2, 4 ) =  4.0 ;
    tB( 3, 0 ) = -4.0 ;
    tB( 3, 2 ) =  2.0 ;
    tB( 3, 3 ) =  7.0 ;
    tB( 4, 1 ) =  8.0 ;
    tB( 4, 4 ) = -5.0 ;

    if( comm_rank() == 0 )
    {
        tB.print( "M" );
    }
    SpMatrix tM( tB, SpMatrixType::CSR );

    Vector< real > tY = { -13., 8., 56., 30., -9. } ;
    Vector< real > tX( 5 );

    Solver tSolver( SolverType::PARDISO ) ;

    tSolver.solve( tM, tX, tY ) ;

    if( comm_rank() == 0 )
    {
        tX.print( "X" );
    }
    tSolver.free() ;

    return  gComm.finalize();
}