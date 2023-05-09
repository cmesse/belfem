//
// Created by Christian Messe on 09.07.20.
//

// thermostructural analysis of rocked engines and spacecraft

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
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
    gComm = Communicator( argc, argv );

    /*uint tN = 4 ;

    Matrix< real > tK( tN, tN, 0.0 );

    for( uint k=0; k<tN; ++k )
    {
        tK( k,k ) = 2.0 ;
        if( k>0 ) tK( k, k-1 ) = -1.0 ;
        if( k<tN-1) tK( k,k+1) = -1.0 ;
    }

    tK.print("K"); */


    Matrix< real > tK( 5, 5, 0.0 );

    /*tK( 0, 0 ) =  1.0 ;
    tK( 0, 1 ) = -1.0 ;
    tK( 0, 3 ) = -3.0 ;
    tK( 1, 0 ) = -2.0 ;
    tK( 1, 1 ) =  5.0 ;
    tK( 2, 2 ) =  4.0 ;
    tK( 2, 3 ) =  6.0 ;
    tK( 2, 4 ) =  4.0 ;
    tK( 3, 0 ) = -4.0 ;
    tK( 3, 2 ) =  2.0 ;
    tK( 3, 3 ) =  7.0 ;
    tK( 4, 1 ) =  8.0 ;
    tK( 4, 4 ) = -5.0 ;*/
    tK(0,0) = 8.2 ;
    tK(0,1) = 0.1 ;
    tK(0,4) = 3.1 ;
    tK(1,0) = 1e-12 ;
    tK(1,2) = -4.8 ;
    tK(2,0) = 6.2 ;
    tK(2,1) = 1.1 ;
    tK(2,3) = 2.6 ;
    tK(3,2) = -1.0 ;
    tK(4,3) = 99.9 ;
    tK(4,4) = 4.0 ;

    uint tN = 5 ;
    if( gComm.rank() == 0 )
    {
        tK.print( "K" );
    }

    SpMatrix tM( tK, SpMatrixType::CSR );
    //tM.print("M");

    Vector< real > tX( tN );
    for( uint k=0; k<tN; ++k )
    {
        tX( k ) = ( real ) ( k + 1 ) ;
    }

    Vector< real > tB( tK * tX );
    if( gComm.rank() == 0 )
    {
        tB.print("B");

    }
    //tX.print("x");
    //tB.print("b");

    //tM.print("M");

    // tX *= 0.9 ;

    Solver tSolver( SolverType::MUMPS );

    //tSolver.set_petsc( Preconditioner::LU, KrylovMethod::PREONLY, 1e-8 );
    //tSolver.set_petsc( Preconditioner::GAMG, KrylovMethod::GMRES, 1e-8 );
    tSolver.solve( tM, tX, tB );

    if( gComm.rank() == 0 )
    {
        tX.print("X");
    }

    //petsc::solve( tM, tX, tB, SymmetryMode::PositiveDefiniteSymmetric );

    tSolver.free() ;

    return gComm.finalize();
}