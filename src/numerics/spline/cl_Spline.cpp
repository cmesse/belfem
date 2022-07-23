//
// Created by Christian Messe on 2019-01-26.
//

#include "assert.hpp"
#include "cl_Communicator.hpp"
#include "commtools.hpp"

#include "cl_Spline.hpp"
#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"
#include "fn_Graph_clear.hpp"

#include "fn_sum.hpp"
#include "fn_norm.hpp"
#include "fn_sort.hpp"

#include "cl_Solver.hpp"
#include "fn_Create_Truss_Poly.hpp"
#include "cl_HDF5.hpp"

namespace belfem
{
    namespace spline
    {
//------------------------------------------------------------------------------

        void
        create_helpmatrix(
                const real & aSize,
                const real & aDeltaX,
                SpMatrix   & aA )
        {
            BELFEM_ASSERT( aSize > 3, "Need at least four datapoints for a spline" );

            // help constant
            index_t tN = aSize - 1;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            // step 1: create Graph
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


            Cell<graph::Vertex *> tGraph( aSize, nullptr );

            for ( index_t k = 0; k < aSize; ++k )
            {
                tGraph( k ) = new graph::Vertex();
                tGraph( k )->set_index( k );
                tGraph( k )->set_id( k + 1 );
            }

            // first row
            tGraph( 0 )->init_vertex_container( 2 );
            tGraph( 0 )->insert_vertex( tGraph( 0 ));
            tGraph( 0 )->insert_vertex( tGraph( 1 ));

            // intermediate rows
            for ( index_t k = 1; k < tN; ++k )
            {
                tGraph( k )->init_vertex_container( 3 );
                tGraph( k )->insert_vertex( tGraph( k - 1 ));
                tGraph( k )->insert_vertex( tGraph( k ) );
                tGraph( k )->insert_vertex( tGraph( k + 1 ));
            }

            // last row
            tGraph( tN )->init_vertex_container( 2 );
            tGraph( tN )->insert_vertex( tGraph( tN - 1 ));
            tGraph( tN )->insert_vertex( tGraph( tN ) );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // step 2: Initialize Matrix
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // temporary matrix
            SpMatrix tA( tGraph, SpMatrixType::CSC );


            // free graph
            graph::clear( tGraph );

            real tC = 1.0/aDeltaX;


            tA( 0, 0 ) = tC;
            tA( 1, 0 ) = tC;

            for ( uint k = 1; k < tN; ++k )
            {
                tA( k - 1, k ) = tC;

                tA( k, k ) = 4.0 * tC;

                tA( k + 1, k ) = tC;
            }
            tA( tN - 1, tN ) = tC;
            tA( tN, tN ) = tC;

            // copy matrix to output
            aA = tA;
        }

//------------------------------------------------------------------------------
    } /* namespace spline */

//------------------------------------------------------------------------------

    // constructor that creates empty container
    Spline::Spline( const index_t & aN, const real & aXmin, const real & aXmax )
    {
        BELFEM_ASSERT( aN > 0, "aN must be positive" );
        BELFEM_ASSERT( aXmax > aXmin, "aXmax must be bigger than aXmin" );

        mNumberOfPoints = aN ;
        mNumberOfIntervals = aN - 1;
        mXmin = aXmin;
        mXmax = aXmax;
        mDeltaX = ( aXmax - aXmin ) / ( ( real ) mNumberOfIntervals );
        mInvDeltaX = 1.0/mDeltaX;

    }

//------------------------------------------------------------------------------

      Spline::Spline(
               const Vector<real> & aX,
               const Vector<real> & aY,
                         SpMatrix & aA,
               const real aXref,
               const real aSref,
               const proc_t aMasterProc )
       {
           this->initialize( aX, aY, aA, aXref, aSref );

           if( aMasterProc < gNoOwner )
           {
               BELFEM_ASSERT( comm_rank() == aMasterProc,
                "this spline constructor must be called by the master proc" );

               this->synchronize( aMasterProc );
           }
       }
//------------------------------------------------------------------------------

    Spline::Spline( const proc_t aMasterProc )
    {
        BELFEM_ASSERT( comm_rank() != aMasterProc,
                      "this spline constructor must not be called by the master proc" );

        this->synchronize( aMasterProc );
    }

//------------------------------------------------------------------------------

       void
       Spline::initialize(
                   const Vector< real > & aX,
                   const Vector< real > & aY,
                         SpMatrix       & aA,
                   const         real     aXref ,
                   const         real     aSref )
        {
            mNumberOfPoints = aX.length();

            // number of intervalls
            mNumberOfIntervals = mNumberOfPoints - 1;

            // mimum X-Value
            mXmin = aX( 0 );

            // maximum X-Value
            mXmax = aX( mNumberOfIntervals );

            // stepsize and check input
            mDeltaX = this->check_input( aX, aY, aA );

            // invert stepsize
            mInvDeltaX = 1.0 / mDeltaX;

            // create the RHS
            Vector<real> tB;
            this->create_rhs( aY, tB );

            // derivatives for spline
            Vector<real> tDYDX( tB.length(), 0 );

            // create a solver
            Solver tSolver( SolverType::UMFPACK ) ;

            // solve system
            tSolver.solve( aA, tDYDX, tB ) ;

            // delete the solver
            tSolver.free();

            // create polynomial coefficients from derivatives
            this->create_coeffs( aX, aY, tDYDX );

            if ( aXref > 0 )
            {
                this->create_entropy( aXref, aSref );
            }

       }

//------------------------------------------------------------------------------

    void  Spline::save(
            const string            & aLabel,
            const string            & aPath,
            const enum FileMode   aMode )
    {
        // create a new file
        HDF5 tFile( aPath, aMode, comm_rank() > 1 );

        // get status from file
        herr_t & tStatus = tFile.status();

        // create a group
        hid_t tGroup = tFile.create_group( aLabel );

        this->save( tGroup, tStatus );

        // close file
        tFile.close();

    }
//------------------------------------------------------------------------------

    void
    Spline::save(         hid_t  & aGroup,
                          herr_t & aStatus )
    {
        // save size
        hdf5::save_scalar_to_file( aGroup, "min",  mXmin, aStatus );
        hdf5::save_scalar_to_file( aGroup, "max",  mXmax, aStatus );
        hdf5::save_scalar_to_file( aGroup, "step", mDeltaX, aStatus );
        hdf5::save_scalar_to_file( aGroup, "n",    mNumberOfPoints, aStatus );
        hdf5::save_matrix_to_file( aGroup, "data", mData, aStatus );
    }

//------------------------------------------------------------------------------

    void
    Spline::update_data(
            SpMatrix             & aHelpMatrix,
            const Vector< real > & aValues,
            const         real     aXref,
            const         real     aSref  )
    {
        if( gComm.rank() == 0 )
        {
            Vector<real> tB;
            this->create_rhs( aValues, tB );

            // derivatives for spline
            Vector<real> tDYDX( tB.length(), 0.0 );

            // create a solver
            Solver tSolver( SolverType::UMFPACK ) ;

            // solve system
            tSolver.solve( aHelpMatrix, tDYDX, tB ) ;

            tSolver.free();

            // create polynomial coefficients from derivatives
            this->create_coeffs( aValues, tDYDX );

            if ( aXref > 0.0 )
            {
                this->create_entropy( aXref, aSref );
            }

        }
    }

//------------------------------------------------------------------------------
// private :
//------------------------------------------------------------------------------

    real
    Spline::check_input(
            const Vector <real> & aX,
            const Vector <real> & aY,
            SpMatrix & aA )
    {




        BELFEM_ASSERT( aX.length() > 3, "Need at least four datapoints for this spline" );

        size_t tN = aX.length() - 1;

#if !defined( NDEBUG ) || defined( DEBUG )
        Vector <real> tXsort( aX );
        belfem::sort( tXsort );
        BELFEM_ASSERT(  tXsort == aX, "X-Vector is not consecutive" );
#endif

        Vector <real> tDeltaX( tN );

        for ( size_t k = 0; k < tN; ++k )
        {
            tDeltaX( k ) = aX( k + 1 ) - aX( k );
        }

        real aDeltaX = sum( tDeltaX ) / tN;

        tDeltaX -= aDeltaX;

        BELFEM_ASSERT( norm( tDeltaX ) < 1e-9,
                    "X-Vector is not equidistant ( %f ) ",
                    norm( tDeltaX ));

        BELFEM_ASSERT( aDeltaX > 0, "Invalid input" );

        BELFEM_ASSERT( aX.length() == aY.length(), "Length of vectors does not match" );

        BELFEM_ASSERT( aX.length() == aA.n_cols(), "Size of help matrix does not match" );
        BELFEM_ASSERT( aA.n_cols() == aA.n_rows(), "Help Matrix must be quadratic" );
        return aDeltaX;
    }

//------------------------------------------------------------------------------

    void
    Spline::create_rhs( const Vector <real> & aY, Vector <real> & aB )
    {
        BELFEM_ASSERT( aY.length() == mNumberOfPoints,
            "Expect length %lu for aY, but is %lu",
                      ( long unsigned int ) mNumberOfPoints,
                      ( long unsigned int ) aY.length() );

        // allocate memory for RHS
        aB.set_size( mNumberOfPoints );

        // first entry
        aB( 0 ) = 2.0*( aY( 1 ) - aY( 0 ) )/3.0;

        // mid entries
        for( uint k=1; k<mNumberOfIntervals; ++k )
        {
            aB( k ) = aY( k+1 ) - aY( k-1 );
        }

        // last entry
        aB( mNumberOfIntervals )
            = 2.0*( aY( mNumberOfIntervals ) - aY( mNumberOfIntervals-1 ) )/3.0;

        // scale rhs
        aB *= 3.0*std::pow( mDeltaX, -2 );
    }

//------------------------------------------------------------------------------

    void
    Spline::create_coeffs(const Vector< real > & aX,
                          const Vector< real > & aY,
                          const Vector< real > & aDYDX )
    {
        // allocate data
        mData.set_size( 4, mNumberOfPoints );

        // current interval
        Vector< real > tX( 2 );

        // funciton values
        Vector< real > tF( 4 );

        // coefficients
        Vector< real > tC( 4 );

        // loop over all intervals
        for( uint k=0; k<mNumberOfIntervals; ++k )
        {
            tX( 0 ) = aX( k );
            tX( 1 ) = aX( k+1 );

            tF( 0 ) = aY( k );
            tF( 1 ) = aDYDX( k );
            tF( 2 ) = aY( k+1 );
            tF( 3 ) = aDYDX( k+1 );

            create_truss_poly( tX, tF, tC );

            mData.set_col( k, tC );
        }

        // last entry
        mData.set_col( mNumberOfIntervals, mData.col( mNumberOfIntervals-1 ) );
    }

//------------------------------------------------------------------------------

    void
    Spline::create_coeffs( const Vector< real > & aY,
                          const Vector< real > & aDYDX )
    {
        uint tN = aY.length() - 1;

        // allocate data
        mData.set_size( 4, tN+1 );

        // current interval
        Vector< real > tX( 2 );

        // funciton values
        Vector< real > tF( 4 );

        // coefficients
        Vector< real > tC( 4 );

        tX( 1 ) = mXmin;
        // loop over all intervals
        for( uint k=0; k<tN; ++k )
        {
            // shift X
            tX( 0 ) = tX( 1 );
            tX( 1 ) += mDeltaX;

            tF( 0 ) = aY( k );
            tF( 1 ) = aDYDX( k );
            tF( 2 ) = aY( k+1 );
            tF( 3 ) = aDYDX( k+1 );

            create_truss_poly( tX, tF, tC );

            mData.set_col( k, tC );
        }

        // last entry
        mData.set_col( tN, mData.col( tN-1 ) );
    }

//------------------------------------------------------------------------------

    void
    Spline::create_entropy( const real & aXref, const real & aSref )
    {
        // create temporary matrix
        Matrix< real > tTemp( mData );

        uint tN = mData.n_cols();

        // reallocate original matrix
        mData.set_size( 5, tN );

        // copy back original coeffs
        for( index_t k=0; k<tN; ++k )
        {
            for( index_t i=0; i<4; ++i )
            {
                mData( i, k ) = tTemp( i, k );
            }
        }

        real tX = mXmin;

        // first entry
        mData( 4, 0 ) = 0;

        // all other entries
        for( index_t k=1; k<tN; ++k )
        {
            tX += mDeltaX;

            // left side
            real tY0 =    ( 1.5 * mData( 0, k-1 )   * tX
                         +  2.0 * mData( 1, k-1 ) ) * tX
                         +        mData( 2, k-1 )   * std::log( tX )
                         +        mData( 4, k-1 );

            real tY1 =    ( 1.5 * mData( 0, k )   * tX
                         +  2.0 * mData( 1, k ) ) * tX
                         +        mData( 2, k )   * std::log( tX );


            mData( 4, k ) = tY0-tY1;
        }

        // shift offset
        real tDeltaS = aSref - this->entropy( aXref );

        for( index_t k=0; k<tN; ++k )
        {
            mData( 4, k ) += tDeltaS;
        }

    }
//------------------------------------------------------------------------------


    void
    Spline::synchronize( const proc_t aMasterProc )
    {
        if( gComm.rank() == aMasterProc )
        {
            // create commtable
            proc_t tCommSize = gComm.size();
            Vector< proc_t > tCommTable( tCommSize );
            for( proc_t p=0; p<tCommSize; ++p )
            {
                tCommTable( p ) = p ;
            }

            // prepare data that is to be sent
            Vector< uint > tIData( 2 );
            tIData( 0 ) = mNumberOfPoints ;
            tIData( 1 ) = mNumberOfIntervals ;

            Vector< real > tRData( 3 );
            tRData( 0 ) = mXmin ;
            tRData( 1 ) = mXmax ;
            tRData( 2 ) = mDeltaX ;

            comm_barrier() ;

            // send data to others
            send_same( tCommTable, tIData );
            send_same( tCommTable, tRData );
            send_same( tCommTable, mData );
        }
        else
        {
            comm_barrier() ;

            Vector< uint > tIData;
            receive( aMasterProc, tIData );

            Vector< real > tRData;
            receive( aMasterProc, tRData );

            receive( aMasterProc, mData );

            mNumberOfPoints    = tIData( 0 );
            mNumberOfIntervals = tIData( 1 );
            mXmin   = tRData( 0 );
            mXmax   = tRData( 1 );
            mDeltaX = tRData( 2 );
            mInvDeltaX = 1.0 / mDeltaX ;
        }
    }

//------------------------------------------------------------------------------
}
