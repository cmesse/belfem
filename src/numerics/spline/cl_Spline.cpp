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
#include "fn_trans.hpp"

namespace belfem
{
    namespace spline
    {
//------------------------------------------------------------------------------


        
        void
        create_helpmatrix(
                const real & aSize,
                const real & aDeltaX,
                SpMatrix   & aA,
                const SplineBC aStartBC,
                const SplineBC aEndBC )
        {
            BELFEM_ASSERT( aSize > 3, "Need at least four datapoints for a spline" );

            // help constant
            index_t n = aSize - 1;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            // step 1: create Graph
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            Graph tGraph( aSize, nullptr );

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
            for ( index_t k = 1; k < n; ++k )
            {
                tGraph( k )->init_vertex_container( 3 );
                tGraph( k )->insert_vertex( tGraph( k - 1 ));
                tGraph( k )->insert_vertex( tGraph( k ) );
                tGraph( k )->insert_vertex( tGraph( k + 1 ));
            }

            // last row
            tGraph( n )->init_vertex_container( 2 );
            tGraph( n )->insert_vertex( tGraph( n - 1 ));
            tGraph( n )->insert_vertex( tGraph( n ) );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // step 2: Initialize Matrix
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // temporary matrix
            SpMatrix tA( tGraph, SpMatrixType::CSC );


            // free graph
            graph::clear( tGraph );

            real tC = 1.0/aDeltaX;

            // inner part of matrix
            for ( uint k = 1; k < n; ++k )
            {
                tA( k - 1, k ) = tC;

                tA( k, k ) = 4.0 * tC;

                tA( k + 1, k ) = tC;
            }

            switch ( aStartBC )
            {
                case SplineBC::NoCurvature:
                    tA( 0, 0 ) = 2.0 * tC;
                    tA( 1, 0 ) = tC;
                    break;
                case SplineBC::Parabolic:
                    tA( 0, 0 ) = tC;
                    tA( 1, 0 ) = tC;
                    break ;
                case SplineBC::Tangent:
                    tA( 0, 0 ) = 1. ;
                    tA( 0, 1 ) = 0. ;
                    tA( 1, 0 ) = tC;
                    break;
            }


            switch ( aEndBC )
            {
                case SplineBC::NoCurvature:
                    tA( n-1, n ) = tC;
                    tA( n, n ) = 2.0 * tC;
                    break;
                case SplineBC::Parabolic:
                    tA( n-1, n ) = tC;
                    tA( n, n ) = tC;
                    break ;
                case SplineBC::Tangent:
                    tA( n-1, n ) = tC;
                    tA( n, n-1 ) = 0. ;
                    tA( n, n ) = 1. ;
                    break;
            }

            // move matrix to output
            aA = std::move(tA);
        }

//------------------------------------------------------------------------------
    } /* namespace spline */

//------------------------------------------------------------------------------

    Spline::Spline(
        const string & aFile,
        const string & aLabel,
        const proc_t aMaster ) :
      mCommRank( comm_rank() )
    {
        if( mCommRank == aMaster )
        {
            HDF5 tFile( aFile, FileMode::OPEN_RDONLY );
            hid_t tGroup = tFile.select_group( aLabel );
            herr_t tStatus = this->load( tGroup );

            BELFEM_ERROR( tStatus == 0, "Error loading spline fom hdf5 file %s : %g ",
                aFile.c_str(), ( int ) tStatus );


        }
        this->synchronize( aMaster );
    }

    Spline::Spline( const hid_t aGroup,  const proc_t aMaster ):
      mCommRank( comm_rank() )
    {
        if( mCommRank == aMaster )
        {
            hid_t tGroup = aGroup ;
            herr_t tStatus = this->load( tGroup );

            BELFEM_ERROR( tStatus == 0, "Error loading spline from hdf5 file: %d", ( int ) tStatus );
        }
        this->synchronize( aMaster );
    }

    // constructor that creates empty container
    Spline::Spline( const index_t & aN, const real aXmin, const real aXmax ) :
      mCommRank( comm_rank() )
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
            const proc_t aMasterProc ):
      mCommRank( comm_rank() )
    {
        this->initialize( aX, aY, aA,
            spline::SplineBC::NoCurvature,
            spline::SplineBC::NoCurvature,
            0.0,
            0.0,
            aXref,
            aSref );

        if( aMasterProc < gNoOwner )
        {
            BELFEM_ASSERT( mCommRank == aMasterProc,
                           "this spline constructor must be called by the master proc" );

            this->synchronize( aMasterProc );
        }
    }

    // constructor with helpmatrix
    Spline::Spline( const Vector< real > & aX,
            const Vector< real > & aY,
                  SpMatrix       & aA,
            const spline::SplineBC aStartBC,
            const spline::SplineBC aEndBC,
            const         real     adYdX0,
            const         real     adYdX1,
            const         real     aXref,
            const         real     aSref,
            const         proc_t   aMasterProc ) :
            mCommRank( comm_rank() )
    {
        this->initialize( aX, aY, aA,
                   aStartBC,
                   aEndBC,
                   adYdX0,
                   adYdX1,
                   aXref,
                   aSref );

        if( aMasterProc < gNoOwner )
        {
            BELFEM_ASSERT( mCommRank == aMasterProc,
                           "this spline constructor must be called by the master proc" );

            this->synchronize( aMasterProc );
        }
    }


//------------------------------------------------------------------------------

    Spline::Spline( const proc_t aMasterProc ):
      mCommRank( comm_rank() )
    {
        BELFEM_ASSERT( mCommRank != aMasterProc,
                       "this spline constructor must not be called by the master proc" );
        this->synchronize( aMasterProc );
    }

//------------------------------------------------------------------------------

    void
    Spline::initialize(
            const Vector< real > & aX,
            const Vector< real > & aY,
            SpMatrix       & aA,
            const spline::SplineBC aStartBC,
            const spline::SplineBC aEndBC,
            const         real     adYdX0,
            const         real     adYdX1,
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
        this->create_rhs( aY, tB, aStartBC, aEndBC, adYdX0, adYdX1 );

        // derivatives for spline
        Vector<real> tDYDX( tB.length(), 0 );

        // create a solver
#ifdef BELFEM_SUPERLU
        Solver tSolver( SolverType::SUPERLU ) ;
#elif BELFEM_UMFPACK
        Solver tSolver( SolverType::UMFPACK ) ;
#else
        Solver tSolver( gDefaultSolver ) ;
#endif

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

        tStatus = this->save( tGroup );

        // close file
        tFile.close();

    }
//------------------------------------------------------------------------------

    herr_t
    Spline::save(         hid_t  & aGroup )
    {
        herr_t aStatus = 0 ;

        hdf5::save_scalar_to_file( aGroup, "xmin",    mXmin, aStatus );
        hdf5::save_scalar_to_file( aGroup, "xmax",    mXmax, aStatus );
        hdf5::save_scalar_to_file( aGroup, "step",    mDeltaX, aStatus );
        hdf5::save_scalar_to_file( aGroup, "npoints", mNumberOfPoints, aStatus );
        hdf5::save_scalar_to_file( aGroup, "extra", static_cast< uint >( mExtraMode ), aStatus );

        // save transposed so that load() can read it back directly
        Matrix< real > tData = trans( mData );
        hdf5::save_matrix_to_file( aGroup, "coeffs", tData, aStatus );

        return aStatus ;
    }

//------------------------------------------------------------------------------

    herr_t
    Spline::load( hid_t & aGroup )
    {
        herr_t aStatus = 0 ;

        hdf5::load_scalar_from_file( aGroup, "xmin",  mXmin, aStatus );
        hdf5::load_scalar_from_file( aGroup, "xmax",  mXmax, aStatus );
        hdf5::load_scalar_from_file( aGroup, "step", mDeltaX, aStatus );
        hdf5::load_scalar_from_file( aGroup, "npoints", mNumberOfPoints, aStatus );

        if ( hdf5::group_exists( aGroup, "extra" ) )
        {
            uint tMode ;
            hdf5::load_scalar_from_file( aGroup, "extra", tMode, aStatus );
            mExtraMode = static_cast< spline::ExtraMode >( tMode );
        }
        else
        {
            mExtraMode = spline::ExtraMode::None; ;
        }

        Matrix< real > tData ;
        hdf5::load_matrix_from_file( aGroup, "coeffs", tData, aStatus );
        mData = trans( tData );

        mNumberOfIntervals = mNumberOfPoints - 1 ;
        mInvDeltaX = 1.0 / mDeltaX ;

        return aStatus ;
    }

//------------------------------------------------------------------------------

    void
    Spline::update_data(
            SpMatrix             & aHelpMatrix,
            const Vector< real > & aValues,
            const spline::SplineBC aStartBC,
            const spline::SplineBC aEndBC,
            const         real     adYdX0,
            const         real     adYdX1,
            const         real     aXref,
            const         real     aSref)
    {
        if( mCommRank == 0 )
        {
            Vector<real> tB;
            this->create_rhs( aValues, tB, aStartBC, aEndBC, adYdX0, adYdX1 );

            // derivatives for spline
            Vector<real> tDYDX( tB.length(), 0.0 );

            // SuperLU, not UMFPACK: SuiteSparse is legacy here and is not
            // BSD-3 clean, so it is off by default and cannot be relied on.
            Solver tSolver( SolverType::SUPERLU ) ;

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

        BELFEM_ASSERT( ( norm( tDeltaX )/( aX( tN-1) - aX( 0 ) )  ) < 1e-9,
                       "X-Vector is not equidistant ( %f ) ",
                       ( double ) norm( tDeltaX ));

        BELFEM_ASSERT( aDeltaX > 0, "Invalid input" );

        BELFEM_ASSERT( aX.length() == aY.length(), "Length of vectors does not match" );

        BELFEM_ASSERT( aX.length() == aA.n_cols(), "Size of help matrix does not match" );
        BELFEM_ASSERT( aA.n_cols() == aA.n_rows(), "Help Matrix must be quadratic" );
        return aDeltaX;
    }

//------------------------------------------------------------------------------

    void
    Spline::create_rhs(
        const Vector <real> & aY,
        Vector <real> & aB,
        const spline::SplineBC aStartBC,
        const spline::SplineBC aEndBC,
        const         real     adYdX0,
        const         real     adYdX1 )
    {
        BELFEM_ASSERT( aY.length() == mNumberOfPoints,
                       "Expect length %lu for aY, but is %lu",
                       ( long unsigned int ) mNumberOfPoints,
                       ( long unsigned int ) aY.length() );

        // allocate memory for RHS
        aB.set_size( mNumberOfPoints );

        // mid entries
        for( uint k=1; k<mNumberOfIntervals; ++k )
        {
            aB( k ) = aY( k+1 ) - aY( k-1 );
        }
        // scale rhs
        aB *= 3.0 / ( mDeltaX * mDeltaX );

        // first entry
        switch (  aStartBC )
        {
            case spline::SplineBC::NoCurvature :
            {
                aB( 0 ) = 3.0*( aY( 1 ) - aY( 0 ) )  / ( mDeltaX * mDeltaX );
                break ;
            }
            case spline::SplineBC::Parabolic :
            {
                aB( 0 ) = 2.0*( aY( 1 ) - aY( 0 ) )  / ( mDeltaX * mDeltaX );
                break ;
            }
            case spline::SplineBC::Tangent :
            {
                aB( 0 ) = adYdX0 ;
                break;
            }
        }


        // last entry
        switch (  aEndBC )
        {
            case spline::SplineBC::NoCurvature :
            {
                aB( mNumberOfIntervals ) = 3.0*( aY( mNumberOfIntervals ) - aY( mNumberOfIntervals-1 ) )
                / ( mDeltaX * mDeltaX );
                break ;
            }
            case spline::SplineBC::Parabolic :
            {
                aB( mNumberOfIntervals ) = 2.0*( aY( mNumberOfIntervals ) - aY( mNumberOfIntervals-1 ) )
                / ( mDeltaX * mDeltaX );
                break ;
            }
            case spline::SplineBC::Tangent :
            {
                aB( mNumberOfIntervals ) = adYdX1 ;
                break;
            }
        }
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
    Spline::create_entropy( const real aXref, const real & aSref )
    {
        this->add_row_to_data();

        uint tN = mData.n_cols();

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

        // set mode before using entropy() for the offset calculation
        mExtraMode = spline::ExtraMode::Entropy;

        // shift offset so that entropy(aXref) == aSref
        real tDeltaS = aSref - this->entropy( aXref );

        for( index_t k=0; k<tN; ++k )
        {
            mData( 4, k ) += tDeltaS;
        }
    }

//------------------------------------------------------------------------------

    void
    Spline::create_integral( const real aXref, const real aYref )
    {
        this->add_row_to_data();

        uint tN = mData.n_cols();

        real tX = mXmin;

        mData( 4, 0 ) = 0.0;

        for( index_t k=1; k<tN; ++k )
        {
            tX += mDeltaX;

            // left side
            real tY0 =    ((( 0.25 * mData( 0, k-1 ) * tX
                            + mData( 1, k-1 )/3.0   ) * tX
                            + 0.5 * mData( 2, k-1 ) ) * tX
                            +       mData( 3, k-1 ) ) * tX
                            +       mData( 4, k-1 );

            real tY1 =   ((( 0.25 * mData( 0, k ) * tX
                            + mData( 1, k )/3.0   ) * tX
                            + 0.5 * mData( 2, k ) ) * tX
                            +       mData( 3, k ) ) * tX ;

            mData( 4, k ) = tY0-tY1;

            mExtraMode = spline::ExtraMode::Integral;
        }

        tX = std::isnan( aXref ) ? mXmin : aXref;
        real tDeltaY = aYref - this->integrate( tX );

        if ( std::abs( tDeltaY ) < 1e-10 ) return;
        for( index_t k=0; k<tN; ++k )
        {
            mData( 4, k )+= tDeltaY;
        }

        mExtraMode = spline::ExtraMode::Integral;
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
            broadcast( tIData, aMasterProc );
            broadcast( tRData, aMasterProc );
            broadcast( mData, aMasterProc );

            uint tMode = static_cast< uint >( mExtraMode );
            broadcast( tMode, aMasterProc );

        }
        else
        {
            comm_barrier() ;

            Vector< uint > tIData;
            broadcast( tIData, aMasterProc );

            Vector< real > tRData;
            broadcast( tRData, aMasterProc );

            broadcast( mData, aMasterProc );

            uint tMode ;
            broadcast( tMode, aMasterProc );

            mNumberOfPoints    = tIData( 0 );
            mNumberOfIntervals = tIData( 1 );
            mXmin   = tRData( 0 );
            mXmax   = tRData( 1 );
            mDeltaX = tRData( 2 );
            mInvDeltaX = 1.0 / mDeltaX ;

            mExtraMode = static_cast< spline::ExtraMode >( tMode );
        }
    }

//------------------------------------------------------------------------------

    void
    Spline::save_to_database( const string & aDatabase, const string & aLabel )
    {
        // check if file exists
        FileMode tMode = file_exists( aDatabase ) ?
                         FileMode::OPEN_RDWR : FileMode::NEW ;

        // open database
        HDF5 tFile( aDatabase, tMode );

        // create a new group
        tFile.create_group( aLabel );

        // write the dimension
        uint tDimension = 1 ;
        tFile.save_data( "dimension", tDimension );

        // write the number of points
        Vector<  uint > tNumPoints( tDimension );
        tNumPoints( 0 ) = mNumberOfPoints ;
        tFile.save_data( "npoints", tNumPoints );

        // write the offset
        Vector< double > tOffset( tDimension );
        tOffset( 0 ) = mXmin ;
        tFile.save_data( "offset", tOffset );

        // write the spline order, zero is a special indicator
        uint tOrder = 0 ;
        tFile.save_data( "order", tOrder );

        // write the stepsize
        Vector< double > tStep( tDimension );
        tStep( 0 ) = mDeltaX ;
        tFile.save_data("step", tStep );


        Matrix< real > tCoeffs = trans( mData );

        // save the data
        tFile.save_data( "coeffs", tCoeffs );

        uint tExtra = static_cast< uint >( mExtraMode );
        tFile.save_data( "extra", tExtra );

        tFile.close_active_group();
        tFile.close();
    }

//------------------------------------------------------------------------------

    void
    Spline::add_row_to_data()
    {
        BELFEM_ASSERT( mData.n_rows() == 4, "Expected 4 rows before adding integral row" );

        uint tN = mData.n_cols();

        // create temporary matrix
        Matrix< real > tTemp = std::move( mData );

        // reallocate original matrix
        mData.set_size( 5, tN, 0.0 );

        // copy back original coeffs
        for( index_t k=0; k<tN; ++k )
        {
            for( index_t i=0; i<4; ++i )
            {
                mData( i, k ) = tTemp( i, k );
            }
        }
    }

//------------------------------------------------------------------------------
}
