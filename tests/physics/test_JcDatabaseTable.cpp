/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

/*
 * The tabulated critical-current functions and the tables they accept.
 *
 * The shipped tables in share/material store log10 jc and log10 n over
 * T, log10( B / T ) from -2 and the unfolded angle theta = 0 .. pi. Below
 * Bmin = 10 mT a table has no data. There JcFunctionDatabase continues n from
 * the table alone, and continues jc towards the measured self-field level
 * stored with a REBCO table; bscco jc has no such curve and keeps the clamp.
 * These cases pin the contracts: the value does not rise with |B| anywhere,
 * jc does not rise with T, the continuation joins the table C1, its level at
 * B = 0 carries no angle where the measured curve dominates, the exact
 * derivatives are the derivatives of the value, n stays a power law, and a
 * table that breaks the reader's assumptions is refused at load.
 * n is monotone in B only; nothing here asserts a temperature trend of n.
 * The tables are HDF5 files, so the file compiles to nothing without HDF5.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <memory>
#include <string>
#include <unistd.h>

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"
#include "cl_Mesh.hpp"
#include "cl_HDF5.hpp"
#include "cl_Database.hpp"
#include "cl_JcFunction_Database.hpp"
#include "cl_MaterialFactory.hpp"

#ifdef BELFEM_HDF5

using belfem::real ;
using belfem::index_t ;
using belfem::material::JcFunction ;

namespace
{
    const real gPi = belfem::constant::pi ;

    // Bmin of every shipped table: 10^-2 T
    const real gBmin = 1e-2 ;

    struct ShippedTable
    {
        const char * mFile ;
        real mTmin ;
        real mTmax ;
    };

    // the first five carry a measured self-field curve; bscco does not
    const ShippedTable gShipped[ 6 ] = {
        { "sp-ap.hdf5",      3.0,  91.0 },
        { "superox.hdf5",    3.0,  91.0 },
        { "fesc.hdf5",       3.0,  91.0 },
        { "fysc.hdf5",       3.0,  91.0 },
        { "sst-1.hdf5",      3.0,  91.0 },
        { "bscco-2223.hdf5", 4.0, 108.0 } };

    const char * gLabels[ 2 ] = { "jc", "n" };

    // the largest rise from one sample to the next on the grids below:
    // 2.1e-14 in B and 1.2e-14 in T, both roundoff
    const real gMaxRelativeRise = 1e-13 ;

    std::unique_ptr< JcFunction >
    shipped( const std::string & aFile, const std::string & aLabel )
    {
        belfem::MaterialFactory tFactory ;
        return std::unique_ptr< JcFunction >( tFactory.create_jc_function( aFile, aLabel ) );
    }

    //! a table ( T, log10 B, theta ) written to a scratch HDF5 file in the
    //! layout JcFunctionDatabase reads; the value at every node is
    //! aValue( T, log10 B, theta ). aExtra may add groups beside the table.
    //! Removed on scope exit
    class ScopedTable
    {
        const std::string mPath ;

    public:

        ScopedTable( const std::string & aName,
                     const std::string & aLabel,
                     const belfem::uint aOrder,
                     const belfem::Vector< index_t > & aNumNodes,
                     const belfem::Vector< real > & aStep,
                     const belfem::Vector< real > & aOrigin,
                     const std::function< real( real, real, real ) > & aValue,
                     const std::function< void( belfem::HDF5 & ) > & aExtra = nullptr ) :
            mPath( ::testing::TempDir() + "belfem_test_jctable_" + aName + "_"
                   + std::to_string( ( long ) getpid() ) + ".hdf5" )
        {
            belfem::Mesh tMesh( aOrder, aNumNodes, aStep, aOrigin );

            belfem::Vector< real > & tField = tMesh.create_field( aLabel );

            const bool t3D = aNumNodes.length() == 3 ;
            for ( index_t i = 0; i < aNumNodes( 0 ); ++i )
            {
                for ( index_t j = 0; j < aNumNodes( 1 ); ++j )
                {
                    for ( index_t k = 0; k < ( t3D ? aNumNodes( 2 ) : 1 ); ++k )
                    {
                        belfem::mesh::Node * tNode = t3D ? tMesh.node( i, j, k ) : tMesh.node( i, j );
                        tField( tNode->index() ) = aValue( tNode->x(), tNode->y(), t3D ? tNode->z() : 0.0 );
                    }
                }
            }

            belfem::Database tDatabase( &tMesh, aLabel, false );

            belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
            tFile.create_group( aLabel );
            tDatabase.save( tFile.active_group() );
            tFile.close_active_group();
            if ( aExtra )
            {
                aExtra( tFile );
            }
            tFile.close();
        }

        ~ScopedTable()
        {
            std::remove( mPath.c_str() );
        }

        const std::string &
        path() const
        {
            return mPath ;
        }
    };

    //! a constant order-1 table over T = 4 .. 24 K and log10 B = -2 .. -1 whose
    //! angular window starts at aAngleOrigin and spans aAngleWidth
    ScopedTable
    window_table( const std::string & aName, const real aAngleOrigin, const real aAngleWidth )
    {
        return ScopedTable( aName, "jc", 1,
            belfem::Vector< index_t >( { 3, 3, 5 } ),
            belfem::Vector< real >( { 10.0, 0.5, aAngleWidth / 4.0 } ),
            belfem::Vector< real >( { 4.0, -2.0, aAngleOrigin } ),
            []( real, real, real ) { return 9.0 ; } );
    }

    //! one HEX27 element in T ( 3 .. 5 K, xi = T - 4 ) and log10 B
    //! ( -2 .. -1 ), over the full window in theta; the value is
    //! 9 + ( log10 B + 2 ) aSlope( xi )
    ScopedTable
    face_table( const std::string & aName, const std::function< real( real ) > & aSlope )
    {
        return ScopedTable( aName, "jc", 2,
            belfem::Vector< index_t >( { 3, 3, 3 } ),
            belfem::Vector< real >( { 1.0, 0.5, 0.5 * gPi } ),
            belfem::Vector< real >( { 3.0, -2.0, 0.0 } ),
            [ aSlope ]( real T, real u, real ) { return 9.0 + ( u + 2.0 ) * aSlope( T - 4.0 ); } );
    }
}

// =============================================================================
// tables the reader refuses
// =============================================================================

TEST( JcDatabaseTable, ShippedTablesAreAccepted )
{
    for ( const ShippedTable & tTable : gShipped )
    {
        for ( const char * tLabel : gLabels )
        {
            EXPECT_NO_THROW( shipped( tTable.mFile, tLabel ) ) << tTable.mFile << " " << tLabel ;
        }
    }
}

TEST( JcDatabaseTable, FoldedAngularWindowIsRefused )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture";
    }

    // width pi, but [ -pi/2, pi/2 ] cannot serve an unfolded theta > pi/2.
    // Mutation: turn the coverage guard into a width test and this goes red.
    const ScopedTable tTable = window_table( "folded", -0.5 * gPi, gPi );

    EXPECT_THROW( belfem::material::JcFunctionDatabase( tTable.path(), "jc" ),
                  std::runtime_error );
}

TEST( JcDatabaseTable, HalfWindowIsRefused )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture";
    }

    // starts at 0 but stops at pi/2: a second refusal, from the upper bound
    const ScopedTable tTable = window_table( "half", 0.0, 0.5 * gPi );

    EXPECT_THROW( belfem::material::JcFunctionDatabase( tTable.path(), "jc" ),
                  std::runtime_error );
}

TEST( JcDatabaseTable, FullWindowIsAccepted )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture";
    }

    // the control for the two refusals: the same fixture over [ 0, pi ]
    const ScopedTable tTable = window_table( "full", 0.0, gPi );

    EXPECT_NO_THROW( belfem::material::JcFunctionDatabase( tTable.path(), "jc" ) );
}

TEST( JcDatabaseTable, TwoAxisTableIsRefused )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture";
    }

    const ScopedTable tTable( "twoaxis", "jc", 1,
        belfem::Vector< index_t >( { 3, 3 } ),
        belfem::Vector< real >( { 10.0, 0.5 } ),
        belfem::Vector< real >( { 4.0, -2.0 } ),
        []( real, real, real ) { return 9.0 ; } );

    // the reader's own refusal, not an accident further down: in a release
    // build nothing else stops a 2-axis table before it is read out of range
    try
    {
        belfem::material::JcFunctionDatabase tJc( tTable.path(), "jc" );
        ADD_FAILURE() << "a 2-axis table was accepted" ;
    }
    catch ( const std::runtime_error & aError )
    {
        EXPECT_NE( std::string( aError.what() ).find( "needs T, log10 B and theta" ), std::string::npos )
            << aError.what() ;
    }
}

TEST( JcDatabaseTable, TableRisingInsideAnElementIsRefused )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture";
    }

    // d log10 jc / d log10 B = 0.1 - ( xi - 0.5 )^2 on the lowest-field face:
    // negative at all three temperature nodes, +0.1 at xi = 0.5. No
    // continuation below Bmin can be monotone there, and a nodal check
    // passes it.
    // Mutation: check the nodal slopes instead of the Bernstein bound and
    // this goes red
    const ScopedTable tTable = face_table( "rising",
        []( const real xi ) { return 0.1 - ( xi - 0.5 ) * ( xi - 0.5 ); } );

    EXPECT_THROW( belfem::material::JcFunctionDatabase( tTable.path(), "jc" ),
                  std::runtime_error );
}

TEST( JcDatabaseTable, TableFallingOnItsFaceIsAccepted )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture";
    }

    // the control: a slope of -0.5 + 0.2 xi, negative on the whole face
    const ScopedTable tTable = face_table( "falling",
        []( const real xi ) { return -0.5 + 0.2 * xi ; } );

    EXPECT_NO_THROW( belfem::material::JcFunctionDatabase( tTable.path(), "jc" ) );
}

TEST( JcDatabaseTable, NBelowTheFloorIsRefused )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture";
    }

    // n = 1.005 at Bmin is no longer a usable power law; 1.02 is
    const auto tTable = [ ]( const std::string & aName, const real aN )
    {
        return ScopedTable( aName, "n", 1,
            belfem::Vector< index_t >( { 3, 3, 5 } ),
            belfem::Vector< real >( { 10.0, 0.5, 0.25 * gPi } ),
            belfem::Vector< real >( { 4.0, -2.0, 0.0 } ),
            [ aN ]( real, real, real ) { return std::log10( aN ); } );
    };

    const ScopedTable tLow  = tTable( "nlow", 1.005 );
    const ScopedTable tHigh = tTable( "nhigh", 1.02 );

    EXPECT_THROW( belfem::material::JcFunctionDatabase( tLow.path(), "n" ), std::runtime_error );
    EXPECT_NO_THROW( belfem::material::JcFunctionDatabase( tHigh.path(), "n" ) );
}

TEST( JcDatabaseTable, NonFiniteTableIsRefused )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture";
    }

    // one NaN node on the lowest-field face: no bound can be certified.
    // Mutation: let the bound reductions skip a NaN coefficient and this
    // goes red
    const auto tTable = [ ]( const std::string & aName, const std::string & aLabel )
    {
        return ScopedTable( aName, aLabel, 1,
            belfem::Vector< index_t >( { 3, 3, 5 } ),
            belfem::Vector< real >( { 10.0, 0.5, 0.25 * gPi } ),
            belfem::Vector< real >( { 4.0, -2.0, 0.0 } ),
            []( real T, real u, real ) { return T == 14.0 && u == -2.0 ? BELFEM_QUIET_NAN : 1.0 ; } );
    };

    const ScopedTable tJc = tTable( "nanjc", "jc" );
    const ScopedTable tN  = tTable( "nann", "n" );

    EXPECT_THROW( belfem::material::JcFunctionDatabase( tJc.path(), "jc" ), std::runtime_error );
    EXPECT_THROW( belfem::material::JcFunctionDatabase( tN.path(), "n" ), std::runtime_error );
}

TEST( JcDatabaseTable, AngleSnappedOntoTheWindowHasNoTangent )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture";
    }

    // a window starting 5e-7 rad above zero passes the coverage guard; an
    // angle below it is held on the window end, so the value does not move
    // with the angle there and neither may the tangent. Both the table and
    // the continuation branch.
    // Mutation: drop the snap test from deval_dbeta and this goes red
    const real tStart = 5e-7 ;
    const ScopedTable tTable( "snap", "n", 1,
        belfem::Vector< index_t >( { 3, 3, 5 } ),
        belfem::Vector< real >( { 10.0, 0.5, 0.25 * gPi } ),
        belfem::Vector< real >( { 4.0, -2.0, tStart } ),
        []( real, real u, real theta ) { return 1.0 + 0.1 * theta - 0.1 * ( u + 2.0 ); } );

    const belfem::material::JcFunctionDatabase tN( tTable.path(), "n" );

    for ( const real tB : { 0.5 * gBmin, 0.1 } )
    {
        EXPECT_EQ( tN.eval( tB, 1e-7, 10.0 ), tN.eval( tB, 4e-7, 10.0 ) ) << "B = " << tB << " T" ;
        EXPECT_EQ( tN.deval_dbeta( tB, 2.5e-7, 10.0 ), 0.0 ) << "B = " << tB << " T" ;
        EXPECT_GT( tN.deval_dbeta( tB, 0.3, 10.0 ), 0.0 ) << "B = " << tB << " T" ;
    }
}

// =============================================================================
// the continuation on synthetic tables, against closed forms
// =============================================================================

TEST( JcDatabaseTable, NContinuesAsTheParabolaOfItsEdge )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture";
    }

    // log10 n = log10 20 - 0.4 ( log10 B + 2 ): n = 20 at Bmin, falling with
    // dn/dB = -0.4 * 20 / Bmin there. The parabola with that value and slope
    // and zero slope at B = 0 gives n( 0 ) = 20 ( 1 + 0.4 / 2 ) = 24 and
    // n( Bmin / 2 ) = 20 + 8 * 0.75 / 2 = 23
    const ScopedTable tTable( "nparabola", "n", 1,
        belfem::Vector< index_t >( { 3, 3, 5 } ),
        belfem::Vector< real >( { 10.0, 0.5, 0.25 * gPi } ),
        belfem::Vector< real >( { 4.0, -2.0, 0.0 } ),
        []( real, real u, real ) { return std::log10( 20.0 ) - 0.4 * ( u + 2.0 ); } );

    const belfem::material::JcFunctionDatabase tN( tTable.path(), "n" );

    for ( const real tTheta : { 0.0, 0.3, 2.0 } )
    {
        EXPECT_NEAR( tN.eval( 0.0, tTheta, 10.0 ), 24.0, 1e-12 );
        EXPECT_NEAR( tN.eval( 0.5 * gBmin, tTheta, 10.0 ), 23.0, 1e-12 );
        EXPECT_NEAR( tN.deval_dB( 0.0, tTheta, 10.0 ), 0.0, 1e-9 );
        EXPECT_NEAR( tN.deval_dB( 0.5 * gBmin, tTheta, 10.0 ), -0.4 * 20.0 / gBmin * 0.5, 1e-9 );
    }
}

TEST( JcDatabaseTable, JcWithoutAMeasuredCurveKeepsTheClamp )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture";
    }

    // no /source rows: below Bmin the value is the edge value and the field
    // slope is zero, as for bscco
    const ScopedTable tTable = face_table( "clamp",
        []( const real xi ) { return -0.5 + 0.2 * xi ; } );

    const belfem::material::JcFunctionDatabase tJc( tTable.path(), "jc" );

    EXPECT_DOUBLE_EQ( tJc.eval( 0.3 * gBmin, 1.0, 4.2 ), tJc.eval( gBmin, 1.0, 4.2 ) );
    EXPECT_EQ( tJc.deval_dB( 0.3 * gBmin, 1.0, 4.2 ), 0.0 );
}

namespace
{
    //! a jc table with jc = 10 at Bmin falling as B^-0.1, T = 4 .. 74 K, and
    //! a measured self-field curve 100, 90, 95, 80, 60 at 10, 20, 30, 45,
    //! 50 K, calibrated to 60 at 77.5 K ( the curve's end value, so the
    //! scale is 1 )
    ScopedTable
    measured_table( const std::string & aName )
    {
        return ScopedTable( aName, "jc", 1,
            belfem::Vector< index_t >( { 8, 3, 5 } ),
            belfem::Vector< real >( { 10.0, 0.5, 0.25 * gPi } ),
            belfem::Vector< real >( { 4.0, -2.0, 0.0 } ),
            []( real, real u, real ) { return 1.0 - 0.1 * ( u + 2.0 ); },
            []( belfem::HDF5 & aFile )
            {
                const real tT[ 5 ] = { 10.0, 20.0, 30.0, 45.0, 50.0 };
                const real tIcw[ 5 ] = { 100.0, 90.0, 95.0, 80.0, 60.0 };

                // two angles per temperature, averaged by the reader, and a
                // row at finite field that it must ignore
                belfem::Matrix< real > tPoints( 11, 4 );
                for ( index_t k = 0; k < 5; ++k )
                {
                    for ( index_t a = 0; a < 2; ++a )
                    {
                        tPoints( 2 * k + a, 0 ) = tT[ k ] ;
                        tPoints( 2 * k + a, 1 ) = 0.0 ;
                        tPoints( 2 * k + a, 2 ) = a * 0.5 * gPi ;
                        tPoints( 2 * k + a, 3 ) = tIcw[ k ] * ( a == 0 ? 0.999 : 1.001 );
                    }
                }
                tPoints( 10, 0 ) = 20.0 ;
                tPoints( 10, 1 ) = 0.5 ;
                tPoints( 10, 2 ) = 0.0 ;
                tPoints( 10, 3 ) = 1000.0 ;

                aFile.create_group( "source" );
                aFile.save_data( "columns",
                    belfem::Cell< std::string >( { "temperature", "field", "angle", "icw" } ) );
                aFile.save_data( "points", tPoints );
                aFile.close_active_group();

                aFile.create_group( "meta" );
                aFile.save_data( "Icw_77p5K_sf_A_per_m", 60.0 );
                aFile.save_data( "t_eff_m", 1.0 );
                aFile.close_active_group();
            } );
    }
}

TEST( JcDatabaseTable, MeasuredCurveIsProjectedNonIncreasing )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture";
    }

    // the raw curve rises from 20 to 30 K; the anchor uses its running
    // maximum from high to low T, 100, 95, 95, 80, 60. The knot slopes are
    // -0.5 and -4 at the ends ( one-sided secants ), 0 at 20 and 30 K
    // ( a zero secant ), and at 45 K, between secants -1 over 15 K and -4
    // over 5 K, the weighted harmonic mean 60 / ( 25 / -1 + 35 / -4 )
    // = -16/9. The Hermite on [ 45, 50 ] K then gives 70 + 5 ( -16/9 ) / 8
    // + 20 / 8 = 71.3889 at 47.5 K, and outside the knots the curve continues
    // linearly. The measured level dominates the requirement by a factor of
    // four or more, so the smooth maximum adds less than 1e-5.
    // Mutation: drop the projection, or swap the two weights, and this goes
    // red
    const ScopedTable tTable = measured_table( "projected" );
    const belfem::material::JcFunctionDatabase tJc( tTable.path(), "jc" );

    const real tExpected[ 9 ][ 2 ] = {
        {  4.0, 103.0 }, { 10.0, 100.0 }, { 20.0, 95.0 }, { 25.0, 95.0 },
        { 30.0,  95.0 }, { 45.0,  80.0 }, { 47.5, 70.0 + 2.5 - 10.0 / 9.0 },
        { 50.0,  60.0 }, { 54.0,  44.0 } };

    for ( const auto & tPoint : tExpected )
    {
        EXPECT_NEAR( tJc.eval( 0.0, 0.7, tPoint[ 0 ] ), tPoint[ 1 ], 1e-4 ) << "T = " << tPoint[ 0 ] << " K" ;
    }

    real tPrevious = tJc.eval( 0.0, 0.7, 4.0 );
    for ( real T = 4.1; T <= 54.0; T += 0.1 )
    {
        const real tValue = tJc.eval( 0.0, 0.7, T );
        EXPECT_LE( tValue, tPrevious * ( 1.0 + gMaxRelativeRise ) ) << "T = " << T << " K" ;
        tPrevious = tValue ;
    }
}

TEST( JcDatabaseTable, AnchorFallsThroughTheZeroOfTheContinuedCurve )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture";
    }

    // beyond 50 K the measured curve continues as 60 - 4 ( T - 50 ), which
    // crosses zero at 65 K; the anchor is then the requirement, and must
    // keep falling smoothly through the crossing.
    // Mutation: scale the smooth maximum's width with |a| + |b| and this
    // goes red
    const ScopedTable tTable = measured_table( "zero" );
    const belfem::material::JcFunctionDatabase tJc( tTable.path(), "jc" );

    real tPrevious = tJc.eval( 0.0, 0.7, 60.0 );
    for ( real T = 60.01; T <= 70.0; T += 0.01 )
    {
        const real tValue = tJc.eval( 0.0, 0.7, T );
        EXPECT_LE( tValue, tPrevious * ( 1.0 + gMaxRelativeRise ) ) << "T = " << T << " K" ;
        EXPECT_LE( tJc.deval_dT( 0.0, 0.7, T ), 0.0 ) << "T = " << T << " K" ;
        tPrevious = tValue ;
    }
}

TEST( JcDatabaseTable, JcContinuesAsTheQuadraticToItsAnchor )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture";
    }

    // at 10 K the anchor is J = 100, the edge V = 10 with S = -1 / Bmin.
    // Delta = ( V - J ) / Bmin = -90 / Bmin and beta = S / Delta = 1/90, so
    // the origin slope is s0 = 2 Delta - S = -179 / Bmin and at B = Bmin / 2
    // the Hermite gives 0.5 J + 0.125 s0 Bmin + 0.5 V - 0.125 S Bmin = 32.75.
    // It is the quadratic, so its slope there is the mean of the end slopes,
    // ( -179 - 1 ) / 2 / Bmin = -90 / Bmin
    const ScopedTable tTable = measured_table( "quadratic" );
    const belfem::material::JcFunctionDatabase tJc( tTable.path(), "jc" );

    EXPECT_NEAR( tJc.eval( 0.5 * gBmin, 0.7, 10.0 ), 32.75, 1e-4 );
    EXPECT_NEAR( tJc.deval_dB( 0.5 * gBmin, 0.7, 10.0 ) * gBmin, -90.0, 1e-3 );
    EXPECT_NEAR( tJc.deval_dB( 0.0, 0.7, 10.0 ) * gBmin, -179.0, 1e-3 );
}

// =============================================================================
// the contracts on the shipped tables
// =============================================================================

TEST( JcDatabaseTable, NothingRisesWithField )
{
    // B = 0, then a log grid from 1 uT through Bmin to 10^1.9 T, for T every
    // 2 K and theta every 5 degrees.
    // Mutation: read the field axis backwards and this goes red
    for ( const ShippedTable & tTable : gShipped )
    {
        for ( const char * tLabel : gLabels )
        {
            const std::unique_ptr< JcFunction > tF = shipped( tTable.mFile, tLabel );

            for ( real T = tTable.mTmin + 0.5; T < tTable.mTmax; T += 2.0 )
            {
                for ( int a = 0; a <= 36; ++a )
                {
                    const real tTheta = gPi * ( real ) a / 36.0 ;
                    real tPrevious = tF->eval( 0.0, tTheta, T );

                    for ( int k = 0; k <= 158; ++k )
                    {
                        const real tB = 1e-6 * std::pow( 10.0, ( real ) k / 20.0 );
                        const real tValue = tF->eval( tB, tTheta, T );

                        EXPECT_LE( ( tValue - tPrevious ) / tPrevious, gMaxRelativeRise )
                                << tTable.mFile << " " << tLabel << ": T = " << T
                                << " K, theta = " << tTheta << ", B = " << tB << " T" ;

                        tPrevious = tValue ;
                    }
                }
            }
        }
    }
}

TEST( JcDatabaseTable, JcDoesNotRiseWithTemperature )
{
    // from below the table's Tmin to beyond its Tmax, at fields on both sides
    // of Bmin. Mutation: read the measured curve's temperature axis
    // backwards and this goes red
    const real tFields[ 7 ] = { 0.0, 0.1 * gBmin, 0.5 * gBmin, 0.99 * gBmin, gBmin, 0.1, 1.0 };

    for ( const ShippedTable & tTable : gShipped )
    {
        const std::unique_ptr< JcFunction > tJc = shipped( tTable.mFile, "jc" );

        for ( const real tB : tFields )
        {
            for ( int a = 0; a <= 36; ++a )
            {
                const real tTheta = gPi * ( real ) a / 36.0 ;
                real tPrevious = tJc->eval( tB, tTheta, 1.9 );

                for ( real T = 2.0; T <= tTable.mTmax + 1.0; T += 0.25 )
                {
                    const real tValue = tJc->eval( tB, tTheta, T );

                    EXPECT_LE( ( tValue - tPrevious ) / tPrevious, gMaxRelativeRise )
                            << tTable.mFile << ": T = " << T << " K, theta = " << tTheta
                            << ", B = " << tB << " T" ;

                    tPrevious = tValue ;
                }
            }
        }
    }
}

TEST( JcDatabaseTable, ContinuationJoinsTheTableC1 )
{
    // value and exact field slope just below Bmin against the table's at
    // Bmin, relative to the natural scale value / Bmin. Secant tests cannot
    // resolve this where the table is flat in B: the slope there is smaller
    // than the rounding noise of a difference quotient.
    // Mutation: drop the edge slope from the continuation's field derivative
    // and this goes red
    const real tBelow = gBmin * ( 1.0 - 1e-9 );

    for ( int s = 0; s < 5; ++s )
    {
        for ( const char * tLabel : gLabels )
        {
            const std::unique_ptr< JcFunction > tF = shipped( gShipped[ s ].mFile, tLabel );

            for ( real T = 3.5; T < 91.0; T += 1.0 )
            {
                for ( int a = 0; a <= 36; ++a )
                {
                    const real tTheta = gPi * ( real ) a / 36.0 ;
                    const real tValue = tF->eval( gBmin, tTheta, T );
                    const real tScale = tValue / gBmin ;

                    EXPECT_NEAR( tF->eval( tBelow, tTheta, T ), tValue, 1e-8 * tValue )
                            << gShipped[ s ].mFile << " " << tLabel << ": T = " << T << " K, theta = " << tTheta ;
                    EXPECT_NEAR( tF->deval_dB( tBelow, tTheta, T ), tF->deval_dB( gBmin, tTheta, T ), 1e-6 * tScale )
                            << gShipped[ s ].mFile << " " << tLabel << ": T = " << T << " K, theta = " << tTheta ;
                }
            }
        }
    }
}

TEST( JcDatabaseTable, SelfFieldLevelIsAngleFreeWhereTheCurveDominates )
{
    // at B = 0 the continuation returns the measured, angle-free level
    // wherever that level clears the table's requirement at every angle by
    // a margin well above the smooth maximum's width. On sp-ap that holds at
    // 83, 84, 88 and 89 K, by at least 4 percent; the smooth maximum then
    // exceeds the measured level by at most ( 1e-3 )^2 / ( 4 * 0.04 ) = 6e-6
    // of it. Elsewhere the level is lifted, and carries the angle.
    // Mutation: anchor on the requirement alone and this goes red
    const std::unique_ptr< JcFunction > tJc = shipped( "sp-ap.hdf5", "jc" );

    for ( const real T : { 83.0, 84.0, 88.0, 89.0 } )
    {
        const real tReference = tJc->eval( 0.0, 0.0, T );

        for ( int a = 1; a <= 36; ++a )
        {
            const real tTheta = gPi * ( real ) a / 36.0 ;
            EXPECT_NEAR( tJc->eval( 0.0, tTheta, T ), tReference, 1e-5 * tReference )
                    << "T = " << T << " K, theta = " << tTheta ;
        }
    }
}

TEST( JcDatabaseTable, TangentsAreTheDerivativesOfTheValue )
{
    // the three exact derivatives against central differences, below Bmin
    // and inside the table, at points off every element face and away from
    // Tmin and Tmax. The temperatures cross many points where the jc anchor
    // switches between the measured level and the requirement, so the smooth
    // maximum's own derivative is exercised. Measured agreement is 2e-8 of
    // the scale; the band is 1e-6.
    // Mutation: drop the mixed-derivative term of the edge slope's T tangent
    // and this goes red
    const real tStep[ 3 ] = { 1e-6 * gBmin, 1e-5, 1e-6 };

    for ( const ShippedTable & tTable : gShipped )
    {
        for ( const char * tLabel : gLabels )
        {
            const std::unique_ptr< JcFunction > tF = shipped( tTable.mFile, tLabel );

            for ( real T = tTable.mTmin + 2.7; T < tTable.mTmax - 1.0; T += 6.0 )
            {
                for ( real tDegrees = 3.3; tDegrees < 180.0; tDegrees += 14.0 )
                {
                    for ( const real x : { 0.13, 0.47, 0.81, 2.3, 37.7 } )
                    {
                        const real tTheta = tDegrees * gPi / 180.0 ;
                        const real tB = x * gBmin ;
                        const real tValue = tF->eval( tB, tTheta, T );

                        const real tDB = ( tF->eval( tB + tStep[ 0 ], tTheta, T )
                                         - tF->eval( tB - tStep[ 0 ], tTheta, T ) ) / ( 2.0 * tStep[ 0 ] );
                        const real tDT = ( tF->eval( tB, tTheta, T + tStep[ 1 ] )
                                         - tF->eval( tB, tTheta, T - tStep[ 1 ] ) ) / ( 2.0 * tStep[ 1 ] );
                        const real tDA = ( tF->eval( tB, tTheta + tStep[ 2 ], T )
                                         - tF->eval( tB, tTheta - tStep[ 2 ], T ) ) / ( 2.0 * tStep[ 2 ] );

                        const std::string tWhere = std::string( tTable.mFile ) + " " + tLabel
                                + ": T = " + std::to_string( T ) + " K, theta = " + std::to_string( tDegrees )
                                + " deg, B / Bmin = " + std::to_string( x );

                        EXPECT_NEAR( tF->deval_dB( tB, tTheta, T ), tDB, 1e-6 * tValue / std::min( tB, gBmin ) ) << tWhere ;
                        EXPECT_NEAR( tF->deval_dT( tB, tTheta, T ), tDT, 1e-6 * tValue ) << tWhere ;
                        EXPECT_NEAR( tF->deval_dbeta( tB, tTheta, T ), tDA, 1e-6 * tValue ) << tWhere ;
                    }
                }
            }
        }
    }
}

TEST( JcDatabaseTable, NStaysAPowerLawAtZeroField )
{
    // n( 0 ) >= n( Bmin ) by construction; the smallest on the shipped tables
    // is 1.29 ( superox near Tc ), well above 1.01
    for ( const ShippedTable & tTable : gShipped )
    {
        const std::unique_ptr< JcFunction > tN = shipped( tTable.mFile, "n" );

        for ( real T = tTable.mTmin - 1.0; T <= tTable.mTmax + 1.0; T += 0.5 )
        {
            for ( int a = 0; a <= 36; ++a )
            {
                const real tTheta = gPi * ( real ) a / 36.0 ;
                EXPECT_GT( tN->eval( 0.0, tTheta, T ), 1.01 )
                        << tTable.mFile << ": T = " << T << " K, theta = " << tTheta ;
            }
        }
    }
}

TEST( JcDatabaseTable, BelowTminTheValueIsHeld )
{
    // the temperature window is clamped: below Tmin the value is the value at
    // Tmin and the T tangent is zero, on both sides of Bmin
    const struct { const char * mFile ; const char * mLabel ; real mTmin ; } tCases[ 2 ] = {
        { "sp-ap.hdf5", "jc", 3.0 },
        { "bscco-2223.hdf5", "n", 4.0 } };

    for ( const auto & tCase : tCases )
    {
        const std::unique_ptr< JcFunction > tF = shipped( tCase.mFile, tCase.mLabel );

        for ( const real tB : { 0.0, 0.5 * gBmin, 0.1 } )
        {
            for ( const real tTheta : { 0.0, 1.0, 2.5 } )
            {
                EXPECT_EQ( tF->eval( tB, tTheta, 1.9 ), tF->eval( tB, tTheta, tCase.mTmin ) )
                        << tCase.mFile << " " << tCase.mLabel << ": B = " << tB << " T" ;
                EXPECT_EQ( tF->deval_dT( tB, tTheta, 1.9 ), 0.0 )
                        << tCase.mFile << " " << tCase.mLabel << ": B = " << tB << " T" ;
            }
        }
    }
}

#endif // BELFEM_HDF5
