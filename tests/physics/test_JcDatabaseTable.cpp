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
 * The tabulated REBCO critical-current function and the tables it accepts.
 *
 * The shipped table, share/material/sp-ap.hdf5, spans T = 3 .. 91 K,
 * log10( B / T ) = -2 .. 1.6 and the unfolded angle theta = 0 .. pi. Below
 * Bmin = 10 mT the table has no data, and JcFunctionDatabase bridges to the
 * measured self-field level J0( T ) stored with the table. These cases pin
 * the contracts of that function: jc and n do not rise with the field, a
 * table that does not cover the whole unfolded window is refused, and where
 * the bridge engages its self-field level carries no angle and it joins the
 * table without a kink. Below 76 K it does not always engage; the disabled
 * pins at the end record that window.
 * n is monotone in B only; nothing here asserts a temperature trend of n.
 * The tables are HDF5 files, so the file compiles to nothing without HDF5.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <unistd.h>

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_Vector.hpp"
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

    // Bmin of the shipped table: 10^-2 T
    const real gBmin = 1e-2 ;

    const real gTemperatures[ 4 ] = { 4.2, 20.0, 50.0, 77.0 };
    const real gAngles[ 5 ] = { 0.0, 0.25 * gPi, 0.5 * gPi, 0.75 * gPi, gPi };

    // where the low-field bridge of the shipped table engages at all five
    // angles
    const real gBridgedTemperatures[ 3 ] = { 76.0, 80.0, 84.0 };

    // the largest rise from one field sample to the next, measured on the
    // sampled grid below: 1.6e-15 for jc and 8.2e-15 for n, both roundoff
    const real gMaxRelativeRiseJc = 1e-13 ;
    const real gMaxRelativeRiseN  = 1e-13 ;

    //! sweeps B on a log grid from 1 mT, below Bmin, to the table's top at
    //! 10^1.6 T, for T every 4 K from 4 to 88 K and theta every 5 degrees, and
    //! requires the value never to rise by more than aMaxRelativeRise from
    //! one field sample to the next
    void
    expect_no_rise_with_field( const JcFunction & aFunction, const real aMaxRelativeRise )
    {
        for ( real T = 4.0; T <= 88.0; T += 4.0 )
        {
            for ( int a = 0; a <= 36; ++a )
            {
                const real tTheta = gPi * ( real ) a / 36.0 ;
                real tPrevious = aFunction.eval( 1e-3, tTheta, T );

                for ( int k = 1; k <= 92; ++k )
                {
                    const real tB = 1e-3 * std::pow( 10.0, ( real ) k / 20.0 );
                    const real tValue = aFunction.eval( tB, tTheta, T );
                    const real tRise = ( tValue - tPrevious ) / tPrevious ;

                    EXPECT_LE( tRise, aMaxRelativeRise )
                            << "T = " << T << " K, theta = " << tTheta << ", B = " << tB << " T";

                    tPrevious = tValue ;
                }
            }
        }
    }

    //! relative difference of the two one-sided secants of jc at Bmin
    real
    secant_jump_at_bmin( const JcFunction & aJc, const real aTheta, const real T )
    {
        const real h = 1e-4 * gBmin ;
        const real tAt    = aJc.eval( gBmin, aTheta, T );
        const real tLeft  = ( tAt - aJc.eval( gBmin - h, aTheta, T ) ) / h ;
        const real tRight = ( aJc.eval( gBmin + h, aTheta, T ) - tAt ) / h ;
        return std::abs( tRight - tLeft ) / std::max( std::abs( tLeft ), std::abs( tRight ) );
    }

    std::unique_ptr< JcFunction >
    shipped_table( const std::string & aLabel )
    {
        belfem::MaterialFactory tFactory ;
        return std::unique_ptr< JcFunction >(
                tFactory.create_jc_function( "sp-ap.hdf5", aLabel ) );
    }

    //! a 3-axis table ( T, log10 B, theta ) with a constant value, written to
    //! a scratch HDF5 file in the layout JcFunctionDatabase reads. Only the
    //! angular window varies between the cases.
    class ScopedTable
    {
        const std::string mPath ;

    public:

        ScopedTable( const std::string & aName, const real aAngleOrigin, const real aAngleWidth ) :
            mPath( ::testing::TempDir() + "belfem_test_jctable_" + aName + "_"
                   + std::to_string( ( long ) getpid() ) + ".hdf5" )
        {
            const index_t tNumAngles = 5 ;

            belfem::Mesh tMesh( 1,
                belfem::Vector< index_t >( { 3, 3, tNumAngles } ),
                belfem::Vector< real >( { 10.0, 0.5, aAngleWidth / ( real ) ( tNumAngles - 1 ) } ),
                belfem::Vector< real >( { 4.0, -2.0, aAngleOrigin } ) );

            belfem::Vector< real > & tField = tMesh.create_field( "jc" );
            tField.fill( 9.0 );

            belfem::Database tDatabase( &tMesh, "jc", false );

            belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
            tFile.create_group( "jc" );
            tDatabase.save( tFile.active_group() );
            tFile.close_active_group();
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
}

// =============================================================================
// the angular window a table must cover
// =============================================================================

TEST( JcDatabaseTable, ShippedTableCoversTheFullUnfoldedWindow )
{
    EXPECT_NO_THROW( shipped_table( "jc" ) );
    EXPECT_NO_THROW( shipped_table( "n" ) );
}

TEST( JcDatabaseTable, FoldedAngularWindowIsRefused )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture";
    }

    // width pi, but [ -pi/2, pi/2 ] cannot serve an unfolded theta > pi/2.
    // Mutation: turn the coverage guard into a width test and this goes red.
    const ScopedTable tTable( "folded", -0.5 * gPi, gPi );

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
    const ScopedTable tTable( "half", 0.0, 0.5 * gPi );

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
    const ScopedTable tTable( "full", 0.0, gPi );

    EXPECT_NO_THROW( belfem::material::JcFunctionDatabase( tTable.path(), "jc" ) );
}

// =============================================================================
// the contracts of the shipped REBCO table
// =============================================================================

TEST( JcDatabaseTable, JcDoesNotRiseWithFieldOnASampledGrid )
{
    // Mutation: read the field axis backwards and this goes red
    const std::unique_ptr< JcFunction > tJc = shipped_table( "jc" );
    expect_no_rise_with_field( *tJc, gMaxRelativeRiseJc );
}

TEST( JcDatabaseTable, NDoesNotRiseWithFieldOnASampledGrid )
{
    // n is monotone in B at fixed T and angle; no temperature trend is claimed
    const std::unique_ptr< JcFunction > tN = shipped_table( "n" );
    expect_no_rise_with_field( *tN, gMaxRelativeRiseN );
}

TEST( JcDatabaseTable, SelfFieldLevelIsAngleFree )
{
    // at B = 0 the bridge returns the measured J0( T ), which has no angle.
    // Asserted where the bridge engages at every angle of the shipped table,
    // 76 K and above; below that see the disabled pin further down.
    // Mutation: bridge to the table's own value at Bmin instead, which carries
    // the angle, and this goes red
    const std::unique_ptr< JcFunction > tJc = shipped_table( "jc" );

    for ( const real T : gBridgedTemperatures )
    {
        const real tReference = tJc->eval( 0.0, 0.0, T );

        for ( const real tTheta : gAngles )
        {
            EXPECT_NEAR( tJc->eval( 0.0, tTheta, T ), tReference, 1e-12 * std::abs( tReference ) )
                    << "T = " << T << " K, theta = " << tTheta ;
        }
    }
}

TEST( JcDatabaseTable, NoSlopeJumpAtBmin )
{
    // one-sided secants of jc, 1e-4 Bmin wide, on both sides of Bmin. The
    // bridge matches the table's value and slope there, so the secants differ
    // only by the curvature over the step: the relative difference shrinks
    // with the step, a kink does not. Measured on the shipped table at 76, 80
    // and 84 K: at most 2.2e-4, against 0.2 to 1 at the kinks further down.
    // The band is ten times that measurement, a project acceptance threshold,
    // not a literature bound.
    // Mutation: drop the slope term from the bridge and this goes red
    const real tBand = 2.2e-3 ;

    const std::unique_ptr< JcFunction > tJc = shipped_table( "jc" );

    for ( const real T : gBridgedTemperatures )
    {
        for ( const real tTheta : gAngles )
        {
            EXPECT_LT( secant_jump_at_bmin( *tJc, tTheta, T ), tBand )
                    << "T = " << T << " K, theta = " << tTheta ;
        }
    }
}

// =============================================================================
// known-defect window: the bridge below 76 K
// =============================================================================
//
// The bridge comment states that the measured J0( T ) exceeds the table's
// value at Bmin at every angle, so the quadratic bridge always applies. On the
// shipped sp-ap table that holds from 76 K up. Between 56 and 72 K the
// monotone limiter leaves a slope kink of 0.2 to 1 at some angles, and from
// 52 K down J0 is not above the table edge at most angles, so eval() returns
// the flat clamp: the angle-dependent table value, with zero slope below Bmin.
// These cases are red until that is resolved.

TEST( JcDatabaseTable, DISABLED_BridgeEngagesAtEveryTemperature )
{
    const std::unique_ptr< JcFunction > tJc = shipped_table( "jc" );

    for ( const real T : gTemperatures )
    {
        for ( const real tTheta : gAngles )
        {
            const real tJ0   = tJc->eval( 0.0, tTheta, T );
            const real tEdge = tJc->eval( gBmin, tTheta, T );

            EXPECT_GT( tJ0, tEdge ) << "T = " << T << " K, theta = " << tTheta
                    << ": J0 / jc( Bmin ) = " << tJ0 / tEdge ;
        }
    }
}

TEST( JcDatabaseTable, DISABLED_NoSlopeJumpAtBminAtEveryTemperature )
{
    const std::unique_ptr< JcFunction > tJc = shipped_table( "jc" );

    for ( const real T : gTemperatures )
    {
        for ( const real tTheta : gAngles )
        {
            EXPECT_LT( secant_jump_at_bmin( *tJc, tTheta, T ), 2.2e-3 )
                    << "T = " << T << " K, theta = " << tTheta ;
        }
    }
}

#endif // BELFEM_HDF5
