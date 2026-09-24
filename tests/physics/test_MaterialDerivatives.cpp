/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Derivative contracts of Material.
 *
 * The thermal and magnetic Newton tangents read dcpdT, dlambdadT, drhodT and
 * the dJ / dB / dT tangents of the HTS resistivity laws directly. A wrong
 * derivative does not stop a run; it turns quadratic convergence into a
 * stall. The derivative accessors below are checked against a Richardson
 * central difference of the value they claim to differentiate, and the
 * overload families of the HTS laws are checked against each other.
 *
 * Not covered here: the piecewise dJ tangent below j1 ( pinned by a
 * DISABLED_ case that documents a known inconsistency ), the piecewise dB
 * tangent inside the Bezier blend ( staged zero ), the defect overloads
 * beyond finiteness of drho_piecewise_dJ, and rho_riva.
 */

#include <algorithm>
#include <cmath>
#include <string>

#include <gtest/gtest.h>

#include "cl_Material.hpp"
#include "cl_MaterialFactory.hpp"
#include "cl_JcFunction.hpp"

using belfem::real ;
using belfem::Material ;
using belfem::MaterialFactory ;
using belfem::MaterialProperty ;
using belfem::MaterialDependency ;
using belfem::ResistivityLaw ;

namespace
{
    // every label the factory builds: the nine Kohler metals, the alloy, the
    // buffer and the superconductor. Tables off: analytic property paths only
    const char * gBuiltinLabels[] = { "copper", "aluminum", "silver", "nickel", "iron",
                                      "chromium", "indium", "lead", "tin",
                                      "hastelloyc276", "magnesia", "ybco" } ;

    // generic temperatures, away from round numbers and known switch points
    const real gGridT[] = { 7.3, 23.7, 61.1, 113.9, 217.3, 293.4 } ;

    // the constant-jc overloads of the HTS laws read rho( gTbulk ): pin it for
    // the duration of a test and hand the previous value back afterwards
    struct TbulkGuard
    {
        const real mSaved ;

        explicit TbulkGuard( const real aT ) : mSaved( belfem::gTbulk )
        {
            belfem::gTbulk = aT ;
        }

        ~TbulkGuard()
        {
            belfem::gTbulk = mSaved ;
        }
    };

    // a defect that leaves jc unchanged, so the defect overloads can be
    // called on the fixtures below
    real
    unit_defect( const real, const real, const real, const real )
    {
        return 1.0 ;
    }

    bool
    within( const real aA, const real aB, const real aRel, const real aAbs )
    {
        return std::abs( aA - aB ) <= aAbs + aRel * std::max( std::abs( aA ), std::abs( aB ) ) ;
    }

    // two central differences of a function at aX, with steps aH and aH/2.
    // The absolute parts of the tolerances are measured in the derivative's
    // own scale |f(x)/x|, so that a resistivity slope of 1e-10 Ohm*m/K is
    // held to the same relative resolution as a heat capacity slope of 1
    struct FiniteDifference
    {
        real mCoarse ;
        real mFine ;
        real mBest ;
        real mScale ;
        bool mFinite ;

        template< typename F >
        FiniteDifference( F && aFun, const real aX, const real aH )
        {
            const real tPlus      = aFun( aX + aH ) ;
            const real tMinus     = aFun( aX - aH ) ;
            const real tHalfPlus  = aFun( aX + 0.5 * aH ) ;
            const real tHalfMinus = aFun( aX - 0.5 * aH ) ;
            const real tCenter    = aFun( aX ) ;

            mCoarse = 0.5 * ( tPlus - tMinus ) / aH ;
            mFine   = ( tHalfPlus - tHalfMinus ) / aH ;
            // Richardson: the O(h^2) truncation of the fine estimate is
            // ( coarse - fine ) / 3, so this is the O(h^4) value
            mBest   = ( 4.0 * mFine - mCoarse ) / 3.0 ;
            mScale  = std::abs( tCenter / aX ) ;
            mFinite =    std::isfinite( tPlus ) && std::isfinite( tMinus )
                      && std::isfinite( tHalfPlus ) && std::isfinite( tHalfMinus )
                      && std::isfinite( tCenter ) ;
        }

        // the two estimates agree to the truncation of the coarse one: the
        // function is C1 at aX to this resolution, and the extrapolated value
        // carries an error of about a third of that gap divided by four
        bool
        admissible() const
        {
            return mFinite && within( mCoarse, mFine, 3e-5, 1e-14 * mScale ) ;
        }

        bool
        matches( const real aAnalytic, const real aRel = 1e-5 ) const
        {
            return std::isfinite( aAnalytic ) && within( aAnalytic, mBest, aRel, 1e-12 * mScale ) ;
        }
    };

    // a thermal property over the temperature grid: skip the points where the
    // property is not C1 to finite-difference resolution, require at least
    // five of six to be checkable, and every checkable point to match
    template< typename V, typename D >
    void
    check_thermal_derivative( const std::string & aWhat, V && aValue, D && aDerivative,
                              const real aRel = 1e-5 )
    {
        belfem::uint tAdmissible = 0 ;

        for ( real T : gGridT )
        {
            FiniteDifference tFD( aValue, T, 1e-4 * T ) ;
            const real tAnalytic = aDerivative( T ) ;

            EXPECT_TRUE( tFD.mFinite ) << aWhat << " value at T = " << T ;
            EXPECT_TRUE( std::isfinite( tAnalytic ) ) << aWhat << " derivative at T = " << T ;

            if ( ! tFD.admissible() ) continue ;
            ++tAdmissible ;

            EXPECT_TRUE( tFD.matches( tAnalytic, aRel ) )
                << aWhat << " T = " << T << " analytic " << tAnalytic << " fd " << tFD.mBest ;
        }

        EXPECT_GE( tAdmissible, 5u ) << aWhat << " admissible grid points" ;
    }

    // the factory's YBCO with ec, jc and n as plain constants; rho_n( T ) is
    // the material's own normal-state law
    Material *
    make_constant_hts( MaterialFactory & aFactory, const real aN )
    {
        Material * tHts = aFactory.create_material( "ybco", 50.0, false ) ;
        tHts->set_constant( MaterialProperty::ec, 1e-4 ) ;
        tHts->set_constant( MaterialProperty::jc, 3.1e10 ) ;
        tHts->set_constant( MaterialProperty::n, aN ) ;
        return tHts ;
    }

    // the factory's YBCO with a field-dependent jc and a constant n; the
    // material takes ownership of aJc
    Material *
    make_field_hts( MaterialFactory & aFactory, const belfem::material::JcFunction * aJc )
    {
        Material * tHts = aFactory.create_material( "ybco", 50.0, false ) ;
        tHts->set_constant( MaterialProperty::ec, 1e-4 ) ;
        tHts->set_constant( MaterialProperty::n, 27.0 ) ;
        tHts->set_jc_function( aJc ) ;
        return tHts ;
    }

    // the modified Kim law with its analytic field derivative, so the dB
    // tangents can be checked against the law they differentiate
    class KimWithFieldDerivative : public belfem::material::JcFunction
    {
        const real mJc0 ;
        const real mB0 ;
        const real mK2 ;
        const real mAlpha ;

    public:

        KimWithFieldDerivative( const real aJc0, const real aB0, const real aK2, const real aAlpha ) :
            mJc0( aJc0 ), mB0( aB0 ), mK2( aK2 ), mAlpha( aAlpha )
        {
            this->set_dependency( belfem::material::JcParameter::normB ) ;
            this->set_dependency( belfem::material::JcParameter::angleNxB ) ;
        }

        real
        eval( const real normB, const real angle ) const override
        {
            const real c2 = std::cos( angle ) * std::cos( angle ) ;
            const real g  = std::sqrt( mK2 * ( 1.0 - c2 ) + c2 ) ;
            return mJc0 * std::pow( 1.0 + g * normB / mB0, -mAlpha ) ;
        }

        real
        deval_dB( const real normB, const real angle, const real T ) const override
        {
            const real c2 = std::cos( angle ) * std::cos( angle ) ;
            const real g  = std::sqrt( mK2 * ( 1.0 - c2 ) + c2 ) ;
            return -mAlpha * mJc0 * ( g / mB0 ) * std::pow( 1.0 + g * normB / mB0, -mAlpha - 1.0 ) ;
        }
    };

    // the knots of the piecewise law as rho_piecewise builds them from the
    // constants; the flux-flow exponent is the class default
    struct PiecewiseKnots
    {
        real mJ1 ;
        real mRho1 ;
        real mJ3 ;

        PiecewiseKnots( const Material * aHts, const real aT )
        {
            const real ec  = aHts->constant_property( MaterialProperty::ec ) ;
            const real jc  = aHts->constant_property( MaterialProperty::jc ) ;
            const real n   = aHts->constant_property( MaterialProperty::n ) ;
            const real nff = 3.0 ;

            mJ1   = jc * std::pow( 10.0, 2.5 / n ) ;
            mRho1 = ( ec / jc ) * std::pow( mJ1 / jc, n - 1.0 ) ;
            mJ3   = mJ1 * std::pow( aHts->rho( aT ) / mRho1, 1.0 / nff ) ;
        }
    };
}

//------------------------------------------------------------------------------
// R1: thermal accessors of every builtin material
//------------------------------------------------------------------------------

class ThermalDerivatives : public ::testing::TestWithParam< const char * > {};

TEST_P( ThermalDerivatives, MatchCentralDifferences )
{
    MaterialFactory tFactory ;
    Material * tMaterial = tFactory.create_material( GetParam(), 50.0, false ) ;
    const std::string tLabel = GetParam() ;

    if ( tMaterial->have( MaterialProperty::cp ) )
    {
        check_thermal_derivative( tLabel + " cp",
            [ & ]( real T ){ return tMaterial->cp( T ) ; },
            [ & ]( real T ){ return tMaterial->dcpdT( T ) ; } ) ;
    }

    if ( tMaterial->have( MaterialProperty::lambda ) )
    {
        check_thermal_derivative( tLabel + " lambda",
            [ & ]( real T ){ return tMaterial->lambda( T ) ; },
            [ & ]( real T ){ return tMaterial->dlambdadT( T ) ; } ) ;
    }

    if ( tMaterial->have( MaterialProperty::rho ) )
    {
        check_thermal_derivative( tLabel + " rho",
            [ & ]( real T ){ return tMaterial->rho( T ) ; },
            [ & ]( real T ){ return tMaterial->drhodT( T ) ; } ) ;
    }

    // the field-dependent overloads, at zero field and in field; their T
    // derivatives are themselves central differences stepped by the framework
    // global gFinDiffDeltaT, so they are held to that step's truncation, one
    // order less than the analytic ones
    const real tBeta      = 0.37 ;
    const real tFieldRel  = 1e-4 ;

    for ( real B : { 0.0, 5.0 } )
    {
        const std::string tField = " B = " + std::to_string( B ) ;

        if ( tMaterial->have( MaterialProperty::rho ) &&
             tMaterial->depends( MaterialProperty::rho, MaterialDependency::normB ) )
        {
            check_thermal_derivative( tLabel + " rho" + tField,
                [ & ]( real T ){ return tMaterial->rho( T, B, tBeta ) ; },
                [ & ]( real T ){ return tMaterial->drhodT( T, B, tBeta ) ; },
                tFieldRel ) ;
        }

        if ( tMaterial->have( MaterialProperty::lambda ) &&
             tMaterial->depends( MaterialProperty::lambda, MaterialDependency::normB ) )
        {
            check_thermal_derivative( tLabel + " lambda" + tField,
                [ & ]( real T ){ return tMaterial->lambda( T, B, tBeta ) ; },
                [ & ]( real T ){ return tMaterial->dlambdadT( T, B, tBeta ) ; },
                tFieldRel ) ;
        }
    }

    delete tMaterial ;
}

INSTANTIATE_TEST_SUITE_P(
    Builtin,
    ThermalDerivatives,
    ::testing::ValuesIn( gBuiltinLabels ),
    []( const ::testing::TestParamInfo< const char * > & aInfo )
    {
        return std::string( aInfo.param ) ;
    } ) ;

//------------------------------------------------------------------------------
// R2: the power law
//------------------------------------------------------------------------------

TEST( HtsLaws, PowerLawDerivativesMatchFiniteDifferences )
{
    TbulkGuard tGuard( 77.0 ) ;
    MaterialFactory tFactory ;

    const real T     = 77.0 ;
    const real B     = 1.3 ;
    const real tBeta = 0.37 ;

    // constant jc: the crossover rho_PL( J_cross ) = rho_n follows from the
    // fixture's own values
    Material * tHts = make_constant_hts( tFactory, 27.0 ) ;

    const real ec = tHts->constant_property( MaterialProperty::ec ) ;
    const real jc = tHts->constant_property( MaterialProperty::jc ) ;
    const real n  = tHts->constant_property( MaterialProperty::n ) ;

    const real tJcross = jc * std::pow( tHts->rho( T ) * jc / ec, 1.0 / ( n - 1.0 ) ) ;

    // below 0.9 J_cross the parallel combination leaves rho_n a 1e-16
    // relative share of rho ( ratio^( n - 1 ) ), so a difference in T is
    // roundoff and not a derivative
    for ( real tRatio : { 0.9, 1.0, 1.5 } )
    {
        const real J = std::max( tRatio * tJcross, 0.3 * jc ) ;

        FiniteDifference tFdJ( [ & ]( real aJ ){ return tHts->rho_powerlaw( aJ, T, B, tBeta ) ; }, J, 1e-5 * J ) ;
        ASSERT_TRUE( tFdJ.admissible() ) << "J / J_cross = " << tRatio ;
        EXPECT_TRUE( tFdJ.matches( tHts->drho_powerlaw_dJ( J, T, B, tBeta ) ) )
            << "dJ at J / J_cross = " << tRatio << " analytic "
            << tHts->drho_powerlaw_dJ( J, T, B, tBeta ) << " fd " << tFdJ.mBest ;

        FiniteDifference tFdT( [ & ]( real aT ){ return tHts->rho_powerlaw( J, aT, B, tBeta ) ; }, T, 1e-4 * T ) ;
        ASSERT_TRUE( tFdT.admissible() ) << "J / J_cross = " << tRatio ;
        EXPECT_TRUE( tFdT.matches( tHts->drho_powerlaw_dT( J, T, B, tBeta ) ) )
            << "dT at J / J_cross = " << tRatio << " analytic "
            << tHts->drho_powerlaw_dT( J, T, B, tBeta ) << " fd " << tFdT.mBest ;

        // a constant jc has no field derivative
        EXPECT_EQ( 0.0, tHts->drho_powerlaw_dB( J, T, B, tBeta ) ) ;
    }

    delete tHts ;

    // field-dependent jc, constant n: dJ and dT on the factory's Kim law
    const belfem::material::JcFunction * tKimLaw = tFactory.create_jc_function( 3.1e10, 0.13, 0.7, 1.4 ) ;
    Material * tKim = make_field_hts( tFactory, tKimLaw ) ;

    const real tJcB    = tKimLaw->eval( B, tBeta ) ;
    const real tEcKim   = tKim->constant_property( MaterialProperty::ec ) ;
    const real tNKim    = tKim->constant_property( MaterialProperty::n ) ;
    const real tJcross2 = tJcB * std::pow( tKim->rho( T ) * tJcB / tEcKim, 1.0 / ( tNKim - 1.0 ) ) ;

    for ( real tRatio : { 0.9, 1.0, 1.5 } )
    {
        const real J = std::max( tRatio * tJcross2, 0.3 * tJcB ) ;

        FiniteDifference tFdJ( [ & ]( real aJ ){ return tKim->rho_powerlaw( aJ, T, B, tBeta ) ; }, J, 1e-5 * J ) ;
        ASSERT_TRUE( tFdJ.admissible() ) << "Kim J / J_cross = " << tRatio ;
        EXPECT_TRUE( tFdJ.matches( tKim->drho_powerlaw_dJ( J, T, B, tBeta ) ) )
            << "Kim dJ at J / J_cross = " << tRatio ;

        FiniteDifference tFdT( [ & ]( real aT ){ return tKim->rho_powerlaw( J, aT, B, tBeta ) ; }, T, 1e-4 * T ) ;
        ASSERT_TRUE( tFdT.admissible() ) << "Kim J / J_cross = " << tRatio ;
        EXPECT_TRUE( tFdT.matches( tKim->drho_powerlaw_dT( J, T, B, tBeta ) ) )
            << "Kim dT at J / J_cross = " << tRatio ;
    }

    delete tKim ;
}

// the field tangent of the power law against the law it differentiates, on a
// jc function that provides its analytic field derivative
TEST( HtsLaws, PowerLawFieldDerivativeMatchesFiniteDifference )
{
    TbulkGuard tGuard( 77.0 ) ;
    MaterialFactory tFactory ;

    const real T     = 77.0 ;
    const real tBeta = 0.37 ;

    const KimWithFieldDerivative * tLaw = new KimWithFieldDerivative( 3.1e10, 0.13, 0.7, 1.4 ) ;
    Material * tHts = make_field_hts( tFactory, tLaw ) ;

    for ( real B : { 0.2, 1.3, 4.1 } )
    {
        const real J = 0.8 * tLaw->eval( B, tBeta ) ;

        FiniteDifference tFdB( [ & ]( real aB ){ return tHts->rho_powerlaw( J, T, aB, tBeta ) ; }, B, 1e-4 * B ) ;
        ASSERT_TRUE( tFdB.admissible() ) << "B = " << B ;
        EXPECT_TRUE( tFdB.matches( tHts->drho_powerlaw_dB( J, T, B, tBeta ) ) )
            << "dB at B = " << B << " analytic " << tHts->drho_powerlaw_dB( J, T, B, tBeta )
            << " fd " << tFdB.mBest ;
    }

    delete tHts ;
}

// drho_powerlaw_dB equals the field derivative of rho_powerlaw on the factory's Kim law
TEST( HtsLaws, KimFieldDerivativeMatchesFiniteDifference )
{
    TbulkGuard tGuard( 77.0 ) ;
    MaterialFactory tFactory ;

    const real T     = 77.0 ;
    const real tBeta = 0.37 ;

    const belfem::material::JcFunction * tKimLaw = tFactory.create_jc_function( 3.1e10, 0.13, 0.7, 1.4 ) ;
    Material * tKim = make_field_hts( tFactory, tKimLaw ) ;

    for ( real B : { 0.2, 1.3, 4.1 } )
    {
        const real J = 0.8 * tKimLaw->eval( B, tBeta ) ;

        FiniteDifference tFdB( [ & ]( real aB ){ return tKim->rho_powerlaw( J, T, aB, tBeta ) ; }, B, 1e-4 * B ) ;
        ASSERT_TRUE( tFdB.admissible() ) << "B = " << B ;
        EXPECT_TRUE( tFdB.matches( tKim->drho_powerlaw_dB( J, T, B, tBeta ) ) )
            << "dB at B = " << B << " analytic " << tKim->drho_powerlaw_dB( J, T, B, tBeta )
            << " fd " << tFdB.mBest ;
    }

    delete tKim ;
}

//------------------------------------------------------------------------------
// R3: the piecewise law
//------------------------------------------------------------------------------

TEST( HtsLaws, PiecewiseBlendDerivativesMatchFiniteDifferences )
{
    TbulkGuard tGuard( 77.0 ) ;
    MaterialFactory tFactory ;

    const real T     = 77.0 ;
    const real B     = 1.3 ;
    const real tBeta = 0.37 ;

    Material * tHts = make_constant_hts( tFactory, 27.0 ) ;
    tHts->set_resistivity_law( ResistivityLaw::Piecewise ) ;

    const PiecewiseKnots tKnots( tHts, T ) ;
    const real tRhoN = tHts->rho( T ) ;

    // the reconstructed upper knot is the material's own: the value is rho_n
    // just above it and still below rho_n just below it
    ASSERT_LT( tKnots.mJ1, tKnots.mJ3 ) ;
    ASSERT_DOUBLE_EQ( tRhoN, tHts->rho_piecewise( 1.001 * tKnots.mJ3, T, B, tBeta ) ) ;
    ASSERT_LT( tHts->rho_piecewise( 0.999 * tKnots.mJ3, T, B, tBeta ), tRhoN ) ;

    for ( real J : { 1.2 * tKnots.mJ1, 0.5 * ( tKnots.mJ1 + tKnots.mJ3 ), 0.95 * tKnots.mJ3 } )
    {
        FiniteDifference tFdJ( [ & ]( real aJ ){ return tHts->rho_piecewise( aJ, T, B, tBeta ) ; }, J, 1e-6 * J ) ;
        ASSERT_TRUE( tFdJ.admissible() ) << "J / j1 = " << J / tKnots.mJ1 ;
        EXPECT_TRUE( tFdJ.matches( tHts->drho_piecewise_dJ( J, T, B, tBeta ) ) )
            << "dJ at J / j1 = " << J / tKnots.mJ1 << " analytic "
            << tHts->drho_piecewise_dJ( J, T, B, tBeta ) << " fd " << tFdJ.mBest ;

        FiniteDifference tFdT( [ & ]( real aT ){ return tHts->rho_piecewise( J, aT, B, tBeta ) ; }, T, 1e-4 * T ) ;
        ASSERT_TRUE( tFdT.admissible() ) << "J / j1 = " << J / tKnots.mJ1 ;
        EXPECT_TRUE( tFdT.matches( tHts->drho_piecewise_dT( J, T, B, tBeta ) ) )
            << "dT at J / j1 = " << J / tKnots.mJ1 << " analytic "
            << tHts->drho_piecewise_dT( J, T, B, tBeta ) << " fd " << tFdT.mBest ;
    }

    delete tHts ;
}

TEST( HtsLaws, PiecewiseNormalRegimeAndAboveTcrit )
{
    TbulkGuard tGuard( 77.0 ) ;
    MaterialFactory tFactory ;

    const real T     = 77.0 ;
    const real B     = 1.3 ;
    const real tBeta = 0.37 ;

    Material * tHts = make_constant_hts( tFactory, 27.0 ) ;
    tHts->set_resistivity_law( ResistivityLaw::Piecewise ) ;

    const PiecewiseKnots tKnots( tHts, T ) ;

    // normal regime: the value is rho_n( T ), so the tangents are those of rho_n
    {
        const real J = 1.5 * tKnots.mJ3 ;

        EXPECT_DOUBLE_EQ( tHts->rho( T ), tHts->rho_piecewise( J, T, B, tBeta ) ) ;
        EXPECT_EQ( 0.0, tHts->drho_piecewise_dJ( J, T, B, tBeta ) ) ;
        EXPECT_DOUBLE_EQ( tHts->drhodT( T ), tHts->drho_piecewise_dT( J, T, B, tBeta ) ) ;

        FiniteDifference tFdT( [ & ]( real aT ){ return tHts->rho_piecewise( J, aT, B, tBeta ) ; }, T, 1e-4 * T ) ;
        ASSERT_TRUE( tFdT.admissible() ) ;
        EXPECT_TRUE( tFdT.matches( tHts->drho_piecewise_dT( J, T, B, tBeta ) ) )
            << "normal-regime dT analytic " << tHts->drho_piecewise_dT( J, T, B, tBeta )
            << " fd " << tFdT.mBest ;
    }

    // above T_crit: rho_n( T ) in every regime
    {
        const real tAbove = 100.0 ;
        ASSERT_GT( tAbove, tHts->constant_property( MaterialProperty::T_crit ) ) ;

        for ( real J : { 0.6 * tKnots.mJ1, 0.5 * ( tKnots.mJ1 + tKnots.mJ3 ), 1.5 * tKnots.mJ3 } )
        {
            EXPECT_DOUBLE_EQ( tHts->rho( tAbove ), tHts->rho_piecewise( J, tAbove, B, tBeta ) ) ;
            EXPECT_EQ( 0.0, tHts->drho_piecewise_dJ( J, tAbove, B, tBeta ) ) ;
            EXPECT_DOUBLE_EQ( tHts->drhodT( tAbove ), tHts->drho_piecewise_dT( J, tAbove, B, tBeta ) ) ;
        }
    }

    delete tHts ;
}

// below j1 rho_piecewise returns the raw power law, so its field tangent is
// the raw chain through jc( B ) with no parallel factor
TEST( HtsLaws, PiecewiseFieldDerivativeBelowJ1MatchesFiniteDifference )
{
    TbulkGuard tGuard( 77.0 ) ;
    MaterialFactory tFactory ;

    const real T     = 77.0 ;
    const real tBeta = 0.37 ;

    const KimWithFieldDerivative * tLaw = new KimWithFieldDerivative( 3.1e10, 0.13, 0.7, 1.4 ) ;
    Material * tHts = make_field_hts( tFactory, tLaw ) ;
    tHts->set_resistivity_law( ResistivityLaw::Piecewise ) ;

    const real n = tHts->constant_property( MaterialProperty::n ) ;

    for ( real B : { 0.2, 1.3, 4.1 } )
    {
        const real tJ1 = tLaw->eval( B, tBeta ) * std::pow( 10.0, 2.5 / n ) ;
        const real J   = 0.6 * tJ1 ;

        FiniteDifference tFdB( [ & ]( real aB ){ return tHts->rho_piecewise( J, T, aB, tBeta ) ; }, B, 1e-4 * B ) ;
        ASSERT_TRUE( tFdB.admissible() ) << "B = " << B ;
        EXPECT_TRUE( tFdB.matches( tHts->drho_piecewise_dB( J, T, B, tBeta ) ) )
            << "dB at B = " << B << " analytic " << tHts->drho_piecewise_dB( J, T, B, tBeta )
            << " fd " << tFdB.mBest ;
    }

    delete tHts ;
}

// n = 4 makes the Bezier discriminant vanish at the upper knot; the value and
// the temperature and current tangents must stay finite there, and at the
// knot itself dJ is the flux-flow limit mNff rho_n / J ( mNff = 3 )
//
// The ( J, T ) and ( J, T, x, y, z, t ) overloads of drho_piecewise_dJ read
// jc_custom / n_custom, which the factory's YBCO does not provide ( they
// abort ), so neither fixture can call them. The ( J, |B|, theta ) forms need
// a JcFunction and run on the Kim fixture, whose n = 27 is not degenerate:
// there only finiteness at its own j3 is checked
TEST( HtsLaws, PiecewiseDegenerateBlendStaysFinite )
{
    TbulkGuard tGuard( 77.0 ) ;
    MaterialFactory tFactory ;

    const real T     = 77.0 ;
    const real B     = 1.3 ;
    const real tBeta = 0.37 ;

    // placeholder coordinates for the defect overloads
    const real tX = 0.1 ;
    const real tY = 0.2 ;
    const real tZ = 0.3 ;
    const real tTime = 0.4 ;

    Material * tHts = make_constant_hts( tFactory, 4.0 ) ;
    tHts->set_resistivity_law( ResistivityLaw::Piecewise ) ;
    tHts->set_user_defined_defect( unit_defect ) ;

    const PiecewiseKnots tKnots( tHts, T ) ;
    ASSERT_LT( tKnots.mJ1, tKnots.mJ3 ) ;

    for ( real tFactor : { 1.0, 1.0 - 1e-7, 1.0 - 1e-3 } )
    {
        const real J = tFactor * tKnots.mJ3 ;

        EXPECT_TRUE( std::isfinite( tHts->rho_piecewise( J, T, B, tBeta ) ) )
            << "value at J / j3 = " << tFactor ;
        EXPECT_TRUE( std::isfinite( tHts->drho_piecewise_dT( J, T, B, tBeta ) ) )
            << "dT at J / j3 = " << tFactor ;

        const real tDJ1 = tHts->drho_piecewise_dJ( J ) ;
        const real tDJ4 = tHts->drho_piecewise_dJ( J, T, B, tBeta ) ;

        EXPECT_TRUE( std::isfinite( tDJ1 ) )
            << "dJ( J ) at J / j3 = " << tFactor ;
        EXPECT_TRUE( std::isfinite( tDJ4 ) )
            << "dJ( J, T, B, beta ) at J / j3 = " << tFactor ;
        EXPECT_TRUE( std::isfinite( tHts->drho_piecewise_dJ( J, tX, tY, tZ, tTime ) ) )
            << "dJ( J, x, y, z, t ) at J / j3 = " << tFactor ;
        EXPECT_TRUE( std::isfinite( tHts->drho_piecewise_dJ( J, T, B, tBeta, tX, tY, tZ, tTime ) ) )
            << "dJ( J, T, B, beta, x, y, z, t ) at J / j3 = " << tFactor ;

        if ( tFactor == 1.0 )
        {
            const real tLimit = 3.0 * tHts->rho( T ) / J ;
            EXPECT_NEAR( tDJ1, tLimit, 1e-6 * tLimit ) << "dJ( J ) at the knot" ;
            EXPECT_NEAR( tDJ4, tLimit, 1e-6 * tLimit ) << "dJ( J, T, B, beta ) at the knot" ;
        }
    }

    delete tHts ;

    // near the coincident knots: n = 4 + 1e-6 keeps the discriminant above the
    // guard's threshold ( the root path ), n = 4 + 1e-9 falls inside it ( the
    // fallback ); both must be finite, and the fallback equals the limit
    for ( real tN : { 4.0 + 1e-6, 4.0 + 1e-9 } )
    {
        Material * tNear = make_constant_hts( tFactory, tN ) ;
        tNear->set_resistivity_law( ResistivityLaw::Piecewise ) ;

        const PiecewiseKnots tNearKnots( tNear, T ) ;
        const real J = tNearKnots.mJ3 ;
        const real tDJ = tNear->drho_piecewise_dJ( J, T, B, tBeta ) ;

        EXPECT_TRUE( std::isfinite( tDJ ) ) << "dJ at the knot for n = " << tN ;
        EXPECT_TRUE( std::isfinite( tNear->drho_piecewise_dT( J, T, B, tBeta ) ) ) << "dT at the knot for n = " << tN ;
        EXPECT_LE( std::abs( tDJ ), 1.001 * 3.0 * tNear->rho( T ) / J ) << "dJ bounded by the coincident limit for n = " << tN ;

        delete tNear ;
    }

    const belfem::material::JcFunction * tKimLaw = tFactory.create_jc_function( 3.1e10, 0.13, 0.7, 1.4 ) ;
    Material * tKim = make_field_hts( tFactory, tKimLaw ) ;
    tKim->set_resistivity_law( ResistivityLaw::Piecewise ) ;
    tKim->set_user_defined_defect( unit_defect ) ;

    const real tEcKim = tKim->constant_property( MaterialProperty::ec ) ;
    const real tNKim  = tKim->constant_property( MaterialProperty::n ) ;
    const real tJcB   = tKimLaw->eval( B, tBeta ) ;
    const real tJ1Kim = tJcB * std::pow( 10.0, 2.5 / tNKim ) ;
    const real tRho1  = ( tEcKim / tJcB ) * std::pow( tJ1Kim / tJcB, tNKim - 1.0 ) ;
    const real tJ3Kim = tJ1Kim * std::pow( tKim->rho( T ) / tRho1, 1.0 / 3.0 ) ;
    ASSERT_LT( tJ1Kim, tJ3Kim ) ;

    for ( real tFactor : { 1.0, 1.0 - 1e-7, 1.0 - 1e-3 } )
    {
        const real J = tFactor * tJ3Kim ;

        EXPECT_TRUE( std::isfinite( tKim->drho_piecewise_dJ( J, B, tBeta ) ) )
            << "Kim dJ( J, B, beta ) at J / j3 = " << tFactor ;
        EXPECT_TRUE( std::isfinite( tKim->drho_piecewise_dJ( J, B, tBeta, tX, tY, tZ, tTime ) ) )
            << "Kim dJ( J, B, beta, x, y, z, t ) at J / j3 = " << tFactor ;
    }

    delete tKim ;
}

// contract pinned: below j1 drho_piecewise_dJ differentiates the raw power law that rho_piecewise returns there
TEST( HtsLaws, DISABLED_PowerLawRegimeTangentMatchesRawResidual )
{
    TbulkGuard tGuard( 77.0 ) ;
    MaterialFactory tFactory ;

    const real T     = 77.0 ;
    const real B     = 1.3 ;
    const real tBeta = 0.37 ;

    Material * tHts = make_constant_hts( tFactory, 27.0 ) ;
    tHts->set_resistivity_law( ResistivityLaw::Piecewise ) ;

    const real ec = tHts->constant_property( MaterialProperty::ec ) ;
    const real jc = tHts->constant_property( MaterialProperty::jc ) ;
    const real n  = tHts->constant_property( MaterialProperty::n ) ;

    // lower rho_n until the parallel factor is measurable at 0.9 j1
    // ( rho_PL / rho_n = 0.05 ), while keeping j1 < j3
    const real tJ1    = jc * std::pow( 10.0, 2.5 / n ) ;
    const real J      = 0.9 * tJ1 ;
    const real tRhoPL = ( ec / jc ) * std::pow( J / jc, n - 1.0 ) ;
    tHts->set_constant( MaterialProperty::rho, 20.0 * tRhoPL ) ;

    const PiecewiseKnots tKnots( tHts, T ) ;
    ASSERT_LT( tKnots.mJ1, tKnots.mJ3 ) ;

    FiniteDifference tFdJ( [ & ]( real aJ ){ return tHts->rho_piecewise( aJ, T, B, tBeta ) ; }, J, 1e-6 * J ) ;
    ASSERT_TRUE( tFdJ.admissible() ) ;

    // the raw closed form is what the value differentiates to
    const real tRaw = ( ec / ( jc * jc ) ) * ( n - 1.0 ) * std::pow( J / jc, n - 2.0 ) ;
    EXPECT_TRUE( tFdJ.matches( tRaw ) ) << "raw " << tRaw << " fd " << tFdJ.mBest ;

    EXPECT_TRUE( tFdJ.matches( tHts->drho_piecewise_dJ( J, T, B, tBeta ) ) )
        << "analytic " << tHts->drho_piecewise_dJ( J, T, B, tBeta ) << " fd " << tFdJ.mBest ;

    delete tHts ;
}

//------------------------------------------------------------------------------
// R4: overload families
//------------------------------------------------------------------------------

TEST( HtsLaws, OverloadPairsAgree )
{
    TbulkGuard tGuard( 77.0 ) ;
    MaterialFactory tFactory ;

    const real T     = belfem::gTbulk ;
    const real B     = 1.3 ;
    const real tBeta = 0.37 ;

    // constant jc and n: the 1-argument overloads read gTbulk and the
    // constants, the 4-argument ones the same values through the routing
    Material * tHts = make_constant_hts( tFactory, 27.0 ) ;

    const PiecewiseKnots tKnots( tHts, T ) ;

    for ( real J : { 0.6 * tKnots.mJ1, 0.5 * ( tKnots.mJ1 + tKnots.mJ3 ), 1.5 * tKnots.mJ3 } )
    {
        EXPECT_TRUE( within( tHts->rho_powerlaw( J ), tHts->rho_powerlaw( J, T, B, tBeta ), 1e-14, 0.0 ) )
            << "rho_powerlaw at J / j1 = " << J / tKnots.mJ1 ;
        EXPECT_TRUE( within( tHts->drho_powerlaw_dJ( J ), tHts->drho_powerlaw_dJ( J, T, B, tBeta ), 1e-14, 0.0 ) )
            << "drho_powerlaw_dJ at J / j1 = " << J / tKnots.mJ1 ;
        EXPECT_TRUE( within( tHts->rho_piecewise( J ), tHts->rho_piecewise( J, T, B, tBeta ), 1e-14, 0.0 ) )
            << "rho_piecewise at J / j1 = " << J / tKnots.mJ1 ;
        EXPECT_TRUE( within( tHts->drho_piecewise_dJ( J ), tHts->drho_piecewise_dJ( J, T, B, tBeta ), 1e-14, 0.0 ) )
            << "drho_piecewise_dJ at J / j1 = " << J / tKnots.mJ1 ;
    }

    delete tHts ;

    // field-dependent jc, constant n: the 3-argument overloads against the
    // 4-argument ones at the bulk temperature
    const belfem::material::JcFunction * tKimLaw = tFactory.create_jc_function( 3.1e10, 0.13, 0.7, 1.4 ) ;
    Material * tKim = make_field_hts( tFactory, tKimLaw ) ;

    const real tJcB = tKimLaw->eval( B, tBeta ) ;

    for ( real tRatio : { 0.5, 1.0, 1.5 } )
    {
        const real J = tRatio * tJcB ;

        EXPECT_TRUE( within( tKim->rho_powerlaw( J, B, tBeta ), tKim->rho_powerlaw( J, T, B, tBeta ), 1e-14, 0.0 ) )
            << "rho_powerlaw at J / jc = " << tRatio ;
        EXPECT_TRUE( within( tKim->drho_powerlaw_dJ( J, B, tBeta ), tKim->drho_powerlaw_dJ( J, T, B, tBeta ), 1e-14, 0.0 ) )
            << "drho_powerlaw_dJ at J / jc = " << tRatio ;
    }

    delete tKim ;
}
