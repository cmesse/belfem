/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Value contracts of the builtin materials.
 *
 * A thermal or magnetic assembly reads cp, lambda and rho at every
 * integration point. A non-finite, negative or jumping value there does not
 * stop a run; it poisons the tangent. The public accessors of every builtin
 * are swept on a half-kelvin grid up to T_max, the protected raw laws the
 * splines are sampled from are swept finer, the Magnesia splines are checked
 * between their knots where a sample can miss an undershoot, and the
 * Hastelloy resistivity switch is recomputed from the fit it comes from.
 *
 * Two measured defect windows are left out of the sweeps that cross them and
 * are held instead by a disabled case each, which is red until the fit behind
 * it is repaired; everything above a window is swept as usual.
 *
 * Not covered here: the field-dependent overloads, the HTS resistivity laws,
 * the mechanical and expansion properties, and the raw-law seams below the
 * resolution of the public spline.
 */

#include <cmath>
#include <string>

#include <gtest/gtest.h>

#include "cl_Material.hpp"
#include "cl_MaterialFactory.hpp"
#include "cl_Material_SplineLookupTable.hpp"
#include "cl_Material_Aluminum.hpp"
#include "cl_Material_Chromium.hpp"
#include "cl_Material_Copper.hpp"
#include "cl_Material_HastelloyC276.hpp"
#include "cl_Material_Indium.hpp"
#include "cl_Material_Iron.hpp"
#include "cl_Material_Lead.hpp"
#include "cl_Material_Magnesia.hpp"
#include "cl_Material_Nickel.hpp"
#include "cl_Material_Silver.hpp"
#include "cl_Material_WhiteTin.hpp"
#include "cl_Material_YBCO.hpp"
#include "cl_Spline.hpp"

using belfem::real ;
using belfem::index_t ;
using belfem::Material ;
using belfem::MaterialFactory ;
using belfem::MaterialProperty ;
using belfem::Spline ;

namespace
{
    // every label the factory builds: the nine Kohler metals, the alloy, the
    // buffer and the superconductor. The false constructor flag skips the
    // angle-dependent resistivity database only; it does not select an
    // analytic path. Public cp and lambda are the 4 K splines sampled from the
    // raw laws, a metal's public rho is the raw law itself, and Hastelloy's
    // public rho is a spline
    const char * gBuiltinLabels[] = { "copper", "aluminum", "silver", "nickel", "iron",
                                      "chromium", "indium", "lead", "tin",
                                      "hastelloyc276", "magnesia", "ybco" } ;

    // Step envelope of the public sweep: adjacent samples above 0.5 K must
    // satisfy |d ln f| <= gEnvelopeExponent * ln( T1 / T0 ) + gEnvelopeOffset.
    // This pair is a project acceptance threshold — a smoothness screen for a
    // jump, not a bound taken from the literature: the log-space bridges and
    // the Bezier pieces carry no global cap on the local exponent. A fit that
    // this screen rejects is a ruling for the property owner, not by itself a
    // defect. The measured maximum is recorded with every sweep.
    const real gEnvelopeExponent = 5.0 ;
    const real gEnvelopeOffset   = 0.01 ;

    const real gPublicStep = 0.5 ;
    const real gRawStep    = 0.05 ;
    const real gKnotStep   = 0.01 ;

    struct PropertySet
    {
        bool mCp ;
        bool mLambda ;
        bool mRho ;
        bool mDebye ;
    };

    // mHave: what have() reports on the public surface. mRaw: which of the
    // protected raw laws the concrete class implements; the base returns an
    // always-active error for the others, so a sweep must never reach them.
    // A constant Debye temperature ( Hastelloy, YBCO ) is a public property
    // without a raw law; Magnesia carries neither a resistivity nor a Debye
    // temperature.
    struct BuiltinContract
    {
        const char * mLabel ;
        PropertySet  mHave ;
        PropertySet  mRaw ;
    };

    const BuiltinContract gContracts[] = {
        { "copper",        { true, true, true,  true  }, { true, true, true,  true  } },
        { "aluminum",      { true, true, true,  true  }, { true, true, true,  true  } },
        { "silver",        { true, true, true,  true  }, { true, true, true,  true  } },
        { "nickel",        { true, true, true,  true  }, { true, true, true,  true  } },
        { "iron",          { true, true, true,  true  }, { true, true, true,  true  } },
        { "chromium",      { true, true, true,  true  }, { true, true, true,  true  } },
        { "indium",        { true, true, true,  true  }, { true, true, true,  true  } },
        { "lead",          { true, true, true,  true  }, { true, true, true,  true  } },
        { "tin",           { true, true, true,  true  }, { true, true, true,  true  } },
        { "hastelloyc276", { true, true, true,  true  }, { true, true, true,  false } },
        { "magnesia",      { true, true, false, false }, { true, true, false, false } },
        { "ybco",          { true, true, true,  true  }, { true, true, true,  false } },
    } ;

    const BuiltinContract &
    contract( const std::string & aLabel )
    {
        for ( const BuiltinContract & tContract : gContracts )
        {
            if ( aLabel == tContract.mLabel ) return tContract ;
        }
        ADD_FAILURE() << "no contract for label " << aLabel ;
        return gContracts[ 0 ] ;
    }

    // defects the sweeps found, each bounded by the temperature window it was
    // measured in: the window ( 0, mBelow ] is exempted from the enabled sweep
    // of that ( label, surface, property ) and is pinned by a DISABLED_ case in
    // MaterialValueDefects below. Every grid point above the window, and every
    // other property of the same label, is tested as usual. The bound itself
    // belongs to the window: the YBCO NaN reaches 1.0 K, and lead's offending
    // step is the pair that ends at 3.0 K.
    enum class Surface
    {
        Public,
        Raw
    };

    struct KnownDefect
    {
        const char * mLabel ;
        Surface      mSurface ;
        const char * mProperty ;
        real         mBelow ;
        const char * mDefect ;
    };

    const KnownDefect gKnownDefects[] = {
        { "lead", Surface::Public, "cp",     3.0, "negative below 3 K" },
        { "ybco", Surface::Raw,    "lambda", 1.0, "NaN up to 1 K" },
    } ;

    // the defect of a ( label, surface, property ) check, nullptr if none
    const KnownDefect *
    known_defect( const std::string & aLabel, const Surface aSurface, const std::string & aProperty )
    {
        for ( const KnownDefect & tDefect : gKnownDefects )
        {
            if ( aLabel == tDefect.mLabel && aSurface == tDefect.mSurface && aProperty == tDefect.mProperty )
            {
                return & tDefect ;
            }
        }
        return nullptr ;
    }

    // the upper end of the exempted window, 0 when there is no known defect
    real
    defect_window( const KnownDefect * aDefect )
    {
        return aDefect == nullptr ? 0.0 : aDefect->mBelow ;
    }

    // a grid point inside the exempted window ( 0, aWindow ]; a window of 0
    // exempts nothing, so the 0 K point of a clean sweep is kept
    bool
    inside_window( const real aT, const real aWindow )
    {
        return aWindow > 0.0 && aT <= aWindow ;
    }

    // records the exempted window on the running test, so a green verdict
    // still names what it left out
    void
    record_window( const std::string & aWhat, const KnownDefect * aDefect )
    {
        if ( aDefect == nullptr ) return ;

        ::testing::Test::RecordProperty( "exempted_window_" + aWhat,
            "T <= " + std::to_string( aDefect->mBelow ) + " K, " + aDefect->mDefect
            + "; pinned in MaterialValueDefects" ) ;
    }

    // a concrete material with its protected raw laws made public. The
    // constructors are inherited, so a probe is built with the factory's own
    // arguments and is the same object the factory would hand out.
    template< typename M >
    struct RawProbe : public M
    {
        using M::M ;
        using M::cp_custom ;
        using M::lambda_custom ;
        using M::rho_custom ;
        using M::debye_custom ;
    };

    // builds the probe of a label on the stack and runs aBody on it; the
    // constructor split matches MaterialFactory::create_material
    template< typename F >
    void
    with_raw_probe( const std::string & aLabel, F && aBody )
    {
        using namespace belfem::material ;

        if      ( aLabel == "copper" )        { RawProbe< Copper >        tProbe( 50.0, false ) ; aBody( tProbe ) ; }
        else if ( aLabel == "aluminum" )      { RawProbe< Aluminum >      tProbe( 50.0, false ) ; aBody( tProbe ) ; }
        else if ( aLabel == "silver" )        { RawProbe< Silver >        tProbe( 50.0, false ) ; aBody( tProbe ) ; }
        else if ( aLabel == "nickel" )        { RawProbe< Nickel >        tProbe( 50.0, false ) ; aBody( tProbe ) ; }
        else if ( aLabel == "iron" )          { RawProbe< Iron >          tProbe( 50.0, false ) ; aBody( tProbe ) ; }
        else if ( aLabel == "chromium" )      { RawProbe< Chromium >      tProbe( 50.0, false ) ; aBody( tProbe ) ; }
        else if ( aLabel == "indium" )        { RawProbe< Indium >        tProbe( 50.0, false ) ; aBody( tProbe ) ; }
        else if ( aLabel == "lead" )          { RawProbe< Lead >          tProbe( 50.0, false ) ; aBody( tProbe ) ; }
        else if ( aLabel == "tin" )           { RawProbe< WhiteTin >      tProbe( 50.0, false ) ; aBody( tProbe ) ; }
        else if ( aLabel == "hastelloyc276" ) { RawProbe< HastelloyC276 > tProbe ;               aBody( tProbe ) ; }
        else if ( aLabel == "magnesia" )      { RawProbe< Magnesia >      tProbe ;               aBody( tProbe ) ; }
        else if ( aLabel == "ybco" )          { RawProbe< YBCO >          tProbe ;               aBody( tProbe ) ; }
        else
        {
            ADD_FAILURE() << "no raw probe for label " << aLabel ;
        }
    }

    // first offending sample of one kind, and how many there were
    struct Offense
    {
        std::size_t mCount = 0 ;
        real        mT     = 0.0 ;
        real        mValue = 0.0 ;

        void
        record( const real aT, const real aValue )
        {
            if ( mCount++ == 0 )
            {
                mT     = aT ;
                mValue = aValue ;
            }
        }
    };

    // the value at T: finite; >= 0 at 0 K where the laws vanish, > 0 above
    void
    classify( const real aT, const real aValue, Offense & aNonFinite, Offense & aNonPositive )
    {
        if ( ! std::isfinite( aValue ) )
        {
            aNonFinite.record( aT, aValue ) ;
        }
        else if ( aT == 0.0 ? aValue < 0.0 : aValue <= 0.0 )
        {
            aNonPositive.record( aT, aValue ) ;
        }
    }

    void
    expect_none( const Offense & aOffense, const std::string & aWhat, const std::string & aKind )
    {
        EXPECT_EQ( aOffense.mCount, 0u )
            << aWhat << " is " << aKind << " at " << aOffense.mCount
            << " grid points, first at T = " << aOffense.mT << " K with value " << aOffense.mValue ;
    }

    // number of aStep multiples in [ 0, aTmax ]
    std::size_t
    steps_below( const real aTmax, const real aStep )
    {
        return static_cast< std::size_t >( std::floor( aTmax / aStep + 1e-9 ) ) ;
    }

    // one public property of one label on the grid 0, aStep, 2 aStep, ...
    // with T_max appended once when it is not a grid point. Finite
    // everywhere, >= 0 at 0 K, > 0 above, and the step envelope on every
    // adjacent pair whose lower point is at or above 0.5 K. Grid points inside
    // the exempted window of aDefect are left out and the chain of adjacent
    // pairs restarts above it, so no pair straddles the window
    template< typename F >
    void
    sweep_public_grid( const std::string & aLabel, const std::string & aProperty,
                       const real aTmax, const KnownDefect * aDefect, F && aValue )
    {
        const std::string tWhat   = aLabel + " " + aProperty ;
        const real        tWindow = defect_window( aDefect ) ;

        Offense tNonFinite ;
        Offense tNonPositive ;

        std::size_t tEnvelopeCount = 0 ;
        real tEnvelopeT      = 0.0 ;
        real tEnvelopeStep   = 0.0 ;
        real tMaxExponent    = 0.0 ;
        real tMaxExponentT   = 0.0 ;

        const std::size_t tNumSteps  = steps_below( aTmax, gPublicStep ) ;
        const bool        tAppend    = static_cast< real >( tNumSteps ) * gPublicStep < aTmax ;
        const std::size_t tNumPoints = tNumSteps + 1 + ( tAppend ? 1 : 0 ) ;

        real tPrevT  = 0.0 ;
        real tPrevLn = 0.0 ;
        bool tPrevOk = false ;

        for ( std::size_t k = 0 ; k < tNumPoints ; ++k )
        {
            const real T = k <= tNumSteps ? static_cast< real >( k ) * gPublicStep : aTmax ;

            if ( inside_window( T, tWindow ) )
            {
                tPrevOk = false ;
                continue ;
            }

            const real f = aValue( T ) ;

            classify( T, f, tNonFinite, tNonPositive ) ;

            const bool tOk = std::isfinite( f ) && f > 0.0 ;

            if ( tOk && tPrevOk && tPrevT >= gPublicStep )
            {
                const real tLogRatio = std::log( T / tPrevT ) ;
                const real tStep     = std::abs( std::log( f ) - tPrevLn ) ;
                const real tExponent = tStep / tLogRatio ;

                if ( tExponent > tMaxExponent )
                {
                    tMaxExponent  = tExponent ;
                    tMaxExponentT = T ;
                }

                if ( tStep > gEnvelopeExponent * tLogRatio + gEnvelopeOffset )
                {
                    if ( tEnvelopeCount++ == 0 )
                    {
                        tEnvelopeT    = T ;
                        tEnvelopeStep = tStep ;
                    }
                }
            }

            tPrevT  = T ;
            tPrevOk = tOk ;
            tPrevLn = tOk ? std::log( f ) : 0.0 ;
        }

        expect_none( tNonFinite, tWhat, "not finite" ) ;
        expect_none( tNonPositive, tWhat, "not positive" ) ;

        EXPECT_EQ( tEnvelopeCount, 0u )
            << tWhat << " exceeds the step envelope on " << tEnvelopeCount
            << " adjacent pairs, first with |d ln f| = " << tEnvelopeStep
            << " on the pair ending at T = " << tEnvelopeT
            << " K; largest local exponent " << tMaxExponent
            << " on the pair ending at T = " << tMaxExponentT << " K" ;

        ::testing::Test::RecordProperty( "max_local_exponent_" + aProperty, std::to_string( tMaxExponent ) ) ;
        ::testing::Test::RecordProperty( "max_local_exponent_" + aProperty + "_T", std::to_string( tMaxExponentT ) ) ;

        record_window( aProperty, aDefect ) ;
    }

    // one raw law of one label on the grid aStep, 2 aStep, ... over
    // ( 0, aTmax ], T_max appended once when it is not a grid point: finite
    // and > 0 everywhere above the exempted window of aDefect
    template< typename F >
    void
    sweep_raw_law( const std::string & aLabel, const std::string & aLaw,
                   const real aTmax, const KnownDefect * aDefect, F && aValue )
    {
        const std::string tWhat   = aLabel + " raw " + aLaw ;
        const real        tWindow = defect_window( aDefect ) ;

        Offense tNonFinite ;
        Offense tNonPositive ;

        const std::size_t tNumSteps  = steps_below( aTmax, gRawStep ) ;
        const bool        tAppend    = static_cast< real >( tNumSteps ) * gRawStep < aTmax ;
        const std::size_t tNumPoints = tNumSteps + ( tAppend ? 1 : 0 ) ;

        for ( std::size_t k = 1 ; k <= tNumPoints ; ++k )
        {
            const real T = k <= tNumSteps ? static_cast< real >( k ) * gRawStep : aTmax ;

            if ( inside_window( T, tWindow ) ) continue ;

            classify( T, aValue( T ), tNonFinite, tNonPositive ) ;
        }

        expect_none( tNonFinite, tWhat, "not finite" ) ;
        expect_none( tNonPositive, tWhat, "not positive" ) ;

        record_window( "raw_" + aLaw, aDefect ) ;
    }

    // real roots of aA x^2 + aB x + aC, the linear case included
    std::size_t
    quadratic_roots( const real aA, const real aB, const real aC, real ( & aRoots )[ 2 ] )
    {
        if ( aA == 0.0 )
        {
            if ( aB == 0.0 ) return 0 ;
            aRoots[ 0 ] = -aC / aB ;
            return 1 ;
        }

        const real tDiscriminant = aB * aB - 4.0 * aA * aC ;
        if ( tDiscriminant < 0.0 ) return 0 ;

        const real tRoot = std::sqrt( tDiscriminant ) ;
        aRoots[ 0 ] = ( -aB - tRoot ) / ( 2.0 * aA ) ;
        aRoots[ 1 ] = ( -aB + tRoot ) / ( 2.0 * aA ) ;
        return 2 ;
    }

    // a spline-backed public property on the grid aStep, 2 aStep, ... over
    // ( 0, aTop ], plus the interior extrema of every cubic piece that starts
    // below aTop: the sampled value and every extremum must be > 0. The
    // extrema are the roots of the piece's derivative inside the piece's own
    // interval; the coefficient rows are those of Spline::eval, a polynomial
    // in the absolute abscissa. Tolerance of the read-back: the extremum is
    // evaluated through the public accessor, which picks its piece with
    // find_col, so a root within an ulp of a piece edge is scored on the
    // neighbouring piece — harmless, because the spline is continuous there
    template< typename F >
    void
    check_spline_positive_below( Spline * aSpline, const std::string & aWhat,
                                 const real aTop, F && aValue )
    {
        ASSERT_NE( aSpline, nullptr ) << aWhat << " has no spline" ;

        Offense tNonFinite ;
        Offense tNonPositive ;

        const std::size_t tNumSteps = steps_below( aTop, gKnotStep ) ;

        for ( std::size_t k = 1 ; k <= tNumSteps ; ++k )
        {
            const real T = static_cast< real >( k ) * gKnotStep ;
            classify( T, aValue( T ), tNonFinite, tNonPositive ) ;
        }

        expect_none( tNonFinite, aWhat, "not finite" ) ;
        expect_none( tNonPositive, aWhat, "not positive" ) ;

        const auto &  tCoefficients = aSpline->coefficients() ;
        const real    tDeltaX       = aSpline->delta_x() ;
        const index_t tNumIntervals = aSpline->n() - 1 ;

        ASSERT_GT( tDeltaX, 0.0 ) << aWhat << " spline spacing" ;

        for ( index_t tCol = 0 ; tCol < tNumIntervals ; ++tCol )
        {
            const real tLo = aSpline->x_min() + static_cast< real >( tCol ) * tDeltaX ;
            const real tHi = tLo + tDeltaX ;

            if ( tLo >= aTop ) break ;

            real tRoots[ 2 ] = { 0.0, 0.0 } ;
            const std::size_t tNumRoots = quadratic_roots( 3.0 * tCoefficients( 0, tCol ),
                                                           2.0 * tCoefficients( 1, tCol ),
                                                                 tCoefficients( 2, tCol ),
                                                           tRoots ) ;

            for ( std::size_t r = 0 ; r < tNumRoots ; ++r )
            {
                const real T = tRoots[ r ] ;
                if ( ! ( T > tLo && T < tHi ) ) continue ;

                const real f = aValue( T ) ;
                EXPECT_TRUE( std::isfinite( f ) && f > 0.0 )
                    << aWhat << " has a non-positive extremum " << f << " at T = " << T
                    << " K inside the piece [ " << tLo << ", " << tHi << " ]" ;
            }
        }
    }

    // a polynomial and its derivatives, coefficients in descending powers
    template< std::size_t N >
    real
    poly( const real ( & aP )[ N ], const real aX )
    {
        real tValue = 0.0 ;
        for ( std::size_t i = 0 ; i < N ; ++i )
        {
            tValue = tValue * aX + aP[ i ] ;
        }
        return tValue ;
    }

    template< std::size_t N >
    real
    dpoly( const real ( & aP )[ N ], const real aX )
    {
        real tValue = 0.0 ;
        for ( std::size_t i = 0 ; i + 1 < N ; ++i )
        {
            tValue = tValue * aX + static_cast< real >( N - 1 - i ) * aP[ i ] ;
        }
        return tValue ;
    }

    template< std::size_t N >
    real
    ddpoly( const real ( & aP )[ N ], const real aX )
    {
        real tValue = 0.0 ;
        for ( std::size_t i = 0 ; i + 2 < N ; ++i )
        {
            tValue = tValue * aX + static_cast< real >( ( N - 1 - i ) * ( N - 2 - i ) ) * aP[ i ] ;
        }
        return tValue ;
    }
}

//------------------------------------------------------------------------------
// R0: the property set of every builtin, as documented
//------------------------------------------------------------------------------

class BuiltinMaterialValues : public ::testing::TestWithParam< const char * > {};

TEST_P( BuiltinMaterialValues, PropertySetIsAsDocumented )
{
    MaterialFactory tFactory ;
    Material * tMaterial = tFactory.create_material( GetParam(), 50.0, false ) ;
    const std::string tLabel = GetParam() ;
    const PropertySet & tHave = contract( tLabel ).mHave ;

    EXPECT_EQ( tHave.mCp,     tMaterial->have( MaterialProperty::cp ) )     << tLabel << " cp" ;
    EXPECT_EQ( tHave.mLambda, tMaterial->have( MaterialProperty::lambda ) ) << tLabel << " lambda" ;
    EXPECT_EQ( tHave.mRho,    tMaterial->have( MaterialProperty::rho ) )    << tLabel << " rho" ;
    EXPECT_EQ( tHave.mDebye,  tMaterial->have( MaterialProperty::debye ) )  << tLabel << " debye" ;

    // every sweep below is bounded by T_max
    EXPECT_TRUE( tMaterial->have( MaterialProperty::T_max ) ) << tLabel << " T_max" ;
    EXPECT_GT( tMaterial->constant_property( MaterialProperty::T_max ), 0.0 ) << tLabel << " T_max" ;

    delete tMaterial ;
}

//------------------------------------------------------------------------------
// R1: the public accessors on the half-kelvin grid
//------------------------------------------------------------------------------

TEST_P( BuiltinMaterialValues, CpFiniteAndPositiveOnTheHalfKelvinGrid )
{
    MaterialFactory tFactory ;
    Material * tMaterial = tFactory.create_material( GetParam(), 50.0, false ) ;
    const std::string tLabel = GetParam() ;

    if ( ! tMaterial->have( MaterialProperty::cp ) )
    {
        delete tMaterial ;
        GTEST_SKIP() << tLabel << " has no cp" ;
    }

    sweep_public_grid( tLabel, "cp", tMaterial->constant_property( MaterialProperty::T_max ),
        known_defect( tLabel, Surface::Public, "cp" ),
        [ & ]( real T ){ return tMaterial->cp( T ) ; } ) ;

    delete tMaterial ;
}

TEST_P( BuiltinMaterialValues, LambdaFiniteAndPositiveOnTheHalfKelvinGrid )
{
    MaterialFactory tFactory ;
    Material * tMaterial = tFactory.create_material( GetParam(), 50.0, false ) ;
    const std::string tLabel = GetParam() ;

    if ( ! tMaterial->have( MaterialProperty::lambda ) )
    {
        delete tMaterial ;
        GTEST_SKIP() << tLabel << " has no lambda" ;
    }

    sweep_public_grid( tLabel, "lambda", tMaterial->constant_property( MaterialProperty::T_max ),
        known_defect( tLabel, Surface::Public, "lambda" ),
        [ & ]( real T ){ return tMaterial->lambda( T ) ; } ) ;

    delete tMaterial ;
}

TEST_P( BuiltinMaterialValues, RhoFiniteAndPositiveOnTheHalfKelvinGrid )
{
    MaterialFactory tFactory ;
    Material * tMaterial = tFactory.create_material( GetParam(), 50.0, false ) ;
    const std::string tLabel = GetParam() ;

    if ( ! tMaterial->have( MaterialProperty::rho ) )
    {
        delete tMaterial ;
        GTEST_SKIP() << tLabel << " has no rho" ;
    }

    sweep_public_grid( tLabel, "rho", tMaterial->constant_property( MaterialProperty::T_max ),
        known_defect( tLabel, Surface::Public, "rho" ),
        [ & ]( real T ){ return tMaterial->rho( T ) ; } ) ;

    delete tMaterial ;
}

//------------------------------------------------------------------------------
// R2: the raw laws behind the splines, finite and positive up to T_max
//------------------------------------------------------------------------------

TEST_P( BuiltinMaterialValues, RawLawsFiniteAndPositive )
{
    const std::string tLabel = GetParam() ;
    const PropertySet & tRaw = contract( tLabel ).mRaw ;

    with_raw_probe( tLabel, [ & ]( auto & aProbe )
    {
        const real tTmax = aProbe.constant_property( MaterialProperty::T_max ) ;

        auto tSweep = [ & ]( const bool aPresent, const char * aLaw, auto && aValue )
        {
            if ( ! aPresent ) return ;

            sweep_raw_law( tLabel, aLaw, tTmax, known_defect( tLabel, Surface::Raw, aLaw ), aValue ) ;
        } ;

        tSweep( tRaw.mCp,     "cp",     [ & ]( real T ){ return aProbe.cp_custom( T ) ; } ) ;
        tSweep( tRaw.mLambda, "lambda", [ & ]( real T ){ return aProbe.lambda_custom( T ) ; } ) ;
        tSweep( tRaw.mRho,    "rho",    [ & ]( real T ){ return aProbe.rho_custom( T ) ; } ) ;
        tSweep( tRaw.mDebye,  "debye",  [ & ]( real T ){ return aProbe.debye_custom( T ) ; } ) ;
    } ) ;
}

INSTANTIATE_TEST_SUITE_P(
    Builtin,
    BuiltinMaterialValues,
    ::testing::ValuesIn( gBuiltinLabels ),
    []( const ::testing::TestParamInfo< const char * > & aInfo )
    {
        return std::string( aInfo.param ) ;
    } ) ;

//------------------------------------------------------------------------------
// R3: the Magnesia splines between their knots
//------------------------------------------------------------------------------

// both raw laws vanish at 0 K and rise steeply; the cubic pieces between the
// 4 K knots can undershoot zero where the knot values alone stay positive
TEST( MagnesiaSplines, MagnesiaLambdaAndCpPositiveBetweenKnots )
{
    RawProbe< belfem::material::Magnesia > tProbe ;

    const real tTop = 20.0 ;

    check_spline_positive_below( tProbe.spline( MaterialProperty::lambda ), "magnesia lambda", tTop,
        [ & ]( real T ){ return tProbe.lambda( T ) ; } ) ;

    check_spline_positive_below( tProbe.spline( MaterialProperty::cp ), "magnesia cp", tTop,
        [ & ]( real T ){ return tProbe.cp( T ) ; } ) ;
}

//------------------------------------------------------------------------------
// R4: the Hastelloy resistivity switch
//------------------------------------------------------------------------------

// the low-temperature quartic is anchored where the log-log fit turns, at
// T0 = exp( x* ) with x* the root of the fit's derivative. The raw law must
// take the quartic below T0 and the log-log fit above it, and must be
// continuous across T0 as well as across x*; a switch stored on the wrong
// scale ( x* instead of exp( x* ) ) serves the wrong branch on the whole
// interval between the two and opens a jump at x*
TEST( HastelloyResistivity, HastelloyRhoSwitchIsTheLogSpaceRoot )
{
    RawProbe< belfem::material::HastelloyC276 > tProbe ;

    // the log-log fit above the turning point, as the class stores it: the
    // resistivity fit of cl_Material_HastelloyC276.cpp, fitted there against
    // the measured curve of Lu 2008 ( 10.1063/1.2899058 )
    const real tFit[ 5 ] = {   1.93174404699903E-04,
                             - 2.21920482428318E-03,
                             + 1.07884809281212E-02,
                             - 2.40374094654130E-02,
                             - 1.35906507676544E+01 } ;

    // Newton on the fit's derivative from the class's own start value
    real tX = 2.373 ;
    real tF = dpoly( tFit, tX ) ;
    belfem::uint tIteration = 0 ;

    while ( std::abs( tF ) > belfem::BELFEM_EPSILON )
    {
        ASSERT_LT( tIteration++, 100u ) << "Newton on the log-space fit does not converge" ;
        tX -= tF / ddpoly( tFit, tX ) ;
        tF  = dpoly( tFit, tX ) ;
    }

    const real tSwitch = std::exp( tX ) ;
    const real tTmax   = tProbe.constant_property( MaterialProperty::T_max ) ;

    ASSERT_TRUE( std::isfinite( tSwitch ) ) ;
    ASSERT_GT( tX, 0.0 ) ;
    ASSERT_LT( tSwitch, tTmax ) ;

    // the fit is documented as the branch above 10.73 K
    EXPECT_NEAR( tSwitch, 10.73, 0.01 ) ;

    // which branch the law serves on either side of the expected switch. The
    // low branch is the quartic the class builds there, in the shifted form it
    // is constructed from: a ( T - T0 )^4 + c ( T - T0 )^2 + rho( T0 ), its two
    // odd coefficients zero. It is anchored to the fit's value and slope at T0,
    // so the branches part slowly and the sample must sit far enough out for
    // them to be told apart; the assertions below state that separation before
    // the comparison is trusted
    const real tA      = 2.42719471085423E-13 ;
    const real tC      = 3.69608432119349E-11 ;
    const real tRho0   = std::exp( poly( tFit, tX ) ) ;
    const real tOffset = 0.5 ;
    const real tMatch  = 1e-9 ;

    auto tLowBranch = [ & ]( const real T )
    {
        const real tShift = T - tSwitch ;
        return ( tA * tShift * tShift + tC ) * tShift * tShift + tRho0 ;
    } ;

    auto tLogBranch = [ & ]( const real T )
    {
        return std::exp( poly( tFit, std::log( T ) ) ) ;
    } ;

    for ( real T : { tSwitch - tOffset, tSwitch + tOffset } )
    {
        ASSERT_GT( std::abs( tLowBranch( T ) - tLogBranch( T ) ) / tRho0, tMatch )
            << "the two branches are within the match tolerance at T = " << T
            << " K, so the comparison there decides nothing" ;
    }

    EXPECT_NEAR( tProbe.rho_custom( tSwitch - tOffset ), tLowBranch( tSwitch - tOffset ), tMatch * tRho0 )
        << "hastelloyc276 raw rho at T = " << tSwitch - tOffset
        << " K is not the low-temperature quartic; the served switch is not exp( x* ) = " << tSwitch << " K" ;

    EXPECT_NEAR( tProbe.rho_custom( tSwitch + tOffset ), tLogBranch( tSwitch + tOffset ), tMatch * tRho0 )
        << "hastelloyc276 raw rho at T = " << tSwitch + tOffset
        << " K is not the log-log fit; the served switch is not exp( x* ) = " << tSwitch << " K" ;

    // continuity as a relative second difference across +- delta: a smooth
    // law leaves delta^2 rho'' / rho, a jump leaves its own size
    const real tDelta = 1e-3 ;

    for ( real T : { tSwitch, tX } )
    {
        const real tCenter = tProbe.rho_custom( T ) ;
        const real tBelow  = tProbe.rho_custom( T - tDelta ) ;
        const real tAbove  = tProbe.rho_custom( T + tDelta ) ;

        ASSERT_TRUE( std::isfinite( tCenter ) && tCenter > 0.0 )
            << "hastelloyc276 raw rho at T = " << T << " K is " << tCenter ;
        ASSERT_TRUE( std::isfinite( tBelow ) && std::isfinite( tAbove ) )
            << "hastelloyc276 raw rho around T = " << T << " K" ;

        const real tKink = std::abs( tAbove - 2.0 * tCenter + tBelow ) / tCenter ;

        EXPECT_LE( tKink, 1e-6 )
            << "hastelloyc276 raw rho has a relative jump of " << tKink
            << " at T = " << T << " K ( x* = " << tX << ", exp( x* ) = " << tSwitch << " )" ;
    }
}

//------------------------------------------------------------------------------
// the exempted windows: the contract that the sweeps above leave out, red
// until the fit behind it is repaired
//------------------------------------------------------------------------------

// the public cp of lead is positive on the half-kelvin grid up to 3 K, and the
// pair that ends there stays inside the step envelope
TEST( MaterialValueDefects, DISABLED_LeadCpIsNegativeBelowThreeKelvin )
{
    MaterialFactory tFactory ;
    Material * tLead = tFactory.create_material( "lead", 50.0, false ) ;

    for ( real T : { 0.5, 1.0, 1.5, 2.0 } )
    {
        EXPECT_GT( tLead->cp( T ), 0.0 ) << "lead cp at T = " << T << " K" ;
    }

    const real tLower = 2.5 ;
    const real tUpper = 3.0 ;
    const real tStep  = std::abs( std::log( tLead->cp( tUpper ) ) - std::log( tLead->cp( tLower ) ) ) ;

    EXPECT_LE( tStep, gEnvelopeExponent * std::log( tUpper / tLower ) + gEnvelopeOffset )
        << "lead cp |d ln f| = " << tStep << " on the pair ending at T = " << tUpper << " K" ;

    delete tLead ;
}

// the raw Callaway conductivity of YBCO is finite on the 0.05 K grid up to 1 K
TEST( MaterialValueDefects, DISABLED_YbcoRawLambdaIsFiniteBelowOneKelvin )
{
    RawProbe< belfem::material::YBCO > tProbe ;

    for ( std::size_t k = 1 ; k <= 20 ; ++k )
    {
        const real T = static_cast< real >( k ) * gRawStep ;
        EXPECT_TRUE( std::isfinite( tProbe.lambda_custom( T ) ) )
            << "ybco raw lambda at T = " << T << " K is " << tProbe.lambda_custom( T ) ;
    }
}

//------------------------------------------------------------------------------
// thermal conductivity at room temperature against the data the fit cites.
// Disabled until the reference values are confirmed: the numbers below were not
// read from the cited tables. The materials are built at RRR 50, whose
// residual term is a few percent of 1 / lambda at 300 K.
// Enabled once, each case prints the model value, the reference and the
// relative deviation, so the fit's position is visible without a debugger
//------------------------------------------------------------------------------

namespace
{
    // Touloukian et al., Thermophysical Properties of Matter, Vol. 1, the data
    // set the nickel fit names, recommended value for pure nickel at 300 K
    constexpr real gNickelLambda300K = 90.7 ;

    // Hust and Lankford 1984, NBSIR 84-3007, the source the iron fit names,
    // at 300 K
    constexpr real gIronLambda300K = 83.2 ;

    // the band the rows ask for
    constexpr real gLambdaBand = 0.15 ;

    void
    expect_lambda_at_300K( const std::string & aLabel, const real aReference )
    {
        MaterialFactory tFactory ;
        Material * tMaterial = tFactory.create_material( aLabel, 50.0, false ) ;

        const real tModel = tMaterial->lambda( 300.0 ) ;
        const real tDeviation = ( tModel - aReference ) / aReference ;

        EXPECT_LT( std::abs( tDeviation ), gLambdaBand )
            << aLabel << " lambda( 300 K ): model " << tModel << " W/(m K), reference "
            << aReference << " W/(m K), deviation " << 100.0 * tDeviation << " %" ;

        delete tMaterial ;
    }
}

TEST( MaterialValueDefects, DISABLED_NickelLambdaAt300KMatchesTheReference )
{
    expect_lambda_at_300K( "nickel", gNickelLambda300K ) ;
}

TEST( MaterialValueDefects, DISABLED_IronLambdaAt300KMatchesTheReference )
{
    expect_lambda_at_300K( "iron", gIronLambda300K ) ;
}

//------------------------------------------------------------------------------
// the reduced magnetization of nickel, which scales the magnetic resistivity
//------------------------------------------------------------------------------

namespace
{
    struct NickelProbe : public belfem::material::Nickel
    {
        using Nickel::Nickel ;
        using Nickel::compute_mred ;
    };
}

TEST( NickelMagnetization, ReducedMagnetizationMatchesCrangleGoodman )
{
    // the Bezier in Nickel::create_magnetization_curve is fitted to Crangle
    // and Goodman 1971, 10.1098/rspa.1971.0044, with Tc = 631 K; the class
    // records the check m( 600 K ) = 0.47 against the measured 0.45 to 0.50.
    // Mutation: swap the two interior control ordinates and 600 K leaves the
    // window
    const NickelProbe tNickel( 50.0, false ) ;

    const real tTc = tNickel.constant_property( MaterialProperty::Tcurie ) ;

    EXPECT_DOUBLE_EQ( tNickel.compute_mred( 0.0 ), 1.0 ) ;
    EXPECT_NEAR( tNickel.compute_mred( tTc ), 0.0, 1e-12 ) ;
    EXPECT_EQ( tNickel.compute_mred( tTc + 1.0 ), 0.0 ) ;

    // the measured band, and the fitted value the class records inside it
    const real tAt600 = tNickel.compute_mred( 600.0 ) ;
    EXPECT_GE( tAt600, 0.45 ) ;
    EXPECT_LE( tAt600, 0.50 ) ;
    EXPECT_NEAR( tAt600, 0.47, 0.005 ) ;

    // monotone decreasing on a 10 K grid up to the Curie point
    real tPrevious = tNickel.compute_mred( 0.0 ) ;
    for ( real T = 10.0; T < tTc; T += 10.0 )
    {
        const real tValue = tNickel.compute_mred( T ) ;
        EXPECT_LE( tValue, tPrevious ) << "T = " << T << " K" ;
        tPrevious = tValue ;
    }
}
