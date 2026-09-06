/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include "cl_Material_SplineLookupTable.hpp"
namespace belfem
{
    namespace material
    {
        SplineLookupTable::SplineLookupTable( const MaterialType aType, const bool aIsIsotropic )
            : Material( aType, aIsIsotropic )
        {
            mSplines.set_size( gNumNonConstantMaterialProperties, nullptr );
        }

        SplineLookupTable::~SplineLookupTable()
        {
            for ( Spline * tSpline : mSplines )
            {
                if ( tSpline != nullptr ) delete tSpline ;
            }
        }


        void
        SplineLookupTable::set_spline( const MaterialProperty aProperty, Spline * aSpline )
        {
            index_t tIndex = static_cast< size_t >( aProperty ) ;
            mConstantProperties( tIndex ) = BELFEM_QUIET_NAN ;
            mPropertyDependencies( tIndex )->reset() ;

            // replace the stored spline; with the default nullptr argument,
            // the existing spline is kept and only the dispatch is refreshed
            if ( aSpline != nullptr && mSplines( tIndex ) != aSpline )
            {
                this->reset_spline( aProperty );
                mSplines( tIndex ) = aSpline ;
            }
            this->set_have( aProperty );
            this->set_dependency(  aProperty, MaterialDependency::T );

            switch ( aProperty )
            {
               case MaterialProperty::cp :
               {
                   mFunctionCp      = & Material::cp_spline ;
                   mFunctiondCpdT   = & Material::dcpdT_spline ;
                   mFunctiond2CpdT2 = & Material::d2cpdT2_spline ;
                   break ;
               }
               case MaterialProperty::lambda :
               {
                   mFunctionLambda    = & Material::lambda_spline ;
                   mFunctiondLambdadT = & Material::dlambdadT_spline ;
                   break ;
               }
               case MaterialProperty::rho :
               {
                   mFunctionRho    = & Material::rho_spline ;
                   mFunctiondRhodT = & Material::drhodT_spline ;
                   break;
               }
               case MaterialProperty::mu :
               {
                   BELFEM_ERROR( false, "Material::set_spline must not be called for mu-property");
                   break;
               }
               case MaterialProperty::E :
               {
                   mFunctionE = & Material::E_spline ;
                   break ;
               }
               case MaterialProperty::nu :
               {
                   mFunctionNu = & Material::nu_spline ;
                   break;
               }
               case MaterialProperty::alpha :
               {
                   mFunctionAlpha = & Material::alpha_spline ;
                   break;
               }
               case MaterialProperty::Rp02 :
               {
                   mFunctionRp02 = & Material::Rp02_spline ;
                   break;
               }
               case MaterialProperty::rho_i :
               {
                   mFunctionRhoI = & Material::rho_i_spline ;
                   break;
               }
               case MaterialProperty::debye :
               {
                   mFunctionDebye = & Material::debye_spline ;
                   break;
               }
               default:
               {
                   BELFEM_ERROR( false, "Unknown material property ");
               }
           }
        }



        void
        SplineLookupTable::reset_spline( const MaterialProperty aProperty )
        {
            index_t tIndex = static_cast< size_t >( aProperty ) ;
            if ( mSplines( tIndex ) != nullptr ) delete mSplines( tIndex ) ;
            mSplines( tIndex ) = nullptr ;
        }

        void
        SplineLookupTable::create_spline( real (Material::*aFunction)(const real aT) const,
            const MaterialProperty aProperty,
            const uint aStartBC,
            const uint aEndBC,
            const real adYdX0,
            const real adYdX1 )
        {
            // we always make steps in 4-K width, we need to be equidistant and with 4 K, Helium is captured properly
            const real dT = 4.0 ;

            size_t tNumPoints = std::ceil( this->constant_property( MaterialProperty::T_max ) / dT ) + 1 ;

            Vector< real > tX( tNumPoints );
            Vector< real > tY( tNumPoints );

            real T = 0.0 ;
            for ( size_t k=0; k<tNumPoints; ++k )
            {
                tX( k ) = T ;
                tY( k ) = ( this->*aFunction )( T ) ;
                T += dT ;

            }

            spline::SplineBC tStartBC = static_cast< spline::SplineBC >( aStartBC );
            spline::SplineBC tEndBC = static_cast< spline::SplineBC >( aEndBC );
            SpMatrix tA ;
            spline::create_helpmatrix( tNumPoints, dT, tA, tStartBC, tEndBC );
            this->set_spline(  aProperty, new Spline( tX, tY, tA, tStartBC, tEndBC, adYdX0, adYdX1 ) );
        }
        void
        SplineLookupTable::create_cryo_expansion(
        const Bezier * aThermalExpansion,
        Vector< real > & aThermalExpansionCryo,
        const real aTSwitch )
        {
            // the cryogenic branch is fitted against cp, so create_cp() must have
            // run first; the spline below then samples the composite alpha_custom
            this->create_low_temperature_alpha( aThermalExpansion, aThermalExpansionCryo, aTSwitch );
            this->finish_cryo_expansion( aThermalExpansionCryo );
        }

        void
        SplineLookupTable::create_cryo_expansion(
        const real alpha, const real dalphadT, const real d2alphadT2,
        Vector< real > & aThermalExpansionCryo )
        {
            this->create_low_temperature_alpha( alpha, dalphadT, d2alphadT2, aThermalExpansionCryo );
            this->finish_cryo_expansion( aThermalExpansionCryo );
        }

        void
        SplineLookupTable::create_cryo_expansion_anchored(
        const real alpha,
        Vector< real > & aThermalExpansionCryo,
        const real aTSwitch )
        {
            BELFEM_ERROR( this->have( MaterialProperty::E ) && this->have( MaterialProperty::nu ),
                "%s: the anchored cryogenic alpha branch needs E and nu ( for K )",
                this->label().c_str() );

            real T = this->set_alpha_switch_temperature( aTSwitch );

            // The Grueneisen relation C = alpha / cp = C* K* / K( T ) supplies the
            // value and the first two derivatives of alpha AT THE SPLIT only:
            //     alpha'  = C ( cp' - cp K'/K )
            //     alpha'' = C ( cp'' - 2 cp' K'/K + cp ( 2 (K'/K)^2 - K''/K ) )
            // Below the split ln( C ) is the cubic ( or quadratic ) that
            // create_low_temperature_alpha() fits to that 2-jet with p'( 0 ) = 0,
            // not 1/K( T ). K' and K'' by central differences on the material's
            // own K( T ), so the split must sit at least h above 0 K.
            real h  = 1.0 ;

            BELFEM_ERROR( T > h,
                "%s: the anchored alpha branch needs a split above %g K, got %g K",
                this->label().c_str(), ( double ) h, ( double ) T );

            real K0 = this->K( T );
            real Kp = ( this->K( T + h ) - this->K( T - h ) ) / ( 2. * h );
            real Kpp = ( this->K( T + h ) - 2. * K0 + this->K( T - h ) ) / ( h * h );
            real r1 = Kp / K0 ;
            real r2 = Kpp / K0 ;

            real cp    = this->cp( T );
            real cp1   = this->dcpdT( T );
            real cp2   = this->d2cpdT2( T );

            BELFEM_ERROR( std::isfinite( K0 ) && K0 > 0.0 && std::isfinite( cp ) && cp > 0.0
                          && std::isfinite( alpha ) && alpha > 0.0,
                "%s: anchored alpha branch at %g K needs K > 0, cp > 0 and alpha > 0 "
                "( K = %g, cp = %g, alpha = %g )",
                this->label().c_str(), ( double ) T, ( double ) K0, ( double ) cp, ( double ) alpha );

            real C     = alpha / cp ;

            real dalphadT   = C * ( cp1 - cp * r1 );
            real d2alphadT2 = C * ( cp2 - 2. * cp1 * r1 + cp * ( 2. * r1 * r1 - r2 ) );

            // The dln(C)/dT > 0 guard is skipped: it rejects fitted expansion
            // curves whose C falls with T, and here dln(C)/dT = -K'/K comes from
            // the material's own moduli. That is load-bearing, not cosmetic -
            // Magnesia's K rises imperceptibly at the split ( -K'/K = -5e-5 ), and
            // the strict check would refuse it; the resulting drift of C over the
            // branch is of order 1 %.
            this->create_low_temperature_alpha( alpha, dalphadT, d2alphadT2, aThermalExpansionCryo, false );
            this->finish_cryo_expansion( aThermalExpansionCryo );
        }

        void
        SplineLookupTable::finish_cryo_expansion( const Vector< real > & aThermalExpansionCryo )
        {
            // alpha = exp( p(T) ) * cp(T) with p'(0) = 0 and cp(0) = 0, hence
            // dalpha/dT( 0 ) = exp( p(0) ) * dcp/dT( 0 ): the Sommerfeld
            // coefficient scaled by C supplies the boundary condition, not the
            // Bezier, whose slope at the origin is a fitting artifact.
            real dAlphadT = std::exp( polyval( aThermalExpansionCryo, 0.0 ) ) * this->dcpdT( 0.0 );

            this->create_spline( MaterialProperty::alpha, dAlphadT );
            this->spline( MaterialProperty::alpha )->create_integral( this->constant_property( MaterialProperty::T_ref_density ), 0.0 );
        }

    real
    SplineLookupTable::set_alpha_switch_temperature( const real aTSwitch )
    {
        if ( std::isnan( aTSwitch ) )
        {
            // default split: 0.618 * theta, capped at gTAlphaSwitchMax
            BELFEM_ERROR( this->have( MaterialProperty::debye0K),
                "%s doesn't have a low temperature debye temperature set", mLabel.c_str() );

            mTAlphaSwitch = std::min(
                this->constant_property( MaterialProperty::debye0K ) * ( constant::phi - 1. ),
                gTAlphaSwitchMax );
        }
        else
        {
            // a material whose expansion curve flattens early ( a plateau, an
            // anomaly ) passes its own, lower value; 0 K means "no branch yet"
            mTAlphaSwitch = aTSwitch ;
        }

        // 0 K is legal and means "no cryogenic branch yet" ( a provisional
        // spline on the plain curve, replaced once cp exists )
        BELFEM_ERROR( mTAlphaSwitch >= 0.0 && mTAlphaSwitch <= gTAlphaSwitchMax,
            "%s: alpha split temperature %g K is outside [ 0, %g ]",
            mLabel.c_str(), ( double ) mTAlphaSwitch, ( double ) gTAlphaSwitchMax );

        return mTAlphaSwitch ;
    }

//------------------------------------------------------------------------------

    void
    SplineLookupTable::create_low_temperature_alpha( const Bezier * aBezier, Vector< real > & aPoly, const real aTSwitch )
    {
        real T = this->set_alpha_switch_temperature( aTSwitch );

        // the Bezier stores dL/L, the relative length is L = 1 + dL/L
        real l      = 1. + aBezier->y( T );
        real dldT   = aBezier->dydx( T );
        real dl2dT2 = aBezier->d2ydx2( T );
        real dl3dT3 = aBezier->d3ydx3( T );

        // with alpha = L'/L follows
        //     dalpha/dT   = L''/L - alpha^2
        //     d2alpha/dT2 = L'''/L - 3*alpha*dalpha/dT - alpha^3
        real alpha      = dldT / l ;
        real dalphadT   = dl2dT2 / l - alpha * alpha ;
        real d2alphadT2 = dl3dT3 / l - 3. * alpha * dalphadT - alpha * alpha * alpha ;

        this->create_low_temperature_alpha( alpha, dalphadT, d2alphadT2, aPoly );
    }

//------------------------------------------------------------------------------

    void
    SplineLookupTable::create_low_temperature_alpha(
            const real alpha, const real dalphadT, const real d2alphadT2, Vector< real > & aPoly,
            const bool aCheckGuard )
    {
        // the branch is fitted against cp and its first two derivatives, so cp
        // must already exist AND already carry its final routing: whatever
        // create_cp() / create_debye() leave behind is what alpha_custom will
        // call at run time, and the fit has to be made against the same curve
        BELFEM_ERROR( this->have( MaterialProperty::cp ),
            "%s: cp must be created before the low temperature alpha branch is fitted",
            mLabel.c_str() );

        // the caller evaluated alpha and its derivatives at the split temperature
        BELFEM_ERROR( ! std::isnan( mTAlphaSwitch ) && mTAlphaSwitch > 0.0,
            "%s: set_alpha_switch_temperature() must run before the cryogenic alpha fit",
            mLabel.c_str() );
        real T = mTAlphaSwitch ;

        real cp      = this->cp( T );
        real dcpdT   = this->dcpdT( T );
        real d2cpdT2 = this->d2cpdT2( T );

        // C = alpha / cp, hence h = ln( C ) = ln( alpha ) - ln( cp )
        real f = std::log( alpha );
        real g = std::log( cp );
        real h = f - g ;

        real dfdT = dalphadT / alpha ;
        real dgdT = dcpdT / cp ;
        real dhdT = dfdT - dgdT ;

        // C must not decrease with temperature. Where it does, the fitted
        // expansion curve has collapsed relative to the heat capacity, and
        // extrapolating it downwards inflates alpha at cryogenic temperatures.
        // The anchored branch derives the slope from K( T ) itself, so there
        // is no fitted curve to reject and the check is skipped there.
        BELFEM_ERROR( ( ! aCheckGuard ) || dhdT > BELFEM_EPSILON,
            "dln(C)/dT = %g is not positive at T = %g K for %s: the thermal expansion "
            "curve is inconsistent with cp at the split temperature",
            ( double ) dhdT, ( double ) T, mLabel.c_str() );

        real d2fdT2 = ( alpha * d2alphadT2 - dalphadT * dalphadT ) / ( alpha * alpha );
        real d2gdT2 = ( cp * d2cpdT2 - dcpdT * dcpdT ) / ( cp * cp );
        real d2hdT2 = d2fdT2 - d2gdT2 ;

        // make sure that there is no turning point or extremum
        if ( T * d2hdT2 < 2. * dhdT + BELFEM_EPSILON )
        {
            real d = -6.*T*T ;

            aPoly = {
                2.*( dhdT - d2hdT2*T ),
                3. * d2hdT2 * T * T - 6. * dhdT * T,
                0.,
                d*h + T*T*T * ( 4. * dhdT - d2hdT2*T )
            };
            aPoly /= d ;
        }
        else
        {
            aPoly = { 0.5*dhdT/T, 0., h - 0.5 * dhdT * T } ;
        }

        // Diagnostic: the Grueneisen parameter implied by this anchor. It is not
        // used by the model - it is an independent cross check on alpha, cp, rho
        // and K at the split temperature, and it reads O(2) for every metal.
        // Logged rather than enforced, because the elastic data of several
        // materials is still provisional; the assertion below only catches a
        // scale or unit error, not imprecision.
        if (    this->have( MaterialProperty::E )
             && this->have( MaterialProperty::nu )
             && this->have( MaterialProperty::ref_density ) )
        {
            real rho0 = this->constant_property( MaterialProperty::ref_density );
            real K    = this->E( T ) / ( 3. * ( 1. - 2. * this->nu( T ) ) );
            real cv   = cp - 9. * alpha * alpha * T * K / rho0 ;
            real grueneisen = 3. * alpha * K / ( rho0 * cv );

            BELFEM_ERROR( grueneisen > 0.2 && grueneisen < 10.,
                "%s: implied Grueneisen parameter %g at %g K is outside any physical "
                "range - check alpha, cp, ref_density and the elastic data",
                mLabel.c_str(), ( double ) grueneisen, ( double ) T );
        }
    }

    }
}
