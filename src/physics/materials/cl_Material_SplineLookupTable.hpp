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

#ifndef BELFEM_CL_MATERIAL_SPLINELOOKUPTABLE_HPP
#define BELFEM_CL_MATERIAL_SPLINELOOKUPTABLE_HPP
#include "cl_Material.hpp"
#include "cl_Spline.hpp"
#include "cl_Vector.hpp"
#include "cl_Bezier.hpp"
#include "fn_polyval.hpp"

namespace belfem
{
    namespace material
    {
        class SplineLookupTable : public Material
        {
        protected:

            // Spline interpolations for tabulated properties
            Cell< Spline * > mSplines ;

        public:

            SplineLookupTable( const MaterialType aType, const bool aIsIsotropic=true );

            ~SplineLookupTable() override ;

            /**
             * @brief Assign a spline for property interpolation
            */
            void
            set_spline( const MaterialProperty aProperty, Spline * aSpline=nullptr ) ;

            /**
            * @brief Get spline for a property
            */
            Spline *
            spline( const MaterialProperty aProperty ) ;

            real
            spline_property( const MaterialProperty aProperty, const real aX ) const override ;

            real
            dspline_property( const MaterialProperty aProperty, const real aX ) const override ;

            real
            ddspline_property( const MaterialProperty aProperty, const real aX ) const override ;

            /**
             * @brief Relative length after thermal expansion, from integrated alpha spline
             */
            real
            l( real T ) const override ;

        protected:

            // un-hide the convenience overload create_spline( aProperty, adYdX0, adYdX1 )
            using Material::create_spline ;

            void
            reset_spline( const MaterialProperty aProperty ) override ;

            void
            create_spline( real (Material::*aFunction)(const real aT) const,
                const MaterialProperty aProperty,
                const uint aStartBC,
                const uint aEndBC,
                const real adYdX0,
                const real adYdX1 ) override ;

            void
            create_cryo_expansion(
                const Bezier * aThermalExpansion,
                Vector< real > & aThermalExpansionCryo,
                const real aTSwitch = BELFEM_QUIET_NAN ) ;

            /**
             * @brief cryogenic branch from alpha and its derivatives at the split
             *        temperature, for materials whose alpha( T ) is not a dL/L Bezier
             */
            void
            create_cryo_expansion(
                const real alpha, const real dalphadT, const real d2alphadT2,
                Vector< real > & aThermalExpansionCryo ) ;

            /**
             * @brief cryogenic branch anchored on the VALUE of alpha at the split
             *        only. The first two derivatives at the split come from the
             *        Grueneisen relation C = alpha / cp = C* K( T* ) / K( T ), with
             *        K from E and nu; below the split ln( C ) is the fitted
             *        polynomial of create_low_temperature_alpha(). For materials
             *        whose expansion curve is data-backed at the split but not
             *        below it ( Magnesia, HastelloyC276, YBCO ).
             */
            void
            create_cryo_expansion_anchored(
                const real alpha,
                Vector< real > & aThermalExpansionCryo,
                const real aTSwitch = BELFEM_QUIET_NAN ) ;


            /**
             * @brief composite thermal expansion for a Material::alpha_custom()
             *        override: the Grueneisen polynomial aPoly below the switch
             *        temperature, the Bezier aBezier above it. Not itself a
             *        Material::alpha_custom overload, hence the distinct name.
             */
            real
            alpha_composite( const Bezier * aBezier, const Vector< real > & aPoly, const real T ) const ;


            /**
             * @brief Fix the split temperature of the cryogenic expansion branch
             * @param aTSwitch explicit value [K], or NaN for min( 0.618 theta_D, gTAlphaSwitchMax )
             * @return the split temperature
             */
            real
            set_alpha_switch_temperature( const real aTSwitch = BELFEM_QUIET_NAN );

            /**
             * @brief Fit the cryogenic branch alpha = C( T ) cp( T ) below the split
             *        temperature from a dL/L Bezier ( the pure metals )
             */
            void
            create_low_temperature_alpha( const Bezier * aBezier , Vector< real > & aPoly,
                                          const real aTSwitch = BELFEM_QUIET_NAN );

            /**
             * @brief Same fit from alpha and its first two derivatives at the split
             *        temperature, which the caller evaluated from any alpha( T )
             *        source after set_alpha_switch_temperature()
             */
            void
            create_low_temperature_alpha( const real alpha, const real dalphadT, const real d2alphadT2,
                                          Vector< real > & aPoly, const bool aCheckGuard = true );

        private:

            void
            finish_cryo_expansion( const Vector< real > & aThermalExpansionCryo );

        };

        inline Spline *
        SplineLookupTable::spline( const MaterialProperty aProperty )
        {
            return mSplines( static_cast< size_t >( aProperty ) ) ;
        }

        inline real
        SplineLookupTable::spline_property( const MaterialProperty aProperty, const real aX ) const
        {
            BELFEM_ASSERT( mSplines( static_cast< size_t >( aProperty ) ) != nullptr,
                "Spline for property %s is not set", to_string( aProperty ).c_str() );
            return mSplines( static_cast< size_t >( aProperty ) )->eval( aX ) ;
        }

        inline real
        SplineLookupTable::dspline_property( const MaterialProperty aProperty, const real aX ) const
        {
            BELFEM_ASSERT( mSplines( static_cast< size_t >( aProperty ) ) != nullptr,
                "Spline for property %s is not set", to_string( aProperty ).c_str() );
            return mSplines( static_cast< size_t >( aProperty ) )->deval( aX ) ;
        }

        inline real
        SplineLookupTable::ddspline_property( const MaterialProperty aProperty, const real aX ) const
        {
            BELFEM_ASSERT( mSplines( static_cast< size_t >( aProperty ) ) != nullptr,
                "Spline for property %s is not set", to_string( aProperty ).c_str() );
            return mSplines( static_cast< size_t >( aProperty ) )->ddeval( aX ) ;
        }

        inline real
        SplineLookupTable::l( const real T ) const
        {
            if ( std::abs(T- this->constant_property( MaterialProperty::T_ref_density )) < BELFEM_EPSILON )
            {
                return 1.0 ;
            }

            return std::exp(
            mSplines( static_cast< uint >( MaterialProperty::alpha ) )->integrate( T ) );
        }

//------------------------------------------------------------------------------

        inline real
        SplineLookupTable::alpha_composite(
            const Bezier * aBezier,
            const Vector< real > & aPoly,
            const real T ) const
        {
            if ( std::abs( T ) < BELFEM_EPSILON ) return 0. ;

            if ( T < this->alpha_switch_temperature() )
            {
                // aPoly holds ln( C ) with C = alpha / cp ( Grueneisen equation of state )
                return std::exp( polyval( aPoly, T ) ) * this->cp( T );
            }
            else
            {
                real xi = aBezier->xi_by_x( T );
                real dLdxi, dTdxi ;
                aBezier->dpoint( xi, dTdxi, dLdxi );
                real L = 1 + aBezier->y_by_xi( xi );
                return dLdxi / ( L * dTdxi );
            }
        }
    }
}

#endif //BELFEM_CL_MATERIAL_SPLINELOOKUPTABLE_HPP
