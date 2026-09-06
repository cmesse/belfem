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

#ifndef BELFEM_CL_MATERIAL_METAL_HPP
#define BELFEM_CL_MATERIAL_METAL_HPP


#include "cl_Material_SplineLookupTable.hpp"
#include "cl_Database.hpp"
#include "cl_BhCurve.hpp"
#include "cl_Bezier.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_ddpolyval.hpp"

namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        /**
         * @brief Base class for metallic materials (pure metals and alloys)
         *
         * Provides common functionality for metals including:
         * - Bloch-Grüneisen resistivity model for intrinsic resistivity (rho_i)
         * - Matthiessen's rule for combining intrinsic and residual resistivity
         * - Debye integrals (J3, J4, J5) for thermodynamic calculations
         * - Thermal conductivity using Wiedemann-Franz law with corrections
         * - Magnetoresistance via Kohler's rule (implemented in derived classes)
         * - Temperature-dependent Debye temperature calculations
         *
         * The class uses precomputed splines for efficient evaluation of Debye
         * integrals J_n(Z) where Z = θ/T (θ = Debye temperature).
         *
         * Reference: Bloch-Grüneisen law, Matthiessen's rule, Kohler's rule
         */
        class Metal : public SplineLookupTable
        {
            // Maximum values and integration limits for Debye integrals
            real mZ3Max = 0. ;
            real mZ4Max = 0. ;
            real mZ5Max = 0. ;
            real mZ6Max = 0. ;
            real mZ7Max = 0. ;
            real mZnMax = 0. ;

            real mJ3Max = 0. ;
            real mJ4Max = 0. ;
            real mJ5Max = 0. ;
            real mJ6Max = 0. ;
            real mJ7Max = 0. ;
            real mJnMax = 0. ;

            // we need to store this to avoid log(0)
            real mInvLog10 = 1.0/std::log( 10. ) ;
            real mDatabaseTmin = BELFEM_QUIET_NAN ;
            real mDatabaseTmax = BELFEM_QUIET_NAN ;
            real mDatabaseBmin = BELFEM_QUIET_NAN ;
            real mDatabaseBmax = BELFEM_QUIET_NAN ;
            Database * mRhoData = nullptr ;
            const BhCurve * mBhCurve = nullptr ;

            // Coefficients for thermal conductivity model
            Vector< real > mLambdaCoefficients ;

            bool mComputeTables = true ;

            real ( Metal::*mFunJ )( const real Z ) const ;

        protected:

            // Splines for Debye integrals J_n(θ/T)
            Spline * mSplineJ3 = nullptr ;  // ∫₀^Z x³eˣ/(eˣ-1)² dx
            Spline * mSplineJ4 = nullptr ;  // ∫₀^Z x⁴eˣ/(eˣ-1)² dx
            Spline * mSplineJ5 = nullptr ;  // ∫₀^Z x⁵eˣ/(eˣ-1)² dx
            Spline * mSplineJn = nullptr ;

            Cell< Vector< real > > mCpPolys ;
            Vector< real > mTCpSwitch ;
            Bezier * mCpBezierLow     = nullptr ;
            Bezier * mCpBezierMedium  = nullptr ;
            Bezier * mCpBezierHigh    = nullptr ;

            Vector< real > mWachtmanYoung ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * @brief Constructor
             * @param aLabel Material name (e.g., "Copper", "Silver")
             * @param aType Material type ( PureMetal, LookupAlloy or HTS )
             * @param aBuildTables Build the property lookup tables on construction
             */
            Metal( const string & aLabel,
                   const MaterialType aType,
                   const bool aBuildTables = true );

            /**
             * @brief Destructor - deletes Debye integral splines
             */
            ~Metal() override;

            /**
             * @brief Set residual resistivity ratio
             * @param RRR Residual resistivity ratio ρ(273.15K)/ρ(0K)
             *
             * Iteratively solves for ρ₀ (residual resistivity) given RRR and
             * the reference intrinsic resistivity ρᵢ(T_ref). Uses Matthiessen's
             * rule: ρ_total = ρᵢ(T) + ρ₀
             */
            void
            set_RRR( const real RRR ) override ;

            /**
             * @brief Compute Debye temperature from measured resistivity
             * @param T Temperature [K]
             * @param rho Measured electrical resistivity [Ω·m]
             * @param theta_guess Initial guess for the Debye temperature [K];
             *        defaults to the material's tabulated value
             * @return Debye temperature [K]
             *
             * Inverse calculation using the Bloch-Grüneisen resistivity model.
             */
            virtual real
            compute_debye_from_rho( const real T, const real rho, const real theta_guess=BELFEM_QUIET_NAN );

            /**
             * @brief Compute Debye temperature from specific heat capacity
             * @param T Temperature [K]
             * @param cp Specific heat capacity [J/(kg·K)] (default: use material's cp(T))
             * @param theta_guess Initial guess for the Debye temperature [K];
             *        defaults to the material's tabulated value
             * @return Debye temperature [K]
             *
             * Inverse calculation from Debye model for heat capacity.
             */
            real
            compute_debye_from_cp( const real T, const real cp=BELFEM_QUIET_NAN, const real theta_guess=BELFEM_QUIET_NAN );

            /**
             * @brief Compute Debye temperature from the specific heat at constant volume
             * @param T Temperature [K] - cryogenic only, see below
             * @param cv Specific heat [J/(kg·K)] (default: use material's cp(T), valid where cp = cv)
             * @param theta_guess Initial guess for the Debye temperature [K]
             * @return Debye temperature [K]
             *
             * Same inversion as compute_debye_from_cp() without the dilation term,
             * so it does not need the elastic data. Refuses ( BELFEM_ERROR ) where
             * the neglected cp - cv exceeds 0.1 % of cp, i.e. above roughly 100 K.
             */
            real
            compute_debye_from_cv( const real T, const real cv=BELFEM_QUIET_NAN, const real theta_guess=BELFEM_QUIET_NAN );

            // Bring base class single-parameter versions into scope
            // (needed because C++ name hiding hides Material::rho(T) when we override rho(T,B,beta))
            using Material::rho;
            using Material::lambda;
            using Material::drhodT;
            using Material::dlambdadT;

            /**
             * @brief Electrical resistivity with magnetoresistance (Kohler's rule)
             * @param T Temperature [K]
             * @param B Magnetic field magnitude [T]
             * @param beta Angle between current and field [rad]
             * @return Electrical resistivity [Ω·m]
             *
             * Applies Kohler's rule: ρ(B,T) = ρ₀(T) · [1 + K(B·S, β)]
             * where S = ρ_ref/ρ₀(T) is the similarity parameter.
             */
            real
            rho( const real T, const real B, const real beta ) const override ;

            real
            drhodT( const real T, const real B, const real beta ) const override ;

            real
            drhodB( const real T, const real B, const real beta ) const override ;

            real
            drhodbeta( const real T, const real B, const real beta ) const override ;

            /**
             * @brief Thermal conductivity with magnetoresistance
             * @param T Temperature [K]
             * @param B Magnetic field magnitude [T]
             * @param beta Angle between current and field [rad]
             * @return Thermal conductivity [W/(m·K)]
             *
             * Uses Wiedemann-Franz relation corrected for magnetoresistance.
             */
            real
            lambda( const real T, const real B, const real beta ) const override ;

            real
            dlambdadT( const real T, const real B, const real beta ) const override ;

            real
            dlambdadB( const real T, const real B, const real beta ) const override ;

            real
            dlambdadbeta( const real T, const real B, const real beta ) const override ;

//------------------------------------------------------------------------------

            /**
            * @brief Relative length after thermal expansion
            * @param T Temperature [K] (default: room temperature)
            * @return Relative length l/l₀ [-]
            */
            real
            l( real T ) const override ;

//------------------------------------------------------------------------------

            void
            set_bh_curve( const BhCurve * aCurve ) override ;

//------------------------------------------------------------------------------
            /*
             * @brief Evaluate property using spline interpolation
             */
            real
            spline_property( const MaterialProperty aProperty, const real aX ) const override ;

//------------------------------------------------------------------------------

            /**
             * @brief Set reference intrinsic resistivity for Bloch-Grüneisen model
             * @param T_ref Reference temperature [K] (typically 273.15 K)
             * @param rho_i_ref Reference intrinsic resistivity [Ω·m]
             * @param theta Debye temperature [K] (default: use material's debye(T_ref))
             *
             * Computes the Bloch-Grüneisen parameter A from the reference point.
             */
            void
            set_rho_i_ref( const real T_ref, const real rho_i_ref, const real theta=BELFEM_QUIET_NAN );

            /**
             * @brief Debye integral Jₙ(Z) = ∫₀^Z xⁿeˣ/(eˣ-1)² dx with n the
             *        Bloch–Grüneisen exponent ( 5 by default, see
             *        set_bloch_gruen_parameter; dispatched to J3/J4/J5/Jn )
             * @param Z Reduced temperature Z = θ/T
             * @return Integral value
             */
            real
            J( const real Z ) const ;


            /**
             * @brief Debye integral J₃(Z) = ∫₀^Z x³eˣ/(eˣ-1)² dx
             * @param Z Reduced temperature Z = θ/T
             * @return Integral value
             */
            real
            J3( const real Z ) const ;

            /**
             * @brief Debye integral J₄(Z) = ∫₀^Z x⁴eˣ/(eˣ-1)² dx
             * @param Z Reduced temperature Z = θ/T
             * @return Integral value (used for specific heat)
             */
            real
            J4( const real Z ) const ;

            /**
             * @brief Debye integral J₅(Z) = ∫₀^Z x⁵eˣ/(eˣ-1)² dx
             * @param Z Reduced temperature Z = θ/T
             * @return Integral value (used for Bloch-Grüneisen resistivity)
             */
            real
            J5( const real Z ) const ;

            /**
                 * @brief Debye integral for arbitrary n
                 * @param Z Reduced temperature Z = θ/T
                 * @return Integral value (used for Bloch-Grüneisen resistivity)
                 */
            real
            Jn( const real Z ) const ;


            /**
             * @brief Total resistivity at zero field (overrides base class)
             * @param T Temperature [K]
             * @return ρ(T) = ρᵢ(T) + ρ₀ [Ω·m]
             *
             * Matthiessen's rule: intrinsic + residual resistivity.
             */
            real
            rho_custom(const real T) const override;

            /**
             * @brief Intrinsic resistivity (Bloch-Grüneisen model)
             * @param T Temperature [K]
             * @return Intrinsic resistivity ρᵢ(T) [Ω·m]
             *
             * Uses ρᵢ(T) = A · J₅(θ/T) / (θ/T)⁵
             */
            real
            rho_i_custom( const real T ) const override ;

            /**
             * @brief Thermal conductivity at zero field
             * @param T Temperature [K]
             * @return Thermal conductivity [W/(m·K)]
             *
             * Uses empirical model with Wiedemann-Franz corrections.
             */
            real
            lambda_custom( const real T ) const  override ;

            /**
             * @brief Initialize resistivity properties
             *
             * Sets up dependencies and creates splines for ρ and ρᵢ.
             */
            void
            create_rho();

            /**
             * @brief Set thermal conductivity model coefficients
             * @param aCoeffs Vector of 8 coefficients for empirical lambda model
             *
             * For pure metals only. Coefficients parameterize the thermal
             * conductivity including electron-phonon scattering effects.
             */
            void
            set_lambda_coefficients( const Vector< real > & aCoeffs );

            /**
             * @brief Phonon group velocity
             * @param T Temperature [K]
             * @return Group velocity [m/s]
             */
            real
            group_velocity( const real T ) const;

            /**
             * @brief Grüneisen parameter
             * @param T Temperature [K]
             * @return Dimensionless Grüneisen parameter γ
             */
            real
            grueneisen( const real T ) const ;

            /**
             * @brief Grüneisen parameter and speed of sound
             * @param T Temperature [K]
             * @param gamma Output: Grüneisen parameter
             * @param c Output: Speed of sound [m/s]
             */
            void
            grueneisen( const real T, real & gamma, real & c ) const;

            /**
             * @brief Specific heat from Debye temperature
             * @param T Temperature [K]
             * @param theta Debye temperature [K]
             * @return Specific heat capacity [J/(kg·K)]
             */
            real
            cp_from_debye( const real T, const real theta ) const ;

            /**
             * @brief Specific heat at constant volume from Debye temperature
             *        ( Sommerfeld + Debye, no dilation term )
             */
            real
            cv_from_debye( const real T, const real theta ) const ;

            /**
             * @brief Kohler function for magnetoresistance
             * @param B Magnetic field magnitude [T]
             * @param S Similarity parameter S = ρ_ref/ρ₀(T)
             * @param beta Angle between field and current [rad]
             * @return Relative magnetoresistance Δρ/ρ₀
             *
             * Must be implemented in derived classes (Copper, Silver, etc.).
             * Returns NaN in base class.
             */
            virtual real
            kohler( const real B, const real S, const real beta ) const ;

            /**
             * @brief Find critical B·S for given relative resistivity increase
             * @param delta_rho_res Relative resistivity increase Δρ/ρ
             * @return Critical B·S value [T]
             */
            real
            kohler_find_bs_crit( const real delta_rho_res );

            void
            populate_rho_database();

//------------------------------------------------------------------------------

            real
            H_bhcurve( const real B ) const override ;

            real
            mu_bhcurve( real H, const real T ) const override ;

            void
            dmudH_bhcurve( const real H, real & mu, real & dmudH ) const override ;

            const Vector< real > &
            lambda_coefficients() const ;

            void
            set_table_flags( const bool aFlag ) override ;

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            /**
          * @brief Create spline for Debye integral J_n(Z)
          * @param aExponent Integral order ( 3, 4, 5, or any real such as Iron's 4.5 )
          * @return Pointer to created spline
          */
            Spline *
            create_J_spline( const real aExponent );

            real
            rho_table( const real T, const real B, const real beta ) const override ;

            real
            drhodT_table( const real T, const real B, const real beta ) const override ;

            real
            drhodB_table( const real T, const real B, const real beta ) const override ;

            real
            drhodbeta_table( const real T, const real B, const real beta ) const override ;

            real
            rho_kohler( const real T, const real B, const real beta ) const override ;

            real
            drhodT_kohler( const real T, const real B, const real beta ) const override ;

            real
            drhodB_kohler( const real T, const real B, const real beta ) const override ;

            real
            drhodbeta_kohler( const real T, const real B, const real beta ) const override ;

            real
            cp_custom( const real T) const override;

            real
            dcpdT_custom( const real T) const override;

            real
            d2cpdT2_custom( const real T) const override;


            void
            set_kohler_dependencies();

            void
            set_bloch_gruen_parameter( const real n );

            void
            create_cp(
                const Vector< real > & Px,
                const Vector< real > & Py,
                const Vector< real > & Qx,
                const Vector< real > & Qy,
                const Vector< real > & Rx = {},
                const Vector< real > & Ry = {} );

            real
            E_custom(const real T) const override;

            real
            dEdT_custom(const real T) const override;

            void
            create_mech( const real E0, const real b, const real T1, const real T2, const real nu2 );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            /**
             * @brief Regula falsi for the Debye temperature: solves
             *        aFunction( T, theta ) = y for theta
             */
            real
            invert_debye( const real T, const real y, const real theta_guess,
                          real ( Metal::*aFunction )( const real, const real ) const );

            void
            save_rho_database( const std::string & aPath );

            void
            load_rho_database( const std::string & aPath );

            /**
             * @brief Objective function for Debye temperature iteration from ρ
             * @param T Temperature [K]
             * @param rho Measured resistivity [Ω·m]
             * @param theta Trial Debye temperature [K]
             * @return Residual
             */
            real
            fun_debye_from_rho( const real T, const real rho, const real theta );

            void
            populate_rho_database_serial();

            void
            populate_rho_database( Mesh * aMesh, Vector< real > & aRho );


        };

        inline real Metal::J3( const real Z ) const
        {
            BELFEM_ASSERT( mSplineJ3 != nullptr, "lookup table for J3 is not initialized" ) ;

            if ( Z > mZ3Max ) return mJ3Max ;
            return mSplineJ3->eval( Z ) ;
        }

        inline real Metal::J4( const real Z ) const
        {
            BELFEM_ASSERT( mSplineJ4 != nullptr, "lookup table for J4 is not initialized" ) ;

            if ( Z > mZ4Max ) return mJ4Max ;
            return mSplineJ4->eval( Z ) ;
        }

        inline real Metal::J5( const real Z ) const
        {
            BELFEM_ASSERT( mSplineJ5 != nullptr, "lookup table for J5 is not initialized" ) ;

            if ( Z > mZ5Max ) return mJ5Max ;
            return mSplineJ5->eval( Z ) ;
        }


        inline real Metal::Jn( const real Z ) const
        {
            BELFEM_ASSERT( mSplineJn != nullptr, "lookup table for Jn is not initialized" ) ;

            if ( Z > mZnMax ) return mJnMax ;
            return mSplineJn->eval( Z ) ;
        }


        inline real
        Metal::J( const real Z ) const
        {
            return ( this->*mFunJ )( Z );
        }

        inline real
        Metal::rho_custom(const real T) const
        {
            return this->rho_i( T ) + this->constant_property( MaterialProperty::rho_0 );
        }

        inline real
        Metal::rho( const real T, const real B, const real beta ) const
        {
            return (this->*mFunctionRhoKohler)( T, B, beta );
        }

        inline real
        Metal::rho_table( const real T, const real B, const real beta ) const
        {
            BELFEM_ASSERT( mRhoData != nullptr, "lookup table for rho is not initialized" ) ;

            real theta = std::clamp( T, mDatabaseTmin, mDatabaseTmax );
            real log10B = std::log(std::clamp(  B, mDatabaseBmin, mDatabaseBmax ) )*mInvLog10 ;
            real angle = beta < 0 ? beta + constant::pi : beta > constant::pi ? beta - constant::pi : beta;

            return std::exp(mRhoData->evaluate( theta, log10B, angle )) ;
        }

        inline real
        Metal::drhodT( const real T, const real B, const real beta ) const
        {
            return (this->*mFunctiondRhoKohlerdT)( T, B, beta );
        }

        inline real
        Metal::drhodB( const real T, const real B, const real beta ) const
        {
            return (this->*mFunctiondRhoKohlerdB)( T, B, beta );
        }

        inline real
        Metal::drhodbeta( const real T, const real B, const real beta ) const
        {
            return (this->*mFunctiondRhoKohlerdbeta)( T, B, beta );
        }

        inline real
        Metal::drhodT_table( const real T, const real B, const real beta ) const
        {
            BELFEM_ASSERT( mRhoData != nullptr, "lookup table for rho is not initialized" ) ;

            real theta = std::clamp( T, mDatabaseTmin, mDatabaseTmax );
            real log10B = std::log(std::clamp(  B, mDatabaseBmin, mDatabaseBmax ) )*mInvLog10 ;
            real angle = beta < 0 ? beta + constant::pi : beta > constant::pi ? beta - constant::pi : beta;

            real y    = mRhoData->evaluate( theta, log10B, angle ) ;
            real dydT = mRhoData->evaluate_derivx( theta, log10B, angle ) ;
            return std::exp( y ) * dydT ;
        }

        inline real
        Metal::drhodB_table( const real T, const real B, const real beta ) const
        {
            BELFEM_ASSERT( mRhoData != nullptr, "lookup table for rho is not initialized" ) ;

            real theta = std::clamp( T, mDatabaseTmin, mDatabaseTmax );
            real Bc = std::clamp(  B, mDatabaseBmin, mDatabaseBmax ) ;

            real log10B = std::log( Bc )*mInvLog10 ;
            real angle = beta < 0 ? beta + constant::pi : beta > constant::pi ? beta - constant::pi : beta;

            real y    = mRhoData->evaluate( theta, log10B, angle ) ;
            real dydlogB = mRhoData->evaluate_derivy( theta, log10B, angle ) ;

            return std::exp( y ) * dydlogB *mInvLog10 / Bc ;

        }

        inline real
        Metal::drhodbeta_table( const real T, const real B, const real beta ) const
        {
            BELFEM_ASSERT( mRhoData != nullptr, "lookup table for rho is not initialized" ) ;

            real theta = std::clamp( T, mDatabaseTmin, mDatabaseTmax );
            real Bc = std::clamp(  B, mDatabaseBmin, mDatabaseBmax ) ;

            real log10B = std::log( Bc )*mInvLog10 ;
            real angle = beta < 0 ? beta + constant::pi : beta > constant::pi ? beta - constant::pi : beta;

            real y    = mRhoData->evaluate( theta, log10B, angle ) ;
            real dydbeta = mRhoData->evaluate_derivz( theta, log10B, angle ) ;

            return std::exp( y ) * dydbeta ;

        }

        inline real
        Metal::lambda( real T, const real B, const real beta ) const
        {
            return  this->lambda( T ) * this->rho( T ) / this->rho( T, B, beta );
        }

        inline real
        Metal::dlambdadT( real T, const real B, const real beta ) const
        {
            real a = this->lambda( T ) ;
            real b = this->rho( T ) ;
            real c = this->rho( T, B, beta );

            real da = this->dlambdadT( T );
            real db = this->drhodT( T );
            real dc = this->drhodT( T, B, beta );

            return ( c * ( da*b + a * db ) - a * b * dc ) /( c * c );
        }

        // lambda( T, B, beta ) = lambda(T) * rho(T) / rho(T,B,beta):
        // only the denominator depends on B and beta
        inline real
        Metal::dlambdadB( const real T, const real B, const real beta ) const
        {
            real a = this->lambda( T ) ;
            real b = this->rho( T ) ;
            real c = this->rho( T, B, beta );

            return - a * b * this->drhodB( T, B, beta ) / ( c * c );
        }

        inline real
        Metal::dlambdadbeta( const real T, const real B, const real beta ) const
        {
            real a = this->lambda( T ) ;
            real b = this->rho( T ) ;
            real c = this->rho( T, B, beta );

            return - a * b * this->drhodbeta( T, B, beta ) / ( c * c );
        }

        inline real
        Metal::H_bhcurve( const real B ) const
        {
            return mBhCurve->H( B );
        }


        inline real
        Metal::mu_bhcurve( real H, const real T ) const
        {
            return mBhCurve->mu( H );
        }


        inline void
        Metal::dmudH_bhcurve( const real H, real & mu, real & dmudH ) const
        {
            mBhCurve->dmudH( H, mu, dmudH );
        }

        inline
        const Vector< real > & Metal::lambda_coefficients() const
        {
            return mLambdaCoefficients;
        }

        inline real
        Metal::cp_custom( const real T ) const
        {
            if ( T < mTCpSwitch( 0 ) )
            {
                return polyval( mCpPolys( 0 ), T );
            }
            if ( T < mTCpSwitch( 1 ) )
            {
                return std::exp( polyval( mCpPolys( 1 ), std::log( T ) ) );
            }
            if ( T < mTCpSwitch( 2 ) )
            {
                return std::exp( mCpBezierLow->y(  std::log( T ) ) );
            }
            if ( T < mTCpSwitch( 3 ) )
            {
                return std::exp( mCpBezierMedium->y(  std::log( T ) ) );
            }
            if ( mCpBezierHigh != nullptr && T < mTCpSwitch( mTCpSwitch.length() - 1 ) )
            {
                return std::exp( mCpBezierHigh->y(  std::log( T ) ) );
            }
            return polyval( mCpPolys( 2 ), T );
        }

        inline real
        Metal::dcpdT_custom( const real T) const
        {
            // Subclasses that manage cp themselves ( YBCO, HastelloyC276 )
            // override cp_custom but inherit this function; without the
            // Bezier machinery, fall back to the finite difference the
            // base class provides
            if ( mCpBezierLow == nullptr )
            {
                return Material::dcpdT_custom( T );
            }

            if ( T < mTCpSwitch( 0 ) )
            {
                return dpolyval( mCpPolys( 0 ), T );
            }
            real x = std::log( T );
            if ( T < mTCpSwitch( 1 ) )
            {
                return std::exp( polyval( mCpPolys( 1 ), x ) ) * dpolyval( mCpPolys( 1 ), x ) / T;
            }
            if ( T < mTCpSwitch( 2 ) )
            {
                return std::exp( mCpBezierLow->y( x ) ) * mCpBezierLow->dydx( x ) / T ;
            }
            if ( T < mTCpSwitch( 3 ) )
            {

                return std::exp( mCpBezierMedium->y( x ) ) * mCpBezierMedium->dydx( x ) / T ;
            }
            if ( mCpBezierHigh != nullptr && T < mTCpSwitch( mTCpSwitch.length() - 1 ) )
            {
                return std::exp( mCpBezierHigh->y( x ) ) * mCpBezierHigh->dydx( x ) / T ;
            }

            // the extrapolation above T3 is linear in T, not in log-log space
            return dpolyval( mCpPolys( 2 ), T );
        }

        inline real
        Metal::d2cpdT2_custom( const real T) const
        {
            if ( mCpBezierLow == nullptr )
            {
                return Material::d2cpdT2_custom( T );
            }
            if ( T < mTCpSwitch( 0 ) )
            {
                return ddpolyval( mCpPolys( 0 ), T );
            }
            real x = std::log( T );
            // with y = ln( cp ) over x = ln( T ):
            // d2cp/dT2 = cp * ( y'' + y'^2 - y' ) / T^2
            if ( T < mTCpSwitch( 1 ) )
            {
                real y  = polyval( mCpPolys( 1 ), x );
                real dy = dpolyval( mCpPolys( 1 ), x );
                return std::exp( y ) * ( ddpolyval( mCpPolys( 1 ), x ) + dy*dy - dy ) / ( T * T );
            }
            if ( T < mTCpSwitch( 2 ) )
            {
                real dy = mCpBezierLow->dydx( x );
                return std::exp( mCpBezierLow->y( x ) )
                    * ( mCpBezierLow->d2ydx2( x ) + dy*dy - dy ) / ( T * T );
            }
            if ( T < mTCpSwitch( 3 ) )
            {
                real dy = mCpBezierMedium->dydx( x );
                return std::exp( mCpBezierMedium->y( x ) )
                    * ( mCpBezierMedium->d2ydx2( x ) + dy*dy - dy ) / ( T * T );
            }
            if ( mCpBezierHigh != nullptr && T < mTCpSwitch( mTCpSwitch.length() - 1 ) )
            {
                real dy = mCpBezierHigh->dydx( x );
                return std::exp( mCpBezierHigh->y( x ) )
                    * ( mCpBezierHigh->d2ydx2( x ) + dy*dy - dy ) / ( T * T );
            }

            return  ddpolyval( mCpPolys( 2 ), T );
        }

        inline real
        Metal::E_custom( const real T ) const
        {
            BELFEM_ASSERT( mWachtmanYoung.length() >= 3,
                "Wachtman coefficients not assigned for %s",
                this->label().c_str() );

            if ( T < BELFEM_EPSILON ) return  mWachtmanYoung( 0 );

            real E0 = mWachtmanYoung( 0 );
            real  b = mWachtmanYoung( 1 );
            real T0 = mWachtmanYoung( 2 );
            return E0 - b * T * std::exp( - T0/T );
        }

        inline real
        Metal::dEdT_custom( const real T ) const
        {
            BELFEM_ASSERT( mWachtmanYoung.length() >= 3,
               "Wachtman coefficients not assigned for %s",
               this->label().c_str() );

            if ( T < BELFEM_EPSILON ) return  0 ;

            real  b = mWachtmanYoung( 1 );
            real T0 = mWachtmanYoung( 2 );
            return - b *  std::exp( - T0/T ) * ( T + T0 ) / T ;
        }

    }
}

#endif // BELFEM_CL_MATERIAL_METAL_HPP