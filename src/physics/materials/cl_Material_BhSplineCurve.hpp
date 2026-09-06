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

#ifndef BELFEM_CL_MATERIAL_BHSPLINECURVE_HPP
#define BELFEM_CL_MATERIAL_BHSPLINECURVE_HPP

#include "cl_BhCurve.hpp"
#include "cl_Spline.hpp"
namespace belfem
{
    namespace material
    {
        class BhSplineCurve : public BhCurve
        {
            Spline * mNuSpline = nullptr ;  //!< ν(B) = H/B cubic spline [A/(T·m)]
            Spline * mMuSpline = nullptr ;  //!< 1/μ(H) = H/B cubic spline inverse

            real mBsat = BELFEM_QUIET_NAN ;  //!< Saturation flux density [T]
            real mHsat = BELFEM_QUIET_NAN ;  //!< Saturation field intensity [A/m]
            real mMsat = BELFEM_QUIET_NAN ;  //!< Saturation magnetization [A/m]

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            BhSplineCurve(const string & aPath, const string & aLabel);

//------------------------------------------------------------------------------

            ~BhSplineCurve() override ;

//------------------------------------------------------------------------------
               /**
             * @brief Reluctivity ν(B) = H/B
             * @param B Magnetic flux density magnitude [T]
             * @return Reluctivity ν [A/(T·m)]
             *
             * For B < Bsat: Uses cubic spline interpolation
             * For B ≥ Bsat: Returns ν₀ - Msat/B (linear + offset)
             */
            real
            nu( const real B ) const override ;
//------------------------------------------------------------------------------

            /**
             * @brief Permeability μ(H) = B/H
             * @param H Magnetic field intensity magnitude [A/m]
             * @return Permeability μ [T·m/A]
             *
             * For H < Hsat: Uses cubic spline interpolation
             * For H ≥ Hsat: Returns μ₀·(1 + Msat/H)
             */
            real
            mu( const real H ) const override ;

//------------------------------------------------------------------------------

            /**
             * @brief Permeability and its derivative with respect to H
             * @param H Magnetic field intensity magnitude [A/m]
             * @param[out] mu Permeability μ(H) [T·m/A]
             * @param[out] dmudH Derivative dμ/dH [T·m/A²]
             *
             * Computes both μ and dμ/dH efficiently for Newton-Raphson
             * iterations in nonlinear magnetic FEM solvers.
             *
             * For H < Hsat: Uses spline evaluation and derivative
             * For H ≥ Hsat: Analytical formulas for linear regime
             */
            void
            dmudH( const real H, real & mu, real & dmudH ) const override ;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            load_data( const string & aPath, const string & aLabel );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline real
        BhSplineCurve::nu( const real B ) const
        {
            if ( B < mBsat )
            {
                return mNuSpline->eval( B ) * constant::nu0 ;
            }
            else
            {
                return ( B * constant::nu0 - mMsat ) / B ;
            }
        }

//------------------------------------------------------------------------------

        inline real
        BhSplineCurve::mu( const real H ) const
        {
            if ( H < mHsat )
            {
                return constant::mu0  / mMuSpline->eval( H )  ;
            }
            else
            {
                return ( mMsat + H ) / H * constant::mu0 ;
            }
        }

//------------------------------------------------------------------------------

        inline void
        BhSplineCurve::dmudH( const real H, real & mu, real & dmudH ) const
        {
            if ( H < mHsat )
            {
                index_t k = mMuSpline->find_col( H );
                real f = mMuSpline->eval( H, k );
                real df = mMuSpline->deval( H, k );

                mu = constant::mu0 / f ;
                dmudH = - mu * df / f;
            }
            else
            {
                mu    = ( mMsat + H ) / H * constant::mu0 ;
                dmudH =  -mMsat * constant::mu0 / ( H * H );
            }
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_MATERIAL_BHSPLINECURVE_HPP
