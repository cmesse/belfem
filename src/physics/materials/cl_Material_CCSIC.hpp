//
// Created by Christian Messe on 14.04.20.
//

#ifndef BELFEM_CL_CCSIC_HPP
#define BELFEM_CL_CCSIC_HPP

#include "typedefs.hpp"
#include "cl_OrthotropicMaterial.hpp"
#include "cl_Spline.hpp"
namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------
        class CCSIC : public OrthotropicMaterial
        {
            // density at room temperature
            const real mRho = 1900.0;

            // elastic parameters
            // source: 10.18419/opus-9381
            const real mE1 = 60.0e9;
            const real mE3 = 20.0e9;
            const real mG12 = 8.0e9;
            const real mG23 = 9.0e9;
            const real mNu13 = 0.2;
            const real mNu12 = 0.03;

            // implicit elastic parameters
            const real mE2 = mE1;
            const real mG13 = mG23;
            const real mNu23 = mNu13;
            const real mNu31 = mNu13 * mE3 / mE1;

            // spline for specific heat
            Spline * mHeat;

            // specific conductivity in 0/90 plane
            Vector <real> mLambda1;

            // specific conductivity in thickness direction
            Vector <real> mLambda3;

            // emittance
            Vector <real> mEpsilon;

            // thermal strain in fiber direction
            Vector <real> mMu1;

            // thermal strain perpendicular to fiber direction
            Vector <real> mMu3;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            CCSIC();

            ~CCSIC();

//----------------------------------------------------------------------------

            /**
              * Elasticity Matrix in 3D
              */
            void
            C( Matrix <real> & aC, const real aT = BELFEM_TREF ) const;

//----------------------------------------------------------------------------

            /**
             * Elasticity matrix in plane stress
             */
            void
            C_ps( Matrix <real> & aC, const real aT = BELFEM_TREF ) const;

//----------------------------------------------------------------------------

            /**
             * Elasticity matrix in rotation symmetry
             */
            void
            C_rot( Matrix <real> & aC, const real aT = BELFEM_TREF ) const;

//----------------------------------------------------------------------------

            // density in kg/m^3
            real
            rho( const real aT ) const;

//----------------------------------------------------------------------------

            // specific heat capacity in J/kg K
            real
            c( const real aT = BELFEM_TREF ) const;

//----------------------------------------------------------------------------

            /**
            * thermal conductivity in 2D
            */
            void
            lambda_p( Matrix <real> & aLambda, const real aT = BELFEM_TREF ) const;

//----------------------------------------------------------------------------

            /**
              * thermal conductivity in 3D
              */
            void
            lambda3d( Matrix <real> & aLambda, const real aT = BELFEM_TREF ) const;

//----------------------------------------------------------------------------

            /**
             * thermal absortivity
             * @param aT
             * @return
             */
            real
            epsilon( const real aT ) const;

//----------------------------------------------------------------------------

            void
            mu_ps( Matrix <real> & aMu, const real aT, const real aTref ) const;

//----------------------------------------------------------------------------

            void
            mu( Matrix <real> & aMu, const real aT, const real aTref ) const;

//----------------------------------------------------------------------------
        private:
//----------------------------------------------------------------------------

            void
            create_heat_spline();

//----------------------------------------------------------------------------
        };
//----------------------------------------------------------------------------
    } /* end namespace material */
}  /* end namespace belfem */
#endif //BELFEM_CL_CCSIC_HPP
