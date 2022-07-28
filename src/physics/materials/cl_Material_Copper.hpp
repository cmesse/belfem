//
// Created by Christian Messe on 29.07.20.
//

#ifndef BELFEM_CL_COPPER_HPP
#define BELFEM_CL_COPPER_HPP


#include "typedefs.hpp"
#include "cl_Vector.hpp"

#include "cl_IsotropicMaterial.hpp"
#include "cl_Spline.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------

        class Copper : public IsotropicMaterial
        {

            Spline * mCpSpline = nullptr ;
            Spline * mRhoSpline = nullptr ;

            const real mSwitchCT0 = 9.;
            const real mSwitchCT1 = 15.;
            const real mSwitchCT2 = 50.;
            const real mSwitchCT3 = 100.;
            const real mSwitchCT4 = 200.;
            const real mSwitchCT5 = 220.;
            const real mSwitchCT6 = 280.;
            const real mSwitchCT7 = 300. ;
                  real mSwitchCT8 = 1125.0 ;


            Vector <real> mSpecificHeatPoly0;
            Vector <real> mSpecificHeatPoly1;
            Vector <real> mSpecificHeatPoly2;
            Vector <real> mSpecificHeatPoly3;
            Vector <real> mSpecificHeatPoly4;
            Vector <real> mSpecificHeatPoly5;
            Vector <real> mSpecificHeatPoly6;
            Vector <real> mSpecificHeatPoly7;
            Vector <real> mSpecificHeatPoly8;
            Vector <real> mSpecificHeatPoly9;

            const real mSwitchRT0 =   4.0 ;
            const real mSwitchRT1 =  10.0 ;
            const real mSwitchRT2 =  15.0 ;
            const real mSwitchRT3 =  50.0 ;
            const real mSwitchRT4 =  75.0 ;
            const real mSwitchRT5 = 200.0 ;
            const real mSwitchRT6 = 250.0 ;

            Vector< real > mResistivityPoly0 ;
            Vector< real > mResistivityPoly1 ;
            Vector< real > mResistivityPoly2 ;
            Vector< real > mResistivityPoly3 ;
            Vector< real > mResistivityPoly4 ;
            Vector< real > mResistivityPoly5 ;
            Vector< real > mResistivityPoly6 ;
            Vector< real > mResistivityPoly7 ;

            Vector <real> mThermalConductivityPoly0;
            Vector <real> mThermalConductivityPoly1;
            Vector <real> mThermalConductivityPoly2;
            Vector <real> mThermalConductivityPoly3;
            Vector <real> mThermalConductivityPoly4;
            Vector <real> mThermalConductivityPoly5;
            Vector <real> mThermalConductivityPoly6;
            Vector <real> mThermalConductivityPoly7;
            Vector <real> mThermalConductivityPoly8;
            Vector <real> mThermalConductivityPoly9;

            real mSwitchLambdaT0;
            real mSwitchLambdaT1;
            real mSwitchLambdaT2;
            real mSwitchLambdaT3;
            real mSwitchLambdaT4;
            real mSwitchLambdaT5;
            real mSwitchLambdaT6;
            real mSwitchLambdaT7;
            real mSwitchLambdaT8;

            Vector <real> mYoungPoly0;
            Vector <real> mYoungPoly1;
            Vector <real> mYoungPoly2;

            Vector <real> mShearPoly0;
            Vector <real> mShearPoly1;
            Vector <real> mShearPoly2;

            const real mSwitchMechT0 = 295.0;
            const real mSwitchMechT1 = 1000.0;

            real mRRR = 100.0 ;
            real mTref = 273.0 ;
            real mRhoRef = BELFEM_QUIET_NAN ;

            // nist data for electric resistivity by means of Kohler's rule
            const Vector < real > mKohlerA = { 0.01827, -0.1839, 0.6229, 0.3168, -2.662 } ;

            // the polynomial polyval(mA, log10(x)) has a minumum around
            // x=0.588672. This is also dumb because the correction funciton
            // at this point goes to 1.002 rather than 1.000
            // instead, we find the value where the correction function is around
            // 1 percent
            const real mKohlerKmin = 0.01;
                  real mKohlerXmin ;

            Vector< real > mKohlerB ; // = { -0.00001983, 0.00031806, 0, 0  };

            // parameter for NIST thermal conductivity
            real mBeta ;
            Vector< real > mPlambda = { 2.443e-8, 1.754e-8, 2.763, 1102., -0.165, 70, 1.756, 0, 0.838, 0.0003, 0.1661 };
            const Vector< real > mPrho = { 1.553e-8, 1.171e-17, 4.49, 3.841e10, 1.14, 50.0,  6.428, 0.4531 };

            // help parameter to find peak value for lambda T
            const Vector< real > mQ = { 2.7015e-04, - 1.8831e-03, 1.2689e-02, 1.6394e-03};
            // value where dpolylva( mA, aY ) = 0

            real
            ( Copper::*mCpFunction )( const real aT ) const ;

            real
            ( Copper::*mLambda0Function )( const real aT ) const ;

            real
            ( Copper::*mRho0Function )( const real aT ) const ;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            Copper();

//----------------------------------------------------------------------------

            ~Copper();

//----------------------------------------------------------------------------

            /**
             * specific heat capacity in J/(kg*K)
             */
            real
            c( const real aT = BELFEM_TREF ) const;

//----------------------------------------------------------------------------

            /**
             * thermal conductivity in W/(m*K)
             */
            real
            lambda( const real aT = BELFEM_TREF ) const;

//----------------------------------------------------------------------------

            /**
             * thermal conductivity in W/(m*K)
             */
            real
            lambda( const real aT, const real aB ) const;

//----------------------------------------------------------------------------

            /**
             * Young's Modulus in GPa
             */
            real
            E( const real aT = BELFEM_TREF ) const;

//----------------------------------------------------------------------------

            /**
             * Shear Modulus in GPa
             */
            real
            G( const real aT = BELFEM_TREF ) const;

//----------------------------------------------------------------------------

            /**
             * set the purity of this material
             */
            void
            set_rrr( const real aRRR ) ;

//----------------------------------------------------------------------------

            // electric resistance
            real
            rho_el ( const real aJ=0, const real aT=0, const real aB=0 ) const ;

//----------------------------------------------------------------------------

            void
            use_splines( const bool aSwitch ) ;

//----------------------------------------------------------------------------
        private:
//----------------------------------------------------------------------------

            void
            create_specific_heat_polys();

//----------------------------------------------------------------------------

            real
            c_poly( const real aT ) const;

//----------------------------------------------------------------------------

            real
            c_spline( const real aT ) const;

//----------------------------------------------------------------------------

            void
            create_specific_heat_spline();

//----------------------------------------------------------------------------

            void
            create_conductivity_polys();

//----------------------------------------------------------------------------

            void
            create_density_poly();

//----------------------------------------------------------------------------

            void
            create_mech_polys();

//----------------------------------------------------------------------------

            void
            create_expansion_polys();

//----------------------------------------------------------------------------

            void
            create_resistivity_polys();

//----------------------------------------------------------------------------

            void
            create_resistivity_spline();

//----------------------------------------------------------------------------

            /**
             * resistivity for zero field
             */
            real
            rho_el0( const real aT ) const;

//----------------------------------------------------------------------------

            real
            lambda0_nist( const real aT ) const;

//----------------------------------------------------------------------------

            real
            lambda_poly( const real aT ) const ;

//----------------------------------------------------------------------------

            /**
             * resistivity for zero field
             */
            real
            rho_el0_nist( const real aT ) const;

//----------------------------------------------------------------------------

            /**
             * polynomial created from nist dataset
             */
            real
            rho_el0_poly( const real aT ) const;

//----------------------------------------------------------------------------

            /**
             * polynomial created from nist dataset
             */
            real
            rho_el0_spline( const real aT ) const;

//----------------------------------------------------------------------------

            /**
             * computes the temperature where the thermal conductivity peaks
             */
            real
            compute_T_lambda_peak() ;

//----------------------------------------------------------------------------
        };
//----------------------------------------------------------------------------

        inline real
        Copper::rho_el0_spline( const real aT ) const
        {
            return mRhoSpline->eval( aT );
        }


//----------------------------------------------------------------------------

        inline real
        Copper::c_spline( const real aT ) const
        {
            return mCpSpline->eval( aT );
        }

//----------------------------------------------------------------------------

        inline real
        Copper::c( const real aT ) const
        {
            return ( this->*mCpFunction )( aT ) ;
        }

//----------------------------------------------------------------------------

        inline real
        Copper::rho_el0( const real aT ) const
        {
            return ( this->*mRho0Function )( aT );
        }

//----------------------------------------------------------------------------

        inline real
        Copper::lambda( const real aT ) const
        {
            return ( this->*mLambda0Function )( aT );
        }

//----------------------------------------------------------------------------
    } /* end namespace material */
}  /* end namespace belfem */

#endif //BELFEM_CL_COPPER_HPP
