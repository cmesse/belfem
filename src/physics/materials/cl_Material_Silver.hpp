//
// Created by christian on 4/22/22.
//

#ifndef BELFEM_CL_MATERIAL_SILVER_HPP
#define BELFEM_CL_MATERIAL_SILVER_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

#include "cl_IsotropicMaterial.hpp"
#include "cl_Spline.hpp"

namespace belfem
{
    namespace material
    {
//----------------------------------------------------------------------------

        class Silver : public IsotropicMaterial
        {
            real mRRR = 206 ;

            real mTref = 273.0 ;

            Spline * mCpSpline = nullptr ;
            Spline * mRhoSpline = nullptr ;

            real mRhoRef ; // density at reference temperature

            const real mSwitchET0 = 300.0;
            const real mSwitchET1 = 600.0;
            Vector< real > mYoungPoly0 ;
            Vector< real > mYoungPoly1 ;
            Vector< real > mYoungPoly2 ;
            Vector< real > mNuPoly ;

            Vector <real> mSpecificHeatPoly0;
            Vector <real> mSpecificHeatPoly1;
            Vector <real> mSpecificHeatPoly2;
            Vector <real> mSpecificHeatPoly3;
            Vector <real> mSpecificHeatPoly4;
            Vector <real> mSpecificHeatPoly5;

            const real mSwitchCT0 = 8.0;
            const real mSwitchCT1 = 50.0;
            const real mSwitchCT2 = 80.0;
            const real mSwitchCT3 = 300.0;
            const real mSwitchCT4 = 500.0;

            const real mSwitchRT0 = 4.0;
            const real mSwitchRT1 = 10.0;
            const real mSwitchRT2 = 15.0;
            const real mSwitchRT3 = 50.0;
            const real mSwitchRT4 = 75.0;

            const real mSwitchLT0 = 50.0;
            const real mSwitchLT1 = 80.0;

            Vector< real > mResistivityPoly0 ;
            Vector< real > mResistivityPoly1 ;
            Vector< real > mResistivityPoly2 ;
            Vector< real > mResistivityPoly3 ;
            Vector< real > mResistivityPoly4 ;
            Vector< real > mResistivityPoly5 ;

            Vector <real> mThermalConductivityPoly1;
            Vector <real> mThermalConductivityPoly2;
            Vector <real> mThermalConductivityPoly3;

            const Vector< real > mKohlerA = {  0.0151959,
                                              -0.0931658,
                                              -0.0966798,
                                               2.31218,
                                              -4.36736 };
            Vector< real > mKohlerB ;
            const real mKohlerKmin = 0.001 ;
                  real mKohlerXmin ;

            const Vector< real > mPrho = { 1.467e-8,
                                           1.18777e-15,
                                           3.56474,
                                           3.18825e10,
                                           1.09259,
                                           49.7615,
                                           3.01754,
                                           -0.557109 };


            real
            ( Silver::*mCpFunction )( const real aT ) const ;

            real
            ( Silver::*mRho0Function )( const real aT ) const ;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            Silver() ;

//----------------------------------------------------------------------------

            ~Silver();

 //----------------------------------------------------------------------------

            /**
             * set the purity of this material
             */
            void
            set_rrr( const real aRRR ) ;

//----------------------------------------------------------------------------

            /**
             * specific heat capacity in J/(kg*K)
             */
            real
            c( const real aT = BELFEM_TREF ) const;

//----------------------------------------------------------------------------

            real
            E( const real aT= BELFEM_TREF ) const ;

//----------------------------------------------------------------------------

            real
            nu( const real aT=BELFEM_TREF ) const;

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
             * Shear Modulus in Pa
             */
            real
            G( const real aT=BELFEM_TREF ) const;

//----------------------------------------------------------------------------

            real
            rho_el ( const real aJ, const real aT, const real aB ) const ;

//----------------------------------------------------------------------------

            void
            use_splines( const bool aSwitch );

//----------------------------------------------------------------------------
        private:
//----------------------------------------------------------------------------


            void
            create_density_poly();

//----------------------------------------------------------------------------
            void
            create_mech_polys();

//----------------------------------------------------------------------------

            void
            create_specific_heat_polys();

//--------------------------------------------------------------------------

            void
            create_specific_heat_spline();

//--------------------------------------------------------------------------

            void
            create_resistivity_polys();

//--------------------------------------------------------------------------

            void
            create_resistivity_spline();

//----------------------------------------------------------------------------

            void
            create_conductivity_polys();

//----------------------------------------------------------------------------

            real
            rho_el0( const real aT ) const ;

//----------------------------------------------------------------------------

            real
            rho_el0_nist( const real aT ) const ;

//----------------------------------------------------------------------------

            real
            rho_el0_poly( const real aT ) const ;

//----------------------------------------------------------------------------

            real
            rho_el0_spline( const real aT ) const ;

//----------------------------------------------------------------------------

            real
            c_poly( const real aT ) const;

//----------------------------------------------------------------------------

            real
            c_spline( const real aT ) const;

//--------------------------------------------------------------------------
        };

//----------------------------------------------------------------------------

        inline real
        Silver::c( const real aT ) const
        {
            return ( this->*mCpFunction )( aT );
        }

//----------------------------------------------------------------------------

        inline real
        Silver::c_spline( const real aT ) const
        {
            return mCpSpline->eval( aT );
        }

//----------------------------------------------------------------------------

        inline real
        Silver::rho_el0_spline( const real aT ) const
        {
            return mRhoSpline->eval( aT );
        }

//----------------------------------------------------------------------------


        inline real
        Silver::rho_el0( const real aT ) const
        {
            return ( this->*mRho0Function )( aT );
        }

//----------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_MATERIAL_SILVER_HPP
