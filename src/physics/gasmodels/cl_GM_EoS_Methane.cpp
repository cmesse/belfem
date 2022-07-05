//
// Created by Christian Messe on 18.08.20.
//


#include "cl_GM_EoS_Methane.hpp"
#include "cl_Gas.hpp"
#include "cl_GM_Statevals.hpp"

namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        EoS_Methane::EoS_Methane( Gas & aParent ) :
                Helmholtz( aParent, "CH4" )
        {
            this->init_tables();
            //this->set_reference_point();
            this->set_reference_point( this->T_vap( 1.01325e5), 1.01325e5 );
        }

//----------------------------------------------------------------------------

        EoS_Methane::~EoS_Methane()
        {
            this->delete_cubic_eos();
        }

//----------------------------------------------------------------------------

        void
        EoS_Methane::init_tables()
        {
            // triple point
            mTtriple = 90.6941;

            // critical data
            mTcrit = 190.564;
            mPcrit = 4.5922e6;
            mRhocrit = 162.66;
            mVcrit = 1.0 / mRhocrit;

            // synch data with data object
            this->set_critical_point_in_data_object();

            // maximum temperature for this model
            mTmax = 625.0;

            // vapor pressure coefficients
            mNvap = { -6.036219, 1.409353, -0.4945199, -1.443048 };
            mKvap = { 1.0, 1.5, 2.0, 4.5 };

            // initial guess for liquid specific volume
            mVliq = { 4.188150e-8, -1.189574e-6, 1.978329e-3 };

            // vapor temperature, initial solution
            this->init_Tvap_poly();

            mA = { 9.91243972,  -6.33270087, 3.0016, 0.008449,
                   4.6942, 3.4865, 1.6572, 1.4115 };

            mB = { 0., 0., 0., -3.4004324, -10.26951575, -20.43932747,
                   -29.93744884, -79.13351945 };

            mN = { 4.367901028e-2, 6.709236199e-1, -1.765577859e0,
                   8.582330241e-1, -1.206513052e0, 5.120467220e-1,
                   -4.000010791e-4, -1.247842423e-2, 3.100269701e-2,
                   1.754748522e-3, -3.171921605e-6, -2.240346840e-6,
                   2.947056156e-7, 1.830487909e-1, 1.511883679e-1,
                   -4.289363877e-1, 6.894002446e-2, -1.408313996e-2,
                   -3.063054830e-2, -2.969906708e-2, -1.932040831e-2,
                   -1.105739959e-1, 9.952548995e-2, 8.548437825e-3,
                   -6.150555662e-2, -4.291792423e-2, -1.813207290e-2,
                   3.445904760e-2, -2.385919450e-3, -1.159094939e-2,
                   6.641693602e-2, -2.371549590e-2, -3.961624905e-2,
                   1.387292044e-2, 3.389489599e-2, -2.927378753e-3,
                   9.324799946e-5, -6.287171518e0, 1.271069467e1,
                   -6.423953466e0 };

            mD = { 1., 1., 1., 2., 2., 2., 2., 3., 4., 4., 8., 9., 10.,
                   1., 1., 1., 2., 4., 5., 6., 1., 2., 3., 4., 4., 3.,
                   5., 5., 8., 2., 3., 4., 4., 4., 5., 6., 2., 0., 0., 0 };

            mT = { -0.5, 0.5, 1., 0.5, 1., 1.5, 4.5, 0., 1., 3., 1., 3., 3.,
                   0., 1., 2., 0., 0., 2., 2., 5., 5., 5., 2., 4., 12., 8.,
                   10., 10., 10., 14., 12., 18., 22., 18., 14., 2., 0., 1., 2. };

            mC = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1.,
                   1., 1., 1., 1., 1., 2., 2., 2., 2., 2., 3., 3., 3., 3., 4.,
                   4., 4., 4., 4., 4., 4., 0., 0., 0., 0 };

            mAlpha = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                       0., 0., 0., 0., 0., 0., 0., 0., -20., -40., -40., -40. };

            mBeta = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                      0., 0., 0., 0., 0., 0., -200., -250., -250., -250. };

            mPsi = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                     0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                     0., 0., 0., 0., 0., 0., 0., 0., 1.07, 1.11, 1.11, 1.11 };

            mGamma = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                       0., 0., 0., 0., 0., 0., 1., 1., 1., 1 };

            mDeltaPowD.set_size( 41, BELFEM_QUIET_NAN );
            mTauPowT.set_size( 41, BELFEM_QUIET_NAN );
            mE.set_size( 8, BELFEM_QUIET_NAN );
            mF.set_size( 42, BELFEM_QUIET_NAN );
            mG.set_size( 42, BELFEM_QUIET_NAN );
        }

//----------------------------------------------------------------------------

        real
        EoS_Methane::compute_phi0()
        {
            this->update_e() ;

            return std::log( mDelta )
                    + mA( 0 ) + mA( 1 ) * mTau
                    + mA( 2 ) * std::log( mTau )
                    + mA( 3 ) * std::log( 1.0 - mE( 3 ) )
                    + mA( 4 ) * std::log( 1.0 - mE( 4 ) )
                    + mA( 5 ) * std::log( 1.0 - mE( 5 ) )
                    + mA( 6 ) * std::log( 1.0 - mE( 6 ) )
                    + mA( 7 ) * std::log( 1.0 - mE( 7 ) ) ;
        }

//----------------------------------------------------------------------------

        real
        EoS_Methane::compute_phi0_t()
        {
            this->update_e() ;


            return  mA( 1 ) + mA( 2 ) / mTau
                 + mA( 3 ) * mB( 3 ) * mE( 3 ) / ( mE( 3 ) - 1.0 )
                 + mA( 4 ) * mB( 4 ) * mE( 4 ) / ( mE( 4 ) - 1.0 )
                 + mA( 5 ) * mB( 5 ) * mE( 5 ) / ( mE( 5 ) - 1.0 )
                 + mA( 6 ) * mB( 6 ) * mE( 6 ) / ( mE( 6 ) - 1.0 )
                 + mA( 7 ) * mB( 7 ) * mE( 7 ) / ( mE( 7 ) - 1.0 );
        }

//----------------------------------------------------------------------------

        real
        EoS_Methane::compute_phi0_tt()
        {
            this->update_e() ;

            real aValue = -mA( 2 ) / ( mTau * mTau );
            for ( uint k = 3; k < 8; ++k )
            {
                aValue -= mA( k ) * mE( k ) * std::pow( mB( k ) / ( mE( k ) - 1.0) , 2 );
            }

            return aValue ;
        }

//----------------------------------------------------------------------------

        real
        EoS_Methane::compute_phir()
        {
            this->update_f() ;
            this->update_g() ;

            real aValue = 0.0 ;

            for( uint k=0; k<13; ++k )
            {
                aValue += mF( k );
            }

            for( uint k=13; k<40; ++k )
            {
                aValue += mF( k ) * mG( k );
            }

            return aValue ;
        }

//----------------------------------------------------------------------------

        real
        EoS_Methane::compute_phir_d()
        {
            this->update_f();
            this->update_g() ;

            real aValue = 0.0 ;

            for( uint k=0; k<13 ; ++k )
            {
                aValue += mF( k ) * mD( k ) ;
            }

            for( uint k=13; k<36; ++k )
            {
                aValue += mF( k ) * mG( k )
                        * ( mD( k ) - mC( k )  * std::pow( mDelta, mC( k ) ) );
            }

            for( uint k=36; k<40; ++k )
            {
                aValue +=   mF( k ) * mG( k )
                            * ( 2.0 * mAlpha( k ) * mDelta * ( mDelta - mPsi( k ) )
                                + mD( k ) ) ;
            }

            return aValue / mDelta ;
        }

//----------------------------------------------------------------------------

        real
        EoS_Methane::compute_phir_dd()
        {
            this->update_f();
            this->update_g() ;

            real aValue = 0.0 ;

            for( uint k=0; k<13 ; ++k )
            {
                aValue += mF( k ) * ( mD( k ) * ( mD( k ) - 1.0 ) ) ;
            }

            real tDeltaPowC ;

            for( uint k=13; k<36 ; ++k )
            {
                tDeltaPowC = std::pow( mDelta, mC( k ) ) ;

                aValue += mF( k ) * mG( k ) *
                     ( mD( k ) * ( mD( k ) - 1.0 )
                      - tDeltaPowC * mC( k ) * ( 2.0 * mD( k )
                      + mC( k ) - mC( k ) * tDeltaPowC - 1.0 ) );

            }

            for( uint k=36; k<40; ++k )
            {
                aValue += mF( k ) * mG( k ) *
                          ( 2.0 * mDelta * mDelta * mAlpha( k ) *
                            ( 1.0 + 2.0 * mAlpha(k)* std::pow( mDelta - mPsi( k ) , 2 ) )
                            + mD( k ) * ( 4.0 * mDelta * mAlpha( k ) * ( mDelta - mPsi( k ) ) - 1.0 )
                            + mD( k ) * mD( k ) ) ;
            }

            return aValue / ( mDelta * mDelta ) ;
        }

//----------------------------------------------------------------------------

        real
        EoS_Methane::compute_phir_t()
        {
            this->update_f() ;
            this->update_g() ;

            real aValue = 0.0 ;

            for( uint k=0; k<13 ; ++k )
            {
                aValue += mF( k ) * mT( k );
            }

            for( uint k=13; k<36 ; ++k )
            {
                aValue += mF( k ) * mT( k ) * mG( k );
            }

            for( uint k=36; k<40; ++k )
            {
                aValue +=   mF( k ) * mG( k )
                            * ( 2.0 * mBeta( k ) * mTau * ( mTau - mGamma( k ) )
                                + mT( k ) ) ;
            }

            return aValue / mTau ;
        }

//----------------------------------------------------------------------------

        real
        EoS_Methane::compute_phir_tt()
        {
            this->update_f() ;
            this->update_g() ;

            real aValue = 0.0 ;

            for( uint k=0; k<13 ; ++k )
            {
                aValue += mF( k ) * mT( k ) * ( mT( k ) - 1.0 );
            }

            for( uint k=13; k<36 ; ++k )
            {
                aValue += mF( k ) * mT( k ) * ( mT( k ) - 1.0 ) * mG( k );
            }

            aValue /= mTau ;

            for( uint k=36; k<40; ++k )
            {
                aValue +=   mF( k ) * mG( k ) * mTau *
                        ( 2.0 * mBeta( k ) * mTau * ( mTau - mGamma ( k ) ) + mT( k ) );
            }

            return aValue / mTau ;
        }

//----------------------------------------------------------------------------

        real
        EoS_Methane::compute_phir_dt()
        {
            this->update_f() ;
            this->update_g() ;

            real aValue = 0.0 ;

            for( uint k=0; k<13 ; ++k )
            {
                aValue += mF( k ) * mD( k ) * mT( k );
            }

            for( uint k=13; k<36 ; ++k )
            {
                aValue += mF( k ) * mG( k ) *  mT( k )
                        * ( mD( k ) - mC( k ) * std::pow( mDelta, mC( k ) ) ) ;
            }

            for( uint k=36; k<40; ++k )
            {
                aValue += mF( k ) * mG( k )
                          * ( 2.0 * mAlpha( k ) * mDelta
                              * ( mDelta - mPsi( k ) ) + mD( k ) )
                          * ( 2.0 * mBeta( k ) * mTau
                              * ( mTau - mGamma ( k ) ) + mT( k ) );
            }

            return aValue / ( mDelta * mTau );
        }

//----------------------------------------------------------------------------
    }
}