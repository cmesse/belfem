//
// Created by Christian Messe on 17.08.20.
//
#include "cl_GM_EoS_Hydrogen.hpp"
#include "assert.hpp"
#include "constants.hpp"
#include "fn_polyfit.hpp"
#include "fn_linspace.hpp"

namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        EoS_Hydrogen::EoS_Hydrogen( Gas & aParent, const HelmholtzModel aModel ) :
            Helmholtz( aParent, "H2" )
        {
            this->select_table( aModel );

            this->set_reference_point() ;
        }

//----------------------------------------------------------------------------

        EoS_Hydrogen::~EoS_Hydrogen()
        {
            this->delete_cubic_eos();
        }

//----------------------------------------------------------------------------

        void
        EoS_Hydrogen::select_table( const HelmholtzModel aModel )
        {
            mD = { 1.0 , 4.0 , 1.0 , 1.0 , 2.0 , 2.0 , 3.0 , 1.0 , 3.0 , 2.0 , 1.0 , 3.0 , 1.0 , 1.0 };
            switch ( aModel )
            {
                case( HelmholtzModel::ParaHydrogen ) :
                {
                    mTtriple = 13.8033 ;

                    mTcrit   = 32.938 ;
                    mPcrit   = 1.2858e6 ;
                    mRhocrit = mM * 15.538e3 ;
                    mVcrit   = 1.0 / mRhocrit ;

                    mA     = { -1.448589113, 1.884521239,  4.30256,  13.0289,  -47.7365,  50.0013,  -18.6261,  0.993973,  0.536078 };

                    mB     = { 0.0,  0.0,  -15.14967515,  -25.09259821,  -29.47355638,  -35.40591414,  -40.72499848,  -163.79258,  -309.2173174 };

                    mN     = { -7.33375,  0.01,  2.60375,  4.66279,  0.68239, -1.47078,
                               0.135801, -1.05327,  0.328239, -0.0577833,  0.0449743,
                               0.0703464, -0.0401766,  0.11951e0};

                    mT     = { 0.6855,  1.0,  1.0,  0.489,  0.774,  1.133,  1.386,
                               1.619,  1.162,  3.96,  5.276,  0.99,  6.791,  3.19e0};

                    mAlpha   = { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                               -1.7437, -0.5516, -0.0634, -2.1341, -1.777e0};

                    mBeta  = { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                               -0.194, -0.2019, -0.0301, -0.2383, -0.3253 };

                    mGamma = { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                               0.8048,  1.5248,  0.6648,  0.6832,  1.493e0};

                    mPsi   = { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                               1.5487,  0.1785,  1.28,  0.6319,  1.7104e0};

                    mNvap = { -4.87767, 1.03359, 0.82668, -0.129412 };
                    mKvap = { 1.0, 1.5, 2.65, 7.4 };

                    break ;
                }
                case( HelmholtzModel::NormalHydrogen ) :
                {
                    mTtriple = 13.957 ;

                    mTcrit = 33.145 ;
                    mPcrit = 1.2964e6;
                    mRhocrit = mM * 15.508e3 ;

                    mVcrit   = 1.0 / mRhocrit ;

                    mA     = { -1.457985648,  1.888076782,  1.616,  -0.4117,  -0.792,  0.758,  1.217 };

                    mB     = { 0.0,  0.0,  -16.02051591,  -22.6580178,  -60.00905114,  -74.94343038,  -206.9392065 };

                    mN     = { -6.93643,  0.01,  2.1101,  4.52059,  0.732564, -1.34086,  0.130985,-0.777414,
                                0.351944, -0.0211716,  0.0226312,  0.032187, -0.0231752,  0.0557346 };

                    mT     = { 0.6844,  1.0,  0.989,  0.489,  0.803,  1.1444,  1.409,
                               1.754,  1.311,  4.187,  5.646,  0.791,  7.249,  2.986e0};

                    mAlpha   = { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                               -1.685, -0.489, -0.103, -2.506, -1.607 };

                    mBeta  = { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                               -0.171, -0.2245, -0.1304, -0.2785, -0.3967 };

                    mGamma = { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                               0.7164,  1.3444,  1.4517,  0.7204,  1.5445 };

                    mPsi   = { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                               1.506,  0.156,  1.736,  0.67,  1.662e0};

                    mNvap  = { -4.89789, 0.988558, 0.349689, 0.499356 };
                    mKvap  = { 1.0, 1.5, 2.0, 2.85 };

                    break ;
                }
                case( HelmholtzModel::OrthoHydrogen ):
                {
                    mTtriple = 14.008 ;

                    mTcrit = 33.22 ;
                    mPcrit = 1.31065e6 ;
                    mRhocrit = mM * 15.445e3 ;
                    mVcrit   = 1.0 / mRhocrit ;

                    mA     = { -1.467544234,  1.884506886,  2.54151,  -2.3661,  1.00365,  1.22447 };

                    mB     = { 0.0,  0.0,  -25.76760987,  -43.46779049,  -66.04455148,  -209.7531607 };

                    mN     = { -6.83148,  0.01,  2.11505,  4.38353,  0.211292, -1.00939,  0.142086,
                               -0.87696,  0.804927, -0.710775,  0.0639688,  0.0710858, -0.087654,  0.647088e0};

                    mT     = { 0.7333,  1.0,  1.1372,  0.5136,  0.5638,  1.6248,  1.829,  2.404,
                               2.105,  4.1,  7.658,  1.259,  7.589,  3.946e0};

                    mAlpha   = { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                               -1.169, -0.894, -0.04, -2.072, -1.306 };

                    mBeta  = { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                               -0.4555, -0.4046, -0.0869, -0.4415, -0.5743e0};

                    mGamma = { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                               1.5444,  0.6627,  0.763,  0.6587,  1.4327e0};

                    mPsi   = { 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                               0.6366,  0.3876,  0.9437,  0.3976,  0.9626 };

                    mNvap  = { -4.88684, 1.05310, 0.856947, -0.185355 };
                    mKvap  = { 1.0, 1.5, 2.7, 6.2 } ;

                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false , "Invalid Model for class Helmholtz_Hydrogen" );
                }
            }

            // synch data with data object
            this->set_critical_point_in_data_object();

            mNab = mA.length() ;

            this->init_Tvap_poly() ;

            // maximum temperature for this model
            mTmax = 1000.0 ;

            // initial guess poly for liquid specific volume
            mVliq = { 7.552060e-6, -8.277210e-5, 1.266180e-2 };

            // init power containers
            mDeltaPowD.set_size( 15, BELFEM_QUIET_NAN );
            mTauPowT.set_size( 15, BELFEM_QUIET_NAN );

            mE.set_size( mNab, BELFEM_QUIET_NAN );
            mF.set_size( 16, BELFEM_QUIET_NAN );
            mG.set_size( 16, BELFEM_QUIET_NAN );

        }

//----------------------------------------------------------------------------

        real
        EoS_Hydrogen::compute_phi0()
        {
            this->update_e() ;

            // equation (31)
            real aValue = std::log( mDelta ) + 1.5*std::log( mTau )
                        + mA( 0 ) + mA( 1 ) * mTau ;

            for( uint k=2; k<mNab; ++k )
            {
                aValue += mA( k ) * std::log( 1.0 - mE( k ) );
            }
            return aValue ;
        }

//----------------------------------------------------------------------------

        real
        EoS_Hydrogen::compute_phir()
        {
            this->update_f();
            this->update_g() ;

            real aValue = 0.0 ;

            for( uint k=0; k<7; ++k )
            {
                aValue += mF( k );
            }
            for( uint k=7; k<9; ++k )
            {
                aValue += mF( k ) * mG( k );
            }
            for( uint k=9; k<14; ++k )
            {
                aValue += mF( k ) * mG( k );
            }
            return aValue ;
        }

//----------------------------------------------------------------------------

        real
        EoS_Hydrogen::compute_phi0_t()
        {
            this->update_e() ;

            // equation (31)
            real aValue = 1.5 / mTau + mA( 1 ) ;

            for( uint k=2; k<mNab; ++k )
            {
                aValue += mA( k ) * mB( k ) * mE( k ) / ( mE( k ) - 1.0 );
            }

            return aValue ;
        }

//----------------------------------------------------------------------------

        real
        EoS_Hydrogen::compute_phi0_tt()
        {
            this->update_e() ;

            // equation (31)
            real aValue = -1.5 / ( mTau * mTau );

            for( uint k=2; k<mNab; ++k )
            {
                aValue -= mA( k ) * mE( k ) * std::pow( mB( k ) / ( mE( k ) - 1 ) , 2 );
            }

            return aValue ;
        }

//----------------------------------------------------------------------------

        real
        EoS_Hydrogen::compute_phir_d()
        {
            this->update_f();
            this->update_g() ;

            real aValue = 0.0 ;

            for( uint k=0; k<7 ; ++k )
            {
                aValue += mF( k ) * mD( k ) ;
            }

            for( uint k=7; k<9; ++k )
            {
                aValue += mF( k ) * mG( k ) * ( mD( k ) - mDelta )  ;
            }

            for( uint k=9; k<14; ++k )
            {
                aValue +=   mF( k ) * mG( k )
                          * ( 2.0 * mAlpha( k ) * mDelta * ( mDelta - mPsi( k ) )
                              + mD( k ) ) ;
            }

            return aValue / mDelta ;
        }

//----------------------------------------------------------------------------

        real
        EoS_Hydrogen::compute_phir_dd()
        {
            this->update_f();
            this->update_g() ;

            real aValue = 0.0 ;

            for( uint k=0; k<7 ; ++k )
            {
                aValue += mF( k ) * ( mD( k ) * ( mD( k ) - 1.0 ) ) ;
            }

            for( uint k=7; k<9 ; ++k )
            {
                aValue += mF( k ) * mG( k ) *
                        ( std::pow(( mD( k ) -mDelta ), 2) - mD( k ) ) ;
            }


            for( uint k=9; k<14; ++k )
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
        EoS_Hydrogen::compute_phir_t()
        {
            this->update_f();
            this->update_g() ;

            real aValue = 0.0 ;

            for( uint k=0; k<7 ; ++k )
            {
                aValue += mF( k ) * mT( k );
            }

            for( uint k=7; k<9 ; ++k )
            {
                aValue += mF( k ) * mT( k ) * mG( k );
            }

            for( uint k=9; k<14; ++k )
            {
                aValue +=   mF( k ) * mG( k )
                            * ( 2.0 * mBeta( k ) * mTau * ( mTau - mGamma( k ) )
                                + mT( k ) ) ;
            }

            return aValue / mTau ;
        }

//----------------------------------------------------------------------------

        real
        EoS_Hydrogen::compute_phir_tt()
        {
            this->update_f();
            this->update_g() ;

            real aValue = 0.0 ;

            for( uint k=0; k<7 ; ++k )
            {
                aValue += mF( k ) * mT( k ) * ( mT( k ) - 1.0 );
            }

            for( uint k=7; k<9 ; ++k )
            {
                aValue += mF( k ) * mT( k ) * ( mT( k ) - 1.0 ) * mG( k );
            }

            real tA ;
            real tB ;
            real tC ;
            for( uint k=9; k<14; ++k )
            {
                tA = mT( k ) * ( mT( k ) - 4. * mBeta( k ) * mGamma( k ) * mTau - 1. );
                tB = 2. * mBeta( k ) * ( mGamma( k ) - mTau ) * ( mGamma( k ) - mTau )  + 1. + 2. * mT( k );
                tC = 2. * mBeta( k ) * mTau * mTau  ;

                aValue += mF( k ) * mG( k ) * ( tA + tB * tC );
            }

            return aValue / ( mTau * mTau ) ;
        }

//----------------------------------------------------------------------------

        real
        EoS_Hydrogen::compute_phir_dt()
        {
            this->update_f();
            this->update_g() ;

            real aValue = 0.0 ;

            for( uint k=0; k<7 ; ++k )
            {
                aValue += mF( k ) * mD( k ) * mT( k );
            }

            for( uint k=7; k<9 ; ++k )
            {
                aValue += mF( k ) * mG( k )
                          *  mT( k ) * ( mD( k ) - mDelta ) ;
            }

            for( uint k=9; k<14; ++k )
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