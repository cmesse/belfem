//
// Created by Christian Messe on 18.08.20.
//

#include "cl_GM_EoS_Oxygen.hpp"
#include "cl_Gas.hpp"
#include "cl_GM_Statevals.hpp"

namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        EoS_Oxygen::EoS_Oxygen( Gas & aParent ) :
            Helmholtz( aParent, "O2" )
        {
            this->init_tables() ;
            this->set_reference_point() ;
        }

        //----------------------------------------------------------------------------

        EoS_Oxygen::~EoS_Oxygen()
        {
            this->delete_cubic_eos();
        }


//----------------------------------------------------------------------------

        void
        EoS_Oxygen::init_tables()
        {
            // triple point
            mTtriple = 54.33 ;

            // critical data
            mTcrit   = 154.581 ;
            mPcrit   = 5.043e6 ;
            mRhocrit = mM * 13.63e3 ;
            mVcrit   = 1.0 / mRhocrit ;

            // synch data with data object
            this->set_critical_point_in_data_object();

            // reference conditions
            mTau0   = mTcrit / mT0 ;
            mDelta0 = mP0 / ( mR * mT0 * mRhocrit ) ;

            mTmax = 300.0 ;

            /* polynomial for vapor pressure
             * based on function from doi: 10.1007/BF02706097
             * but adapted to current critical data
             */
            mNvap = { -2.410514, 1.815051, - 5.812213 } ;

            // initial guess for liquid specific volume
            mVliq = { 1.764476e-08, 5.325537e-7, 6.845122e-4 } ;

            // vapor temperature, initial solution
            this->init_Tvap_poly() ;


            // values for ideal gas function
            mK = { -7.40775e-4, -6.64930e-5, 2.50042e0, -2.14487e1,
                    1.01258e0, -9.44365e-1, 1.45066e1, 7.49148e1,
                    4.14817e0};

            // values for real gas function
            mD = { 1., 1., 1., 2., 2., 2., 3., 3., 3., 6., 7., 7., 8., 1., 1.,
                   2., 2., 3., 3., 5., 6., 7., 8., 10., 2., 3., 3., 4., 4.,
                   5., 5., 5 };

            mT = { 0., 1.5, 2.5, -0.5, 1.5, 2., 0., 1., 2.5, 0., 2., 5., 2., 5.,
                   6., 3.5, 5.5, 3., 7., 6., 8.5, 4., 6.5, 5.5, 22., 11., 18.,
                   11., 23., 17., 18., 23.0 };

            mN = {  3.983768749e-1, -1.846157454e0,   4.183473197e-1,  2.370620711e-2,
                    9.771730573e-2,  3.017891294e-2,  2.273353212e-2,  1.357254086e-2,
                   -4.052698943e-2,  5.454628515e-4,  5.113182277e-4,  2.953466883e-7,
                   -8.687645072e-5, -2.127082589e-1,  8.735941958e-2,  1.275509190e-1,
                   -9.067701064e-2, -3.540084206e-2, -3.623278059e-2,  1.327699290e-2,
                   -3.254111865e-4, -8.313582932e-3,  2.124570559e-3, -8.325206232e-4,
                   -2.626173276e-5,  2.599581482e-3,  9.984649663e-3,  2.199923153e-3,
                   -2.591350486e-2, -1.259630848e-1,  1.478355637e-1, -1.011251078e-2 };


            // allocate help containers
            mE.set_size( 4, BELFEM_REAL_MAX );

            mDeltaPowD.set_size( 33 );
            mTauPowT.set_size( 33 );
            mF.set_size( 34, BELFEM_QUIET_NAN );
            mG.set_size( 3, BELFEM_QUIET_NAN );

        }

//----------------------------------------------------------------------------

        real
        EoS_Oxygen::compute_phi0()
        {
            this->update_e() ;

            return  mK( 0 ) * mTau * mE( 0 )
                  + mK( 1 ) / ( mTau * mTau )
                  + mK( 2 ) * std::log( mTau )
                  + mK( 3 ) * mTau
                  + mK( 4 ) * std::log( mE( 1 ) -1.0 )
                  //+ mK( 5 ) * std::log( mE( 2 ) + 1.0)
                  + mK( 8 )
                  //+ ( mH0 * mTau / mTcrit - mS0 ) / mR
                  + std::log( mDelta / mDelta0 );
        }

//----------------------------------------------------------------------------

        real
        EoS_Oxygen::compute_phi0_t()
        {
            this->update_e() ;

            return       1.5 * mK( 0 ) * mE( 0 )
                     - ( 2.0 * mK( 1 ) / ( mTau * mTau )
                     -         mK( 2 ) ) / mTau
                     +         mK( 3 )
                     +         mK( 4 ) * mK( 6 ) * mE( 1 ) / ( mE( 1 ) - 1.0 )
                     +         mK( 5 ) * mK( 7 ) * mE( 2 ) / ( mE( 2 ) + 1.0 );
        }

//----------------------------------------------------------------------------

        real
        EoS_Oxygen::compute_phi0_tt()
        {
            this->update_e() ;

            return  0.75 * mK( 0 ) / mE( 0 )
                +  ( 6.0 * mK( 1 ) / mTau
                -          mK( 2 ) ) / ( mTau * mTau )
                -          mK( 4 ) * mE( 1 ) * std::pow( mK( 6 ) / ( mE( 1 ) - 1.0 ), 2 )
                +          mK( 5 ) * mE( 2 ) * std::pow( mK( 7 ) / ( mE( 2 ) + 1.0 ), 2 ) ;

        }

//----------------------------------------------------------------------------

        real
        EoS_Oxygen::compute_phir()
        {
            this->update_f() ;
            this->update_g() ;

            real tValue0 = 0.0 ;
            real tValue1 = 0.0 ;
            real tValue2 = 0.0 ;

            for( uint k=0; k<13; ++k )
            {
                tValue0 += mF( k );
            }

            for( uint k=13; k<24; ++k )
            {
                tValue1 += mF( k );
            }

            for( uint k=24; k<32; ++k )
            {
                tValue2 += mF( k );
            }

            return tValue0 + tValue1 * mG( 0 ) + tValue2 * mG( 1 ) ;
        }

//----------------------------------------------------------------------------

        real
        EoS_Oxygen::compute_phir_d()
        {
            this->update_f() ;
            this->update_g() ;

            real tValue0 = 0.0 ;
            real tValue1 = 0.0 ;
            real tValue2 = 0.0 ;

            real tDelta2 = mDelta * mDelta ;
            real tDelta4 = tDelta2 * tDelta2 ;

            for( uint k=0; k<13; ++k )
            {
                tValue0 += mF( k ) * mD( k );
            }

            for( uint k=13; k<24; ++k )
            {
                  tValue1 += mF( k ) * ( mD( k ) - 2.0 * tDelta2 );
            }

            for( uint k=24; k<32; ++k )
            {
                tValue2 += mF( k ) * ( mD( k ) - 4.0 * tDelta4 );
            }

            return ( tValue0  + tValue1 * mG( 0 ) + tValue2 * mG( 1 ) ) / mDelta ;
        }

//----------------------------------------------------------------------------

        real
        EoS_Oxygen::compute_phir_dd()
        {
            this->update_f() ;
            this->update_g() ;

            real tValue0 = 0.0 ;
            real tValue1 = 0.0 ;
            real tValue2 = 0.0 ;

            real tDelta2 = mDelta  * mDelta ;
            real tDelta4 = tDelta2 * tDelta2 ;

            for( uint k=0; k<13; ++k )
            {
                tValue0 += mF( k ) * mD( k ) * ( mD( k ) - 1.0 );
            }

            for( uint k=13; k<24; ++k )
            {
                tValue1 += mF( k ) * ( mD( k ) * ( mD( k ) - 1.0 )
                        + tDelta2 * ( 4.0 * ( tDelta2 - mD( k ) ) - 2.0 ) );
            }

            for( uint k=24; k<32; ++k )
            {
                tValue2 += mF( k ) * ( mD( k ) * ( mD( k ) - 1.0 )
                            + tDelta4 * ( 16.0 * tDelta4 - 8.0 * mD( k ) - 12.0 ) );
            }

            return ( tValue0 + tValue1 * mG( 0 ) + tValue2 * mG( 1 ) ) / tDelta2;
        }

//----------------------------------------------------------------------------

        real
        EoS_Oxygen::compute_phir_t()
        {
            this->update_f();
            this->update_g();


            real tValue0 = 0.0 ;
            real tValue1 = 0.0 ;
            real tValue2 = 0.0 ;

            for( uint k=0; k<13 ; ++k )
            {
                tValue0 += mF( k ) * mT( k );
            }

            for( uint k=13; k<24 ; ++k )
            {
                tValue1 += mF( k ) * mT( k ) ;
            }

            for( uint k=24; k<32 ; ++k )
            {
                tValue2 += mF( k ) * mT( k ) ;
            }

            return ( tValue0 + mG( 0 ) * tValue1 + mG( 1 ) * tValue2 ) / mTau ;
        }

//----------------------------------------------------------------------------

        real
        EoS_Oxygen::compute_phir_tt()
        {
            this->update_f();
            this->update_g();

            real tValue0 = 0.0 ;
            real tValue1 = 0.0 ;
            real tValue2 = 0.0 ;

            for( uint k=0; k<13 ; ++k )
            {
                tValue0 += mF( k ) * mT( k ) * ( mT( k ) - 1.0 );
            }

            for( uint k=13; k<24 ; ++k )
            {
                tValue1 += mF( k ) * mT( k ) * ( mT( k ) - 1.0 );
            }

            for( uint k=24; k<32 ; ++k )
            {
                tValue2 += mF( k ) * mT( k ) * ( mT( k ) - 1.0 );
            }

            return ( tValue0 + mG( 0 ) * tValue1 + mG( 1 ) * tValue2 )
                / ( mTau * mTau );
        }

//----------------------------------------------------------------------------

        real
        EoS_Oxygen::compute_phir_dt()
        {
            this->update_f();
            this->update_g();

            real tValue0 = 0.0 ;
            real tValue1 = 0.0 ;
            real tValue2 = 0.0 ;

            real tDelta2 = mDelta * mDelta ;
            real tDelta4 = tDelta2 * tDelta2 ;

            for( uint k=0; k<13 ; ++k )
            {
                tValue0 += mF( k ) * mD( k ) * mT( k );
            }


           for( uint k=13; k<24 ; ++k )
           {
               tValue1 += mF( k ) * mT( k ) * ( mD( k ) - 2.0 * tDelta2 ) ;
           }

           for( uint k=24; k<32 ; ++k )
           {
               tValue2 += mF( k ) * mT( k ) * ( mD( k ) - 4.0 * tDelta4 ) ;
           }

            return ( tValue0+ tValue1 * mG( 0 ) + tValue2 * mG( 1 ) )
                   / ( mTau * mDelta );
        }

//----------------------------------------------------------------------------

    }
}