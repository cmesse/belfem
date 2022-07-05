//
// Created by Christian Messe on 24.08.20.
//

#include "fn_GM_Helmholz_DerivTest.hpp"
#include "fn_linspace.hpp"
#include "fn_r2.hpp"
namespace belfem
{
    namespace gasmodels
    {
//----------------------------------------------------------------------------

        void
        deriv_test( Helmholtz & aGas, Vector< real > & aR2 )
        {
            aR2.set_size( 7 );

            // number of points
            uint tNa = 50 ;
            uint tNb = 150 ;
            uint tN = tNa + tNb ;

            // create temperature array
            Vector< real > tTa( tNa ) ;
            Vector< real > tTb( tNb ) ;
            Vector< real > tT( tN );

            // reference pressure for test
            real tP = 1e6 ;


            const real tXi1 = 1.0001 ;
            const real tXi2 = 0.9999 ;

            // create temperature table
            linspace( 1.01 * aGas.T_min(), 0.99 * aGas.T_vap( tP ), tNa, tTa );
            linspace( 1.01 * aGas.T_vap( tP ), 0.99 * aGas.T_max(), tNb, tTb );

            for( uint k=0; k<tNa; ++k )
            {
                tT( k ) = tTa( k );
            }
            uint tCount = tNa ;
            for( uint k=0; k<tNb; ++k )
            {
                tT( tCount++ ) = tTb( k );
            }

            // vector with reference result
            Vector< real > tY0( tN );

            // vector with computed result
            Vector< real > tY( tN ) ;

            real tY1 ;
            real tY2 ;
            real tV ;
            real tX1 ;
            real tX2 ;

            const real tTcrit = aGas.T_crit() ;
            const real tVcrit = aGas.v_crit() ;

            // - - - -  - - - - - - - - - - - - - - - - - - -

            // check phi0_t
            for( uint k=0; k<tN; ++k )
            {
                tV = aGas.v( tT( k ), tP );
                tY( k ) = aGas.phi0_t() ;

                aGas.update_Tv( tT( k )*tXi1, tV ) ;
                tX1 = tTcrit / ( tXi1 * tT( k ) );
                tY1 = aGas.phi0() ;

                aGas.update_Tv( tT( k )*tXi2, tV ) ;
                tX2 = tTcrit / ( tXi2 * tT( k ) );
                tY2 = aGas.phi0() ;

                tY0( k ) = ( tY2 - tY1 ) / ( tX2 - tX1 );

            }
            aR2( 0 ) = r2( tY, tY0 );

            // - - - -  - - - - - - - - - - - - - - - - - - -

            // check phi0_tt
            for( uint k=0; k<tN; ++k )
            {
                tV = aGas.v( tT( k ), tP );
                tY( k ) = aGas.phi0_tt() ;

                aGas.update_Tv( tT( k )*tXi1, tV ) ;
                tX1 = tTcrit / ( tXi1 * tT( k ) );
                tY1 = aGas.phi0_t() ;

                aGas.update_Tv( tT( k )*tXi2, tV) ;
                tX2 = tTcrit / ( tXi2 * tT( k ) );
                tY2 = aGas.phi0_t() ;
                tY0( k ) = ( tY2 - tY1 ) / ( tX2 - tX1 );

            }
            aR2( 1 ) = r2( tY, tY0 );

            // - - - -  - - - - - - - - - - - - - - - - - - -

            // check phir_t
            for( uint k=0; k<tN; ++k )
            {
                tV = aGas.v( tT( k ), tP );
                tY( k ) = aGas.phir_t() ;

                aGas.update_Tv( tT( k )*tXi1, tV ) ;
                tX1 = tTcrit / ( tXi1 * tT( k ) );
                tY1 = aGas.phir() ;

                aGas.update_Tv( tT( k )*tXi2, tV ) ;
                tX2 = tTcrit / ( tXi2 * tT( k ) );
                tY2 = aGas.phir() ;

                tY0( k ) = ( tY2 - tY1 ) / ( tX2 - tX1 );

            }
            aR2( 2 ) = r2( tY, tY0 );

            // - - - -  - - - - - - - - - - - - - - - - - - -

            // check phir_d
            for( uint k=0; k<tN; ++k )
            {
                tV = aGas.v( tT( k ), tP );
                tY( k ) = aGas.phir_d() ;

                aGas.update_Tv( tT( k ), tXi1 * tV );
                tX1 = tVcrit / ( tXi1 * tV );
                tY1 = aGas.phir() ;

                aGas.update_Tv( tT( k ), tXi2 * tV );
                tX2 = tVcrit / ( tXi2 * tV );
                tY2 = aGas.phir() ;

                tY0( k ) = ( tY2 - tY1 ) / ( tX2 - tX1 );

            }
            aR2( 3 ) = r2( tY, tY0 );

            // - - - -  - - - - - - - - - - - - - - - - - - -

            // check phir_tt
            for( uint k=0; k<tN; ++k )
            {
                tV = aGas.v( tT( k ), tP );
                tY( k ) = aGas.phir_tt() ;

                aGas.update_Tv( tT( k )*tXi1, tV ) ;
                tX1 = tTcrit / ( tXi1 * tT( k ) );
                tY1 = aGas.phir_t() ;

                aGas.update_Tv( tT( k )*tXi2, tV ) ;
                tX2 = tTcrit / ( tXi2 * tT( k ) );
                tY2 = aGas.phir_t() ;

                tY0( k ) = ( tY2 - tY1 ) / ( tX2 - tX1 );
            }
            aR2( 4 ) = r2( tY, tY0 );

            // - - - -  - - - - - - - - - - - - - - - - - - -
            // check phir_dt
            for( uint k=0; k<tN; ++k )
            {
                tV = aGas.v( tT( k ), tP );
                tY( k ) = aGas.phir_dt() ;

                aGas.update_Tv( tT( k ), tXi1 * tV );
                tX1 = tVcrit / ( tXi1 * tV );
                tY1 = aGas.phir_t() ;

                aGas.update_Tv( tT( k ), tXi2 * tV );
                tX2 = tVcrit / ( tXi2 * tV );
                tY2 = aGas.phir_t() ;

                tY0( k ) = ( tY2 - tY1 ) / ( tX2 - tX1 );

            }
            aR2( 5 ) = r2( tY, tY0 );

            // - - - -  - - - - - - - - - - - - - - - - - - -
            // check phir_dd
            for( uint k=0; k<tN; ++k )
            {
                tV      = aGas.v( tT( k ), tP );
                tY( k ) = aGas.phir_dd() ;


                aGas.update_Tv( tT( k ), tXi1 * tV );
                tX1 = tVcrit / ( tXi1 * tV );
                tY1 = aGas.phir_d() ;

                aGas.update_Tv( tT( k ), tXi2 * tV );
                tX2 = tVcrit / ( tXi2 * tV );
                tY2 = aGas.phir_d() ;

                tY0( k ) = ( tY2 - tY1 ) / ( tX2 - tX1 );

            }
            aR2( 6 ) = r2( tY, tY0 );


        }

//----------------------------------------------------------------------------
    }
}