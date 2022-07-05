//
// Created by Christian Messe on 18.08.20.
//

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Vector.hpp"
#include "fn_r2.hpp"
#include "cl_Gas.hpp"

#define private public
#define protected public
#include "cl_GM_EoS_Hydrogen.hpp"
#undef protected
#undef private

using namespace belfem;
using namespace belfem::gastables;
using namespace belfem::gasmodels;

TEST( GASMODELS, Hydrogen_Vapor )
{
    // crate reference mixture
    Gas tRef( "H2" );

    for( uint g=0; g<2; ++g )
    {
        HelmholtzModel tModel = ( HelmholtzModel ) g;
        // create EoS
        EoS_Hydrogen tH2( tRef, tModel );

//----------------------------------------------------------------------------
//  Vapor Pressure
//----------------------------------------------------------------------------

        /*
         * Table 11 A, Chapter 11 "Technology and
         * Uses of Liquid Hydrogen", 381, Pergamon Press
         */
        Vector< real > tTvap;
        Vector< real > tPvap;

        switch ( tModel )
        {
            case ( HelmholtzModel::ParaHydrogen ) :
            {
                tTvap = { 13.803, 14., 15., 16., 17., 18., 19., 20., 20.268, 21.,
                          22., 23., 24., 25., 26., 27., 28., 29., 30., 31., 32.,
                          32.976 };
                tPvap = { 0.069, 0.077, 0.132, 0.212, 0.325, 0.475, 0.672, 0.922,
                          1.000, 1.233, 1.612, 2.069, 2.610, 3.245, 3.981, 4.828,
                          5.793, 6.886, 8.117, 9.500, 11.051, 12.759 };
                break;
            }
            case ( HelmholtzModel::NormalHydrogen ) :
            {
                tTvap = { 13.947, 14., 15., 16., 17., 18., 19., 20., 20.38,
                          21., 22., 23., 24., 25., 26., 27., 28., 29., 30.,
                          31., 32., 33., 33.18 };

                tPvap = { 0.071, 0.073, 0.125, 0.202, 0.310, 0.456, 0.648, 0.891,
                          1.000, 1.196, 1.569, 2.018, 2.551, 3.178, 3.906, 4.746,
                          5.705, 6.794, 8.023, 9.410, 10.94, 12.65, 12.98 };
                break;
            }
            case ( HelmholtzModel::OrthoHydrogen ) :
            {
                // note, these data are not from the literature, but precomputed from model
                tTvap = { 14.008, 14., 15., 16., 17., 18., 19., 20., 20.38,
                          21., 22., 23., 24., 25., 26., 27., 28., 29., 30.,
                          31., 32., 33., 33.21 };

                tPvap = { 0.075, 0.074, 0.127, 0.204, 0.313, 0.458, 0.649, 0.892, 1.000,
                          1.195, 1.566, 2.012, 2.542, 3.163, 3.884, 4.714, 5.661, 6.734,
                          7.944, 9.301, 10.821, 12.528, 12.916 };
                break;
            }
            default:
            {
                break;
            }
        }
        // convert Units
        tPvap *= 101325;


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Test function p_vap( T_vap )
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        uint tNumPoints = tTvap.length();

        Vector< real > tY( tNumPoints );
        Vector< real > tY0( tNumPoints );

        for ( uint k = 0; k < tNumPoints; ++k )
        {
            tY0( k ) = std::log( tPvap( k ));
            tY( k ) = std::log( tH2.p_vap( tTvap( k )));
        }

        EXPECT_NEAR( r2( tY, tY0 ), 1.0, 1e-4 );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Test function dp_vap/dT( T_vap )
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // test derivative
        for ( uint k = 0; k < tNumPoints; ++k )
        {
            tY0( k ) = ( tH2.p_vap( tTvap( k ) + 0.01 ) - tH2.p_vap( tTvap( k ) - 0.01 )) / 0.02;
            tY( k ) = tH2.dpvap_dT( tTvap( k ), tH2.p_vap( tTvap( k ) ) , tH2.pi_vap( tTvap( k ) ) );
        }
        EXPECT_NEAR( r2( tY, tY0 ), 1.0, 1e-5 );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Test function T_vap(p_vap) )
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // test derivative
        for ( uint k = 0; k < tNumPoints; ++k )
        {
            tY( k ) = tH2.T_vap( tH2.p_vap( tTvap( k )));
        }
        EXPECT_NEAR( r2( tY, tTvap ), 1.0, 1e-5 );

    }
}