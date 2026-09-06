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

// DOP853: Dormand-Prince 8(5,3) explicit Runge-Kutta method.
// Coefficients from Hairer, Norsett, Wanner: Solving Ordinary
// Differential Equations I, 2nd ed., Springer (1993), Section II.6.

#include "fn_ODE_DOP853.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace ode
    {
//------------------------------------------------------------------------------

        void
        DOP853_init( ODE & aODE, Cell< Vector< real > > & aWork )
        {
            // number of entries in Y-vector
            index_t tN = aODE.dimension();

            // k0-k11 (12 stages), ytemp, ya
            aWork.set_size( 14, {} );

            for( index_t k=0; k<14; ++k )
            {
                aWork( k ).set_size( tN, 0.0 );
            }
        }

//------------------------------------------------------------------------------

        Status
        DOP853(  ODE                    & aODE,
                 real                   & aT,
                 Vector< real >         & aY,
                 real                   & aStep,
                 Cell< Vector< real > > & aWork,
                 const real              aEpsilon,
                 const uint              aMaxIterations,
                 const real              aTmax,
                 const bool              aAutoTimestep )
        {
            // test length of work array
            BELFEM_ASSERT( aWork.size() >= 14, "Size of working cell does not match" );

            // stage derivatives (12 stages)
            Vector< real > & k0  = aWork(  0 );
            Vector< real > & k1  = aWork(  1 );
            Vector< real > & k2  = aWork(  2 );
            Vector< real > & k3  = aWork(  3 );
            Vector< real > & k4  = aWork(  4 );
            Vector< real > & k5  = aWork(  5 );
            Vector< real > & k6  = aWork(  6 );
            Vector< real > & k7  = aWork(  7 );
            Vector< real > & k8  = aWork(  8 );
            Vector< real > & k9  = aWork(  9 );
            Vector< real > & k10 = aWork( 10 );
            Vector< real > & k11 = aWork( 11 );

            // scratch buffer for intermediate Y values
            Vector< real > & yt  = aWork( 12 );

            // 8th order solution
            Vector< real > & ya  = aWork( 13 );

            Vector< real > & y0  = aY;

            // -- time nodes (c_i) --
            const real c2  = 0.526001519587677318785587544488e-01;
            const real c3  = 0.789002279381515978178381316732e-01;
            const real c4  = 0.118350341907227396726757197510e+00;
            const real c5  = 0.281649658092772603273242802490e+00;
            const real c6  = 0.333333333333333333333333333333e+00;
            const real c7  = 0.25e+00;
            const real c8  = 0.307692307692307692307692307692e+00;
            const real c9  = 0.651282051282051282051282051282e+00;
            const real c10 = 0.6e+00;
            const real c11 = 0.857142857142857142857142857142e+00;

            // -- Butcher tableau a_ij coefficients --
            const real a21   =  5.26001519587677318785587544488e-02;

            const real a31   =  1.97250569845378994544595329183e-02;
            const real a32   =  5.91751709536136983633785987549e-02;

            const real a41   =  2.95875854768068491816892993775e-02;
            const real a43   =  8.87627564304205475450678981324e-02;

            const real a51   =  2.41365134159266685502369798665e-01;
            const real a53   = -8.84549479328286085344864962717e-01;
            const real a54   =  9.24834003261792003115737966543e-01;

            const real a61   =  3.7037037037037037037037037037e-02;
            const real a64   =  1.70828608729473871279604482173e-01;
            const real a65   =  1.25467687566822425016691814123e-01;

            const real a71   =  3.7109375e-02;
            const real a74   =  1.70252211019544039314978060272e-01;
            const real a75   =  6.02165389804559606850219397283e-02;
            const real a76   = -1.7578125e-02;

            const real a81   =  3.70920001185047927108779319836e-02;
            const real a84   =  1.70383925712239993810214054705e-01;
            const real a85   =  1.07262030446373284651809199168e-01;
            const real a86   = -1.53194377486244017527936158236e-02;
            const real a87   =  8.27378916381402288758473766002e-03;

            const real a91   =  6.24110958716075717114429577812e-01;
            const real a94   = -3.36089262944694129406857109825e+00;
            const real a95   = -8.68219346841726006818189891453e-01;
            const real a96   =  2.75920996994467083049415600797e+01;
            const real a97   =  2.01540675504778934086186788979e+01;
            const real a98   = -4.34898841810699588477366255144e+01;

            const real a101  =  4.77662536438264365890433908527e-01;
            const real a104  = -2.48811461997166764192642586468e+00;
            const real a105  = -5.90290826836842996371446475743e-01;
            const real a106  =  2.12300514481811942347288949897e+01;
            const real a107  =  1.52792336328824235832596922938e+01;
            const real a108  = -3.32882109689848629194453265587e+01;
            const real a109  = -2.03312017085086261358222928593e-02;

            const real a111  = -9.3714243008598732571704021658e-01;
            const real a114  =  5.18637242884406370830023853209e+00;
            const real a115  =  1.09143734899672957818500254654e+00;
            const real a116  = -8.14978701074692612513997267357e+00;
            const real a117  = -1.85200656599969598641566180701e+01;
            const real a118  =  2.27394870993505042818970056734e+01;
            const real a119  =  2.49360555267965238987089396762e+00;
            const real a1110 = -3.0467644718982195003823669022e+00;

            const real a121  =  2.27331014751653820792359768449e+00;
            const real a124  = -1.05344954667372501984066689879e+01;
            const real a125  = -2.00087205822486249909675718444e+00;
            const real a126  = -1.79589318631187989172765950534e+01;
            const real a127  =  2.79488845294199600508499808837e+01;
            const real a128  = -2.85899827713502369474065508674e+00;
            const real a129  = -8.87285693353062954433549289258e+00;
            const real a1210 =  1.23605671757943030647266201528e+01;
            const real a1211 =  6.43392746015763530355970484046e-01;

            // -- 8th order weights (b_i) --
            const real b1  =  5.42937341165687622380535766363e-02;
            const real b6  =  4.45031289275240888144113950566e+00;
            const real b7  =  1.89151789931450038304281599044e+00;
            const real b8  = -5.8012039600105847814672114227e+00;
            const real b9  =  3.1116436695781989440891606237e-01;
            const real b10 = -1.52160949662516078556178806805e-01;
            const real b11 =  2.01365400804030348374776537501e-01;
            const real b12 =  4.47106157277725905176885569043e-02;

            // -- 5th order error coefficients (er_i) --
            const real er1  =  0.1312004499419488073250102996e-01;
            const real er6  = -0.1225156446376204440720569753e+01;
            const real er7  = -0.4957589496572501915214079952e+00;
            const real er8  =  0.1664377182454986536961530415e+01;
            const real er9  = -0.3503288487499736816886487290e+00;
            const real er10 =  0.3341791187130174790297318841e+00;
            const real er11 =  0.8192320648511571246570742613e-01;
            const real er12 = -0.2235530786388629525884427845e-01;

            real h = aStep;

            real t0 = aT;
            real t1;

            // initial derivative
            aODE.compute( t0, y0, k0 );

            real tolb;
            Status aStatus = Status::OK;

            for( uint i=0; i<aMaxIterations; ++i )
            {
                // reset status (may have been TRAPPED on a previous rejected step)
                aStatus = Status::OK ;

                // trap h
                if( t0 + h >= aTmax )
                {
                    h = aTmax - t0;
                    aStatus = Status::TRAPPED;
                }

                // capture the end time before h is adapted
                t1 = t0 + h;

                // -- Stage 2 --
                yt = y0 + h * a21 * k0;
                aODE.compute( t0 + c2 * h, yt, k1 );

                // -- Stage 3 --
                yt = y0 + h * ( a31 * k0 + a32 * k1 );
                aODE.compute( t0 + c3 * h, yt, k2 );

                // -- Stage 4 --
                yt = y0 + h * ( a41 * k0 + a43 * k2 );
                aODE.compute( t0 + c4 * h, yt, k3 );

                // -- Stage 5 --
                yt = y0 + h * ( a51 * k0 + a53 * k2 + a54 * k3 );
                aODE.compute( t0 + c5 * h, yt, k4 );

                // -- Stage 6 --
                yt = y0 + h * ( a61 * k0 + a64 * k3 + a65 * k4 );
                aODE.compute( t0 + c6 * h, yt, k5 );

                // -- Stage 7 --
                yt = y0 + h * ( a71 * k0 + a74 * k3 + a75 * k4
                              + a76 * k5 );
                aODE.compute( t0 + c7 * h, yt, k6 );

                // -- Stage 8 --
                yt = y0 + h * ( a81 * k0 + a84 * k3 + a85 * k4
                              + a86 * k5 + a87 * k6 );
                aODE.compute( t0 + c8 * h, yt, k7 );

                // -- Stage 9 --
                yt = y0 + h * ( a91 * k0 + a94 * k3 + a95 * k4
                              + a96 * k5 + a97 * k6 + a98 * k7 );
                aODE.compute( t0 + c9 * h, yt, k8 );

                // -- Stage 10 --
                yt = y0 + h * ( a101 * k0 + a104 * k3 + a105 * k4
                              + a106 * k5 + a107 * k6 + a108 * k7
                              + a109 * k8 );
                aODE.compute( t0 + c10 * h, yt, k9 );

                // -- Stage 11 --
                yt = y0 + h * ( a111 * k0 + a114 * k3 + a115 * k4
                              + a116 * k5 + a117 * k6 + a118 * k7
                              + a119 * k8 + a1110 * k9 );
                aODE.compute( t0 + c11 * h, yt, k10 );

                // -- Stage 12 --
                yt = y0 + h * ( a121 * k0 + a124 * k3 + a125 * k4
                              + a126 * k5 + a127 * k6 + a128 * k7
                              + a129 * k8 + a1210 * k9 + a1211 * k10 );
                aODE.compute( t1, yt, k11 );

                // -- 8th order solution --
                ya = y0 + h * ( b1 * k0  + b6  * k5  + b7  * k6
                              + b8 * k7  + b9  * k8  + b10 * k9
                              + b11 * k10 + b12 * k11 );

                // -- error estimation (5th order embedded) --
                tolb = 0.0;

                for( uint k=0; k<ya.length(); ++k )
                {
                    real err = h * ( er1  * k0( k )  + er6  * k5( k )
                                   + er7  * k6( k )  + er8  * k7( k )
                                   + er9  * k8( k )  + er10 * k9( k )
                                   + er11 * k10( k ) + er12 * k11( k ) );
                    tolb += err * err;
                }

                tolb = std::sqrt( tolb ) / h;

                if( aAutoTimestep )
                {
                    // adapt timestep (exponent 1/8 for 8th order method)
                    if ( tolb > 0 )
                    {
                        h *= 0.9 * std::pow( aEpsilon / tolb, 0.125 );
                    }
                    else
                    {
                        h *= 1.05;
                    }
                }

                if( tolb < aEpsilon )
                {
                    aT = t1;
                    aY.vector_data() = ya.vector_data();

                    if( aStatus == Status::OK && aAutoTimestep )
                    {
                        aStep = h;
                    }
                    return aStatus;
                }
            }

            return Status::MAXIT;
        }

//------------------------------------------------------------------------------
    }
}
