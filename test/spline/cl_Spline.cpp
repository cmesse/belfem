//
// Created by Christian Messe on 2019-08-16.
//


#include <iostream>

#include <gtest/gtest.h>


#include "typedefs.hpp"
#include "cl_Communicator.hpp"

#include "cl_Vector.hpp"
#include "cl_Spline.hpp"
#include "cl_SpMatrix.hpp"

using namespace belfem;


const real Rm = 8.3144598;
const real M = 28.949e-3;
const real R = Rm / M;

// values for air
const Vector< real > gA = {
        1.009950160e+04,
        -1.968275610e+02,
        5.009155110e+00,
        -5.761013730e-03,
        1.066859930e-05,
        -7.940297970e-09,
        2.185231910e-12,
        -1.767967310e+02,
        -3.921504225e+00 };

// specific heat formula NASA RP-1311 Eq. ( 4.9 )
real
cp0( const real & aT )
{
    return    (  (     gA( 0 )/aT
              +        gA( 1 ) )/aT
              +        gA( 2 )
              + aT * ( gA( 3 )
              + aT * ( gA( 4 )
              + aT * ( gA( 5 )
              + aT *   gA( 6 ) ))))*R;
}

// specific enthalpy formula NASA RP-1311 Eq. ( 4.10 )
real
h0( const real aT )
{
    return  ( -        gA( 0 )/aT
              +        gA( 1 ) * std::log( aT )
              + aT * ( gA( 2 )
              + aT * ( gA( 3 )*0.5
              + aT * ( gA( 4 )/3.0
              + aT * ( gA( 5 )*0.25
              + aT *   gA( 6 )*0.2 ) ) ) )
              + gA( 7 ) ) * R;
}

// specific entropy formula  NASA RP-1311 Eq. ( 4.11 )
real
s0( const real aT )
{
    return (
                   -      ( gA( 0 ) * 0.5 / aT
                   +        gA( 1 ) ) / aT
                   +        gA( 2 ) * std::log( aT )
                   + aT * ( gA( 3 )
                   + aT * ( gA( 4 ) * 0.5
                   + aT * ( gA( 5 ) / 3.0
                   + aT *   gA( 6 ) * 0.25 ) ) )
                   + gA( 8 ) ) * R;
}

// derivative of cp
real
dcp0dT( const real & aT )
{
    return ( - ( 2.0 * gA( 0 ) + gA( 1 ) * aT ) / std::pow( aT, 3 )
             + gA( 3 ) + aT * ( 2.0 * gA( 4 ) + aT * ( 3.0 * gA( 5 )
             + aT * 4.0 * gA( 6 ) ) ) ) * R;
}

// entropy derivative
real
ds0dT( const real & aT )
{
    return ( ( gA( 0 ) + aT * ( gA( 1 ) + aT * gA( 2 ) ) ) / std::pow( aT, 3 )
             + gA( 3 ) + aT * ( gA( 4 ) + aT * ( gA( 5 ) + aT * gA( 6 ) ) ) ) * R;
}



TEST( Spline, spline )
{
    // create reference data
    Vector< real > tX = { 195.,  202.,  209.,  216.,  223.,  230.,  237.,  244.,
                          251.,  258.,  265.,  272.,  279.,  286.,  293.,  300.,
                          307.,  314.,  321.,  328.,  335.,  342.,  349.,  356.,
                          363.,  370.,  377.,  384.,  391.,  398.,  405.,  412.,
                          419.,  426.,  433.,  440.,  447.,  454.,  461.,  468.,
                          475.,  482.,  489.,  496.,  503.,  510.,  517.,  524.,
                          531.,  538.,  545.,  552.,  559.,  566.,  573.,  580.,
                          587.,  594.,  601.,  608.,  615.,  622.,  629.,  636.,
                          643.,  650.,  657.,  664.,  671.,  678.,  685.,  692.,
                          699.,  706.,  713.,  720.,  727.,  734.,  741.,  748.,
                          755.,  762.,  769.,  776.,  783.,  790.,  797.,  804.,
                          811.,  818.,  825.,  832.,  839.,  846.,  853.,  860.,
                          867.,  874.,  881.,  888.,  895.,  902.,  909.,  916.,
                          923.,  930.,  937.,  944.,  951.,  958.,  965.,  972.,
                          979.,  986.,  993.,  1000 };

    uint tN = tX.length();

    Vector< real > tY( tN );

    for ( uint k=0; k<tN; ++k )
    {
        tY( k ) = h0( tX( k ) );
    }

    // create help matrix
    SpMatrix tHelpMatrix;

    // initialize help matrix
    spline::create_helpmatrix( tN, 7.0, tHelpMatrix );

    // create the spline object, last two values are optional
    Spline tSpline( tX, tY, tHelpMatrix, 273.15, s0( 273.25) );

    // create new grid
    Vector< real > tT = { 200, 225, 250, 275, 300, 325, 350, 375, 400,
                          425, 450, 475, 500, 525, 550, 575, 600, 625,
                          650, 675, 700, 725, 750, 775, 800, 825, 850,
                          875, 900, 925, 950, 975, 1000};


    // calculate interpolated values
    tN = tT.length();

    // Errors
    real tR2H = 0.0;
    real tR2Cp = 0.0;
    real tR2dCpdT = 0.0;
    real tR2S = 0.0;
    real tR2dSdT = 0.0;
    for ( uint k=0; k<tN; ++k )
    {
        // enthalpy test
        real tH = h0( tT( k ) );
        tR2H += std::pow( ( tH - tSpline.eval( tT( k ) ) ) / tH, 2 );

        // specific heat test
        real tCp = cp0( tT( k ) );
        tR2Cp +=  std::pow( ( tCp - tSpline.deval( tT( k ) ) ) / tCp, 2 );

        // dcp test
        real tdCpdT = dcp0dT( tT( k ) );
        tR2dCpdT += std::pow( ( tdCpdT - tSpline.ddeval( tT( k ) ) ) / tdCpdT, 2 );

        // entropy test
        real tS = s0( tT( k ) );
        tR2S += std::pow( ( tS - tSpline.entropy( tT( k ) ) ) / tS, 2 );

        // entropy derovative test
        real tdSdT = ds0dT( tT( k ) );
        tR2dSdT += std::pow( ( tdSdT - tSpline.dentropy( tT( k ) ) ) / tdSdT, 2 );
    }

    EXPECT_TRUE(  tR2H     < 1e-12 ); // 1.1215e-15
    EXPECT_TRUE(  tR2Cp    < 1e-9 );  // 4.63261e-12
    EXPECT_TRUE(  tR2dCpdT < 2e-3 );  // 0.00157442
    EXPECT_TRUE(  tR2S     < 1e-6 );  // 8.03886e-08
    EXPECT_TRUE(  tR2dSdT  < 1e-9 );  // 4.63261e-12
}