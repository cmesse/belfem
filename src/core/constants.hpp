//
// Created by Christian Messe on 2019-01-10.
//

#ifndef BELFEM_CONSTANTS_HPP
#define BELFEM_CONSTANTS_HPP
#include <cmath>
#include "typedefs.hpp"

#ifdef BELFEM_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#define BELFEM_TREF 298.15
#define BELFEM_PREF 100000.0

namespace belfem
{
    namespace constant
    {
//------------------------------------------------------------------------------
// MATH
//------------------------------------------------------------------------------

        /**
         * circle number
         */
        const real pi = 3.141592653589793238462643383279502884;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * golden number
         */
        const real phi = 0.5 * ( 1.0 + std::sqrt( 5.0 ) );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * degree
         */
        const real deg = pi / 180.0;

//------------------------------------------------------------------------------
// GENERAL PHYSICS
//------------------------------------------------------------------------------

        /**
         * speed of light in m/s ( exact )
         * http://physics.nist.gov/cgi-bin/cuu/Value?c
         */
        const real c = 299792458.0 ;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * magnetic constant in V*s/(A*m) ( recommended value )
         * http://physics.nist.gov/cgi-bin/cuu/Value?mu0
         */
        const real mu0 = 1.25663706212E-6;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * inverse of magnetic constant in (A*m)/(V*s)
         */
        const real nu0 = 1./mu0 ;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * electric constant in A*s/(V*m)
         */
        const real epsilon0 = 1.0/( mu0*std::pow( c, 2 ) );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * gravity constant in N*m^2 / kg^2
         *
         * https://physics.nist.gov/cgi-bin/cuu/Value?bg|search_for=gravity
         */
        const real G = 6.67430e-11;

//------------------------------------------------------------------------------
// THERMODYNAMICS
//------------------------------------------------------------------------------

        /**
         * thermal calorie in SI units J/cal
         */
        const real calTh = 4.184;

//------------------------------------------------------------------------------

        /**
         * Boltzman constant ( exact )
         * https://physics.nist.gov/cgi-bin/cuu/Value?k|search_for=boltzman
         */
        const real kB = 1.380649e-23;


//------------------------------------------------------------------------------

        /**
         * Avogadro constant ( exact )
         * https://physics.nist.gov/cgi-bin/cuu/Value?na|search_for=avogadro
         */
         const real NA = 6.02214076e23;

//------------------------------------------------------------------------------

        /**
         * atomic mass unit-kilogram relationship
         * https://physics.nist.gov/cgi-bin/cuu/Value?Rukg|search_for=unit
         */
        const real u = 0.001 / NA ;

//------------------------------------------------------------------------------

        /**
         * Gas constant in J/(K*mol) ( exact )
         * http://physics.nist.gov/cgi-bin/cuu/Value?r
         *
         * using the new international definition
         */
        const real Rm = kB * NA; // 8.3144598;

//------------------------------------------------------------------------------

        /**
         * Gas Constant in cal/(K*mol)
         */
        const real Rm_cal = Rm / calTh;

//------------------------------------------------------------------------------

        /**
         * Planck constant in J / Hz ( exact )
         * http://physics.nist.gov/cgi-bin/cuu/Value?h
         *
         * using the new international definition
         */
         const real h = 6.62607015e-34;

//------------------------------------------------------------------------------

        /**
         * Stephan Boltzman constant
         */
         const real sigma = 2.0 * std::pow( pi, 5 )
                            * std::pow( kB, 4 ) /
                ( 15.0 * std::pow( h, 3 ) * std::pow( c, 2 ) );

//------------------------------------------------------------------------------
// EARTH
//------------------------------------------------------------------------------

        /**
         * earth gravitational acceleration at sea level in m/s^2
         * https://physics.nist.gov/cgi-bin/cuu/Value?gn
         */
        const real g0 = 9.80665;

//------------------------------------------------------------------------------

        /**
         * gravity parameter Earth ( G * M_earth )
         */
        const real mu_earth = 3.986004418e14;

//------------------------------------------------------------------------------

        /**
         * reference radius for Earth
         */
        const real R_earth = std::sqrt( mu_earth / g0 );

//------------------------------------------------------------------------------

    }
}

#ifdef BELFEM_GCC
#pragma GCC diagnostic pop
#endif

#endif //BELFEM_CONSTANTS_HPP
