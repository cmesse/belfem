//
// Created by Christian Messe on 25.08.20.
//

#ifndef BELFEM_UNITS_HPP
#define BELFEM_UNITS_HPP


#include "typedefs.hpp"
#include "constants.hpp"

#ifdef BELFEM_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#define BELFEM_TREF 298.15
#define BELFEM_PREF 100000.0

//------------------------------------------------------------------------------

namespace belfem
{
    namespace constant
    {
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // Metric Units
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        const real l    = 0.001 ;                 // liter
        const real km   = 1000  ;                 // kilometers
        const real bar  = 1.0e5 ;
        const real atm  = 1.01325e5 ;
        const real rpm  = constant::pi / 30.0 ; // rpm in rad/s

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // Freedom Units
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        const real in  = 0.0254 ;               // inch
        const real ft  = 12.0 * in ;            // foot
        const real mi  = 5280.0 * ft ;          // mile
        const real gal = 231.0 * in * in * in ; // us gallon
        const real lb  =  0.45359237 ;          // pound
        const real lbf = lb * constant::g0;     // pound force
        const real psi = lbf / ( in * in ) ;    // pound per square inch

        const real oz = gal / 128.0 ;           // us fluid ounce
    }
}

//------------------------------------------------------------------------------
#ifdef BELFEM_GCC
#pragma GCC diagnostic pop
#endif
#endif //BELFEM_UNITS_HPP
