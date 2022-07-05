//
// Created by Christian Messe on 13.11.19.
//

#ifndef BELFEM_EN_MATERIALS_HPP
#define BELFEM_EN_MATERIALS_HPP
#include "typedefs.hpp"
#include "assert.hpp"
#include "stringtools.hpp"

namespace belfem
{
    enum class MaterialType
    {
        Unity       =  0, // all material properties are 1
        Simple      =  1, // a simple material for testing, based on aluminum
        Aluminum    =  2,
        Copper      =  3,
        Inconel718  =  4,
        TI6AL4V     =  5,
        CCSIC       =  6,
        Altramat80  =  7,
        Rohacell51  =  8,
        Air         =  9, // wrapper around an air object
        CuCrZr      = 10,
        Zirconia    = 11,
        Inconel750X = 12,
        Hastelloy   = 13,
        Silver      = 14,
        Maxwell     = 15, // special type that follows the power law, can use spline as temperature function
        UNDEFINED   = 16
    };

    enum class ResistivityLaw
    {
        Constant       = 0,
        DependT        = 1 ,
        DependB        = 2 ,
        DependBT       = 3 ,
        DependJ        = 4 ,
        DependJT       = 5 ,
        DependJB       = 6 ,
        DependJBT      = 7 ,
        UNDEFINED      = 8
    };

    enum class PermeabilityLaw
    {
        Constant  = 0,
        Spline    = 1,
        UNDEFINED = 2
    };


    string
    to_string( const MaterialType & aMaterialType );

    MaterialType
    string_to_material_type( const string & aString, const bool aNoErrorIfUnknown = false );

}
#endif //BELFEM_EN_MATERIALS_HPP
