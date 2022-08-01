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
        Unity          =  0, // all material properties are 1
        Simple         =  1, // a simple material for testing, based on aluminum
        Maxwell        =  2, // special type that follows the power law, can use spline as temperature function
        Air            =  3, // wrapper around an air object
        Aluminum       =  4,
        Copper         =  5,
        Inconel718     =  6,
        TI6AL4V        =  7,
        CCSIC          =  8,
        Altramat80     =  9,
        Rohacell51     = 10,
        CuCrZr         = 11,
        Zirconia       = 12,
        Inconel750X    = 13,
        HastelloyC276  = 14,
        Silver         = 15,
        YBCO           = 16,
        Pb40Sn60       = 17,
        SAE301         = 18,
        SAE316         = 19,
        SAE347         = 20,
        SAE1010        = 21,
        UNDEFINED      = 22
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
