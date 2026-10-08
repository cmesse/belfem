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

/*
 * User plugin for the tape_with_contact deck: the REBCO/Ag contact layer
 * ( material label "rint" ) and the Jc defect that forces current through it
 * ( defect label "Defect" ).
 *
 * A contact resistance is modeled as a thin resistive layer of the tape stack.
 * Its resistivity is rho = R_ct / h, so that the layer's resistance to current
 * crossing it, rho * h per unit area, equals the measured specific contact
 * resistance R_ct.
 */

#include <belfem_user_api>

#include "user_table.hpp"
using namespace belfem;

// Specific contact resistance of the REBCO/Ag interface, in Ohm*m^2.
//   2.5e-12   as-supplied commercial tape ( ~25 nOhm*cm^2 )
//   1e-11     typical, used here           ( 100 nOhm*cm^2 )
//   1e-10     degraded or aged interface   ( 1 uOhm*cm^2 )
// See src/fem/maxwell/doc/contact_impedance_theory.md, section 10.
const real rint_Rct = 1.0e-11 ;

// Thickness of the rint layer in input.conf, in m. Both must agree.
const real rint_h = 0.1e-6 ;

real rint_resistivity( const Material * mat, real T )
{
    return rint_Rct / rint_h ;
}

real rint_cp( const Material * mat, real T )
{
    static usermat::Table tTable( "cp_HTS.txt" ) ;
    return tTable( T ) ;
}

real rint_k( const Material * mat, real T )
{
    static usermat::Table tTable( "k_HTS.txt" ) ;
    return tTable( T ) ;
}

extern "C" void rint_init( Material * mat )
{
    mat->set_constant( MaterialProperty::ref_density, 6390.0 ) ;   // kg/m^3
    mat->set_constant( MaterialProperty::T_ref_density, 293.15 ) ; // K

    mat->set_user_defined_function( MaterialProperty::rho,
                                    MaterialDependency::T,
                                    &rint_resistivity ) ;

    mat->set_user_defined_function( MaterialProperty::cp,
                                    MaterialDependency::T,
                                    &rint_cp ) ;

    mat->set_user_defined_function( MaterialProperty::lambda,
                                    MaterialDependency::T,
                                    &rint_k ) ;
}

// Jc multiplier of a permanent defect, coordinates in m.
//
// Axes from tape.geo: z runs along the tape ( 10 mm ), x across it ( 4 mm ),
// y is the tape normal. The defect sits mid-length, inside the refined box.
//
// The defect spans the full width, a local loss of REBCO thickness as in
// Wozniak et al. 2025, Sec. III-E. With no superconducting path around it,
// every ampere above its Ic must cross the contact twice.
//
// Sizing: sp-ap gives Icw = 38 kA/m at 77.5 K, self field, for a 0.87 um
// film; with the 1.6 um film of the deck that is ~280 A over 4 mm before the
// self-field reduction. Depth 0.35 leaves Ic_d ~ 100 A, so at 150 A
// transport about 50 A runs in the stabilizer past the defect.
real my_defect( const real x, const real y, const real z, const real t )
{
    const real depth = 0.35 ;
    const real az    = 1.0e-3 ;    // 2 mm along the tape
    const real w     = 0.2e-3 ;    // erf edge width
    const real z0    = 5.0e-3 ;

    // smoothed top hat: ~1 inside, ~0 outside, exactly 0.5 at +-az
    const real bz = 0.5 * ( std::erf( ( z - z0 + az ) / w )
                          - std::erf( ( z - z0 - az ) / w ) ) ;

    return 1.0 - ( 1.0 - depth ) * bz ;
}

// Alternative: a crack entering from the edge at x = -2 mm, the typical shape
// of slitting or handling damage. Its 0.4 mm length needs a mesh resolution
// of 20 to 50 um around it; the default tape.geo does not provide that.
//
// real my_defect( const real x, const real y, const real z, const real t )
// {
//     const real depth = 0.0 ;
//     const real az    = 0.2e-3 ;    // 0.4 mm along the tape
//     const real xTip  = -0.5e-3 ;    // reaches 1.5 mm into the width
//     const real w     = 0.05e-3 ;
//     const real z0    = 5.0e-3 ;
//
//     const real bz = 0.5 * ( std::erf( ( z - z0 + az ) / w )
//                           - std::erf( ( z - z0 - az ) / w ) ) ;
//     const real bx = 0.5 * ( 1.0 - std::erf( ( x - xTip ) / w ) ) ;
//
//     return 1.0 - ( 1.0 - depth ) * bz * bx ;
// }

extern "C" void Defect_init( Material * mat )
{
    mat->set_user_defined_defect( &my_defect ) ;
}
