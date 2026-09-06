/*
 * Example User-Defined Defect for BELFEM
 *
 * This file demonstrates how to create a custom defect that can be
 * dynamically loaded by BELFEM. Copy and modify this template for your
 * own defects.
 *
 * Compilation:
 *   1. Copy UserMaterialTemplate.cmake to your directory as CMakeLists.txt
 *   2. Edit CMakeLists.txt to set BELFEM_DIR and material name
 *   3. mkdir build && cd build && cmake .. && make
 *
 */

#include <belfem_user_api>

using namespace belfem;

// =============================================================================
// User defined defect functions
// =============================================================================
// Define your defect function here. The first parameter must

// Axis convention, from tape.geo: x runs along the tape ( domainLength
// 12.5 mm ), y runs across it ( tapeWidth 4 mm ). The refined "blast"
// ellipse is blastLength = 1.618 mm in x by blastWidth = 1.8 mm in y, at
// tapeResolutionQuench = 0.0438 mm.
//
// Sizing, at 77 K with Jc = 3.375e10 A/m^2 over the 1 um x 4 mm layer:
// Ic = 135 A against a transport current of 137 A, so the tape sits at
// I/Ic = 1.016 before the defect does anything. The defect blocks a
// fraction f of the width at a jc multiplier of depth, leaving an
// effective Ic of Ic*( 1 - f + f*depth ) = 80 A, hence I/Ic_eff = 1.7.
//
// Length matters as much as depth. The minimum propagating zone here is
// ~0.95 mm ( k_eff = 119 W/m/K from data/k_*.txt, stabiliser 0.047 Ohm/m ),
// so a 1 mm defect sits exactly on the stability threshold and dissipates
// forever without running away. 2 mm clears it by better than 2x.
real my_defect(const real x, const real y, const real z, const real t)
{
    const real depth = 0.1 ;
    const real ax    = 1.0e-3 ;   // 2.0 mm along the tape, > 2x the MPZ
    const real ay    = 0.9e-3 ;   // blocks 1.8 mm of the 4 mm width
    const real w     = 0.2e-3 ;   // erf edge 0.47 mm ~ 10.6 elements

    const real x0 = 0.0 ;
    const real y0 = 0.0 ;

    // Switch-on. The transport current in data/I_vs_t_regular_smooth.txt
    // is only on plateau from 5 ms to 155 ms, then decays to zero by
    // 250 ms, so the quench has to nucleate and propagate inside that
    // window. Holding off to 20 ms lets the current distribution settle
    // first ( the previous run reached steady state by then ) and still
    // leaves 130 ms of plateau afterwards -- ample, since a fully normal
    // band of this size heats at order 1e3 K/s.
    //
    // The blend is the quintic smoothstep, which is C2. The cubic is only
    // C1, and its curvature jump at either end of the ramp costs BDF5 an
    // order and makes the controller cut the step twice for no reason.
    const real tOn   = 20.0e-3 ;
    const real tRamp =  5.0e-3 ;

    real s = ( t - tOn ) / tRamp ;
    s = s < 0.0 ? 0.0 : ( s > 1.0 ? 1.0 : s ) ;
    s = s * s * s * ( s * ( 6.0 * s - 15.0 ) + 10.0 ) ;

    // smoothed top hat: ~1 inside, ~0 outside, exactly 0.5 at +-a
    const real bx = 0.5 * ( std::erf( ( x - x0 + ax ) / w )
                          - std::erf( ( x - x0 - ax ) / w ) ) ;
    const real by = 0.5 * ( std::erf( ( y - y0 + ay ) / w )
                          - std::erf( ( y - y0 - ay ) / w ) ) ;

    return 1.0 - s * ( 1.0 - depth ) * bx * by ;

}

// =============================================================================
// DEFECT INITIALIZATION FUNCTION
// =============================================================================
// This function is called when the defect is loaded. The function name
// must be: extern "C" void <DefectName>_init(Material* mat)
// where <DefectName> matches the second argument in read_defect()

extern "C" void MyDefect_init(Material* mat)
{

    mat->set_user_defined_defect(&my_defect);

}
