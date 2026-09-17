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
 * User-defined Ic defect for the tapestack3d quench run.
 *
 * Built by src/CMakeLists.txt into src/build/userdefect.so; the deck loads it
 * through
 *
 *     materials { ybco { defect { file : src/build/userdefect.so ;
 *                                 label : MyDefect ; } } }
 *
 * The function returns the jc multiplier at ( x, y, z, t ): 1 = intact,
 * gDepth = fully degraded.
 */

#include <belfem_user_api>
#include "ramps.hpp"

using namespace belfem;

// =============================================================================
// Geometry of the stack this defect is written for ( tapestack3d.geo )
// =============================================================================
// Eight 4 mm wide tapes in the planes y = -0.35, -0.25, ..., +0.35 mm
// ( tapeDistance 0.1 mm ), x in [ -2, 2 ] mm, z in [ 0, 10 ] mm, solder in
// between. Only the TOPMOST tape ( y = +0.35 mm ) is damaged: the y gate
// below sits halfway between it and its neighbour at +0.25 mm.
//
// The defect is a cylinder along y, centred on x = 0, z = 5 mm ( mid-length,
// away from the periodic end planes ), radius gRadius with an erf-smoothed
// wall of width gEdge. On the tape it is a disk that leaves both tape edges
// ( |x| = 2 mm ) at more than 99 % jc.
// =============================================================================

namespace
{
    constexpr real gX0     = 0.0e-3 ;
    constexpr real gZ0     = 5.0e-3 ;

    constexpr real gYGate  = 0.3e-3 ;   // tapes with y > gYGate carry the defect

    constexpr real gRadius = 1.5e-3 ;
    constexpr real gEdge   = 0.3e-3 ;

    constexpr real gDepth  = 0.1 ;      // jc multiplier inside the defect
}

real my_defect( const real x, const real y, const real z, const real t )
{
    if( y < gYGate ) return 1.0 ;       // all tapes below the top one are intact

    const real dx = x - gX0 ;
    const real dz = z - gZ0 ;

    const real r  = std::sqrt( dx * dx + dz * dz );

    // smoothed disk indicator: 0.5 on r = gRadius, ~1 inside, ~0 outside
    const real disk = 0.5 * ( 1.0 - std::erf( ( r - gRadius ) / gEdge ) );

    return 1.0  - ramps::defect_switch( t ) * ( 1.0 - gDepth ) * disk ;
}

// =============================================================================
// DEFECT INITIALIZATION FUNCTION
// =============================================================================
// Called when the defect is loaded. The name must be <label>_init, where
// <label> is the `label` key of the deck's defect block.

extern "C" void MyDefect_init( Material * mat )
{
    mat->set_user_defined_defect( &my_defect );
}
