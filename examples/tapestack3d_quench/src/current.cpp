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
 * User-defined transport current for the tapestack3d model: the linear
 * ramp from ramps.hpp, the same file the Ic defect reads its timing from.
 *
 * Deck:
 *     boundary conditions { current { type : userdefined ; units : A ;
 *                                     file : src/build/userdefect.so ;
 *                                     label : MyCurrent ; } }
 */

#include <belfem_user_api>
#include "ramps.hpp"

using namespace belfem;

real my_current( const real t )
{
    return ramps::current( t );
}

// The name must be <label>_init, where <label> is the `label` key of the
// deck's current block.
extern "C" void MyCurrent_init( SourceFunction * source )
{
    source->set_user_defined( &my_current );
}
