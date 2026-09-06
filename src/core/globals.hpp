/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_GLOBALS_HPP
#define BELFEM_GLOBALS_HPP

#include "typedefs.hpp"


#ifdef BELFEM_INITIALIZE_GLOBALS // only set from within communicator
#define greal real
#define gstring string
#else
#define greal extern real
#define gstring extern string
#endif

/**
 *
 * USER GUIDES:
 *
 *  -  if one adds new parameters, they must also be initialized in
 *    Communicator::set_globals()
 *
 *    recommendation: use BELFEM_QUIET_NAN if there is no clear default value
 *
 *  - if the value is not set on all procs, set value on master (aka root proc )
 *    and synchronize using
 *
 *    broadcast( T & value );
 *
 */
namespace belfem
{
    //! temperature in K of the bulk material when no thermal kernel is chosen;
    //! NaN until the executable or the deck sets it
    greal gTbulk ;

    //! minimim resistrivity in Ohm*m, default: 0
    greal gRhoMin ;

    //! maximim resistrivity in Ohm*m, default: 1e10
    greal gRhoMax ;

    //! path to belfem data files, set by environment variable $BELFEM_DATA
    gstring gBelfemDataPath ;
}

#undef greal

#endif //BELFEM_GLOBALS_HPP