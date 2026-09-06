/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Unit test for belfem::available_memory(), the probe behind the MUMPS
 * ICNTL(23) memory budget.
 */

#include <gtest/gtest.h>
#include <cstddef>

#include "fn_available_memory.hpp"

// On the two supported platforms the probe must find SOMETHING, and the
// number must be a plausible byte count: above one megabyte ( a machine
// that cannot spare that much cannot run the test binary ) and below
// 64 TiB ( a kB figure mistaken for bytes, or bytes for kB, lands outside
// this window on any current host ). Elsewhere the contract is 0
TEST( AvailableMemory, ProbeReturnsPlausibleBytes )
{
    const std::size_t tBytes = belfem::available_memory();

#if defined( __linux__ ) || defined( __APPLE__ )
    EXPECT_GT( tBytes, ( std::size_t ) 1 << 20 );
    EXPECT_LT( tBytes, ( std::size_t ) 1 << 46 );
#else
    EXPECT_EQ( tBytes, ( std::size_t ) 0 );
#endif
}

// two consecutive reads must agree to well within the machine: the probe
// is a measurement, not a random number
TEST( AvailableMemory, ProbeIsStable )
{
    const std::size_t tFirst  = belfem::available_memory();
    const std::size_t tSecond = belfem::available_memory();

    const std::size_t tLarge = tFirst > tSecond ? tFirst : tSecond ;
    const std::size_t tSmall = tFirst > tSecond ? tSecond : tFirst ;

    // half of the larger reading is a generous bound for cache churn
    // between two calls a few microseconds apart
    EXPECT_GE( tSmall * 2, tLarge );
}
