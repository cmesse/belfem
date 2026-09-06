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

//------------------------------------------------------------------------------
// Test fixture for the source-function plugin search path.
//
// This is the suite's only dlopen-loaded artifact. It exists so that
// SourceFunction::read_user_defined() can be tested against a library that is
// reachable ONLY through the shared data directory -- not from the run
// directory and not through the platform loader path. The build places it at
//
//     <build>/tests/circuit/plugindata/material/libtest_source_plugin.so
//
// and the test points gBelfemDataPath at that plugindata/, so the deck-side
// name "libtest_source_plugin.so" carries no slash and no working path. Before
// the search-path fix that name reached dlopen unresolved and failed to load.
//
// It links nothing: set_user_defined() is resolved out of the test executable,
// which is built -rdynamic ( config/compiler/config_gcc.cmake:104 ). That is
// the same contract UserLibraryTemplate.cmake gives a real plugin author.
//------------------------------------------------------------------------------

#include "cl_SourceFunction.hpp"

using belfem::real;
using belfem::SourceFunction;

namespace
{
    //! a waveform no builtin type produces, so a passing test cannot be
    //! satisfied by the object's default state
    real
    test_plugin_waveform( const real aTime )
    {
        return 7.0 * aTime + 3.0 ;
    }
}

//------------------------------------------------------------------------------

extern "C" void
TestPluginSource_init( SourceFunction * aSource )
{
    aSource->set_user_defined( &test_plugin_waveform );
}

//------------------------------------------------------------------------------
