/*
 * Example User-Defined Source Function for BELFEM
 *
 * This file demonstrates how to create a custom source that can be
 * dynamically loaded by BELFEM. Copy and modify this template for your
 * own sources.
 *
 * Compilation:
 *   1. Copy UserLibraryTemplate.cmake to your directory as CMakeLists.txt
 *   2. Edit CMakeLists.txt to set BELFEM_DIR and the library name
 *   3. mkdir build && cd build && cmake .. && make
 *
 */

#include <belfem_user_api>

#include "user_table.hpp"

using namespace belfem;

// =============================================================================
// User defined source functions
// =============================================================================
// Define your source function here.

/**
 * @brief Custom current source with interpolation from data file
 * @param t time (s)
 * @return interpolated current value (A)
 */
real my_current(const real t)
{
    // Transport current as measured, tabulated against time in seconds.
    static usermat::Table tTable( "I_vs_t_regular_smooth.txt" ) ;

    return tTable( t ) ;
}

// =============================================================================
// SOURCE INITIALIZATION FUNCTION
// =============================================================================
// This function is called when the source is loaded. The function name
// must be: extern "C" void <SourceName>_init(SourceFunction* source)
// where <SourceName> is the deck's `label :`, which BELFEM passes as the
// second argument of SourceFunction::read_user_defined().
//
// Note the argument type: a source plugin receives a SourceFunction*, NOT the
// Material* a material or defect plugin receives.

extern "C" void MyCurrent_init(SourceFunction* source)
{

    source->set_user_defined(&my_current);

}
