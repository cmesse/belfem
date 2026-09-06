/*
 * Example User-Defined Source Function for BELFEM
 *
 * This file demonstrates how to create a custom source that can be
 * dynamically loaded by BELFEM. Copy and modify this template for your
 * own sources.
 *
 * A source function maps time to a scalar excitation -- a transport current
 * in A, or an applied flux density in T, depending on what the deck declares.
 * The signature is fixed: real ( const real t ).
 *
 * Compilation:
 *   1. Copy UserLibraryTemplate.cmake to your directory as CMakeLists.txt
 *   2. Edit CMakeLists.txt to set BELFEM_DIR and the library name
 *   3. mkdir build && cd build && cmake .. && make
 *
 * Two complete, running source plugins ship with BELFEM:
 *   examples/disk_pulse/src/bgpulse.cpp          -- analytic waveform
 *   examples/tape_quench_usermat/src/current.cpp -- measured table
 */

#include <belfem_user_api>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace belfem;

// =============================================================================
// User defined source functions
// =============================================================================
// Define your source function here.

// -----------------------------------------------------------------------------
// WHERE THE DATA FILE LIVES
// -----------------------------------------------------------------------------
// The directory holding the table is compiled in from CMake, so the plugin
// works regardless of the directory the solver is launched from.
// UserLibraryTemplate.cmake defines it; see USER_DATA_DIR there.
//
// Do NOT locate the data by deriving a path from __FILE__. That bakes the
// SOURCE directory into the binary, so the plugin breaks the moment the
// sources are moved or the plugin is built somewhere else. This template used
// to do exactly that, and it broke exactly that way.
#ifndef BELFEM_USER_DATA_DIR
#error "BELFEM_USER_DATA_DIR is not set -- see target_compile_definitions in UserLibraryTemplate.cmake"
#endif

/**
 * @brief Custom current source with interpolation from data file
 * @param t time (s)
 * @return interpolated current value (A)
 *
 * The file is a two-column ( time, value ) text table. It is read on the first
 * call and kept for the life of the process. Arguments outside the sampled
 * range are clamped to the end values.
 */
real my_current(const real t)
{
    // Static state, loaded once. A source function is called on every
    // Newton iteration of every timestep, so re-reading the file here would
    // dominate the run.
    static std::vector<real> time_data;
    static std::vector<real> current_data;
    static bool data_loaded = false;

    if (!data_loaded)
    {
        const std::string filename =
                std::string(BELFEM_USER_DATA_DIR) + "my_source.txt";

        std::ifstream file(filename);

        BELFEM_ERROR(file.is_open(),
                     "user source: cannot open %s", filename.c_str());

        std::string line;
        real time_val, current_val;

        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            if (iss >> time_val >> current_val)
            {
                time_data.push_back(time_val);
                current_data.push_back(current_val);
            }
        }

        file.close();
        data_loaded = true;

        BELFEM_ERROR(!time_data.empty(),
                     "user source: no data loaded from %s", filename.c_str());
    }

    // Handle boundary cases
    if (t <= time_data.front())
    {
        return current_data.front();
    }
    if (t >= time_data.back())
    {
        return current_data.back();
    }

    // Find the interval containing t using binary search
    auto it = std::lower_bound(time_data.begin(), time_data.end(), t);
    size_t idx = std::distance(time_data.begin(), it);

    // Linear interpolation
    real t0 = time_data[idx - 1];
    real t1 = time_data[idx];
    real I0 = current_data[idx - 1];
    real I1 = current_data[idx];

    real alpha = (t - t0) / (t1 - t0);
    return I0 + alpha * (I1 - I0);
}

// =============================================================================
// SOURCE INITIALIZATION FUNCTION
// =============================================================================
// This function is called when the source is loaded. The function name
// must be: extern "C" void <SourceName>_init(SourceFunction* source)
// where <SourceName> is the `label :` in the deck, which BELFEM passes as the
// second argument of SourceFunction::read_user_defined().
//
// Note the argument type: a source plugin receives a SourceFunction*, NOT the
// Material* a material or defect plugin receives.
//
// The library path from the deck's `file :` is resolved in this order:
//   1. the path as written, relative to the run directory, or absolute
//   2. the same relative path below <data root>/material
//   3. for a path CONTAINING a slash, the file name alone below that directory
//   4. otherwise the name is handed to dlopen unchanged
// The root is gBelfemDataPath, not $BELFEM_DATA directly. dlopen searches
// $LD_LIBRARY_PATH only for a slashless name, so adding a directory component
// DISABLES the loader search rather than enabling it.

extern "C" void MyCurrent_init(SourceFunction* source)
{

    source->set_user_defined(&my_current);

}
