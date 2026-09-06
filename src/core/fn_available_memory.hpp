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

#ifndef BELFEM_FN_AVAILABLE_MEMORY_HPP
#define BELFEM_FN_AVAILABLE_MEMORY_HPP

#include <cstddef>

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * Bytes a new allocation by this process can claim without swapping,
     * or 0 when that cannot be established.
     *
     * Linux: the kernel's MemAvailable estimate ( /proc/meminfo, kB ), NOT
     * sysinfo().freeram -- "free" excludes the reclaimable page cache, and
     * on a workstation with a warm cache the two differ by an order of
     * magnitude. Capped by the tightest cgroup memory limit above this
     * process ( v2 memory.max / memory.current, v1 memory.limit_in_bytes /
     * memory.usage_in_bytes ), because a batch scheduler enforces that
     * limit and the host figure would overstate what the job may take.
     * Darwin: Mach host statistics, free plus inactive pages.
     * Other platforms: 0.
     *
     * Setup-path code: the file reads are the one-off allocation
     * doc/coding_philosophy.md allows there. Never aborts -- a consumer
     * treats 0 as "unknown" and falls back to whatever it did before.
     */
    std::size_t
    available_memory();

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_AVAILABLE_MEMORY_HPP
