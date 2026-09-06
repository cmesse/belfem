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

#ifndef BELFEM_BANNER_HPP
#define BELFEM_BANNER_HPP

#include <string>

namespace belfem
{
    const std::string gLongName = "BELFEM -- The Berkeley Lab Finite Element Framework";
    const std::string gURL      = "http://belfem.lbl.gov";

//------------------------------------------------------------------------------

    std::string
    exec( const std::string & aCommand );

//------------------------------------------------------------------------------

    /**
     * returns the Unix version ( Linux or Darwin )
     */
    std::string
    uname();

//------------------------------------------------------------------------------

    /**
     * grabs the cpu info for the banner
     */
    std::string
    cpu_info();

//------------------------------------------------------------------------------

    /**
     * returns the BELFEM semantic version as "major.minor.patch" (e.g. "0.1.0"),
     * assigned via project( belfem VERSION ... ) in the top-level CMakeLists.txt
     */
    std::string
    version();

//------------------------------------------------------------------------------

    /**
     * true if the build was produced from a git checkout (i.e. git commit
     * metadata is available); false when built from a source tarball with no
     * .git directory, in which case the git_* accessors return "unknown"
     */
    bool
    is_built_from_git();

//------------------------------------------------------------------------------

    /**
     * returns the full git commit hash of the build, or "unknown"
     * if built outside a git repository (e.g. a source tarball)
     */
    std::string
    git_commit_hash();

//------------------------------------------------------------------------------

    /**
     * returns the abbreviated (short) git commit hash of the build
     */
    std::string
    git_commit_hash_short();

//------------------------------------------------------------------------------

    /**
     * returns the git branch name at build time, or "unknown"
     */
    std::string
    git_branch();

//------------------------------------------------------------------------------

    /**
     * true if the working tree had uncommitted changes at build time
     */
    bool
    git_is_dirty();

//------------------------------------------------------------------------------

    /**
     * prints the banner
     */
    void
    print_banner( const std::string aExecName = "" );

//------------------------------------------------------------------------------

    namespace banners
    {
        bool
        print_default();

        bool
        print_easter();

        bool
        print_stpatrick();

        bool
        print_usa();

        bool
        print_canada();

        bool
        print_thanksgiving();

        bool
        print_christmas();
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_BANNER_HPP