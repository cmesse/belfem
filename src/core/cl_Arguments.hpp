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

#ifndef BELFEM_CL_ARGUMENTS_HPP
#define BELFEM_CL_ARGUMENTS_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
namespace belfem
{
    // a basic argument class to be extended by any other executable
    /**
     * @brief Command-line argument parsing.
     *
     * @ingroup grp_core
     * @see @ref core_core_usage_guide
     */
    class Arguments
    {
//------------------------------------------------------------------------------
    protected:
//------------------------------------------------------------------------------

        // list of arguments that were passed to the executable
        Cell< string > mArguments;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
         * The constructor also handles the flags all BELFEM executables share,
         * GNU style: -v [N], -vN, --verbose [N] and --verbose=N set the info
         * level of the global logger; without a number, the level defaults to
         * InfoLevel::Everything. Subclasses must not reuse -v or --verbose
         * ( -V / --version is the conventional pair for a version flag ).
         */
        Arguments( int & argc, char * argv[] );

//------------------------------------------------------------------------------

        virtual ~Arguments() = default;

//------------------------------------------------------------------------------

        /**
         * returns a cell with arguments as strings
         */
         const Cell< string > &
         data() const ;

//------------------------------------------------------------------------------

        /**
         * access argument by index
         */
        const string &
        data( const index_t aIndex ) const ;

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

        /**
         * scan for -v / --verbose and pass the level to the global logger
         */
        void
        set_verbosity_from_arguments();

//------------------------------------------------------------------------------
    };
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
}
#endif //BELFEM_CL_ARGUMENTS_HPP
