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

#ifndef BELFEM_CL_ASCII_HPP
#define BELFEM_CL_ASCII_HPP



#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "filetools.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * @brief Line-based ASCII file interface.
     *
     * @ingroup grp_io
     * @see @ref io_io_usage_guide
     */
    class Ascii
    {
//------------------------------------------------------------------------------
    protected:
//------------------------------------------------------------------------------

              string         mPath;
        const FileMode       mMode;
        Cell< string >       mBuffer;

        bool                 mChangedSinceLastSave = false;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        Ascii( const string        & aPath,
               const enum FileMode & aMode );

//------------------------------------------------------------------------------

        virtual ~Ascii();

//------------------------------------------------------------------------------

        // save the buffer to the file
        bool
        save();

//------------------------------------------------------------------------------

        void
        print( const std::string & aLine );

//------------------------------------------------------------------------------

        /**
         * return the number of lines
         */
        index_t
        length() const;

//------------------------------------------------------------------------------

        string &
        line( const index_t aLineNumber );

//------------------------------------------------------------------------------

        const string &
        line( const index_t aLineNumber ) const;
//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

        void
        load_buffer( const bool aParallelMode );

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
}
#endif //BELFEM_CL_ASCII_HPP
