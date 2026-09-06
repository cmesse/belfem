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

#ifndef BELFEM_CL_INPUTFILE_HPP
#define BELFEM_CL_INPUTFILE_HPP

#include "cl_Ascii.hpp"
#include "cl_Cell.hpp"
#include "cl_Input_Section.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * @brief Parser for the input.conf configuration format.
     *
     * @ingroup grp_io
     * @see @ref io_io_usage_guide
     */
    class InputFile : public Ascii
    {
        input::Section * mData = nullptr ;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        InputFile( const string & aPath );

        ~InputFile() override;

        void
        print();

//------------------------------------------------------------------------------

        /**
         * return a top-level section by type string
         * ( a labelled section is keyed as "type:label" )
         */
        const input::Section *
        section( const string & aSection  ) const ;

//------------------------------------------------------------------------------

        /**
         * access a section by index
         */
        const input::Section *
        section( const index_t aIndex ) const ;

//------------------------------------------------------------------------------
        /**
          * tell if a section exists
          */
        bool
        section_exists( const string & aSection ) const ;

//------------------------------------------------------------------------------

        /**
         * number of sections
         */
        index_t
        num_sections() const ;

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

        void
        remove_comments();

//------------------------------------------------------------------------------

        void
        tidy_up();

//------------------------------------------------------------------------------
    };

//----------------------------------------------------------------------------
}
#endif //BELFEM_CL_INPUTFILE_HPP
