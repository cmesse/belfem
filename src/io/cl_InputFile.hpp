//
// Created by christian on 9/17/21.
//

#ifndef BELFEM_CL_INPUTFILE_HPP
#define BELFEM_CL_INPUTFILE_HPP

#include "cl_Ascii.hpp"
#include "cl_Cell.hpp"
#include "cl_Input_Section.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    class InputFile : public Ascii
    {
        input::Section * mData = nullptr ;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        InputFile( const string & aPath );

        ~InputFile();

        void
        print();

//------------------------------------------------------------------------------

        /**
         * return a subsection by label
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
