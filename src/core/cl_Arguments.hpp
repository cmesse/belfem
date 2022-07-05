//
// Created by Christian Messe on 31.08.19.
//

#ifndef BELFEM_CL_ARGUMENTS_HPP
#define BELFEM_CL_ARGUMENTS_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
namespace belfem
{
    // a basic argument class to be extended by any other executable
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
    };
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
}
#endif //BELFEM_CL_ARGUMENTS_HPP
