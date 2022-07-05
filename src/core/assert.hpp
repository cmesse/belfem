//
// Created by Christian Messe on 18.11.18.
//

#ifndef BELFEM_ASSERT_HPP
#define BELFEM_ASSERT_HPP

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "typedefs.hpp"
#include "stringtools.hpp"

#ifdef BELFEM_MPI
#include "cl_Communicator.hpp"
extern belfem::Communicator gComm;
#endif

//------------------------------------------------------------------------------

namespace belfem
{
    namespace assert
    {
//------------------------------------------------------------------------------

        template< typename Exception >
        void
        error(
                const std::string    & aLocation,
                const std::string    & aTask,
                const std::string    & aCheck,
                const Exception & aException = Exception()
        )
        {
            std::cerr << "*** ---------------------------------------------------------------------------\n";
            std::cerr << "***\n";
            std::cerr << "*** " << "Error :  Unable to " << aTask << ".\n";
            std::cerr << "***" << std::endl;
            std::cerr << "*** " << "Reason:  " << aCheck<<"\n";
            std::cerr << "***" << std::endl;
            std::istringstream tExceptionMessage( aException.what() );
            std::string tExceptionLine;

            while ( std::getline( tExceptionMessage, tExceptionLine) )
            {
                std::cerr << "***          " << tExceptionLine << std::endl;
                std::cerr << "***" << std::endl;
            }
#ifdef BELFEM_MPI
            if( gComm.size() > 1 )
            {
                std::cerr << "*** " << "Proc:    This error occured on proc " << gComm.rank() << "." << std::endl;
                std::cerr << "***" << std::endl;
            }
#endif
            std::cerr << "*** " << "Where:   This error was encountered inside " << aLocation << "." << std::endl;
            std::cerr << "***" << std::endl;
            std::cerr << "*** ---------------------------------------------------------------------------\n";
            std::cerr << "\n";

#if !defined( NDEBUG ) || defined( DEBUG )
             throw aException;
#else
            exit( -1 );
#endif
        }

//------------------------------------------------------------------------------

        template < typename ... Args >
        inline void
        belfem_assert(
                const std::string & aFile,
                const std::size_t & aLine,
                const std::string & aFunction,
                const std::string & aCheck,
                const Args ...      aArgs )
        {
            std::stringstream tLocation;
            tLocation << belfem::basename( aFile ) << " (line " << aLine << ")";

            std::stringstream tTask;
            tTask << "complete call to function " << aFunction << "()";

            std::stringstream tReason;
            tReason << "Assertion " << aCheck << " failed.";

            // format output message
            std::string tMessage = belfem::sprint( aArgs ... );

            belfem::assert::error(
                    tLocation.str(),
                    tTask.str(),
                    tReason.str(),
                    std::runtime_error( tMessage.c_str() ) );
        }

//------------------------------------------------------------------------------
    } /* namespace assert */
} /* namespace belfem */

//------------------------------------------------------------------------------

#if !defined( NDEBUG ) || defined( DEBUG )
#define BELFEM_ASSERT( aCheck, ... ) \
    do \
    { \
        if ( ! ( aCheck ) ) \
        { \
            belfem::assert::belfem_assert(  \
                __FILE__, \
                __LINE__, \
                __FUNCTION__, \
                #aCheck, \
                __VA_ARGS__ \
                ); \
        } \
    } while ( false )
#else
#define BELFEM_ASSERT( aCheck, ... )
#endif
//---------------------------------------------------------------------------

#define BELFEM_ERROR( aCheck, ... ) \
    do \
    { \
        if ( ! ( aCheck ) ) \
        { \
            belfem::assert::belfem_assert(  \
                __FILE__, \
                __LINE__, \
                __FUNCTION__, \
                #aCheck, \
                __VA_ARGS__ \
                ); \
        } \
    } while ( false )

//------------------------------------------------------------------------------

#endif //BELFEM_ASSERT_HPP