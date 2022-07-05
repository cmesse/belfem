//
// Created by Christian Messe on 2018-12-21.
//

#ifndef BELFEM_CL_COMMUNICATOR_HPP
#define BELFEM_CL_COMMUNICATOR_HPP

#include <string>
#ifdef BELFEM_MPI
#include <mpi.h>
#define COMM_TYPE MPI_Comm
#else
#define COMM_TYPE int
#endif
#include "typedefs.hpp"

#include "cl_Cell.hpp"


namespace  belfem
{
//------------------------------------------------------------------------------
    class Communicator
    {
        std::string mExecutablePath;
        std::string mWorkDir;
        std::string mArguments = "";

        Cell<COMM_TYPE>   mComms;

        proc_t mMyRank = gNoOwner ;
        proc_t mSize = gNoOwner ;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        Communicator(){};

//------------------------------------------------------------------------------

        Communicator(
                int    *argc,
                char ***argv );

//------------------------------------------------------------------------------

        ~Communicator();

//-----------------------------------------------------------------------------

        int
        finalize();

//-----------------------------------------------------------------------------

        COMM_TYPE &
        world();

//-----------------------------------------------------------------------------

        const std::string &
        exec_path();

//-----------------------------------------------------------------------------

        const std::string &
        workdir();

//-----------------------------------------------------------------------------

        const std::string &
        arguments();

//-----------------------------------------------------------------------------

        /**
         * called by arguments object
         * @param aArguments
         */
        void
        set_arguments( const string & aArguments );

//-----------------------------------------------------------------------------

        proc_t
        rank() const;

//-----------------------------------------------------------------------------

        proc_t
        size() const;

//-----------------------------------------------------------------------------
    };
//------------------------------------------------------------------------------

    inline proc_t
    Communicator::rank() const
    {
        return mMyRank ;
    }

//-----------------------------------------------------------------------------

    inline proc_t
    Communicator::size() const
    {
        return mSize ;
    }

//-----------------------------------------------------------------------------
}

// Externally Defined Global Communicator
extern belfem::Communicator gComm;

#endif //BELFEM_CL_COMMUNICATOR_HPP
