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

#ifndef BELFEM_CL_COMMUNICATOR_HPP
#define BELFEM_CL_COMMUNICATOR_HPP

#include <string>
#ifdef BELFEM_MPI
#include <mpi.h>
#define COMM_TYPE MPI_Comm
#else
#define COMM_TYPE int
#endif

#include <random>
#include "typedefs.hpp"

#include "cl_Cell.hpp"


namespace  belfem
{
//------------------------------------------------------------------------------

    class CommunicationObject
    {
        index_t mIndex ;

    public:

        CommunicationObject();

        ~CommunicationObject();

        virtual void
        free();

        index_t
        index() const ;
    };

    class CommunicationObject ;

//------------------------------------------------------------------------------
    /**
     * @brief Global MPI communicator manager.
     *
     * @ingroup grp_comm
     * @see @ref comm_comm_usage_guide
     */
    class Communicator
    {
        std::string mExecutablePath;
        std::string mWorkDir;
        std::string mArgumentString = "";

        Cell<COMM_TYPE>   mComms;

        proc_t mCommRank = gNoOwner ;
        proc_t mSize = gNoOwner ;
        proc_t mNodeSize = 1 ;
        Cell< string > mArguments ;

        int mMaxTag ;

        std::mt19937 mRandom;

        Cell< CommunicationObject * > mObjects ;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        Communicator(){};

//------------------------------------------------------------------------------

        Communicator(
                int argc, char ** argv );

//------------------------------------------------------------------------------

        ~Communicator();

//------------------------------------------------------------------------------

        Communicator( const Communicator & ) = delete;
        Communicator & operator=( const Communicator & ) = delete;
        Communicator( Communicator && ) = delete;
        Communicator & operator=( Communicator && ) = delete;

//------------------------------------------------------------------------------

        void
        init( int & argc, char ** & argv );

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
        argument_string();
//-----------------------------------------------------------------------------

        Cell< std::string > &
        arguments();

//-----------------------------------------------------------------------------

        /**
         * overrides the argument string assembled by init(); no framework
         * caller today ( the Arguments-object call was retired )
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

        /**
         * returns the number of MPI ranks on this node
         */
        proc_t
        node_size() const;

//-----------------------------------------------------------------------------

        inline std::mt19937 &
        random()
        {
             return mRandom ;
        }

        inline int
        max_tag() const
        {
            return mMaxTag ;
        }

        inline
        Cell< CommunicationObject * > &
        objects();


    private:

        int
        init_thread(  int & argc, char ** & argv );

        void
        set_globals();
    };
//------------------------------------------------------------------------------

    inline proc_t
    Communicator::rank() const
    {
        return mCommRank ;
    }

//-----------------------------------------------------------------------------

    inline proc_t
    Communicator::size() const
    {
        return mSize ;
    }

//-----------------------------------------------------------------------------

    inline proc_t
    Communicator::node_size() const
    {
        return mNodeSize ;
    }

    inline Cell< CommunicationObject * > &
    Communicator::objects()
    {
        return mObjects ;
    }

//-----------------------------------------------------------------------------
}

// Externally Defined Global Communicator
extern belfem::Communicator gComm;

//------------------------------------------------------------------------------
// comm_abort is declared HERE rather than in commtools.hpp so that core's
// assert.cpp can reach it: commtools.hpp pulls cl_Vector/cl_Matrix from
// linalg, which is not on the core module's include path, while this header
// already is ( assert.cpp includes it today ). commtools.hpp includes this
// header, so every commtools consumer sees the declaration as before.
// Definition: commtools.cpp.
//------------------------------------------------------------------------------

namespace belfem
{
    /**
     * \brief Kill the whole job. Never returns.
     *
     * This is the error-path twin of the collectives in commtools.hpp, and it obeys a
     * stricter contract than they do ( see error_abort in assert.cpp, whose
     * body this wrapper absorbed on 2026-08-30 ):
     *
     * - MPI_COMM_WORLD, never gComm.world(): world() indexes mComms( 0 ),
     *   which is EMPTY between MPI_Init and the push in Communicator::init.
     *   On an error path that indexing is a nested assert in debug and UB in
     *   release - and the nested assert would re-enter this function.
     * - guarded by MPI_Initialized / MPI_Finalized, because MPI_Abort is
     *   erroneous before init and after finalize.
     * - falls through to std::abort() unconditionally: MPI_Abort is a "best
     *   attempt" and is not required to return, and in a serial build the
     *   caller ( assert::error ) relies on this function not returning -
     *   returning would resume execution after a failed check.
     * - never routed through comm_check or BELFEM_ERROR: both can land back
     *   here, which is recursion on the error path.
     *
     * Defined in commtools.cpp, like comm_barrier.
     */
    [[noreturn]] void
    comm_abort( const int aErrorCode = 1 );
}


#endif //BELFEM_CL_COMMUNICATOR_HPP
