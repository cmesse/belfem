//
// Created by Christian Messe on 2018-12-21.
//
#ifdef BELFEM_PETSC
#include <type_traits>  // for type checking is_same
#include <petscsys.h> // for PetscInt and PetscReal
#endif
#include "cl_Communicator.hpp"

#include "cl_Vector.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    static int gCommunicatorCounter = 0;

//------------------------------------------------------------------------------

    Communicator::Communicator(
            int    *argc,
            char ***argv )
    {
        if( gCommunicatorCounter == 0 )
        {
#ifdef BELFEM_PETSC
            PetscErrorCode tErr = PetscInitialize( argc, argv, NULL, NULL );
            BELFEM_ERROR( ! tErr, "Error while tying to initialize MPI/PETSC: %d", ( int ) tErr );

            // check data types
            BELFEM_ERROR( ( std::is_same< PetscInt, int >::value ) ,
                "<PetscInt> is not identical to datatype <int>" );

            BELFEM_ERROR( ( std::is_same< PetscReal, real >::value ),
                "<PetscReal> is not identical to datatype <real>" );
#elif BELFEM_MPI
            int tErr = MPI_Init( argc, argv );
            BELFEM_ERROR( ! tErr, "Error while tying to initialize MPI: %d", ( int ) tErr );
#endif
        }

        // increment counter
        ++gCommunicatorCounter;

        // create communicator
#ifdef BELFEM_MPI
        mComms.push( MPI_COMM_WORLD );
#else
        mComms.push( 0 );
#endif
        // save executable path
        mExecutablePath = std::string( *argv[ 0 ] );

        // save workdir
        mWorkDir = std::getenv("PWD");

        // assemble argument string
        mArguments = "<no belfem::Arguments object initialized>";

#ifdef BELFEM_MPI
        int tVal;
        MPI_Comm_rank( this->world(), &tVal );
        mMyRank = tVal ;

        MPI_Comm_size( this->world(), &tVal );
        mSize = tVal ;
#else
        mMyRank = 0 ;
        mSize = 1 ;
#endif
    }
//------------------------------------------------------------------------------

    Communicator::~Communicator()

    {
    }

//------------------------------------------------------------------------------


    COMM_TYPE &
    Communicator::world()
    {
        return mComms( 0 );
    }

//-----------------------------------------------------------------------------

    const std::string &
    Communicator::exec_path()
    {
        return mExecutablePath;
    }
//-----------------------------------------------------------------------------

    const std::string &
    Communicator::workdir()
    {
        return mWorkDir;
    }

//-----------------------------------------------------------------------------

    const std::string &
    Communicator::arguments()
    {
        return mArguments;
    }

//-----------------------------------------------------------------------------

    void
    Communicator::set_arguments( const string & aArguments )
    {
        mArguments = aArguments ;
    }

//------------------------------------------------------------------------------

    int
    Communicator::finalize()
    {
        // decrement counter
        --gCommunicatorCounter;

#ifdef BELFEM_PETSC
        MPI_Barrier( MPI_COMM_WORLD );

        if( gCommunicatorCounter == 0 )
        {
             PetscFinalize();
             exit( 0 );
        }
        return gCommunicatorCounter;
#elif BELFEM_MPI

        MPI_Barrier( MPI_COMM_WORLD );

        if( gCommunicatorCounter == 0 )
        {
            MPI_Finalize();
            exit( 0 );
        }
        return gCommunicatorCounter;
#else
        return 0;
#endif
    }

//-----------------------------------------------------------------------------
}