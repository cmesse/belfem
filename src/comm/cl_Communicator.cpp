//
// Created by Christian Messe on 2018-12-21.
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

#include <cstdlib>
#include <iostream>
#include <sys/stat.h>

#include <type_traits>  // for type checking is_same

#ifdef BELFEM_PETSC
#include <petscsys.h> // for PetscInt and PetscReal
#endif
#include "cl_Communicator.hpp"

#define BELFEM_INITIALIZE_GLOBALS
#include "globals.hpp"
#undef BELFEM_INITIALIZE_GLOBALS

#include "cl_Vector.hpp"

#ifdef BELFEM_MPI
// Constructor runs before main(), ensuring MPI sees the environment variable
static void __attribute__((constructor)) setup_mpi_binding() {
    setenv("OMPI_MCA_hwloc_base_binding_policy", "none", 0);
    setenv("PRTE_MCA_hwloc_default_binding_policy", "none", 0);
    setenv("HYDRA_BINDING", "none", 0);
}
#endif


namespace belfem
{
//------------------------------------------------------------------------------

    static int gCommunicatorCounter = 0;

//------------------------------------------------------------------------------

    Communicator::Communicator(
            int argc, char ** argv )
    {
        this->init( argc, argv );
    }

//------------------------------------------------------------------------------

    void
    Communicator::set_globals()
    {
        gTbulk = BELFEM_QUIET_NAN ;

        // numerical guards for the resistivity clamp in
        // calculator::MaxwellData::compute_rho(). The defaults make the clamp
        // a no-op; executables may narrow the window
        gRhoMin = 0 ;
        gRhoMax = 1e10 ;

        // Location of the data files, taken from $BELFEM_DATA and pointing at the
        // share directory. Leaving it unset is not an error: a code that needs a
        // data file may set this path itself, for instance from a config file, and
        // the readers fall back to searching relative to the working directory.
        // The error belongs where a file is actually missing, not here.
        const char * tDataPath = std::getenv( "BELFEM_DATA" );

        if( tDataPath != nullptr && *tDataPath != '\0' )
        {
            gBelfemDataPath = tDataPath ;
        }
#ifdef BELFEM_INSTALL_DATADIR
        // the share directory of an installed tree, compiled in from
        // CMAKE_INSTALL_PREFIX. Checked after $BELFEM_DATA, so an explicit
        // setting always wins, and only used when it exists, so a build-tree
        // binary keeps its working-directory fallbacks when nothing is installed
        else
        {
            struct stat tStat ;
            if( ::stat( BELFEM_INSTALL_DATADIR, &tStat ) == 0 && S_ISDIR( tStat.st_mode ) )
            {
                gBelfemDataPath = BELFEM_INSTALL_DATADIR ;
            }
        }
#endif
    }

//------------------------------------------------------------------------------
    void Communicator::init( int & argc, char ** & argv )
    {

        for( int c=0; c<argc; ++c )
        {
            mArguments.push( std::string( argv[ c ] ) );
        }

        for( int c=1; c<argc ; ++c )
        {
            if( c == argc-1)
            {
                mArgumentString += mArguments( c )  ;
            }
            else
            {
                mArgumentString += mArguments( c ) + " ";
            }
        }

        if( gCommunicatorCounter == 0 )
        {
#ifdef BELFEM_PETSC
#ifdef BELFEM_STRUMPACK
            // Initialize MPI with thread support before PetscInitialize
            // PetscInitialize will detect that MPI is already initialized
            this->init_thread( argc, argv );
#else
            int tMpiInit = MPI_Init( &argc, &argv );

            // NOT BELFEM_ERROR. The error reaction ends in comm_abort(), which
            // asks MPI_Initialized / MPI_Finalized before calling MPI_Abort,
            // and after a FAILED MPI_Init the standard leaves MPI in an
            // undefined state -- so that query could be the second fault.
            // Report and die locally instead
            if ( tMpiInit != MPI_SUCCESS )
            {
                std::cerr << "Error while trying to initialize MPI: "
                          << tMpiInit << std::endl ;
                std::abort();
            }
#endif
            PetscErrorCode tStatus = PetscOptionsSetValue(NULL, "-options_left", "0");
            BELFEM_ERROR( tStatus == 0, "Error while tying to set PETSc option '-options_left 0' : %d", ( int ) tStatus );

            tStatus = PetscInitialize( &argc, &argv, NULL, NULL );
            BELFEM_ERROR( tStatus == 0, "Error while tying to initialize MPI/PETSC: %d", ( int ) tStatus );

            // check data types
            BELFEM_ERROR( ( std::is_same< PetscInt, int >::value ) ,
                "<PetscInt> is not identical to datatype <int>" );

            BELFEM_ERROR( ( std::is_same< PetscReal, real >::value ),
                "<PetscReal> is not identical to datatype <real>" );
#elif BELFEM_MPI
#ifdef BELFEM_STRUMPACK
            int tStatus = this->init_thread( argc, argv );
#else
            int tStatus = MPI_Init( &argc, &argv );
#endif
            // NOT BELFEM_ERROR, for the same reason as the PETSc arm above:
            // the error reaction queries MPI_Initialized before MPI_Abort, and
            // a failed initialization leaves MPI undefined, so that query
            // would be the second fault. ( init_thread already terminates on
            // its own failure; this covers the plain MPI_Init branch and any
            // non-success it forwards. )
            if ( tStatus != MPI_SUCCESS )
            {
                std::cerr << "Error while trying to initialize MPI: "
                          << tStatus << std::endl ;
                std::abort();
            }

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
        mExecutablePath = std::string( argv[ 0 ] );

        // save workdir
        const char * tPwd = std::getenv("PWD");
        mWorkDir = tPwd ? tPwd : ".";

#ifdef BELFEM_MPI
        int tVal;
        MPI_Comm_rank( this->world(), &tVal );
        mCommRank = tVal ;

        MPI_Comm_size( this->world(), &tVal );
        mSize = tVal ;

        // determine the number of ranks on this node
        MPI_Comm tNodeComm ;
        MPI_Comm_split_type( this->world(),
                             MPI_COMM_TYPE_SHARED,
                             0,
                             MPI_INFO_NULL,
                             &tNodeComm );
        MPI_Comm_size( tNodeComm, &tVal );
        mNodeSize = tVal ;
        MPI_Comm_free( &tNodeComm );

        int * tMaxTagPtr ;
        MPI_Comm_get_attr( this->world(), MPI_TAG_UB, &tMaxTagPtr, &tVal );
        mMaxTag = *tMaxTagPtr ;

        if ( mCommRank == 0 )
        {
            // Force first-touch on correct node
            /*cpu_set_t tCpuSet ;
            CPU_ZERO( &tCpuSet );
            CPU_SET( 0, &tCpuSet );
            sched_setaffinity( 0, sizeof( cpu_set_t ), &tCpuSet );*/

            // make sure that we don't exceed the proc limit
            proc_t tMaxRank = static_cast<proc_t>( 0.5 * ( std::sqrt( 2.0 * mMaxTag + 1.0 ) + 1.0 ));

            BELFEM_ERROR( mSize <= tMaxRank ,
    "Due to how the comm_tag() function works, BELFEM doesn't support more than %u procs using the current MPI implementation",
    ( unsigned int ) tMaxRank  );
        }

#else
        mCommRank = 0 ;
        mSize = 1 ;
        mMaxTag = 0 ;
#endif

        // Seed random with rank + time
        mRandom.seed( std::random_device{}() + mCommRank );

        this->set_globals();
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
    Cell< std::string > &
    Communicator::arguments()
    {
        return mArguments ;
    }

//-----------------------------------------------------------------------------

    const std::string &
    Communicator::argument_string()
    {
        return mArgumentString;
    }

//-----------------------------------------------------------------------------

    void
    Communicator::set_arguments( const string & aArguments )
    {
        mArgumentString = aArguments ;
    }

//------------------------------------------------------------------------------

    int
    Communicator::finalize()
    {
        // decrement counter
        --gCommunicatorCounter;

        if( gCommunicatorCounter == 0 )
        {
            for ( CommunicationObject * tObj : mObjects )
            {
                if ( tObj != nullptr )
                {
                    tObj->free() ;
                }
            }
        }

#ifdef BELFEM_PETSC

        MPI_Barrier( MPI_COMM_WORLD );

        if( gCommunicatorCounter == 0 )
        {
             PetscFinalize();
             MPI_Finalize();

            // restore the pre-init sentinels so rank() and size() tell the
            // truth after teardown. The error reaction no longer depends on
            // this -- it asks MPI directly -- but the rank line in the error
            // box does, and a public accessor that keeps reporting a rank
            // count for a communicator that no longer exists is its own trap
            mSize     = gNoOwner ;
            mCommRank = gNoOwner ;
        }
        return gCommunicatorCounter;
#elif BELFEM_MPI

        MPI_Barrier( MPI_COMM_WORLD );

        if( gCommunicatorCounter == 0 )
        {
            MPI_Finalize();

            // restore the pre-init sentinels so rank() and size() tell the
            // truth after teardown. The error reaction no longer depends on
            // this -- it asks MPI directly -- but the rank line in the error
            // box does, and a public accessor that keeps reporting a rank
            // count for a communicator that no longer exists is its own trap
            mSize     = gNoOwner ;
            mCommRank = gNoOwner ;
        }
        return gCommunicatorCounter;
#else
        return 0;
#endif
    }

    int
    Communicator::init_thread( int & argc, char ** & argv )
    {
        // Initialize MPI with thread support for STRUMPACK/SLATE compatibility
        int tProvided;
        int aErr = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &tProvided);

        // see the MPI_Init arm: a failed initialization must not be reported
        // through a reaction path that queries MPI
        if ( aErr )
        {
            std::cerr << "Error while trying to initialize MPI with threads: "
                      << aErr << std::endl ;
            std::abort();
        }

        if (tProvided < MPI_THREAD_MULTIPLE) {
            fprintf(stderr, "Warning: MPI_THREAD_MULTIPLE not available (provided level: %d)\n", tProvided);
        }
        return aErr;
    }

    CommunicationObject::CommunicationObject()
    {
        mIndex = gComm.objects().size() ;
        gComm.objects().push( this ) ;
    }

    CommunicationObject::~CommunicationObject()
    {

    }

    void
    CommunicationObject::free()
    {
        gComm.objects()( mIndex ) = nullptr ;
    }

    index_t
    CommunicationObject::index() const
    {
        return mIndex ;
    }

//-----------------------------------------------------------------------------
}