//
// Created by Christian Messe on 2019-01-18.
//

#include "commtools.hpp"
#include "cl_StringList.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    proc_t
    comm_size()
    {
        return gComm.size();
    }

//------------------------------------------------------------------------------

    proc_t
    comm_rank()
    {
        return gComm.rank();
    }

//------------------------------------------------------------------------------

    void
    comm_barrier()
    {
#ifdef BELFEM_MPI
        MPI_Barrier( gComm.world() );
#endif
    }

//------------------------------------------------------------------------------

    void
    create_master_commlist( Vector< proc_t > & aCommList )
    {
        // number of procs
        proc_t tNumProcs = comm_size() ;
        proc_t tMyRank = comm_rank() ;

        if( tNumProcs > 1 )
        {
            // initialize counter
            uint tCount = 0 ;

            // allocate memory
            aCommList.set_size( tNumProcs-1 );

            // loop over all procs
            for( proc_t p = 0; p<tNumProcs; ++p )
            {
                // we do not communicate with ourselves
                if( p != tMyRank )
                {
                    aCommList( tCount++ ) = p ;
                }
            }
        }
    }

//------------------------------------------------------------------------------

    int
    comm_tag( const proc_t aSource, const proc_t aTarget )
    {
        proc_t tMin = ( aSource < aTarget ) ? ( aSource ) : ( aTarget );
        proc_t tMax = ( aSource > aTarget ) ? ( aSource ) : ( aTarget );

        if ( tMin == aTarget )
        {
            return 1024 * ( tMax * gComm.size() + tMin ) + 1;
        }
        else
        {
            return 1024 * ( tMax * gComm.size() + tMin ) + 513;
        }
    }

//------------------------------------------------------------------------------

    Vector< int >
    split_message( const index_t aMessageLength )
    {
        // number of full packages
        index_t tFull = aMessageLength / gMaxArrayLength;

        Vector< int > aSteps( tFull + 1 );

        for ( index_t k = 0; k < tFull; ++k )
        {
            aSteps( k ) = gMaxArrayLength;
        }

        // length of final message
        aSteps( tFull ) = aMessageLength - tFull * gMaxArrayLength;

        return aSteps;
    }

//------------------------------------------------------------------------------

    /*
     * send a string to other procs
     */
    void
    send( const Vector< proc_t >    & aTargets,
                Cell< std::string > & aData )
    {
#ifdef BELFEM_MPI
        std::size_t tNumLines = aData.size() ;

        // vector with lengths
        Vector< std::size_t > tLengths( tNumLines );
        for( std::size_t k=0; k<tNumLines; ++k )
        {
            tLengths( k ) = aData( k ).size() ;
        }

        // send lengths to other procs
        send_same( aTargets, tLengths );

        uint tListSize = aTargets.length() ;

        if( tListSize > 0 )
        {
            // get total number of procs
            proc_t tCommSize = gComm.size();

            // get my id
            proc_t tMyRank = gComm.rank();

            // Allocate memory for status/request vector
            MPI_Status*  tStatus  = ( MPI_Status*  ) alloca( sizeof( MPI_Status  ) * tListSize * tNumLines );
            MPI_Request* tRequest = ( MPI_Request* ) alloca( sizeof( MPI_Request ) * tListSize * tNumLines );


            std::size_t tCount = 0 ;

            // loop over all procs
            for ( uint p = 0; p < tListSize; ++p )
            {
                // get target
                proc_t tTarget = aTargets( p );

                if ( tTarget < tCommSize && tTarget != tMyRank )
                {
                    for( std::size_t l=0 ; l<tNumLines; ++l )
                    {

                        // create communication tag
                        int tCommTag = tTarget * tNumLines + l ;

                        // send data
                        MPI_Isend( aData( l ).c_str(),
                                   tLengths( l ),
                                   MPI_CHAR,
                                   tTarget,
                                   tCommTag,
                                   gComm.world(),
                                   &tRequest[ tCount ] );

                        // wait until send is complete
                        MPI_Wait( &tRequest[ tCount ], &tStatus[ tCount ] );

                        ++tCount ;
                    }
                }
            }
        }
#endif
    }

//------------------------------------------------------------------------------

    /*
     * receive a list of strings from other proc
     */
    void
    receive( const proc_t          aSource,
             Cell< std::string > & aData )
    {
#ifdef BELFEM_MPI
        Vector< std::size_t > tLengths ;

        // get my id
        proc_t tMyRank = gComm.rank();

        // receive lengths from master
        receive( aSource, tLengths );

        // number of strings
        std::size_t tNumLines = tLengths.length() ;

        // allocate temporary container
        char ** tData ;
        tData = ( char** ) malloc( sizeof( *tData ) * ( tNumLines + 1 ) );

        for( std::size_t l=0; l<tNumLines; ++l )
        {
            tData[ l ] = ( char * ) malloc( sizeof( char ) * ( tLengths( l ) + 1 ) );
        }

        // receive data
        MPI_Status*  tStatus  = ( MPI_Status*  ) alloca( sizeof( MPI_Status  ) * tNumLines );
        MPI_Request* tRequest = ( MPI_Request* ) alloca( sizeof( MPI_Request ) * tNumLines );

        for( std::size_t l=0; l<tNumLines; ++l )
        {
            // create communication tag
            int tCommTag = tMyRank * tNumLines + l ;



            MPI_Irecv( tData[ l ],
                       tLengths( l ),
                       MPI_CHAR,
                       aSource,
                       tCommTag,
                       gComm.world(),
                       & tRequest[ l ] );

            // wait until send is complete
            MPI_Wait( &tRequest[ l ], &tStatus[ l ] );
        }

        // convert data to string
        aData.set_size( tNumLines, "" );
        for( std::size_t l=0; l<tNumLines; ++l )
        {
            aData( l ).assign( tData[ l ], tData[ l ] + tLengths( l ) );
        }

        // free container
        for( std::size_t l=0; l<tNumLines; ++l )
        {
            free( tData[ l ] );
        }

        free( tData );
#endif
    }
}