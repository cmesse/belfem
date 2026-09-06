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

#include <cstring>
#include <cstdlib>
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
    comm_check( const int aErrorCode )
    {
#ifdef BELFEM_MPI
        // early exit if no error
        if ( aErrorCode == MPI_SUCCESS ) return;

        // Get the implementation-specific error message
        char tErrorString[ MPI_MAX_ERROR_STRING ];
        int tStringLength;
        MPI_Error_string( aErrorCode, tErrorString, &tStringLength );
        string tMessage( tErrorString, tStringLength );

        BELFEM_ERROR( aErrorCode == MPI_SUCCESS, tMessage.c_str() );
#endif
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

    [[noreturn]] void
    comm_abort( const int aErrorCode )
    {
#ifdef BELFEM_MPI
        // the full contract is documented at the declaration in cl_Communicator.hpp:
        // MPI_COMM_WORLD ( never gComm.world() ), lifecycle-guarded, no
        // comm_check, and an unconditional std::abort() fall-through
        int tInitialized = 0 ;
        int tFinalized   = 0 ;

        MPI_Initialized( &tInitialized );
        MPI_Finalized( &tFinalized );

        if ( tInitialized && ! tFinalized )
        {
            MPI_Abort( MPI_COMM_WORLD, aErrorCode );
        }
#endif
        std::abort();
    }

//------------------------------------------------------------------------------

    // THE ORDERING CONTRACT OF THE POINT-TO-POINT FABRIC. Every exchange
    // between a rank pair shares exactly TWO tags: this base tag ( sizes and
    // scalars ) and base+1 ( payloads ). MPI matches per ( source, tag ) in
    // FIFO order, so the ONLY thing that keeps a size from being consumed by
    // the wrong collect is that both ranks execute their exchanges in the
    // same order. Consequences: a send whose matching recv is skipped on the
    // other side does not fail at the send -- it silently poisons the NEXT
    // recv on that tag ( MPI_ERR_TRUNCATE if bigger, garbage if it fits,
    // or a hang if nothing ever matches it ); and any
    // conditional early-out must reproduce the communication pattern of the
    // path it replaces, message for message ( see
    // Postprocessor::recover_fields for a worked example )
    int
    comm_tag( const proc_t aSource, const proc_t aTarget )
    {
        if( gComm.max_tag() == 0 ) return 0;

        proc_t tMin = ( aSource < aTarget ) ? ( aSource ) : ( aTarget );
        proc_t tMax = ( aSource > aTarget ) ? ( aSource ) : ( aTarget );

        return 2 * ( tMax * gComm.size() + tMin )  % gComm.max_tag() ;
    }

//------------------------------------------------------------------------------

    // the contract is documented at the declaration in commtools.hpp; this is
    // the ad-hoc drain probe made permanent
    void
    comm_drain_check( const char * aLabel )
    {
#if defined( BELFEM_MPI ) && BELFEM_ASSERTIONS_ACTIVE
        BELFEM_ASSERT( aLabel != nullptr,
            "comm_drain_check: aLabel must not be null" );

        // no fabric to check in serial, and nothing sane to do before
        // gComm.init() has populated rank and size
        if ( gComm.size() < 2 || gComm.size() == gNoOwner ) return ;

        const proc_t tMyRank   = comm_rank();
        const proc_t tCommSize = comm_size();

        // close the exchange this boundary ends: nobody probes until every
        // rank has finished its own sends and recvs
        comm_barrier();

        // local sweep, stopping at the first stray
        int    tHit       = 0 ;
        proc_t tSource    = 0 ;
        int    tTag       = 0 ;
        int    tIsPayload = 0 ;
        int    tBytes     = 0 ;

        for ( proc_t p = 0; p < tCommSize && tHit == 0; ++p )
        {
            if ( p == tMyRank ) continue ;

            const int tBaseTag = comm_tag( tMyRank, p );

            for ( int tKind = 0; tKind < 2 && tHit == 0; ++tKind )
            {
                int tFlag = 0 ;
                MPI_Status tStatus ;

                comm_check( MPI_Iprobe( p, tBaseTag + tKind, gComm.world(),
                        &tFlag, &tStatus ) );

                if ( tFlag )
                {
                    tHit       = 1 ;
                    tSource    = p ;
                    tTag       = tBaseTag + tKind ;
                    tIsPayload = tKind ;

                    comm_check( MPI_Get_count( &tStatus, MPI_BYTE, &tBytes ) );

                    // MPI_Get_count reports MPI_UNDEFINED when the byte
                    // count does not fit an int
                    if ( tBytes == MPI_UNDEFINED ) tBytes = -1 ;
                }
            }
        }

        // fold the verdict: worst rank wins. This is also the fence that
        // keeps a fast rank from sending the next exchange's first message
        // into a peer's still-running sweep, and it puts EVERY rank on the
        // same error path, so the debug throw policy ( assert.cpp ) holds
        // and no peer is stranded in a barrier
        int tGlobalHit = 0 ;
        allreduce( &tHit, &tGlobalHit, 1 );

        if ( tGlobalHit )
        {
            if ( tHit )
            {
                BELFEM_ERROR( false,
                    "comm_drain_check( %s ): stray message on the fabric: source %d -> rank %d, tag %d ( %s ), %d bytes ( -1 = unknown ). A send in the exchange before this boundary has no matching recv -- see the ordering contract on comm_tag().",
                    aLabel,
                    ( int ) tSource,
                    ( int ) tMyRank,
                    tTag,
                    tIsPayload ? "payload" : "size/scalar",
                    tBytes );
            }
            else
            {
                BELFEM_ERROR( false,
                    "comm_drain_check( %s ): a stray message was detected on another rank -- see its error box.",
                    aLabel );
            }
        }
#else
        ( void ) aLabel ;
#endif
    }

//------------------------------------------------------------------------------

    Cell< int >
    comm_split( const index_t aLength )
    {
        // number of full packages
        index_t tDiv = aLength / gMaxCommChunkLength;

        // remaining size
        index_t tMod = aLength % gMaxCommChunkLength;

        index_t tNumMessages = tDiv + ( tMod > 0 ? 1 : 0 );

        Cell< int > aSteps( tNumMessages, gMaxCommChunkLength );

        // length of the final message
        if ( tMod > 0 ) aSteps( tDiv ) = tMod ;

        return aSteps ;
    }

//------------------------------------------------------------------------------

    index_t
    comm_splitcount( const Vector< index_t > & aLengths, const proc_t aRoot   )
    {
        index_t aCount = 0 ;
        proc_t tCommSize = aLengths.length() ;

        for ( proc_t p=0; p<tCommSize; ++p )
        {
            if ( p == aRoot ) continue;

            // we have to split the message using the same logic
            // as in comm_split to count the number of messages we need

            // number of full packages
            index_t tDiv =  aLengths( p ) / gMaxCommChunkLength;

            // remaining size
            index_t tMod =  aLengths( p ) % gMaxCommChunkLength;

            aCount += tDiv + ( tMod > 0 ? 1 : 0 );
        }
        return aCount ;
    }

    index_t
    comm_splitcount( index_t aLength  )
    {
        index_t aCount = 0 ;

        // we have to split the message using the same logic
        // as in comm_split to count the number of messages we need

        // number of full packages
        index_t tDiv =  aLength / gMaxCommChunkLength;

        // remaining size
        index_t tMod =  aLength % gMaxCommChunkLength;

        aCount += tDiv + ( tMod > 0 ? 1 : 0 );

        return aCount * ( comm_size() - 1 );
    }


//==============================================================================
// STRINGS
//==============================================================================

    void
    broadcast( Cell< string > & aData, const proc_t aRoot )
    {
        // get my id
        proc_t tMyRank = gComm.rank();

        // first, we have to compute the string length
        Vector< index_t > tLengths ;

        if ( tMyRank == aRoot && aData.size() > 0 )
        {
            tLengths.set_size( aData.size() );
            for ( index_t p=0; p<aData.size(); ++p )
            {
                tLengths( p ) = aData( p ).length();
            }
        }

        // distribute lengths of the strings
        broadcast( tLengths, aRoot );

        Vector< char > tBuffer ;

        if ( tMyRank == aRoot )
        {
            // compute number of characters
            index_t tCount = 0 ;
            for ( index_t l : tLengths )
            {
                tCount += l ;
            }

            // assemble the buffer vector
            tBuffer.set_size( tCount );
            tCount = 0 ;
            for ( index_t s=0;s<aData.size(); ++s )
            {
                // copy the string
                std::memcpy( tBuffer.data() + tCount, aData( s ).c_str(), aData( s ).length() );
                tCount += tLengths( s );
            }
        }

        // distribute the buffer to all other procs
        broadcast( tBuffer, aRoot );

        if ( tMyRank != aRoot )
        {
            // disassemble the vector
            index_t tCount = 0 ;

            index_t n = tLengths.length() ;
            aData.set_size( n, "" );

            // build the strings
            for ( index_t s=0;s<n; ++s )
            {
                aData( s ) = string( tBuffer.data() + tCount, tLengths( s ) );
                tCount += tLengths( s );
            }
        }
    }


    void
    send( const string & aMessage, const proc_t aTarget )
    {
        Vector< char > tBuffer( aMessage.size() );
        std::memcpy( tBuffer.data(), aMessage.c_str(), aMessage.length() );
        send( tBuffer, aTarget );
    }

    void
    receive( string & aMessage, const proc_t aSource )
    {
        Vector< char > tBuffer;
        receive( tBuffer, aSource );
        aMessage = string( tBuffer.data(), tBuffer.length() );
    }

}
