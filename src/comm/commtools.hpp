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

#ifndef COMMTOOLS_HPP
#define COMMTOOLS_HPP

#include <limits>
#include <type_traits>
#include <string>

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Communicator.hpp"
#include "commtypes.hpp"

extern belfem::Communicator gComm;

namespace belfem
{

    constexpr int gMaxCommChunkLength = 64 * 1024;

//==============================================================================
// HELPERS TO AVOID AMBIGUITIES
//==============================================================================

    template<typename T>
    struct is_scalar : std::true_type {};

    template<typename T>
    struct is_scalar<Vector<T>> : std::false_type {};

    template<typename T>
    struct is_scalar<Matrix<T>> : std::false_type {};

//==============================================================================
// UTILITIES
//==============================================================================

    /**
     * \brief Returns the number of processes in the communicator.
     * \return The number of processes.
     */
    proc_t
    comm_size();

//------------------------------------------------------------------------------

    /**
     * \brief Returns the rank of the current process in the communicator.
     * \return The rank of the current process.
     */
    proc_t
    comm_rank();


//------------------------------------------------------------------------------

    /**
     * \brief Checks the MPI error code and reports any errors.
     * \param aErrorCode The MPI error code to check.
     */
    void
    comm_check( const int aErrorCode );

//------------------------------------------------------------------------------

    /**
     * \brief Synchronizes all processes in the communicator.
     */
    void
    comm_barrier();

//------------------------------------------------------------------------------

    /**
     * \brief Debug-build tripwire for the two-tags-per-pair ordering contract
     *        of the point-to-point fabric ( see comm_tag() in commtools.cpp ).
     *
     * Call at an exchange boundary: a point every rank passes with no exchange
     * in flight. The check synchronizes, probes both fabric tags of every rank
     * pair this rank belongs to ( comm_tag( rank, p ) and +1, never
     * MPI_ANY_TAG -- solver libraries share MPI_COMM_WORLD ), folds the
     * verdict through allreduce, and on a hit raises BELFEM_ERROR on EVERY
     * rank -- detecting ranks name boundary/source/tag/bytes, the others point
     * at the peer's error box. The collective verdict keeps the debug throw
     * policy intact ( assert.cpp: debug throws, on every rank alike ) without
     * stranding peers in a barrier, and doubles as the fence that stops a fast
     * rank from sending the next exchange's first message into a peer's
     * still-running sweep. So a stray dies loudly at the boundary where it was
     * born instead of poisoning an innocent recv an arbitrary distance later.
     *
     * aLabel names the boundary in the diagnostic and must not be null.
     *
     * Collective wherever assertions are active; a call site reachable by a
     * subset of ranks is itself a deadlock, same contract as comm_barrier.
     * Compiles to a no-op in release builds and does nothing in serial runs.
     * Best-effort by nature: a barrier does not flush unmatched messages, so a
     * stray still in flight when the sweep runs can escape it -- absence of an
     * abort is not proof of a clean fabric -- and a stray carrying a tag
     * outside the pair's two fabric tags is invisible to it.
     */
    void
    comm_drain_check( const char * aLabel );

//------------------------------------------------------------------------------

    /**
     * \brief Generates a unique tag for MPI communication between two processes.
     * \param aSource The source process rank.
     * \param aTarget The target process rank.
     * \return A unique tag for communication.
     */
    int
    comm_tag( const proc_t aSource, const proc_t aTarget );

//------------------------------------------------------------------------------

    /**
     * \brief Splits a message into chunks for MPI communication.
     * \param aLength The length of the message to split.
     * \return A cell containing the sizes of each chunk.
     */
    Cell< int >
    comm_split( const index_t aLength );

//------------------------------------------------------------------------------

    /**
     * \brief Calculates the total number of chunks needed for a set of message lengths.
     * \param aLengths Vector of message lengths.
     * \param aRoot The root process rank.
     * \return The total number of chunks.
     */
    index_t
    comm_splitcount( const Vector< index_t > & aLengths, const proc_t aRoot );

    index_t
    comm_splitcount( index_t aLength  );

//------------------------------------------------------------------------------

    /**
     * \brief Broadcasts a scalar value to all processes.
     * \tparam T The type of the scalar.
     * \param aMessage The scalar to broadcast.
     * \param aRoot The root process rank (default is 0).
     */
    template< typename T >
    void
    broadcast( T & aMessage, const proc_t aRoot=0,
        typename std::enable_if<is_scalar<T>::value>::type* = nullptr )
    {
#ifdef BELFEM_MPI
        if( gComm.size() > 1 )
        {
            // Ensure T is a scalar type.
            // Note: std::complex<T> is not supported here; use
            // send/receive for complex scalar communication.
            BELFEM_ASSERT(  std::is_arithmetic<T>::value,
                "Can only broadcast arithmetic types." );

            // call the bcast command
            comm_t tCommType = comm_type<T>();
            comm_check( MPI_Bcast(
                    &aMessage,
                    1,
                    tCommType,
                    aRoot,
                    gComm.world() ) );
        }
#endif
    }

    // for raw array with known length
    template< typename T >
    void
    broadcast( T * aMessage, const proc_t aRoot, const proc_t aLength )
    {
#ifdef BELFEM_MPI
        if( gComm.size() > 1 )
        {
            // Ensure T is a scalar type.
            // Note: std::complex<T> is not supported here; use
            // send/receive for complex array communication.
            BELFEM_ASSERT(  std::is_arithmetic<T>::value,
                "Can only broadcast arithmetic types." );

            // call the bcast command
            comm_t tCommType = comm_type<T>();

            comm_check( MPI_Bcast(
                    aMessage,
                    aLength,
                    tCommType,
                    aRoot,
                    gComm.world() ) );
        }
#endif
    }

//------------------------------------------------------------------------------

    /**
     * \brief Collective MAX-reduction visible on every rank.
     *
     * The reduction operation is MPI_MAX and deliberately not a parameter:
     * the wrapper exists to keep vendor tokens out of the solver files, and
     * the sites it serves ( solver verdicts, test exit codes ) all fold
     * "worst rank wins". A SUM consumer extends this file rather than passing an
     * MPI_Op through it; the MIN sibling is allreduce_min() below.
     *
     * The count is MPI's own int, not int_t: the MPI-3 C binding takes
     * int, and under BELFEM_INT64 an int_t count would narrow silently
     * above 2^31-1. Payloads here are a few elements - bulk data goes
     * through the chunked share / receive pair, as everywhere else.
     *
     * Collective over gComm.world(); every rank must call it, after
     * gComm.init() like every other collective in this file.
     */
    template < typename T >
    void
    allreduce( const T * aSend, T * aRecv, const int aCount )
    {
        // MPI_MAX is meaningless for the complex types comm_type<> also
        // serves - refuse them at compile time
        static_assert( std::is_arithmetic< T >::value,
            "allreduce: MAX-reduction is defined for arithmetic types only" );

        BELFEM_ASSERT( aCount >= 0, "allreduce: negative count %d", aCount );

#ifdef BELFEM_MPI
        comm_check( MPI_Allreduce(
                aSend,
                aRecv,
                aCount,
                comm_type< T >(),
                MPI_MAX,
                gComm.world() ) );
#else
        // serial semantic of a one-rank allreduce: the identity copy
        if ( aSend != aRecv )
        {
            for ( int i = 0; i < aCount; ++i )
            {
                aRecv[ i ] = aSend[ i ];
            }
        }
#endif
    }

//------------------------------------------------------------------------------

    /**
     * \brief Collective MIN-reduction visible on every rank.
     *
     * The MIN sibling of allreduce(), for the sites that fold "the most
     * constrained rank decides" -- a per-process memory budget, where a
     * rank that could not measure ( 0 ) must pull every rank to 0 rather
     * than be outvoted. Same contract as allreduce(): arithmetic types,
     * MPI's own int for the count, collective over gComm.world().
     */
    template < typename T >
    void
    allreduce_min( const T * aSend, T * aRecv, const int aCount )
    {
        static_assert( std::is_arithmetic< T >::value,
            "allreduce_min: MIN-reduction is defined for arithmetic types only" );

        BELFEM_ASSERT( aCount >= 0, "allreduce_min: negative count %d", aCount );

#ifdef BELFEM_MPI
        comm_check( MPI_Allreduce(
                aSend,
                aRecv,
                aCount,
                comm_type< T >(),
                MPI_MIN,
                gComm.world() ) );
#else
        if ( aSend != aRecv )
        {
            for ( int i = 0; i < aCount; ++i )
            {
                aRecv[ i ] = aSend[ i ];
            }
        }
#endif
    }


//==============================================================================
// SCALARS
//==============================================================================

    /**
     * \brief Sends a scalar value to a specific target process.
     * \tparam T The type of the scalar.
     * \param aData The scalar to send.
     * \param aTarget The target process rank.
     */
    template< typename T >
    void
    send(  const T aData, const proc_t aTarget=0, typename std::enable_if<is_scalar<T>::value>::type* = nullptr  )
    {
#ifdef BELFEM_MPI
        // get my id
        proc_t tMyRank = gComm.rank();

        if( aTarget < comm_size() && tMyRank != aTarget )
        {
            // status/request handlers
            MPI_Status  tStatus;
            MPI_Request tRequest;

            comm_check( MPI_Isend( &aData,
                       1,
                       comm_type< T >(),
                       aTarget,
                       comm_tag( aTarget, tMyRank ),
                       gComm.world(),
                       & tRequest ) );


            // wait until receive is complete
            comm_check(  MPI_Wait( &tRequest, &tStatus ) );
        }
#endif
    }

//------------------------------------------------------------------------------

    /**
     * \brief Receives a scalar value from a specific source process.
     * \tparam T The type of the scalar.
     * \param aData The scalar to receive.
     * \param aSource The source process rank (default is 0).
     */
    template< typename T >
    void
    receive(  T & aData, const proc_t aSource=0, typename std::enable_if<is_scalar<T>::value>::type* = nullptr )
    {
#ifdef BELFEM_MPI
        // get my id
        proc_t tMyRank = gComm.rank();

        // get the number of procs
        proc_t tCommSize = gComm.size();


        if( aSource < tCommSize && tMyRank != aSource )
        {
            // status/request handlers
            MPI_Status  tStatus;
            MPI_Request tRequest;

            // call MPI command
            comm_check(  MPI_Irecv( &aData,
                       1,
                       comm_type< T >(),
                       aSource,
                       comm_tag( aSource, tMyRank ),
                       gComm.world(),
                       & tRequest ) );

            // wait until receive is complete
            comm_check( MPI_Wait( &tRequest, &tStatus ) );
        }
#endif
    }

//==============================================================================
// RAW ARRAYS
//==============================================================================

    /**
     * \brief Sends a raw array to a specified target process using MPI.
     *
     * This function sends a pre-allocated raw array \p aData of length \p aLength to the process
     * with rank \p aTarget using non-blocking MPI operations. The array is split into chunks
     * for efficient transmission.
     *
     * \tparam T The type of the array elements.
     * \param aData Pointer to the pre-allocated array to send.
     * \param aLength The number of elements in the array.
     * \param aTarget The rank of the target process in the MPI communicator.
     * \note Requires \p aData to be allocated with at least \p aLength elements.
     */
    template< typename T >
    void
    send( T * aData, const index_t aLength, const proc_t aTarget )
    {
#ifdef BELFEM_MPI
        // get the total number of procs
        proc_t tCommSize = gComm.size();

        // get my rank
        proc_t tMyRank   = gComm.rank();

        if ( aTarget < tCommSize && aTarget != tMyRank )
        {
            // status/request handlers
            MPI_Status  tSizeStatus;
            MPI_Request tSizeRequest;

            // compute the commtag
            int tCommTag = comm_tag( tMyRank, aTarget );

            // send the length
            comm_check(  MPI_Isend( &aLength,
                       1,
                       comm_type< index_t >(),
                       aTarget,
                       tCommTag++,
                       gComm.world(),
                       & tSizeRequest ) );

            // wait until receive is complete
            comm_check( MPI_Wait( &tSizeRequest, &tSizeStatus ) );

            // nothing to do if the vector is empty
            if ( aLength == 0 ) return ;

            // compute the chunks for this message
            Cell< int > tChunkSizes = comm_split( aLength );

            // Allocate memory for status/request vector
            MPI_Status*  tStatus  = ( MPI_Status*  ) malloc( sizeof( MPI_Status  ) * tChunkSizes.size() );
            MPI_Request* tRequest = ( MPI_Request* ) malloc( sizeof( MPI_Request ) * tChunkSizes.size() );

            // initialitze the counter
            index_t tCount = 0 ;

            // offset in data container
            index_t tOffset = 0 ;

            // get the communication type
            comm_t tCommType = comm_type< T >();

            for ( index_t c : tChunkSizes )
            {
                // send data
                comm_check( MPI_Isend(
                        &aData[ tOffset ],
                        c,
                        tCommType,
                        aTarget,
                        tCommTag,
                        gComm.world(),
                        & tRequest[ tCount++] ) );

                // increment offset
                tOffset+= c;
            }

            // wait until send is complete
            comm_check( MPI_Waitall( tCount, tRequest, tStatus ) );

            free( tStatus );
            free( tRequest );
        }
#endif
    }

    /**
     * \brief Receives a raw array from a specified source process using MPI.
     *
     * This function receives a raw array into \p aData from the process with rank \p aSource
     * using non-blocking MPI operations. The received length is stored in \p aLength, and the
     * array is processed in chunks. The array must be pre-allocated with sufficient space.
     *
     * \tparam T The type of the array elements.
     * \param aData Pointer to the pre-allocated array to receive data into.
     * \param aLength Reference to the variable that will be updated with the received length.
     * \param aSource The rank of the source process in the MPI communicator.
     * \note Requires \p aData to be allocated with at least the initial \p aLength elements.
     * \warning Raises BELFEM_ERROR ( active in release too ) if the received length exceeds the allocated space.
     */
    template< typename T >
    void
    receive( T * aData, index_t & aLength, const proc_t aSource )
    {
#ifdef BELFEM_MPI

        // get the total number of procs
        proc_t tCommSize = gComm.size();

        // get my rank
        proc_t tMyRank   = gComm.rank();

        if ( aSource < tCommSize && aSource != tMyRank )
        {
            // status/request handlers
            MPI_Status  tSizeStatus;
            MPI_Request tSizeRequest;

            // compute the commtag
            int tCommTag = comm_tag( aSource, tMyRank );

            // receive the length
            index_t tSize = 0;

            comm_check(  MPI_Irecv( &tSize,
                       1,
                       comm_type< index_t >(),
                       aSource,
                       tCommTag++,
                       gComm.world(),
                       & tSizeRequest ) );

            // wait until receive is complete
            comm_check( MPI_Wait( &tSizeRequest, &tSizeStatus ) );

            // make sure that length is properly defined
            BELFEM_ERROR( tSize <= aLength, "Datastream too large ( length %u but expect <=%u ).",
                ( unsigned int ) tSize, ( unsigned int ) aLength  );

            // return new value
            aLength = tSize;

            // nothing to do if the vector is empty
            if ( aLength == 0 ) return ;

            // compute the chunks for this message
            Cell< int > tChunkSizes = comm_split( tSize );

            // Allocate memory for status/request vector
            MPI_Status*  tStatus  = ( MPI_Status*  ) malloc( sizeof( MPI_Status  ) * tChunkSizes.size() );
            MPI_Request* tRequest = ( MPI_Request* ) malloc( sizeof( MPI_Request ) * tChunkSizes.size() );

            // initialitze the counter
            index_t tCount = 0 ;

            // offset in data container
            index_t tOffset = 0 ;

            // get the communication type
            comm_t tCommType = comm_type< T >();

            for ( index_t c : tChunkSizes )
            {
                // send data
                comm_check( MPI_Irecv(
                        &aData[ tOffset ],
                        c,
                        tCommType,
                        aSource,
                        tCommTag,
                        gComm.world(),
                        & tRequest[ tCount++] ) );

                // increment offset
                tOffset+= c;
            }

            // wait until send is complete
            comm_check( MPI_Waitall( tCount, tRequest, tStatus ) );

            free( tStatus );
            free( tRequest );
        }
#endif
    }

//==============================================================================
// CELLS
//==============================================================================

    /**
     * \brief Sends a cell to a specific target process.
     * \tparam T The type of the cell elements.
     * \param aData The cell to send.
     * \param aTarget The target process rank.
     */
    template< typename T >
    void
    send( Cell< T > & aData, const proc_t aTarget=0 )
    {
#ifdef BELFEM_MPI
        // get my id
        proc_t tMyRank = gComm.rank();

        // nothing to do here
        if ( tMyRank == aTarget ) return ;

        // get the length of the cell
        index_t tSize = aData.size();

        // status/request handlers
        MPI_Status  tSizeStatus;
        MPI_Request tSizeRequest;

        // compute the commtag
        int tCommTag = comm_tag( tMyRank, aTarget );

        // send the length
        comm_check(  MPI_Isend( &tSize,
                   1,
                   comm_type< index_t >(),
                   aTarget,
                   tCommTag++,
                   gComm.world(),
                   & tSizeRequest ) );

        // wait until receive is complete
        comm_check( MPI_Wait( &tSizeRequest, &tSizeStatus ) );

        // nothing to do if the cell is empty
        if ( tSize == 0 ) return ;

        // get the communication type
        comm_t tCommType = comm_type< T >();

        // compute the chunks for this message
        Cell< int > tChunkSizes = comm_split( aData.size() );

        // offset in data container
        index_t tOffset = 0 ;

        // Allocate memory for status/request cell
        MPI_Status*  tStatus  = ( MPI_Status*  ) malloc( sizeof( MPI_Status  ) * tChunkSizes.size() );
        MPI_Request* tRequest = ( MPI_Request* ) malloc( sizeof( MPI_Request ) * tChunkSizes.size() );

        // initialitze the counter
        index_t tCount = 0 ;

        // get the raw pointer of the cell
        const T * tData = aData.data();

        for ( index_t c : tChunkSizes )
        {
            // send data
            comm_check( MPI_Isend( &tData[ tOffset ],
                        c,
                        tCommType,
                        aTarget,
                        tCommTag,
                        gComm.world(),
                        &tRequest[ tCount++ ] ) );

            tOffset+= c;
        }

        // wait until send is complete
        comm_check( MPI_Waitall( tCount, tRequest, tStatus ) );

        free( tStatus );
        free( tRequest );
#endif
    }

//------------------------------------------------------------------------------

    /**
     * \brief Receives a cell from a specific source process.
     * \tparam T The type of the cell elements.
     * \param aData The cell to receive.
     * \param aSource The source process rank (default is 0).
     */
    template< typename T >
    void
    receive( Cell< T > & aData, const proc_t aSource=0 )
    {
#ifdef BELFEM_MPI

        // get my id
        proc_t tMyRank = gComm.rank();

        // nothing to do here
        if ( tMyRank == aSource ) return ;

        // obtain the length of the cell
        index_t tSize = 0 ;
        receive( tSize, aSource );

        // allocate memory for the cell
        aData.set_size( tSize, 0 );

        // nothing more to do if the cell is empty
        if ( tSize == 0 ) return ;

        // get the communication type
        comm_t tCommType = comm_type< T >();

        // compute the chunks for this message
        Cell< int > tChunkSizes = comm_split( aData.size() );

        // offset in data container
        index_t tOffset = 0 ;

        // Inside the length distrubution, we use comm_tag too.
        // For safety, the tag for the messages is incremented.
        int tCommTag = comm_tag( aSource, tMyRank ) + 1 ;

        // Allocate memory for status/request cell
        MPI_Status*  tStatus  = ( MPI_Status*  ) malloc( sizeof( MPI_Status  ) * tChunkSizes.size() );
        MPI_Request* tRequest = ( MPI_Request* ) malloc( sizeof( MPI_Request ) * tChunkSizes.size() );

        // initialitze the counter
        index_t tCount = 0 ;

        // get the raw pointer of the cell
        T * tData = aData.data();

        for ( index_t c : tChunkSizes )
        {
            // send data
            comm_check( MPI_Irecv(
                        &tData[ tOffset ],
                        c,
                        tCommType,
                        aSource,
                        tCommTag,
                        gComm.world(),
                        &tRequest[ tCount++ ] ) );

            tOffset+= c;
        }

        // wait until send is complete
        comm_check( MPI_Waitall( tCount, tRequest, tStatus ) );

        free( tStatus );
        free( tRequest );
#endif
    }

//------------------------------------------------------------------------------

    /**
     * \brief Broadcasts a cell to all processes.
     * \tparam T The type of the cell elements.
     * \param aData The cell to broadcast.
     * \param aRoot The root process rank (default is 0).
     */
    template< typename T >
    void
    broadcast( Cell< T >  & aData, const proc_t aRoot=0 )
    {
#ifdef BELFEM_MPI

        // get my id
        proc_t tMyRank = gComm.rank();

        // Get length of cell
        index_t tSize = tMyRank == aRoot ? aData.size() : 0 ;

        // communicate the cell size
        MPI_Request tSizeRequest ;
        comm_check( MPI_Ibcast(
                            & tSize,
                            1,
                            comm_type< index_t >(),
                            aRoot,
                            gComm.world(),
                            & tSizeRequest ) );

        // wait until the size is received
        MPI_Status tSizeStatus ;
        comm_check( MPI_Wait( &tSizeRequest, &tSizeStatus ) );

        // allocate memory for the cell
        if ( tMyRank != aRoot ) aData.set_size( tSize );

        if ( tSize == 0 ) return ;

        // broadcast the cell data
        MPI_Request tDataRequest ;
        comm_check( MPI_Ibcast(
            aData.data(),
            tSize,
            comm_type< T >(),
            aRoot,
            gComm.world(),
            & tDataRequest ) );

        // wait until data is received
        MPI_Status tDataStatus ;
        comm_check( MPI_Wait( &tDataRequest, &tDataStatus ) );

#endif
    }


//------------------------------------------------------------------------------

    /**
    * \brief Distributes elements of a cell to other processes.
    * \tparam T The type of the cell elements.
    * \param aData The cell to distribute.
    */
    template< typename T >
    void
    distribute( Cell< T > & aData )
    {
#ifdef BELFEM_MPI
        // get my id
        proc_t tMyRank = gComm.rank();

        // get the number of procs
        proc_t tCommSize = gComm.size();

        BELFEM_ASSERT( static_cast< proc_t >( aData.size() )== tCommSize,
        "Length of cell does not match ( is %u, expect commsize %u ).",
        ( unsigned int ) aData.size(), ( unsigned int ) tCommSize );

        // Allocate memory for status/request cell
        MPI_Status*  tStatus  = ( MPI_Status*  ) malloc( sizeof( MPI_Status  ) * tCommSize );
        MPI_Request* tRequest = ( MPI_Request* ) malloc( sizeof( MPI_Request ) * tCommSize );

        for( proc_t p=0; p<tCommSize; ++p )
        {
            if( p == tMyRank )
            {
                tRequest[ p ] = MPI_REQUEST_NULL;
                continue ;
            }

            // send data
            comm_check( MPI_Isend( &aData( p ),
                       1,
                       comm_type< T >(),
                       p,
                       comm_tag( tMyRank, p ),
                       gComm.world(),
                       &tRequest[ p ] ) );

        }

        // wait until send is complete
        comm_check( MPI_Waitall( tCommSize, tRequest, tStatus ) );

        // tidy up memory
        free( tStatus );
        free( tRequest );
#endif
    }

//------------------------------------------------------------------------------

    /**
     * \brief Collects scalar values from all processes into a cell.
     * \tparam T The type of the scalar.
     * \param aData The cell to store collected values.
     * \param aMyValue The value from the current process (default is 0).
     */
    template< typename T >
    void
    collect( Cell< T > & aData, const T aMyValue = 0 )
    {
#ifdef BELFEM_MPI
        // get my id
        proc_t tMyRank = gComm.rank();

        // get the number of procs
        proc_t tCommSize = gComm.size();

        aData.set_size( tCommSize, 0 );

        // Allocate memory for status/request cell
        MPI_Status*  tStatus  = ( MPI_Status*  ) malloc( sizeof( MPI_Status  ) * tCommSize );
        MPI_Request* tRequest = ( MPI_Request* ) malloc( sizeof( MPI_Request ) * tCommSize );

        for( proc_t p=0; p<tCommSize; ++p )
        {
            if( p == tMyRank )
            {
                tRequest[ p ] = MPI_REQUEST_NULL;
                aData( p ) = aMyValue ;
                continue ;
            }

            // send data
            comm_check( MPI_Irecv( &aData( p ),
                       1,
                       comm_type< T >(),
                       p,
                       comm_tag( tMyRank, p ),
                       gComm.world(),
                       &tRequest[ p ] ) );


        }
        // wait until send is complete
        comm_check( MPI_Waitall( tCommSize, tRequest, tStatus ) );

        // tidy up memory
        free( tStatus );
        free( tRequest );
#endif
    }

//==============================================================================
// VECTORS
//==============================================================================

    /**
     * \brief Sends a vector to a specific target process.
     * \tparam T The type of the vector elements.
     * \param aData The vector to send.
     * \param aTarget The target process rank.
     */
    template< typename T >
    void
    send( Vector< T > & aData, const proc_t aTarget=0 )
    {
#ifdef BELFEM_MPI
        // get my id
        proc_t tMyRank = gComm.rank();

        // nothing to do here
        if ( tMyRank == aTarget ) return ;

        // get the length of the vector
        index_t tSize = aData.length();

        // status/request handlers
        MPI_Status  tSizeStatus;
        MPI_Request tSizeRequest;

        // compute the commtag
        int tCommTag = comm_tag( tMyRank, aTarget );

        // send the length
        comm_check(  MPI_Isend( &tSize,
                   1,
                   comm_type< index_t >(),
                   aTarget,
                   tCommTag++,
                   gComm.world(),
                   & tSizeRequest ) );

        // wait until receive is complete
        comm_check( MPI_Wait( &tSizeRequest, &tSizeStatus ) );

        // nothing to do if the vector is empty
        if ( tSize == 0 ) return ;

        // get the communication type
        comm_t tCommType = comm_type< T >();

        // compute the chunks for this message
        Cell< int > tChunkSizes = comm_split( aData.length() );

        // offset in data container
        index_t tOffset = 0 ;

        // Allocate memory for status/request vector
        MPI_Status*  tStatus  = ( MPI_Status*  ) malloc( sizeof( MPI_Status  ) * tChunkSizes.size() );
        MPI_Request* tRequest = ( MPI_Request* ) malloc( sizeof( MPI_Request ) * tChunkSizes.size() );

        // initialitze the counter
        index_t tCount = 0 ;

        // get the raw pointer of the vector
        const T * tData = aData.data();

        for ( index_t c : tChunkSizes )
        {
            // send data
            comm_check( MPI_Isend( &tData[ tOffset ],
                        c,
                        tCommType,
                        aTarget,
                        tCommTag,
                        gComm.world(),
                        &tRequest[ tCount++ ] ) );

            tOffset+= c;
        }

        // wait until send is complete
        comm_check( MPI_Waitall( tCount, tRequest, tStatus ) );

        free( tStatus );
        free( tRequest );
#endif
    }

//------------------------------------------------------------------------------

    /**
     * \brief Receives a vector from a specific source process.
     * \tparam T The type of the vector elements.
     * \param aData The vector to receive.
     * \param aSource The source process rank (default is 0).
     */
    template< typename T >
    void
    receive( Vector< T > & aData, const proc_t aSource=0 )
    {
#ifdef BELFEM_MPI

        // get my id
        proc_t tMyRank = gComm.rank();

        // nothing to do here
        if ( tMyRank == aSource ) return ;

        // obtain the length of the vector
        index_t tSize = 0 ;
        receive( tSize, aSource );

        // allocate memory for the vector
        aData.set_size( tSize, 0 );

        // nothing more to do if the vector is empty
        if ( tSize == 0 ) return ;

        // get the communication type
        comm_t tCommType = comm_type< T >();

        // compute the chunks for this message
        Cell< int > tChunkSizes = comm_split( aData.length() );

        // offset in data container
        index_t tOffset = 0 ;

        // Inside the length distrubution, we use comm_tag too.
        // For safety, the tag for the messages is incremented.
        int tCommTag = comm_tag( aSource, tMyRank ) + 1 ;

        // Allocate memory for status/request vector
        MPI_Status*  tStatus  = ( MPI_Status*  ) malloc( sizeof( MPI_Status  ) * tChunkSizes.size() );
        MPI_Request* tRequest = ( MPI_Request* ) malloc( sizeof( MPI_Request ) * tChunkSizes.size() );

        // initialitze the counter
        index_t tCount = 0 ;

        // get the raw pointer of the vector
        T * tData = aData.data();

        for ( index_t c : tChunkSizes )
        {
            // send data
            comm_check( MPI_Irecv(
                        &tData[ tOffset ],
                        c,
                        tCommType,
                        aSource,
                        tCommTag,
                        gComm.world(),
                        &tRequest[ tCount++ ] ) );

            tOffset+= c;
        }

        // wait until send is complete
        comm_check( MPI_Waitall( tCount, tRequest, tStatus ) );

        free( tStatus );
        free( tRequest );
#endif
    }

//------------------------------------------------------------------------------

    /**
     * \brief Broadcasts a vector to all processes.
     * \tparam T The type of the vector elements.
     * \param aData The vector to broadcast.
     * \param aRoot The root process rank (default is 0).
     */
    template< typename T >
    void
    broadcast( Vector< T >  & aData, const proc_t aRoot=0 )
    {
#ifdef BELFEM_MPI

        // get my id
        proc_t tMyRank = gComm.rank();

        // Get length of vector
        index_t tSize = tMyRank == aRoot ? aData.length() : 0 ;

        // communicate the vector size
        MPI_Request tSizeRequest ;
        comm_check( MPI_Ibcast(
                            & tSize,
                            1,
                            comm_type< index_t >(),
                            aRoot,
                            gComm.world(),
                            & tSizeRequest ) );

        // wait until the size is received
        MPI_Status tSizeStatus ;
        comm_check( MPI_Wait( &tSizeRequest, &tSizeStatus ) );

        // allocate memory for the vector
        if ( tMyRank != aRoot ) aData.set_size( tSize );

        if ( tSize == 0 ) return ;

        // broadcast the vector data
        MPI_Request tDataRequest ;
        comm_check( MPI_Ibcast(
            aData.data(),
            tSize,
            comm_type< T >(),
            aRoot,
            gComm.world(),
            & tDataRequest ) );

        // wait until data is received
        MPI_Status tDataStatus ;
        comm_check( MPI_Wait( &tDataRequest, &tDataStatus ) );

#endif
    }

//------------------------------------------------------------------------------

    /**
     * \brief Distributes elements of a vector to other processes.
     * \tparam T The type of the vector elements.
     * \param aData The vector to distribute.
     */
    template< typename T >
    void
    distribute( Vector< T > & aData )
    {
#ifdef BELFEM_MPI
        // get my id
        proc_t tMyRank = gComm.rank();

        // get the number of procs
        proc_t tCommSize = gComm.size();

        BELFEM_ASSERT( static_cast< proc_t>( aData.length() )== tCommSize,
        "Length of vector does not match ( is %u, expect commsize %u ).",
        ( unsigned int ) aData.length(), ( unsigned int ) tCommSize );

        // Allocate memory for status/request vector
        MPI_Status*  tStatus  = ( MPI_Status*  ) malloc( sizeof( MPI_Status  ) * tCommSize );
        MPI_Request* tRequest = ( MPI_Request* ) malloc( sizeof( MPI_Request ) * tCommSize );

        for( proc_t p=0; p<tCommSize; ++p )
        {
            if( p == tMyRank )
            {
                tRequest[ p ] = MPI_REQUEST_NULL;
                continue ;
            }

            // send data
            comm_check( MPI_Isend( &aData( p ),
                       1,
                       comm_type< T >(),
                       p,
                       comm_tag( tMyRank, p ),
                       gComm.world(),
                       &tRequest[ p ] ) );

        }

        // wait until send is complete
        comm_check( MPI_Waitall( tCommSize, tRequest, tStatus ) );

        // tidy up memory
        free( tStatus );
        free( tRequest );
#endif
    }

//------------------------------------------------------------------------------

    /**
     * \brief Collects scalar values from all processes into a vector.
     * \tparam T The type of the scalar.
     * \param aData The vector to store collected values.
     * \param aMyValue The value from the current process (default is 0).
     */
    template< typename T >
    void
    collect( Vector< T > & aData, const T aMyValue = 0 )
    {
#ifdef BELFEM_MPI
        // get my id
        proc_t tMyRank = gComm.rank();

        // get the number of procs
        proc_t tCommSize = gComm.size();

        aData.set_size( tCommSize, 0 );

        // Allocate memory for status/request vector
        MPI_Status*  tStatus  = ( MPI_Status*  ) malloc( sizeof( MPI_Status  ) * tCommSize );
        MPI_Request* tRequest = ( MPI_Request* ) malloc( sizeof( MPI_Request ) * tCommSize );

        for( proc_t p=0; p<tCommSize; ++p )
        {
            if( p == tMyRank )
            {
                tRequest[ p ] = MPI_REQUEST_NULL;
                aData( p ) = aMyValue ;
                continue ;
            }

            // send data
            comm_check( MPI_Irecv( &aData( p ),
                       1,
                       comm_type< T >(),
                       p,
                       comm_tag( tMyRank, p ),
                       gComm.world(),
                       &tRequest[ p ] ) );


        }
        // wait until send is complete
        comm_check( MPI_Waitall( tCommSize, tRequest, tStatus ) );

        // tidy up memory
        free( tStatus );
        free( tRequest );
#endif
    }

//------------------------------------------------------------------------------

    /**
     * \brief Distributes a collection of vectors to other processes.
     * \tparam T The type of the vector elements.
     * \param aData The collection of vectors to distribute.
     */
    template< typename T >
    void
    distribute( Cell< Vector< T > > & aData )
    {
#ifdef BELFEM_MPI

        // get my id
        proc_t tMyRank = gComm.rank();

        // get the number of procs
        proc_t tCommSize = gComm.size();

        BELFEM_ASSERT( static_cast< proc_t>( aData.size() )== tCommSize,
            "Length of data container does not match ( is %u, expect commsize %u ).",
            ( unsigned int ) aData.size(), ( unsigned int ) tCommSize );

        // populate the vector lengths
        Vector< index_t > tSizes( tCommSize, 0 );
        for ( proc_t p=0; p<tCommSize; ++p )
        {
            tSizes( p ) = aData( p ).length();
        }
        index_t tCount = comm_splitcount( tSizes, tMyRank );

        // send the vector length to the other procs
        distribute( tSizes );

        // Allocate memory for request/status vector
        MPI_Request* tRequest = ( MPI_Request * ) malloc( sizeof( MPI_Request ) * tCount );
        MPI_Status*  tStatus  = ( MPI_Status *  ) malloc( sizeof( MPI_Status  ) * tCount );


        // resetting the counter
        tCount = 0 ;

        // get the communication type
        comm_t tCommType = comm_type< T >();

        for( proc_t p=0; p<tCommSize; ++p )
        {
            // no message if self or there is nothing to send
            if ( p == tMyRank || aData( p ).length() == 0 ) continue ;

            // get the raw pointer of the container
            const T * tData = aData( p ).data();

            // compute the chunks for this message
            Cell< int > tChunkSizes = comm_split( aData( p ).length() );

            // offset in data container
            index_t tOffset = 0 ;

            // Inside the length distrubution, we use comm_tag too.
            // For safety, the tag for the messages is incremented.
            int tCommTag = comm_tag( tMyRank, p ) + 1 ;

            for ( index_t c : tChunkSizes )
            {
                // send data
                comm_check( MPI_Isend( &tData[ tOffset ],
                            c,
                            tCommType,
                            p,
                            tCommTag,
                            gComm.world(),
                            &tRequest[ tCount++ ] ) );

                tOffset+= c;
            }
        }

        // wait until send is complete
        comm_check( MPI_Waitall( tCount, tRequest, tStatus ) );

        // tidy up memory
        free( tStatus );
        free( tRequest );

#endif
    }

      /**
     * \brief Distributes a collection of cells to other processes.
     * \tparam T The type of the cell elements.
     * \param aData The collection of cells to distribute.
     */
    template< typename T >
    void
    distribute( Cell< Cell< T > > & aData )
    {
#ifdef BELFEM_MPI

        // get my id
        proc_t tMyRank = gComm.rank();

        // get the number of procs
        proc_t tCommSize = gComm.size();

        BELFEM_ASSERT( static_cast< proc_t>( aData.size() )== tCommSize,
            "Length of data container does not match ( is %u, expect commsize %u ).",
            ( unsigned int ) aData.size(), ( unsigned int ) tCommSize );

        // populate the vector lengths
        Vector< index_t > tSizes( tCommSize, 0 );
        for ( proc_t p=0; p<tCommSize; ++p )
        {
            tSizes( p ) = aData( p ).size();
        }
        index_t tCount = comm_splitcount( tSizes, tMyRank );

        // send the vector length to the other procs
        distribute( tSizes );

        // Allocate memory for request/status vector
        MPI_Request* tRequest = ( MPI_Request * ) malloc( sizeof( MPI_Request ) * tCount );
        MPI_Status*  tStatus  = ( MPI_Status *  ) malloc( sizeof( MPI_Status  ) * tCount );


        // resetting the counter
        tCount = 0 ;

        // get the communication type
        comm_t tCommType = comm_type< T >();

        for( proc_t p=0; p<tCommSize; ++p )
        {
            // no message if self or there is nothing to send
            if ( p == tMyRank || aData( p ).size() == 0 ) continue ;

            // get the raw pointer of the container
            const T * tData = aData( p ).data();

            // compute the chunks for this message
            Cell< int > tChunkSizes = comm_split( aData( p ).size() );

            // offset in data container
            index_t tOffset = 0 ;

            // Inside the length distrubution, we use comm_tag too.
            // For safety, the tag for the messages is incremented.
            int tCommTag = comm_tag( tMyRank, p ) + 1 ;

            for ( index_t c : tChunkSizes )
            {
                // send data
                comm_check( MPI_Isend( &tData[ tOffset ],
                            c,
                            tCommType,
                            p,
                            tCommTag,
                            gComm.world(),
                            &tRequest[ tCount++ ] ) );

                tOffset+= c;
            }
        }

        // wait until send is complete
        comm_check( MPI_Waitall( tCount, tRequest, tStatus ) );

        // tidy up memory
        free( tStatus );
        free( tRequest );

#endif
    }

//------------------------------------------------------------------------------

    /**
     * \brief Distributes a raw vector based on lengths and offsets
     * \tparam T The type of the array elements.
     * \param aData Pointer to the contiguous array to distribute.
     * \param aOffsets offsets per proc
     */
    template< typename T, typename U >
    void
    distribute( const T * aData, const Vector< U > & aOffsets )
    {
#ifdef BELFEM_MPI

        // get my id
        proc_t tMyRank = gComm.rank();

        // get the number of procs
        proc_t tCommSize = gComm.size();

        BELFEM_ASSERT( static_cast< proc_t>( aOffsets.length() ) == tCommSize + 1,
            "Length of offset container does not match ( is %u, expect commsize %u ).",
            ( unsigned int ) aOffsets.length(), ( unsigned int ) tCommSize + 1 );


        // populate the vector lengths
        Vector< index_t > tSizes( tCommSize, 0 );

        for ( proc_t p=0; p<tCommSize; ++p )
        {
            tSizes( p ) = aOffsets( p+1 ) - aOffsets( p );
        }
        index_t tCount = comm_splitcount( tSizes, tMyRank );

        // send the vector length to the other procs
        distribute( tSizes );

        // Allocate memory for request/status vector
        MPI_Request* tRequest = ( MPI_Request * ) malloc( sizeof( MPI_Request ) * tCount );
        MPI_Status*  tStatus  = ( MPI_Status *  ) malloc( sizeof( MPI_Status  ) * tCount );


        // resetting the counter
        tCount = 0 ;

        // get the communication type (strip const qualifier if present)
        comm_t tCommType = comm_type< typename std::remove_const<T>::type >();

        for( proc_t p=0; p<tCommSize; ++p )
        {
            // no message if self or there is nothing to send
            if ( p == tMyRank || tSizes( p ) == 0 ) continue ;

            // compute the chunks for this message
            Cell< int > tChunkSizes = comm_split(  tSizes( p ) );

            // offset in data container
            index_t tOffset = aOffsets( p ) ;

            // Inside the length distrubution, we use comm_tag too.
            // For safety, the tag for the messages is incremented.
            int tCommTag = comm_tag( tMyRank, p ) + 1 ;

            for ( index_t c : tChunkSizes )
            {
                // send data
                comm_check( MPI_Isend( aData + tOffset,
                            c,
                            tCommType,
                            p,
                            tCommTag,
                            gComm.world(),
                            &tRequest[ tCount++ ] ) );

                tOffset+= c;
            }
        }

        // wait until send is complete
        comm_check( MPI_Waitall( tCount, tRequest, tStatus ) );

        // tidy up memory
        free( tStatus );
        free( tRequest );

#endif
    }

    template< typename T, typename U >
    void
    collect( T * aData, const Vector< U > & aOffsets )
    {
#ifdef BELFEM_MPI
         // get my id
        proc_t tMyRank = gComm.rank();

        // get the number of procs
        proc_t tCommSize = gComm.size();

        // obtain the vector lengths from non-root ranks
        Vector< index_t > tSizes ;
        collect( tSizes );
        // rank 0's data is already in aData; set its size from offsets for consistency
        tSizes( 0 ) = aOffsets( 0 );

        index_t tCount = 0 ;

        for ( proc_t p=1; p<tCommSize; ++p )
        {
            tCount += comm_split( tSizes( p ) ).size();
        }

        // Allocate memory for status/request vector
        MPI_Status*  tStatus  = ( MPI_Status*  ) malloc( sizeof( MPI_Status  ) * tCount );
        MPI_Request* tRequest = ( MPI_Request* ) malloc( sizeof( MPI_Request ) * tCount );

        // resetting the counter
        tCount = 0 ;

        // get the communication type
        comm_t tCommType = comm_type< T >();

        for( proc_t p=1; p<tCommSize; ++p )
        {
            // get the length
            index_t tSize = tSizes( p );

            if ( tSize == 0 ) continue;

            // compute the chunks for this message
            Cell< int > tChunkSizes = comm_split( tSize );

            // offset in data container
            index_t tOffset = aOffsets( p ) ;

            // Inside the length distribution, we use comm_tag too.
            // For safety, the tag for the messages is incremented.
            int tCommTag = comm_tag( tMyRank, p ) + 1 ;

            for ( index_t c : tChunkSizes )
            {
                // receive data
                comm_check( MPI_Irecv( &aData[ tOffset ],
                            c,
                            tCommType,
                            p,
                            tCommTag,
                            gComm.world(),
                            &tRequest[ tCount++ ] ) );

                tOffset+= c;
            }
        }

        // wait until send is complete
        comm_check( MPI_Waitall( tCount, tRequest, tStatus ) );

        // tidy up memory
        free( tStatus );
        free( tRequest );
#endif
    }
//------------------------------------------------------------------------------


    /**
     * \brief Collects vectors from all processes into a cell.
     * \tparam T The type of the vector elements.
     * \param aData The cell to store collected vectors.
     * \param aMyData The vector contributed by this process.
     */
    template< typename T >
    void
    collect( Cell< Vector< T > > & aData, Vector< T > aMyData={} )
    {
#ifdef BELFEM_MPI

        // get my id
        proc_t tMyRank = gComm.rank();

        // get the number of procs
        proc_t tCommSize = gComm.size();

        // allocate memory
        aData.set_size(  tCommSize, {} );

        // obtain the vector lengths
        Vector< index_t > tSizes ;
        collect( tSizes );

        index_t tCount = 0 ;
        for ( proc_t p=0; p<tCommSize; ++p )
        {
            if ( p == tMyRank )
            {
				if( aMyData.length() > 0 )
				{
					aData( p ).vector_data() = std::move( aMyData.vector_data() ) ;
				}
				continue ;
			}

            // allocate memory of vector
            aData( p ).set_size( tSizes( p ), 0 );

            tCount += comm_split( tSizes( p ) ).size();
        }

        // Allocate memory for status/request vector
        MPI_Status*  tStatus  = ( MPI_Status*  ) malloc( sizeof( MPI_Status  ) * tCount );
        MPI_Request* tRequest = ( MPI_Request* ) malloc( sizeof( MPI_Request ) * tCount );

        // resetting the counter
        tCount = 0 ;

        // get the communication type
        comm_t tCommType = comm_type< T >();

        for( proc_t p=0; p<tCommSize; ++p )
        {
            // get the length
            index_t tSize = tSizes( p );

            if ( p == tMyRank || tSize == 0 ) continue;

            // get the raw data container
            T * tData = aData( p ).data() ;

            // compute the chunks for this message
            Cell< int > tChunkSizes = comm_split( tSize );

            // offset in data container
            index_t tOffset = 0 ;

            // Inside the length distribution, we use comm_tag too.
            // For safety, the tag for the messages is incremented.
            int tCommTag = comm_tag( tMyRank, p ) + 1 ;

            for ( index_t c : tChunkSizes )
            {
                // receive data
                comm_check( MPI_Irecv( &tData[ tOffset ],
                            c,
                            tCommType,
                            p,
                            tCommTag,
                            gComm.world(),
                            &tRequest[ tCount++ ] ) );

                tOffset+= c;
            }
        }

        // wait until send is complete
        comm_check( MPI_Waitall( tCount, tRequest, tStatus ) );

        // tidy up memory
        free( tStatus );
        free( tRequest );

#endif
    }

//------------------------------------------------------------------------------

    template< typename T >
    void
    share( Vector< T > & aData )
    {
#ifdef BELFEM_MPI
        // get my id
        proc_t tMyRank = gComm.rank();

        proc_t tCommSize = gComm.size();

        index_t tSize = aData.length();

        MPI_Request* tRequest = ( MPI_Request * ) malloc( sizeof( MPI_Request ) * ( tCommSize-1 ) );
        MPI_Status*  tStatus  = ( MPI_Status *  ) malloc( sizeof( MPI_Status  ) * ( tCommSize-1 ) );

        // get the communication type
        comm_t tIndex_t = comm_type< index_t >();

        // get the communication type
        comm_t tCommType = comm_type< T >();

        index_t tCount = 0 ;

        for ( proc_t p=0; p<tCommSize; ++p )
        {
            if ( p == tMyRank ) continue;

            int tCommTag = comm_tag( tMyRank, p );

            comm_check( MPI_Isend( &tSize,
                        1,
                        tIndex_t,
                        p,
                        tCommTag,
                        gComm.world(),
                        &tRequest[ tCount ++ ] ) );
        }

        // wait until send is complete
        comm_check( MPI_Waitall( tCount, tRequest, tStatus ) );

        // tidy up memory
        free( tStatus );
        free( tRequest );

        tCount = comm_splitcount( tSize );

        tRequest = ( MPI_Request * ) malloc( sizeof( MPI_Request ) * tCount );
        tStatus  = ( MPI_Status *  ) malloc( sizeof( MPI_Status  ) * tCount );

        if ( aData.length() > 0 )
        {
            tCount = 0 ;
            for( proc_t p=0; p<tCommSize; ++p )
            {
                // no message if self or there is nothing to send
                if ( p == tMyRank ) continue ;

                // get the raw pointer of the container
                const T * tData = aData.data();

                // compute the chunks for this message
                Cell< int > tChunkSizes = comm_split( aData.length() );

                // offset in data container
                index_t tOffset = 0 ;

                // Inside the length distrubution, we use comm_tag too.
                // For safety, the tag for the messages is incremented.
                int tCommTag = comm_tag( tMyRank, p ) + 1 ;

                for ( index_t c : tChunkSizes )
                {
                    // send data
                    comm_check( MPI_Isend( &tData[ tOffset ],
                                c,
                                tCommType,
                                p,
                                tCommTag,
                                gComm.world(),
                                &tRequest[ tCount++ ] ) );

                    tOffset+= c;
                }
            }

            // wait until send is complete
            comm_check( MPI_Waitall( tCount, tRequest, tStatus ) );
        }



        // tidy up memory
        free( tStatus );
        free( tRequest );
#endif
    }

//------------------------------------------------------------------------------

    template< typename T >
    void
    share( Cell< T > & aData )
    {
#ifdef BELFEM_MPI
        // get my id
        proc_t tMyRank = gComm.rank();

        proc_t tCommSize = gComm.size();

        index_t tSize = aData.size();

        MPI_Request* tRequest = ( MPI_Request * ) malloc( sizeof( MPI_Request ) * ( tCommSize-1 ) );
        MPI_Status*  tStatus  = ( MPI_Status *  ) malloc( sizeof( MPI_Status  ) * ( tCommSize-1 ) );

        // get the communication type
        comm_t tIndex_t = comm_type< index_t >();

        // get the communication type
        comm_t tCommType = comm_type< T >();

        index_t tCount = 0 ;

        for ( proc_t p=0; p<tCommSize; ++p )
        {
            if ( p == tMyRank ) continue;

            int tCommTag = comm_tag( tMyRank, p );

            comm_check( MPI_Isend( &tSize,
                        1,
                        tIndex_t,
                        p,
                        tCommTag,
                        gComm.world(),
                        &tRequest[ tCount ++ ] ) );
        }

        // wait until send is complete
        comm_check( MPI_Waitall( tCount, tRequest, tStatus ) );

        // tidy up memory
        free( tStatus );
        free( tRequest );

        tCount = comm_splitcount( tSize );

        tRequest = ( MPI_Request * ) malloc( sizeof( MPI_Request ) * tCount );
        tStatus  = ( MPI_Status *  ) malloc( sizeof( MPI_Status  ) * tCount );

        if ( aData.size() > 0 )
        {
            tCount = 0 ;
            for( proc_t p=0; p<tCommSize; ++p )
            {
                // no message if self or there is nothing to send
                if ( p == tMyRank ) continue ;

                // get the raw pointer of the container
                const T * tData = aData.data();

                // compute the chunks for this message
                Cell< int > tChunkSizes = comm_split( aData.size() );

                // offset in data container
                index_t tOffset = 0 ;

                // Inside the length distrubution, we use comm_tag too.
                // For safety, the tag for the messages is incremented.
                int tCommTag = comm_tag( tMyRank, p ) + 1 ;

                for ( index_t c : tChunkSizes )
                {
                    // send data
                    comm_check( MPI_Isend( &tData[ tOffset ],
                                c,
                                tCommType,
                                p,
                                tCommTag,
                                gComm.world(),
                                &tRequest[ tCount++ ] ) );

                    tOffset+= c;
                }
            }

            // wait until send is complete
            comm_check( MPI_Waitall( tCount, tRequest, tStatus ) );
        }



        // tidy up memory
        free( tStatus );
        free( tRequest );
#endif
    }

//==============================================================================
// MATRICES
//==============================================================================

    /**
    * \brief Broadcasts a matrix to all processes.
    * \tparam T The type of the matrix elements.
    * \param aData The matrix to broadcast.
    * \param aRoot The root process rank (default is 0).
    */
    template< typename T >
    void
    broadcast( Matrix< T >  & aData, const proc_t aRoot=0 )
    {
#ifdef BELFEM_MPI

        // get my id
        proc_t tMyRank = gComm.rank();

        // row, columns and transfer length of matrix
        // We transmit spacing()*n_cols rather than n_rows*n_cols so that the
        // raw contiguous buffer (including backend padding, e.g. Blaze
        // SIMD alignment) is transferred as-is. This avoids disassembling
        // and reassembling the matrix. Both ranks use the same binary, so
        // padding layout is identical for a given (rows, cols) pair.
        // Never use capacity() here: a matrix that shrank keeps its old,
        // larger allocation, and transmitting that count overflows the
        // exact-fit buffer on the receiving side.
        index_t tSize[ 3 ];

        tSize[ 0 ] = aData.n_rows();
        tSize[ 1 ] = aData.n_cols();
        tSize[ 2 ] = aData.spacing() * aData.n_cols();

        // broadcast the matrix dimensions and transfer length
        MPI_Request tSizeRequest ;
        comm_check( MPI_Ibcast(
            & tSize,
            3,
            comm_type< index_t >(),
            aRoot,
            gComm.world(),
            & tSizeRequest ) );

        // wait until size is received
        MPI_Status tSizeStatus ;
        comm_check( MPI_Wait( &tSizeRequest, &tSizeStatus ) );

        // allocate memory for the matrix
        if ( tMyRank != aRoot ) aData.set_size( tSize[ 0 ], tSize[ 1 ] );

        // nothing to do if the matrix is empty
        if( tSize[ 0 ] == 0 || tSize[ 1 ] == 0 ) return ;

        BELFEM_ERROR( aData.capacity() >= tSize[ 2 ],
            "broadcast( Matrix ) : transfer length %lu exceeds local buffer capacity %lu",
            ( long unsigned int ) tSize[ 2 ],
            ( long unsigned int ) aData.capacity() );

        // broadcast is unchunked: the count must fit into MPI's int
        BELFEM_ERROR( tSize[ 2 ] <= ( index_t ) std::numeric_limits< int >::max(),
            "broadcast( Matrix ) : matrix too large for unchunked broadcast, use send/receive" );

        // broadcast the matrix data
        MPI_Request tDataRequest ;
        comm_check( MPI_Ibcast(
            aData.data(),
            tSize[ 2 ],
            comm_type< T >(),
            aRoot,
            gComm.world(),
            & tDataRequest ) );

        // wait until data is received
        MPI_Status tDataStatus ;
        comm_check( MPI_Wait( &tDataRequest, &tDataStatus ) );

#endif
    }
    
    template< typename T >
    void
    send( Matrix< T > & aData, const proc_t aTarget=0 )
    {
#ifdef BELFEM_MPI

        // get my id
        proc_t tMyRank = gComm.rank();

        // nothing to do here
        if ( tMyRank == aTarget ) return ;

        // row, columns and transfer length of matrix.
        // We transmit spacing()*n_cols rather than n_rows*n_cols so that the
        // raw contiguous buffer (including backend padding, e.g. Blaze
        // SIMD alignment) is transferred as-is. This avoids disassembling
        // and reassembling the matrix. Both ranks use the same binary, so
        // padding layout is identical for a given (rows, cols) pair.
        // Never use capacity() here: a matrix that shrank keeps its old,
        // larger allocation, and transmitting that count overflows the
        // exact-fit buffer on the receiving side.
        index_t tSize[ 3 ];

        tSize[ 0 ] = aData.n_rows();
        tSize[ 1 ] = aData.n_cols();
        tSize[ 2 ] = aData.spacing() * aData.n_cols();

        // status/request handlers
        MPI_Status  tSizeStatus;
        MPI_Request tSizeRequest;

        // compute the commtag
        int tCommTag = comm_tag( tMyRank, aTarget );

        // send the length
        comm_check(  MPI_Isend( &tSize,
                   3,
                   comm_type< index_t >(),
                   aTarget,
                   tCommTag++,
                   gComm.world(),
                   & tSizeRequest ) );

        // wait until receive is complete
        comm_check( MPI_Wait( &tSizeRequest, &tSizeStatus ) );

        // nothing to do if the matrix is empty
        if( tSize[ 0 ] == 0 || tSize[ 1 ] == 0 ) return ;

        // get the communication type
        comm_t tCommType = comm_type< T >();

        // compute the chunks for this message
        Cell< int > tChunkSizes = comm_split( tSize[ 2 ] );

        // offset in data container
        index_t tOffset = 0 ;

        // Allocate memory for status/request vector
        MPI_Status*  tStatus  = ( MPI_Status*  ) malloc( sizeof( MPI_Status  ) * tChunkSizes.size() );
        MPI_Request* tRequest = ( MPI_Request* ) malloc( sizeof( MPI_Request ) * tChunkSizes.size() );

        // initialitze the counter
        index_t tCount = 0 ;

        // get raw pointer of matrix
        const T * tData =  aData.data();

        for ( index_t c : tChunkSizes )
        {
            // send data
            comm_check( MPI_Isend( &tData[ tOffset ],
                        c,
                        tCommType,
                        aTarget,
                        tCommTag,
                        gComm.world(),
                        &tRequest[ tCount++ ] ) );

            tOffset+= c;
        }

        // wait until send is complete
        comm_check( MPI_Waitall( tCount, tRequest, tStatus ) );

        free( tStatus );
        free( tRequest );
#endif
    }

//------------------------------------------------------------------------------

    /**
     * \brief Receives a matrix from a specific source process.
     * \tparam T The type of the matrix elements.
     * \param aData The matrix to receive.
     * \param aSource The source process rank (default is 0).
     */
    template< typename T >
    void
    receive( Matrix< T > & aData, const proc_t aSource=0 )
    {
#ifdef BELFEM_MPI

        // get my id
        proc_t tMyRank = gComm.rank();

        // nothing to do here
        if ( tMyRank == aSource ) return ;

        // status/request handlers
        MPI_Status  tSizeStatus;
        MPI_Request tSizeRequest;

        // compute the commtag
        int tCommTag = comm_tag( tMyRank, aSource );

        // tSize[2] is the sender's spacing()*n_cols, not n_rows*n_cols. We
        // transfer the raw contiguous buffer including backend padding
        // (e.g. Blaze SIMD alignment) as-is, avoiding matrix disassembly/
        // reassembly. Both ranks use the same binary, so padding layout is
        // identical for a given (rows, cols) pair, and after set_size below
        // our own buffer holds at least spacing()*n_cols elements.
        index_t tSize[ 3 ];

        // receive the dimensions
        comm_check(  MPI_Irecv( tSize,
                   3,
                   comm_type< index_t >(),
                   aSource,
                   tCommTag++,
                   gComm.world(),
                   & tSizeRequest ) );

        // wait until receive is complete
        comm_check( MPI_Wait( &tSizeRequest, &tSizeStatus ) );

        aData.set_size( tSize[ 0 ], tSize[ 1 ] );

        // nothing more to do if the matrix is empty
        if ( tSize[ 0 ] == 0 || tSize[ 1 ] == 0 ) return ;

        BELFEM_ERROR( aData.capacity() >= tSize[ 2 ],
            "receive( Matrix ) : transfer length %lu exceeds local buffer capacity %lu",
            ( long unsigned int ) tSize[ 2 ],
            ( long unsigned int ) aData.capacity() );

        // get the communication type
        comm_t tCommType = comm_type< T >();

        // compute the chunks for this message
        Cell< int > tChunkSizes = comm_split( tSize[ 2 ] );

        // offset in data container
        index_t tOffset = 0 ;

        // Allocate memory for request/status vector
        MPI_Request* tRequest = ( MPI_Request* ) malloc( sizeof( MPI_Request ) * tChunkSizes.size() );
        MPI_Status*  tStatus  = ( MPI_Status*  ) malloc( sizeof( MPI_Status  ) * tChunkSizes.size() );

        // initialitze the counter
        index_t tCount = 0 ;

        // get the raw pointer of the vector
        T * tData = aData.data();

        for ( index_t c : tChunkSizes )
        {
            // send data
            comm_check( MPI_Irecv( &tData[ tOffset ],
                        c,
                        tCommType,
                        aSource,
                        tCommTag,
                        gComm.world(),
                        &tRequest[ tCount++ ] ) );

            tOffset+= c;
        }

        // wait until send is complete
        comm_check( MPI_Waitall( tCount, tRequest, tStatus ) );

        free( tStatus );
        free( tRequest );
#endif
    }

//------------------------------------------------------------------------------

    /**
     * \brief Distributes a collection of matrices to other processes.
     * \tparam T The type of the matrix elements.
     * \param aData The collection of matrices to distribute.
     */
    template< typename T >
    void
    distribute( Cell< Matrix< T > > & aData )
    {
#ifdef BELFEM_MPI

        // get my id
        proc_t tMyRank = gComm.rank();

        // get the number of procs
        proc_t tCommSize = gComm.size();

        BELFEM_ASSERT( static_cast< proc_t>( aData.size() )== tCommSize,
            "Length of data container does not match ( is %u, expect commsize %u ).",
            ( unsigned int ) aData.size(), ( unsigned int ) tCommSize );

        // populate the matrix lengths
        index_t * tSizes = ( index_t *  ) malloc( sizeof( index_t ) * tCommSize * 3 );
        index_t tCount = 0 ;
        for ( proc_t p=0; p<tCommSize; ++p )
        {
            tSizes[ tCount++ ] = aData( p ).n_rows();
            tSizes[ tCount++ ] = aData( p ).n_cols();

            // transfer length: padded footprint of the current shape,
            // NOT capacity() ( which can be stale-large after a shrink )
            tSizes[ tCount++ ] = aData( p ).spacing() * aData( p ).n_cols();
        }

        MPI_Request* tSizeRequest = ( MPI_Request* ) malloc( sizeof( MPI_Request ) * tCommSize );
        MPI_Status*  tSizeStatus  = ( MPI_Status*  ) malloc( sizeof( MPI_Status  ) * tCommSize );

        // send the data
        index_t tOffset = 0 ;
        for ( proc_t p=0; p<tCommSize; ++p )
        {
            if ( p == tMyRank )
            {
                tSizeRequest[ p ] = MPI_REQUEST_NULL;
                tOffset += 3 ;
                continue;
            }

            // send data
            comm_check( MPI_Isend( &tSizes[ tOffset ],
                        3,
                        comm_type< index_t >(),
                        p,
                        comm_tag( tMyRank, p ),
                        gComm.world(),
                        &tSizeRequest[ p ] ) );

            tOffset += 3 ;
        }

        // wait until send is complete
        comm_check( MPI_Waitall( tCommSize, tSizeRequest, tSizeStatus ) );
        free( tSizeRequest );
        free( tSizeStatus );

        // count the chunk sizes
        tOffset = 2 ;
        tCount = 0 ;
        for ( proc_t p=0; p<tCommSize; ++p )
        {
            // no message if self or there is nothing to send
            if ( p == tMyRank || aData( p ).n_rows() == 0 || aData( p ).n_cols() == 0 )
            {
                tOffset += 3 ;
                continue ;
            }
            tCount += comm_split( tSizes[ tOffset ] ).size();
            tOffset += 3 ;
        }

        // clear memory
        free( tSizes );


        // Allocate memory for status/request vector
        MPI_Status*  tStatus  = ( MPI_Status *  ) malloc( sizeof( MPI_Status  ) * tCount );
        MPI_Request* tRequest = ( MPI_Request * ) malloc( sizeof( MPI_Request ) * tCount );

        // resetting the counter
        tCount = 0 ;

        // get the communication type
        comm_t tCommType = comm_type< T >();

        for( proc_t p=0; p<tCommSize; ++p )
        {
            // no message to self or if there is nothing to send
            if ( p == tMyRank || aData( p ).n_rows() == 0 || aData( p ).n_cols() == 0 )
            {
                continue;
            }

            // get the raw pointer of the contaoner
            const T * tData = aData( p ).data();

            // compute the chunks for this message ( same transfer length
            // as announced in tSizes above — never capacity() )
            Cell< int > tChunkSizes = comm_split( aData( p ).spacing() * aData( p ).n_cols() );

            // offset in data container
            index_t tDataOffset = 0 ;
            // Inside the length distrubution, we use comm_tag too.
            // For safety, the tag for the messages is incremented.
            int tCommTag = comm_tag( tMyRank, p ) + 1 ;

            for ( index_t c : tChunkSizes )
            {
                // send data
                comm_check( MPI_Isend(
                            &tData[ tDataOffset ],
                            c,
                            tCommType,
                            p,
                            tCommTag,
                            gComm.world(),
                            &tRequest[ tCount++ ] ) );

                tDataOffset+= c;
            }
        }

        // wait until send is complete
        comm_check( MPI_Waitall( tCount, tRequest, tStatus ) );

        // tidy up memory
        free( tStatus );
        free( tRequest );

#endif
    }

//------------------------------------------------------------------------------

    /**
     * \brief Collects matrices from all processes into a cell.
     * \tparam T The type of the matrix elements.
     * \param aData The cell to store collected matrices.
     */
    template< typename T >
    void
    collect( Cell< Matrix< T > > & aData )
    {
#ifdef BELFEM_MPI

        // get my id
        proc_t tMyRank = gComm.rank();

        // get the number of procs
        proc_t tCommSize = gComm.size();

        // allocate memory
        aData.set_size(  tCommSize, {} );

        // obtain the matrix lengths
        index_t * tSizes = ( index_t *  ) malloc( sizeof( index_t ) * tCommSize * 3 );


        MPI_Request* tSizeRequest = ( MPI_Request* ) malloc( sizeof( MPI_Request ) * tCommSize );
        MPI_Status*  tSizeStatus  = ( MPI_Status*  ) malloc( sizeof( MPI_Status  ) * tCommSize );

        index_t tOffset = 0 ;

        for ( proc_t p=0; p<tCommSize; ++p )
        {
            if ( p == tMyRank )
            {
                tSizeRequest[ p ] = MPI_REQUEST_NULL;
                tSizes[ tOffset++ ] = 0 ;
                tSizes[ tOffset++ ] = 0 ;
                tSizes[ tOffset++ ] = 0 ;
                continue ;
            }

            // receive data
            comm_check( MPI_Irecv( &tSizes[ tOffset ],
                        3,
                        comm_type< index_t >(),
                        p,
                        comm_tag( tMyRank, p ),
                        gComm.world(),
                        &tSizeRequest[ p ] ) );

            tOffset += 3 ;
        }

        // wait until send is complete
        comm_check( MPI_Waitall( tCommSize, tSizeRequest, tSizeStatus ) );

        free( tSizeStatus );
        free( tSizeRequest );

        index_t tCount = 0 ;
        tOffset = 0 ;

        for ( proc_t p=0; p<tCommSize; ++p )
        {
            // skip if self
            if ( p == tMyRank )
            {
                tOffset += 3 ;
                continue ;
            }

            aData( p ).set_size( tSizes[ tOffset ], tSizes[ tOffset+1 ] );

            // skip if empty
            if ( tSizes[ tOffset ] == 0 || tSizes[ tOffset+1 ] == 0 )
            {
                tOffset += 3 ;
                continue ;
            }

            // count memory needs
            tCount += comm_split( tSizes[ tOffset + 2  ] ).size();
            tOffset += 3 ;
        }

        // Allocate memory for status/request vector
        MPI_Status*  tStatus  = ( MPI_Status*  ) malloc( sizeof( MPI_Status  ) * tCount );
        MPI_Request* tRequest = ( MPI_Request* ) malloc( sizeof( MPI_Request ) * tCount );


        // get the communication type
        comm_t tCommType = comm_type< T >();

        // resetting the counter
        tCount = 0 ;

        tOffset = 0 ;

        for( proc_t p=0; p<tCommSize; ++p )
        {
            // no message to self or if there is nothing to receive
            if ( p == tMyRank )
            {
                tOffset += 3 ;
                continue;
            }
            if (  aData( p ).n_rows() == 0 || aData( p ).n_cols() == 0 )
            {
                tOffset += 3 ;
                continue ;
            }

            BELFEM_ERROR( aData( p ).capacity() >= tSizes[ tOffset + 2 ],
                "collect( Matrix ) : transfer length %lu from proc %u exceeds local buffer capacity %lu",
                ( long unsigned int ) tSizes[ tOffset + 2 ],
                ( unsigned int ) p,
                ( long unsigned int ) aData( p ).capacity() );

            // compute the chunks for this message
            Cell< int > tChunkSizes = comm_split( tSizes[ tOffset + 2 ] );

            // get the raw data container
            T * tData = aData( p ).data() ;

            // offset in data container
            index_t tDataOffset = 0 ;

            // Inside the length distribution, we use comm_tag too.
            // For safety, the tag for the messages is incremented.
            int tCommTag = comm_tag( tMyRank, p ) + 1 ;

            for ( index_t c : tChunkSizes )
            {
                // receive data
                comm_check( MPI_Irecv( &tData[ tDataOffset ],
                            c,
                            tCommType,
                            p,
                            tCommTag,
                            gComm.world(),
                            &tRequest[ tCount++ ] ) );

                tDataOffset+= c;
            }

            tOffset += 3 ;
        }

        // wait until send is complete
        comm_check( MPI_Waitall( tCount, tRequest, tStatus ) );

        // tidy up memory
        free( tStatus );
        free( tRequest );
        free( tSizes );


#endif
    }

//==============================================================================
// STRINGS
//==============================================================================

    void
    broadcast( Cell< string > & aData, const proc_t aRoot=0 );

    void
    send( const string & aMessage, const proc_t aTarget=0 );

    void
    receive( string & aMessage, const proc_t aSource=0 );

//------------------------------------------------------------------------------
}
#endif //COMMTOOLS_HPP
