//
// Created by Christian Messe on 2019-01-18.
//

#ifndef BELFEM_COMMTOOLS_HPP
#define BELFEM_COMMTOOLS_HPP
#include <limits>
#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Communicator.hpp"
#include "comm_typedefs.hpp"


extern belfem::Communicator gComm;

namespace belfem
{
    const int gMaxArrayLength = 64 * 1024;

//==============================================================================
// UTILITIES
//==============================================================================

    /**
     * return the number of procs running
     */
    proc_t
    comm_size();

//------------------------------------------------------------------------------

    /**
     * return the id of this proc
     */
    proc_t
    comm_rank();

//------------------------------------------------------------------------------

    void
    comm_barrier();

//------------------------------------------------------------------------------

    // default communication list to send information to all other procs
    void
    create_master_commlist( Vector< proc_t > & aCommList );

//------------------------------------------------------------------------------

    /**
     * creates a tag for safer MPI communication
     *
     * @param[ in ] aSource : sending proc
     * @param[ in ] aTarget : receiving proc
     */
    int
    comm_tag( const proc_t aSource, const proc_t aTarget );

//------------------------------------------------------------------------------

    template< typename T >
    T
    comm_nan( const T & aSample )
    {
        return std::numeric_limits<T>::max();
    }

    template< typename T >
    void
    comm_nan( T & aSample )
    {
        aSample = std::numeric_limits<T>::max();
    }

//------------------------------------------------------------------------------

    template< typename T >
    int
    comm_length( const Vector< T > & aMessage )
    {
        long long unsigned int tLength = aMessage.length();
        long long unsigned int tMaxLength = std::numeric_limits<int>::max();

        BELFEM_ERROR( tLength < tMaxLength,
                     "Vector is too long to be sent with MPI" );

        return ( int ) tLength;
    }

//==============================================================================
// BROADCAST
//==============================================================================

    /**
     * broadcast a value over all procs
     *
     * param[ in ] aRoot : sending processor
     * aMessage[ inoit ] aMessage: Data to be broadcasted
     */
    template< typename T >
    void
    broadcast( const proc_t aRoot, T & aMessage )
    {
#ifdef BELFEM_MPI
        if( gComm.size() > 1 )
        {
            // get the data type
            comm_data_t tDataType = get_comm_datatype(( T ) 0 );

            // we don't use MPI_Bcast, it is not reliable
            MPI_Bcast(
                    &aMessage,
                    1,
                    tDataType,
                    aRoot,
                    gComm.world());
        }
#endif
    }
//==============================================================================
// SCALAR from root to others
//==============================================================================

/* send a specific scalar to other procs
    *
    *  @param[ in ]  aCommunicationList : comm ranks to receive data
    *  @param[ in ]  aDistribute        : Scalar with data to be sent
    *
    */
    template< typename T >
    void
    send( const Vector< proc_t > & aCommunicationList,
          T                      & aData )
    {
#ifdef BELFEM_MPI

        // get total number of procs
        proc_t tCommSize = gComm.size();

        // get size of communication list
        uint tListSize = aCommunicationList.length();

        if( tListSize > 0 )
        {
            // get my id
            proc_t tMyRank = gComm.rank();

            // Allocate memory for status/request vector
            MPI_Status*  tStatus  = ( MPI_Status*  ) alloca( sizeof( MPI_Status  ) * tListSize );
            MPI_Request* tRequest = ( MPI_Request* ) alloca( sizeof( MPI_Request ) * tListSize );

            // get the data type
            comm_data_t tDataType = get_comm_datatype ( ( T ) 0 );

            // loop over all procs
            for( uint p = 0; p < tListSize; ++p )
            {
                // get target
                proc_t tTarget = aCommunicationList( p );

                if( tTarget < tCommSize || tTarget != tMyRank )
                {
                    // create communication tag
                    int tCommTag = comm_tag( tMyRank, tTarget );

                    // send data
                    MPI_Isend( &aData,
                               1,
                               tDataType,
                               tTarget,
                               tCommTag,
                               gComm.world(),
                               &tRequest[ p ] );

                    // wait until send is complete
                    MPI_Wait( &tRequest[ p ], &tStatus[ p ] );
                }
            }
        }
#endif
    }

//------------------------------------------------------------------------------

    /* send a specific scalar to other procs
    *
    *  @param[ in ]  aCommunicationList : comm ranks to receive data
    *  @param[ in ]  aDistribute        : Vector with data to be sent
    *
    */
    template< typename T >
    void
    send( const Vector< proc_t > & aCommunicationList,
          Vector< T >            & aData )
    {
#ifdef BELFEM_MPI

        // get total number of procs
        proc_t tCommSize = gComm.size();

        // get size of communication list
        uint tListSize = aCommunicationList.length();

        BELFEM_ASSERT( aCommunicationList.length() == aData.length(),
                "Length of communication list and data vector do not match" );

        if( tListSize > 0 )
        {
            // get my id
            proc_t tMyRank = gComm.rank();

            // Allocate memory for status/request vector
            MPI_Status*  tStatus  = ( MPI_Status*  ) alloca( sizeof( MPI_Status  ) * tListSize );
            MPI_Request* tRequest = ( MPI_Request* ) alloca( sizeof( MPI_Request ) * tListSize );

            // get the data type
            comm_data_t tDataType = get_comm_datatype ( ( T ) 0 );

             // loop over all procs
            for( uint p = 0; p < tListSize; ++p )
            {
                // get target
                proc_t tTarget = aCommunicationList( p );

                if( tTarget < tCommSize || tTarget != tMyRank )
                {
                    // create communication tag
                    int tCommTag = comm_tag( tMyRank, tTarget );

                    // send data
                    MPI_Isend( &aData( p ),
                               1,
                               tDataType,
                               tTarget,
                               tCommTag,
                               gComm.world(),
                               &tRequest[ p ] );

                    // wait until send is complete
                    MPI_Wait( &tRequest[ p ], &tStatus[ p ] );
                }
            }
        }
#endif
    }

//------------------------------------------------------------------------------

    /* send a specific scalar to other procs
       *
       *  @param[ in ]  aCommunicationList : comm ranks to receive data
       *  @param[ in ]  aDistribute        : Vector with data to be sent
       *
       */
    template< typename T >
    void
    send( const Vector< proc_t > & aCommunicationList,
          Vector< T >            & aData,
                  T              & aMyData )
    {
#ifdef BELFEM_MPI

        // get total number of procs
        proc_t tCommSize = gComm.size();

        // get size of communication list
        uint tListSize = aCommunicationList.length();

        BELFEM_ASSERT( aCommunicationList.length() == aData.length(),
                "Length of communication list and data vector do not match" );

        if( tListSize > 0 )
        {
            // get my id
            proc_t tMyRank = gComm.rank();

            // Allocate memory for status/request vector
            MPI_Status*  tStatus  = ( MPI_Status*  ) alloca( sizeof( MPI_Status  ) * tListSize );
            MPI_Request* tRequest = ( MPI_Request* ) alloca( sizeof( MPI_Request ) * tListSize );

            // get the data type
            comm_data_t tDataType = get_comm_datatype ( ( T ) 0 );

            // Init My data
            comm_nan( aMyData );

            // loop over all procs
            for( uint p = 0; p < tListSize; ++p )
            {
                // get target
                proc_t tTarget = aCommunicationList( p );

                if( tTarget == tMyRank )
                {
                    aMyData = aData( p );
                }
                else if( tTarget < tCommSize )
                {
                    // create communication tag
                    int tCommTag = comm_tag( tMyRank, tTarget );

                    // send data
                    MPI_Isend( &aData( p ),
                               1,
                               tDataType,
                               tTarget,
                               tCommTag,
                               gComm.world(),
                               &tRequest[ p ] );

                    // wait until send is complete
                    MPI_Wait( &tRequest[ p ], &tStatus[ p ] );
                }
            }
        }
#endif
    }

//------------------------------------------------------------------------------

    template< typename T >
    void
    receive(  const proc_t aSource, T & aData )
    {
#ifdef BELFEM_MPI
        // get my id
        proc_t tMyRank = gComm.rank();

        if( aSource < comm_size() && tMyRank != aSource )
        {
            // status/request handlers
            MPI_Status  tStatus;
            MPI_Request tRequest;

            // get the data type
            comm_data_t tDataType = get_comm_datatype ( ( T ) 0 );

            // create communication tag
            int tCommTag = comm_tag( aSource, tMyRank );

            MPI_Irecv( &aData,
                       1,
                       tDataType,
                       aSource,
                       tCommTag,
                       gComm.world(),
                       & tRequest );

            // wait until receive is complete
             MPI_Wait( &tRequest, &tStatus );
        }
#endif
    }


//==============================================================================
// SCALAR from others to root
//==============================================================================

    /**
     * send a scalar to the root
     * @param aTarget      : Rank of target
     * @param aMyData
     */
    template< typename T >
    void
    send( const proc_t aTarget, T & aData )
    {
#ifdef BELFEM_MPI

        // get my id
        proc_t tMyRank = gComm.rank();

        if( aTarget < comm_size() && tMyRank != aTarget )
        {
            //  status/request handlers
            MPI_Request tRequest;
            MPI_Status  tStatus;

            // create communication tag
            int tCommTag = comm_tag( tMyRank, aTarget );

            // get the data type
            comm_data_t tDataType = get_comm_datatype ( ( T ) 0 );

            // receive the data
            MPI_Isend( &aData,
                       1,
                       tDataType,
                       aTarget,
                       tCommTag,
                       gComm.world(),
                       &tRequest );


            // wait until receive is complete
            MPI_Wait( &tRequest, &tStatus );
        }
#endif
    }

//------------------------------------------------------------------------------

    /**
     * receive a sclar from other procs
     *  @param[ in ]  aCommunicationList : romm ranks to receive data
     *  @param[ in ]  aDistribute        : Vector with data to be sent
     *  @param[ in ]  aMyData          : Own data, if exist
     */
    template< typename T >
    void
    receive( const Vector< proc_t > & aCommunicationList,
             Vector< T >            & aData,
             T                        aMyData=std::numeric_limits<T>::quiet_NaN() )
    {
#ifdef BELFEM_MPI
        // get total number of procs
        proc_t tCommSize = gComm.size();

        if ( tCommSize > 1 )
        {
            // get size of communication list
            uint tListSize = aCommunicationList.length();

            // allocate the data container
            aData.set_size( tListSize );

            // get the data type
            comm_data_t tDataType = get_comm_datatype(( T ) 0 );

            // get my id
            proc_t tMyRank = gComm.rank();

            // Allocate memory for status/request vector
            MPI_Request * tRequest = ( MPI_Request * ) alloca( sizeof( MPI_Request ) * tListSize );
            MPI_Status  * tStatus = ( MPI_Status * ) alloca( sizeof( MPI_Status ) * tListSize );

            // loop over all procs
            for ( uint p = 0; p < tListSize; ++p )
            {
                // get source
                proc_t tSource = aCommunicationList( p );

                if ( tSource == tMyRank )
                {
                    aData( p ) = aMyData;
                }
                else if ( tSource < tCommSize )
                {
                    // create communication tag
                    int tCommTag = comm_tag( tSource, tMyRank );

                    // receive the data
                    MPI_Irecv( &aData( p ),
                               1,
                               tDataType,
                               tSource,
                               tCommTag,
                               gComm.world(),
                               &tRequest[ p ] );

                    // wait until receive is complete
                    MPI_Wait( &tRequest[ p ], &tStatus[ p ] );
                }
                else
                {
                    comm_nan( aData( p ));
                }
            }
        }
#endif
    }

//==============================================================================
// VECTOR from root to others
//==============================================================================

    Vector< int >
    split_message( const index_t aMessageLength );

//------------------------------------------------------------------------------

    template< typename T >
    void
    send( const Vector< proc_t > & aCommunicationList,
          Cell< Vector< T >    > & aData )
    {
#ifdef BELFEM_MPI

        // get total number of procs
        proc_t tCommSize = gComm.size();

        // get size of communication list
        uint tListSize = aCommunicationList.length();

        BELFEM_ASSERT( aCommunicationList.length() == aData.size(),
                "Length of communication list and data cell do not match ( %u and %u )",
                ( unsigned int ) aCommunicationList.length(),
                ( unsigned int ) aData.size() );

		if( tListSize > 0 )
        {
            // get my id
            proc_t tMyRank = gComm.rank();

            // Allocate memory for status/request vector
            MPI_Status*  tStatus  = ( MPI_Status*  ) alloca( sizeof( MPI_Status  ) * tListSize );
            MPI_Request* tRequest = ( MPI_Request* ) alloca( sizeof( MPI_Request ) * tListSize );
        
        	// get the data types
            comm_data_t tLengthType = get_comm_datatype ( ( index_t ) 0 );
            comm_data_t tDataType = get_comm_datatype ( ( T ) 0 );
                        
            // loop over all procs
            for( uint p = 0; p < tListSize; ++p )
            {
            	// get proc ID of target
            	proc_t tTarget = aCommunicationList( p );
            	
            	// if target exists and is not myself
            	if( tTarget < tCommSize && tTarget != tMyRank )
                {
                	// create communication tag
                    int tCommTag = comm_tag( tMyRank, tTarget );
                    
                	// get length of vector
                	index_t tLength = aData( p ).length();

                	// send length to target
                	MPI_Isend( &tLength,
                               1,
                               tLengthType,
                               tTarget,
                               tCommTag,
                               gComm.world(),
                               &tRequest[ p ] );
                               
                    // wait until send is complete
                    MPI_Wait( &tRequest[ p ], &tStatus[ p ] );
                    
                    if( tLength > 0 )
                    {
                    	// get vector
                    	Vector< T > & tData = aData( p );
                    	
                    	// calculate number of individual messages to be sent
                    	Vector< int > tLengths = split_message( tData.length() );
                    	index_t tNumMessages = tLengths.length();
                    	
                    	// Allocate status and request containers
                    	MPI_Status*  tStatus2  = ( MPI_Status*  ) alloca( tNumMessages * sizeof( MPI_Status  ) );
            			MPI_Request* tRequest2 = ( MPI_Request* ) alloca( tNumMessages * sizeof( MPI_Request ) );
                    
                    	// offset in array
                    	index_t tOffset = 0;
                    	
                    	// loop over all messages
                    	for( index_t m=0; m<tNumMessages; ++m )
                    	{
                    		// increment comm tag
                    		tCommTag++;

                    		// get raw pointer of vector
                    		T * tData = aData( p ).data();

                    		// send data
                        	MPI_Isend( &tData[ tOffset ],
                               tLengths( m ),
                               tDataType,
                               tTarget,
                               tCommTag,
                               gComm.world(),
                               &tRequest2[ m ] );
                               
                            // increment offset
                            tOffset += tLengths( m );

 							// wait until send is complete
                        	MPI_Wait( &tRequest2[ m ], &tStatus2[ m ] );
                    	}
                    }
                }
            }    
        }
#endif
    }

//------------------------------------------------------------------------------

    // send the same vector to all procs on communication lust
    template< typename T >
    void
    send_same(
            const Vector< proc_t > & aCommunicationList,
            const Vector< T >   & aData )
    {
#ifdef BELFEM_MPI

        // get total number of procs
        proc_t tCommSize = gComm.size();

        // get size of communication list
        uint tListSize = aCommunicationList.length();

		if( tListSize > 0 )
        {
            // get my id
            proc_t tMyRank = gComm.rank();

            // Allocate memory for status/request vector
            MPI_Status*  tStatus  = ( MPI_Status*  ) alloca( sizeof( MPI_Status  ) * tListSize );
            MPI_Request* tRequest = ( MPI_Request* ) alloca( sizeof( MPI_Request ) * tListSize );

        	// get the data types
            comm_data_t tLengthType = get_comm_datatype ( ( index_t ) 0 );
            comm_data_t tDataType = get_comm_datatype ( ( T ) 0 );

            // loop over all procs
            for( uint p = 0; p < tListSize; ++p )
            {
            	// get proc ID of target
            	proc_t tTarget = aCommunicationList( p );

            	// if target exists and is not myself
            	if( tTarget < tCommSize && tTarget != tMyRank )
                {
                	// create communication tag
                    int tCommTag = comm_tag( tMyRank, tTarget );

                	// get length of vector
                	index_t tLength = aData.length();

                	// send length to target
                	MPI_Isend( &tLength,
                               1,
                               tLengthType,
                               tTarget,
                               tCommTag,
                               gComm.world(),
                               &tRequest[ p ] );

                    // wait until send is complete
                    MPI_Wait( &tRequest[ p ], &tStatus[ p ] );

                    if( tLength > 0 )
                    {
                    	// calculate number of individual messages to be sent
                    	Vector< int > tLengths = split_message( aData.length() );
                    	index_t tNumMessages = tLengths.length();

                    	// Allocate status and request containers
                    	MPI_Status*  tStatus2  = ( MPI_Status*  ) alloca( tNumMessages * sizeof( MPI_Status  ) );
            			MPI_Request* tRequest2 = ( MPI_Request* ) alloca( tNumMessages * sizeof( MPI_Request ) );

                    	// offset in array
                    	index_t tOffset = 0;

                    	// loop over all messages
                    	for( index_t m=0; m<tNumMessages; ++m )
                    	{
                    		// increment comm tag
                    		tCommTag++;

                    		// get raw pointer of vector
                    		const T * tData = aData.data();

                    		// send data
                        	MPI_Isend( &tData[ tOffset ],
                               tLengths( m ),
                               tDataType,
                               tTarget,
                               tCommTag,
                               gComm.world(),
                               &tRequest2[ m ] );

                            // increment offset
                            tOffset += tLengths( m );

 							// wait until send is complete
                        	MPI_Wait( &tRequest2[ m ], &tStatus2[ m ] );
                    	}
                    }
                }
            }
        }
#endif
    }
//------------------------------------------------------------------------------

    // send the same matrix to all procs on communication lust
    template< typename T >
    void
    send_same(
            const Vector< proc_t > & aCommunicationList,
            const Matrix< T >      & aData )
    {
#ifdef BELFEM_MPI

        // get total number of procs
        proc_t tCommSize = gComm.size();

        // get size of communication list
        uint tListSize = aCommunicationList.length();

		if( tListSize > 0 )
        {
            // get my id
            proc_t tMyRank = gComm.rank();

            // Allocate memory for status/request vector
            MPI_Status*  tStatus  = ( MPI_Status*  ) alloca( sizeof( MPI_Status  ) * tListSize );
            MPI_Request* tRequest = ( MPI_Request* ) alloca( sizeof( MPI_Request ) * tListSize );

        	// get the data types
            comm_data_t tLengthType = get_comm_datatype ( ( index_t ) 0 );
            comm_data_t tDataType = get_comm_datatype ( ( T ) 0 );

            // row, columns and length of array
            index_t tLength[ 3 ];
            tLength[ 0 ] = aData.n_rows();
            tLength[ 1 ] = aData.n_cols();
            tLength[ 2 ] = aData.capacity();

            // loop over all procs
            for( uint p = 0; p < tListSize; ++p )
            {
            	// get proc ID of target
            	proc_t tTarget = aCommunicationList( p );

            	// if target exists and is not myself
            	if( tTarget < tCommSize && tTarget != tMyRank )
                {
                	// create communication tag
                    int tCommTag = comm_tag( tMyRank, tTarget );

                	// send length to target
                	MPI_Isend( &tLength,
                               3,
                               tLengthType,
                               tTarget,
                               tCommTag,
                               gComm.world(),
                               &tRequest[ p ] );

                    // wait until send is complete
                    MPI_Wait( &tRequest[ p ], &tStatus[ p ] );

                    if( tLength[ 0 ] > 0 && tLength[ 1 ] > 0 )
                    {
                    	// calculate number of individual messages to be sent
                    	Vector< int > tLengths = split_message( tLength[ 2 ] );
                    	index_t tNumMessages = tLengths.length();

                    	// Allocate status and request containers
                    	MPI_Status*  tStatus2  = ( MPI_Status*  ) alloca( tNumMessages * sizeof( MPI_Status  ) );
            			MPI_Request* tRequest2 = ( MPI_Request* ) alloca( tNumMessages * sizeof( MPI_Request ) );

                    	// offset in array
                    	index_t tOffset = 0;

                    	// loop over all messages
                    	for( index_t m=0; m<tNumMessages; ++m )
                    	{
                    		// increment comm tag
                    		tCommTag++;

                    		// get raw pointer of vector
                    		const T * tData = aData.data();

                    		// send data
                        	MPI_Isend( &tData[ tOffset ],
                               tLengths( m ),
                               tDataType,
                               tTarget,
                               tCommTag,
                               gComm.world(),
                               &tRequest2[ m ] );

                            // increment offset
                            tOffset += tLengths( m );

 							// wait until send is complete
                        	MPI_Wait( &tRequest2[ m ], &tStatus2[ m ] );
                    	}
                    }
                }
            }
        }
#endif
    }


//------------------------------------------------------------------------------

    template< typename T >
    void
    receive( const Vector< proc_t > & aCommunicationList,
          Cell< Vector< T >    > & aData )
    {
#ifdef BELFEM_MPI

        // get total number of procs
        proc_t tCommSize = gComm.size();

        // get size of communication list
        uint tListSize = aCommunicationList.length();

		if( tListSize > 0 )
        {
		    bool tResetMyData = false;

		    if( ( proc_t ) aData.size() != tCommSize )
		    {
		        // populate data vector
		        Vector< T > tEmpty;
		        aData.set_size( tListSize, tEmpty );

		        tResetMyData = true;
		    }

            // get my id
            proc_t tMyRank = gComm.rank();

            // Allocate memory for status/request vector
            MPI_Status*  tStatus  = ( MPI_Status*  ) alloca( sizeof( MPI_Status  ) * tListSize );
            MPI_Request* tRequest = ( MPI_Request* ) alloca( sizeof( MPI_Request ) * tListSize );
        
        	// get the data types
            comm_data_t tLengthType = get_comm_datatype ( ( index_t ) 0 );
            comm_data_t tDataType = get_comm_datatype ( ( T ) 0 );
                        
            // loop over all procs
            for( uint p = 0; p < tListSize; ++p )
            {
            	// get proc ID of target
            	proc_t tSource = aCommunicationList( p );
            	
            	// if target exists and is not myself
            	if( tSource < tCommSize && tSource != tMyRank )
                {
                	// create communication tag
                    int tCommTag = comm_tag( tSource, tMyRank );
                    
                	// get length of vector
                	index_t tLength = 0;
                	
                	// send length to target
                	MPI_Irecv( &tLength,
                               1,
                               tLengthType,
                               tSource,
                               tCommTag,
                               gComm.world(),
                               &tRequest[ p ] );

                    // wait until send is complete
                    MPI_Wait( &tRequest[ p ], &tStatus[ p ] );
                    
                    if( tLength > 0 )
                    {
                    	// get vector
                    	Vector< T > & tData = aData( p );
                    	
                    	// set size of vector
                    	tData.set_size( tLength );
                    	
                    	// calculate number of individual messages to be sent
                    	Vector< int > tLengths = split_message( tLength );
                    	index_t tNumMessages = tLengths.length();
                    	
                    	// Allocate status and request containers
                    	MPI_Status*  tStatus2  = ( MPI_Status*  ) alloca( tNumMessages * sizeof( MPI_Status  ) );
            			MPI_Request* tRequest2 = ( MPI_Request* ) alloca( tNumMessages * sizeof( MPI_Request ) );
                    
                    	// offset in array
                    	index_t tOffset = 0;

                    	// loop over all messages
                    	for( index_t m=0; m<tNumMessages; ++m )
                    	{
                    		// increment comm tag
                    		tCommTag++;

                    		// get raw pointer of data
                    	    T * tData = aData( p ).data();

                    		// send data
                        	MPI_Irecv( &tData[ tOffset ],
                               tLengths( m ),
                               tDataType,
                               tSource,
                               tCommTag,
                               gComm.world(),
                               &tRequest2[ m ] );
                               
                            // increment offset
                            tOffset += tLengths( m );

 							// wait until send is complete
                        	MPI_Wait( &tRequest2[ m ], &tStatus2[ m ] );
                    	}
                    }
                    else
                    {
                    	aData( p ) = {};
                    }
                }
                else if( tSource >= tCommSize || tResetMyData )
                {
                    	aData( p ) = {};
                }
            }    
        }
#endif
    }

//------------------------------------------------------------------------------

    template< typename T >
    void
    receive( const Vector< proc_t > & aCommunicationList,
             Cell< Vector< T >    > & aData,
             const Vector< T >      & aMyData )
    {
        receive( aCommunicationList, aData );

        proc_t tMyRank = gComm.rank();

        for( uint p=0; p<aCommunicationList.length(); ++p )
        {
            if( aCommunicationList( p ) == tMyRank )
            {
                aData( p ).vector_data() = aMyData.vector_data();
                break;
            }
        }
    }

//==============================================================================
// VECTOR from others to root
//==============================================================================

    template< typename T >
    void
    send( const proc_t aTarget,
          Vector< T >  & aData )
    {
#ifdef BELFEM_MPI

		// get total number of procs
        proc_t tCommSize = gComm.size();
        
        // get my rank
        proc_t tMyRank   = gComm.rank();
        
        if( aTarget < tCommSize && aTarget != tMyRank )
        {
        	// Allocate memory for status/request vector
            MPI_Status  tStatus;
            MPI_Request tRequest;
            
            // get the data types
            comm_data_t tLengthType = get_comm_datatype ( ( index_t ) 0 );
            comm_data_t tDataType = get_comm_datatype ( ( T ) 0 );
            
        	// get length of vector
            index_t tLength = aData.length();
                	
            // create communication tag
            int tCommTag = comm_tag( tMyRank, aTarget );
                    
            // send length to target
            MPI_Isend( &tLength,
            	1,
                tLengthType,
                aTarget,
                tCommTag,
                gComm.world(),
                &tRequest );
                               
            // wait until send is complete
            MPI_Wait( &tRequest, &tStatus );
            
            if( tLength > 0 )
            {
                // calculate number of individual messages to be sent
                Vector< int > tLengths = split_message( aData.length() );
                index_t tNumMessages = tLengths.length();
                    	
                // Allocate status and request containers
                MPI_Status*  tStatus2  = ( MPI_Status*  ) alloca( tNumMessages * sizeof( MPI_Status  ) );
            	MPI_Request* tRequest2 = ( MPI_Request* ) alloca( tNumMessages * sizeof( MPI_Request ) );
                    
                // offset in array
                index_t tOffset = 0;

                // loop over all messages
                for( index_t m=0; m<tNumMessages; ++m )
                {
                    // increment comm tag
                    tCommTag++;

                    // get raw pointer of data
                    T * tData =  aData.data();

                    // send data
                    MPI_Isend( &tData[ tOffset ],
                        tLengths( m ),
                        tDataType,
                        aTarget,
                        tCommTag,
                        gComm.world(),
                        &tRequest2[ m ] );

                    // increment offset
                    tOffset += tLengths( m );

 					// wait until send is complete
                    MPI_Wait( &tRequest2[ m ], &tStatus2[ m ] );

                }
            }
        }
#endif
    }
 
//------------------------------------------------------------------------------

    template< typename T >
    void
    receive( const proc_t aSource,
             Vector< T >  & aData )
    {
#ifdef BELFEM_MPI

		// get total number of procs
        proc_t tCommSize = gComm.size();
        
        // get my rank
        proc_t tMyRank   = gComm.rank();
        
        if( aSource < tCommSize && aSource != tMyRank )
        {
        	// Allocate memory for status/request vector
            MPI_Status  tStatus;
            MPI_Request tRequest;
            
            // get the data types
            comm_data_t tLengthType = get_comm_datatype ( ( index_t ) 0 );
            comm_data_t tDataType = get_comm_datatype ( ( T ) 0 );
            
        	// get length of vector
            index_t tLength = 0;
                	
            // create communication tag
            int tCommTag = comm_tag( aSource, tMyRank );
                    
            // receive length from target
            MPI_Irecv( &tLength,
            	1,
                tLengthType,
                aSource,
                tCommTag,
                gComm.world(),
                &tRequest );

            // wait until send is complete
            MPI_Wait( &tRequest, &tStatus );

            if( tLength > 0 )
            {
            	// set size of data
            	aData.set_size( tLength );
            		
                // calculate number of individual messages to be sent
                Vector< int > tLengths = split_message( tLength );
                index_t tNumMessages = tLengths.length();
                    	
                // Allocate status and request containers
                MPI_Status*  tStatus2  = ( MPI_Status*  ) alloca( tNumMessages * sizeof( MPI_Status  ) );
            	MPI_Request* tRequest2 = ( MPI_Request* ) alloca( tNumMessages * sizeof( MPI_Request ) );
                    
                // offset in array
                index_t tOffset = 0;

                // loop over all messages
                for( index_t m=0; m<tNumMessages; ++m )
                {
                    // increment comm tag
                    tCommTag++;

                    // get raw pointer of data
                    T * tData =  aData.data();

                    // receive data
                    MPI_Irecv( &tData[ tOffset ],
                        tLengths( m ),
                        tDataType,
                        aSource,
                        tCommTag,
                        gComm.world(),
                        &tRequest2[ m ] );
                               
                    // increment offset
                    tOffset += tLengths( m );

 					// wait until receive is complete
                    MPI_Wait( &tRequest2[ m ], &tStatus2[ m ] );

                }
            }
            else
            {
            	aData = {};
            }
        }
        else
        {
        	aData = {};
        }
#endif
    }

//------------------------------------------------------------------------------

    /**
     *  send raw data to another proc ( in this usecase, the proc receives
     *  a vector
     * @tparam T
     * @param aTarget
     * @param aData
     */
    template< typename T >
    void
    send(
            const proc_t   aTarget,
            const index_t  aNumberOfSamples,
            T             * aData )
    {
#ifdef BELFEM_MPI

        // get total number of procs
        proc_t tCommSize = gComm.size();

        // get my rank
        proc_t tMyRank   = gComm.rank();

        if( aTarget < tCommSize && aTarget != tMyRank )
        {
        	// Allocate memory for status/request vector
            MPI_Status  tStatus;
            MPI_Request tRequest;

            // get the data types
            comm_data_t tLengthType = get_comm_datatype ( ( index_t ) 0 );
            comm_data_t tDataType = get_comm_datatype ( ( T ) 0 );

            // create communication tag
            int tCommTag = comm_tag( tMyRank, aTarget );

            index_t tLength = aNumberOfSamples ;

            // send length to target
            MPI_Isend( &tLength,
            	1,
                tLengthType,
                aTarget,
                tCommTag,
                gComm.world(),
                &tRequest );

            // wait until send is complete
            MPI_Wait( &tRequest, &tStatus );

            if( aNumberOfSamples > 0 )
            {
                // calculate number of individual messages to be sent
                Vector< int > tLengths = split_message( tLength );
                index_t tNumMessages = tLengths.length();

                // Allocate status and request containers
                MPI_Status*  tStatus2  = ( MPI_Status*  ) alloca( tNumMessages * sizeof( MPI_Status  ) );
            	MPI_Request* tRequest2 = ( MPI_Request* ) alloca( tNumMessages * sizeof( MPI_Request ) );

                // offset in array
                index_t tOffset = 0;

                // loop over all messages
                for( index_t m=0; m<tNumMessages; ++m )
                {
                    // increment comm tag
                    tCommTag++;

                    // send data
                    MPI_Isend( &aData[ tOffset ],
                        tLengths( m ),
                        tDataType,
                        aTarget,
                        tCommTag,
                        gComm.world(),
                        &tRequest2[ m ] );

                    // increment offset
                    tOffset += tLengths( m );

 					// wait until send is complete
                    MPI_Wait( &tRequest2[ m ], &tStatus2[ m ] );

                }
            }
        }
#endif
    }

// -------------------------------------------------------------------------

    template< typename T >
    index_t
    receive( const proc_t aSource,
             T             * aData )
    {
#ifdef BELFEM_MPI

        // get total number of procs
        proc_t tCommSize = gComm.size();

        // get my rank
        proc_t tMyRank   = gComm.rank();

        if( aSource < tCommSize && aSource != tMyRank )
        {
            // Allocate memory for status/request vector
            MPI_Status  tStatus;
            MPI_Request tRequest;

            // get the data types
            comm_data_t tLengthType = get_comm_datatype ( ( index_t ) 0 );
            comm_data_t tDataType = get_comm_datatype ( ( T ) 0 );

            // get length of vector
            index_t aNumberOfSamples = 0;

            // create communication tag
            int tCommTag = comm_tag( aSource, tMyRank );

            // receive length from target
            MPI_Irecv( & aNumberOfSamples,
                       1,
                       tLengthType,
                       aSource,
                       tCommTag,
                       gComm.world(),
                       &tRequest );

            // wait until send is complete
            MPI_Wait( &tRequest, &tStatus );

            if( aNumberOfSamples > 0 )
            {

                // calculate number of individual messages to be sent
                Vector< int > tLengths = split_message( aNumberOfSamples );
                index_t tNumMessages = tLengths.length();

                // Allocate status and request containers
                MPI_Status*  tStatus2  = ( MPI_Status*  ) alloca( tNumMessages * sizeof( MPI_Status  ) );
                MPI_Request* tRequest2 = ( MPI_Request* ) alloca( tNumMessages * sizeof( MPI_Request ) );

                // offset in array
                index_t tOffset = 0;

                // loop over all messages
                for( index_t m=0; m<tNumMessages; ++m )
                {
                    // increment comm tag
                    tCommTag++;

                    // receive data
                    MPI_Irecv( aData + tOffset,
                               tLengths( m ),
                               tDataType,
                               aSource,
                               tCommTag,
                               gComm.world(),
                               &tRequest2[ m ] );

                    // increment offset
                    tOffset += tLengths( m );

                    // wait until receive is complete
                    MPI_Wait( &tRequest2[ m ], &tStatus2[ m ] );

                }
            }

            return aNumberOfSamples ;
        }
#endif
        return 0 ;
    }

//==============================================================================
// Matrix from root to others
//==============================================================================

    template< typename T >
    void
    send( const proc_t aTarget,
          Matrix< T >  & aData )
    {
#ifdef BELFEM_MPI

        // get total number of procs
        proc_t tCommSize = gComm.size();

        // get my rank
        proc_t tMyRank   = gComm.rank();

        if( aTarget < tCommSize && aTarget != tMyRank )
        {
        	// Allocate memory for status/request vector
            MPI_Status  tStatus;
            MPI_Request tRequest;

            // get the data types
            comm_data_t tLengthType = get_comm_datatype ( ( index_t ) 0 );
            comm_data_t tDataType = get_comm_datatype ( ( T ) 0 );

        	// row, columns and length of array
            index_t tLength[ 3 ];
            tLength[ 0 ] = aData.n_rows();
            tLength[ 1 ] = aData.n_cols();
            tLength[ 2 ] = aData.capacity();

            // create communication tag
            int tCommTag = comm_tag( tMyRank, aTarget );

            // send length to target
            MPI_Isend( &tLength,
            	3,
                tLengthType,
                aTarget,
                tCommTag,
                gComm.world(),
                &tRequest );

            // wait until send is complete
            MPI_Wait( &tRequest, &tStatus );

            if( tLength[ 0 ] > 0 && tLength[ 1 ] > 0 )
            {
                // calculate number of individual messages to be sent
                Vector< int > tLengths = split_message( tLength[ 2 ] );
                index_t tNumMessages = tLengths.length();

                // Allocate status and request containers
                MPI_Status*  tStatus2  = ( MPI_Status*  ) alloca( tNumMessages * sizeof( MPI_Status  ) );
            	MPI_Request* tRequest2 = ( MPI_Request* ) alloca( tNumMessages * sizeof( MPI_Request ) );

                // offset in array
                index_t tOffset = 0;

                // loop over all messages
                for( index_t m=0; m<tNumMessages; ++m )
                {
                    // increment comm tag
                    tCommTag++;

                    // get raw pointer of matrix
                    T * tData =  aData.data();

                    // send data
                    MPI_Isend( &tData[ tOffset ],
                        tLengths( m ),
                        tDataType,
                        aTarget,
                        tCommTag,
                        gComm.world(),
                        &tRequest2[ m ] );

                    // increment offset
                    tOffset += tLengths( m );

 					// wait until send is complete
                    MPI_Wait( &tRequest2[ m ], &tStatus2[ m ] );

                }
            }
        }
#endif
    }

//------------------------------------------------------------------------------

    template< typename T >
    void
    receive( const proc_t aSource,
             Matrix< T >  & aData )
    {
#ifdef BELFEM_MPI

        // get total number of procs
        proc_t tCommSize = gComm.size();

        // get my rank
        proc_t tMyRank   = gComm.rank();

        if( aSource < tCommSize && aSource != tMyRank )
        {
        	// Allocate memory for status/request vector
            MPI_Status  tStatus;
            MPI_Request tRequest;

            // get the data types
            comm_data_t tLengthType = get_comm_datatype ( ( index_t ) 0 );
            comm_data_t tDataType = get_comm_datatype ( ( T ) 0 );

        	// get length of matrix
            index_t tLength[ 3 ] = { 0, 0, 0 };

            // create communication tag
            int tCommTag = comm_tag( aSource, tMyRank );

            // send length to target
            MPI_Irecv( &tLength,
            	3,
                tLengthType,
                aSource,
                tCommTag,
                gComm.world(),
                &tRequest );

            // wait until send is complete
            MPI_Wait( &tRequest, &tStatus );

            if ( tLength[ 0 ] > 0 && tLength[ 1 ] > 0 )
            {
            	// set size of data
            	aData.set_size( tLength[ 0 ], tLength[ 1 ] );

                // calculate number of individual messages to be sent
                Vector< int > tLengths = split_message( tLength[ 2 ] );
                index_t tNumMessages = tLengths.length();

                // Allocate status and request containers
                MPI_Status*  tStatus2  = ( MPI_Status*  ) alloca( tNumMessages * sizeof( MPI_Status  ) );
            	MPI_Request* tRequest2 = ( MPI_Request* ) alloca( tNumMessages * sizeof( MPI_Request ) );

                // offset in array
                index_t tOffset = 0;

                // loop over all messages
                for( index_t m=0; m<tNumMessages; ++m )
                {
                    // increment comm tag
                    tCommTag++;

                    // get raw pointer of matrix
                    T * tData =  aData.data();

                    // receive data
                    MPI_Irecv( &tData[ tOffset ],
                        tLengths( m ),
                        tDataType,
                        aSource,
                        tCommTag,
                        gComm.world(),
                        &tRequest2[ m ] );

                    // increment offset
                    tOffset += tLengths( m );

 					// wait until send is complete
                    MPI_Wait( &tRequest2[ m ], &tStatus2[ m ] );

                }
            }
            else
            {
            	aData = {{}};
            }
        }
        else
        {
        	aData = {{}};
        }
#endif
    }
//------------------------------------------------------------------------------

    template< typename T >
    void
    receive( const Vector< proc_t > & aCommunicationList,
             Cell< Matrix< T >    > & aData )
    {
#ifdef BELFEM_MPI

        // get total number of procs
        proc_t tCommSize = gComm.size();

        // get size of communication list
        uint tListSize = aCommunicationList.length();

		if( tListSize > 0 )
        {
		    bool tResetMyData = false;

		    if( ( proc_t ) aData.size() != tCommSize )
		    {
		        // populate data vector
		        Matrix< T > tEmpty;
		        aData.set_size( tListSize, tEmpty );

		        tResetMyData = true;
		    }

            // get my id
            proc_t tMyRank = gComm.rank();

            // Allocate memory for status/request vector
            MPI_Status*  tStatus  = ( MPI_Status*  ) alloca( sizeof( MPI_Status  ) * tListSize );
            MPI_Request* tRequest = ( MPI_Request* ) alloca( sizeof( MPI_Request ) * tListSize );

        	// get the data types
            comm_data_t tLengthType = get_comm_datatype ( ( index_t ) 0 );
            comm_data_t tDataType = get_comm_datatype ( ( T ) 0 );

            // loop over all procs
            for( uint p = 0; p < tListSize; ++p )
            {
            	// get proc ID of target
            	proc_t tSource = aCommunicationList( p );

            	// if target exists and is not myself
            	if( tSource < tCommSize && tSource != tMyRank )
                {
                	// create communication tag
                    int tCommTag = comm_tag( tSource, tMyRank );

                	// get length of vector
                	index_t tLength[ 3 ] = { 0, 0, 0 };

                	// send length to target
                	MPI_Irecv( &tLength,
                               3,
                               tLengthType,
                               tSource,
                               tCommTag,
                               gComm.world(),
                               &tRequest[ p ] );

                    // wait until send is complete
                    MPI_Wait( &tRequest[ p ], &tStatus[ p ] );

                    if( tLength[ 0 ] > 0 && tLength[ 1 ] > 0 )
                    {
                    	// get vector
                    	Matrix< T > & tData = aData( p );

                    	// set size of vector
                    	tData.set_size( tLength[ 0 ], tLength[ 1 ] );

                    	// calculate number of individual messages to be sent
                    	Vector< int > tLengths = split_message( tLength[ 2 ] );
                    	index_t tNumMessages = tLengths.length();

                    	// Allocate status and request containers
                    	MPI_Status*  tStatus2  = ( MPI_Status*  ) alloca( tNumMessages * sizeof( MPI_Status  ) );
            			MPI_Request* tRequest2 = ( MPI_Request* ) alloca( tNumMessages * sizeof( MPI_Request ) );

                    	// offset in array
                    	index_t tOffset = 0;

                    	// loop over all messages
                    	for( index_t m=0; m<tNumMessages; ++m )
                    	{
                    		// increment comm tag
                    		tCommTag++;

                    		// get raw pointer of matrix
                    		T * tData = aData( p ).data();

                    		// send data
                        	MPI_Irecv( &tData[ tOffset ],
                               tLengths( m ),
                               tDataType,
                               tSource,
                               tCommTag,
                               gComm.world(),
                               &tRequest2[ m ] );

                            // increment offset
                            tOffset += tLengths( m );

 							// wait until send is complete
                        	MPI_Wait( &tRequest2[ m ], &tStatus2[ m ] );
                    	}
                    }
                    else
                    {
                    	aData( p ) = {{}};
                    }
                }
                else if( tSource >= tCommSize || tResetMyData )
                {
                    	aData( p ) = {{}};
                }
            }
        }
#endif
    }

    /**
     * this is a simple test routine to check the MPI functionality.
     * @param aMessage
     */
    inline uint
    pingpong( const uint aMessage )
    {
        BELFEM_ASSERT( comm_size() == 2, "pingpong can only be played by two players" );
        comm_barrier() ;

        uint aResult = 0 ;
        if( comm_rank() == 0 )
        {
            std::cout << "Ping " << aMessage << std::endl ;
            send( 1, aMessage );
            comm_barrier() ;
            aResult = aMessage ;
        }
        else
        {
            receive( 0, aResult );
            comm_barrier() ;
            std::cout << "Pong " << aResult << std::endl << std::endl ;
        }
        comm_barrier() ;
        return aResult ;

    }

} /* namespace belfem */
#endif //BELFEM_COMMTOOLS_HPP
