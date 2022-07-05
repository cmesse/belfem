//
// Created by Christian Messe on 13.07.20.
//

#include "cl_SolverPetscDistributor.hpp"
#include "commtools.hpp"
#include "fn_sum.hpp"

namespace belfem
{
    namespace solver
    {
//------------------------------------------------------------------------------

        PetscDistributor::PetscDistributor( PetscData & aData, SpMatrix & aMatrix ) :
                mMyRank( comm_rank() ),
                mCommSize( comm_size() ),
#ifdef BELFEM_PETSC
                mData( aData ),
                mMatrix( aMatrix )
#else
                mData( aData )
#endif
        {
#ifdef BELFEM_PETSC
            // split the indices of the matrix among the other procs
            this->create_communication_list() ;

            // make sure that matrix is cpp based and has both indices
            // also enforces cpp indexing (zero-based)
            mMatrix.create_coo_indices();

            this->split_ownership();

            // compute how the memory is distributed
            this->compute_memory() ;

            // create the matrix on the data stricture
            this->initialize_matrix() ;

            this->initialize_vectors() ;

#endif
        }

//------------------------------------------------------------------------------

        PetscDistributor::~PetscDistributor()
        {
#ifdef BELFEM_PETSC
            /*if( mMyRank != 0 )
            {
                // deallocate matrix data container
                free( mMatrixData );
            }*/
#endif
        }

//------------------------------------------------------------------------------

        void
        PetscDistributor::synchronize(
                Vector <PetscReal> & aLHS,
                Vector <PetscReal> & aRHS,
                bool aHaveLHS )
        {
#ifdef BELFEM_PETSC

            this->update_matrix();

            if( aHaveLHS )
            {
                this->distribute_vector( aLHS, mLHS ) ;

                // write values into petsc
                petsctools_set_vector( mLHS , mData.mVectorIndices, mData.mLHS );

                comm_barrier() ;
            }
            else
            {
                mLHS.fill( 0.0 );
            }
            this->distribute_vector( aRHS, mRHS );
            petsctools_set_vector( mRHS , mData.mVectorIndices, mData.mRHS );

            // wait for other procs
            comm_barrier();
#endif
        }


//------------------------------------------------------------------------------

        void
        PetscDistributor::create_communication_list()
        {
#ifdef BELFEM_PETSC
            if( mMyRank == 0 )
            {
                mCommunicationList.set_size( mCommSize ) ;

                for( proc_t k=0; k<mCommSize; ++k )
                {
                    mCommunicationList( k ) = k ;
                }
            }
#endif
        }
//------------------------------------------------------------------------------

        void
        PetscDistributor::split_ownership()
        {
#ifdef BELFEM_PETSC
            mData.mMyNumRows = PETSC_DECIDE ;
            mData.mMyNumCols = PETSC_DECIDE ;

            // communicate size of matrix
            if( mMyRank == 0 )
            {
                mData.mNumRows = mMatrix.n_rows() ;
                mData.mNumCols = mMatrix.n_cols() ;

                // make sure that this is a square matrix
                BELFEM_ASSERT( mData.mNumRows == mData.mNumCols, "Matrix must be square" );

                mNumRowsPerProc.set_size( mCommSize, mData.mNumRows );
                send( mCommunicationList, mNumRowsPerProc );
                PetscSplitOwnership( mData.mComm, &mData.mMyNumRows, &mData.mNumRows );
                receive( mCommunicationList, mNumRowsPerProc );
                mNumRowsPerProc( 0 ) = mData.mMyNumRows ;
            }
            else
            {
                receive( 0, mData.mNumRows );
                mData.mNumCols = mData.mNumRows ;
                PetscSplitOwnership( mData.mComm, & mData.mMyNumRows, &mData.mNumRows );
                send( 0 , mData.mMyNumRows );
            }
#endif
        }

//------------------------------------------------------------------------------

        // to be called from master only during distribute_indices
        void
        PetscDistributor::compute_memory()
        {
#ifdef BELFEM_PETSC


            if( mMyRank == 0 )
            {

                // temporary vector for row offset, we need two of these,  they only differ
                // by the fact that the first one has one more entry
                Vector< PetscInt > tRowOffsetsPerProcA( mCommSize + 1, 0 );
                Vector< PetscInt > tRowOffsetsPerProcB( mCommSize, 0 );

                // compute the row offsets per proc
                PetscInt tRowOffset = 0 ;
                for( PetscInt p=0; p<mCommSize; ++p )
                {
                    tRowOffsetsPerProcA( p ) += tRowOffset ;
                    tRowOffsetsPerProcB( p ) += tRowOffset ;
                    tRowOffset += mNumRowsPerProc( p );
                }
                tRowOffsetsPerProcA( mCommSize ) = tRowOffset ;

                // data containers for MPI distribution
                Cell< Vector< PetscInt > > tAllPointers( mCommSize, Vector< PetscInt >() );
                Cell< Vector< PetscInt > > tAllRows( mCommSize, Vector< PetscInt >() );
                Cell< Vector< PetscInt > > tAllColumns( mCommSize, Vector< PetscInt >() );

                // grab pointers from matrix
                int * tMatPointers = mMatrix.pointers();
                int * tMatCols     = mMatrix.cols() ;
                int * tMatRows     = mMatrix.rows() ;

                // reset offset
                tRowOffset = 0 ;

                // initialize also offset for columns
                PetscInt tColOffset = 0 ;

                // allocate memory for nonzeros
                mNumNnzPerProc.set_size( mCommSize );

                // now we can compute pointers and indices
                for( PetscInt p=0; p<mCommSize; ++p )
                {
                    // grab vector with pointers
                    Vector< PetscInt > & tPointers = p==0 ? mData.mMyPointers : tAllPointers( p );
                    Vector< PetscInt > & tColumns = p==0 ? mData.mMyColumns : tAllColumns( p );
                    Vector< PetscInt > & tRows = p==0 ? mData.mMyRows : tAllRows( p );

                    // compute NNZ count
                    mNumNnzPerProc( p ) =  tMatPointers[ tRowOffsetsPerProcA( p + 1 ) ]
                            - tMatPointers[ tRowOffsetsPerProcA( p ) ] ;

                    // counter for pointers
                    PetscInt tCount = mNumRowsPerProc( p ) + 1 ;

                    // allocate memory for pointers
                    tPointers.set_size( tCount );

                    // populate pointers
                    for( PetscInt i=0; i<tCount; ++i )
                    {
                        tPointers( i ) = tMatPointers[ tRowOffsetsPerProcB( p ) + i ] - tRowOffset ;
                    }
                    tRowOffset += tPointers(  mNumRowsPerProc( p ) );

                    // allocate memory for row and column indices
                    tCount = mNumNnzPerProc( p );
                    tColumns.set_size( tCount );
                    tRows.set_size( tCount );

                    // populate row and column indices
                    for( PetscInt k=0; k<tCount; ++k )
                    {
                        tColumns( k ) = tMatCols[ tColOffset + k ];
                        tRows( k ) = tMatRows[ tColOffset + k ];
                    }

                    // increment offset
                    tColOffset += tCount ;
                }

                // copy data for my number of rows
                mData.mMyNumRows = mNumRowsPerProc( 0 );
                mData.mMyNumCols = PETSC_DECIDE ;

                // copy data for my number of nonzeros
                mData.mMyNumNnz = mNumNnzPerProc( 0 );

                // send number of nonzeros to other procs
                send( mCommunicationList, mNumNnzPerProc );

                // send row offset to other procs
                mData.mMyRowOffset = 0 ;
                send( mCommunicationList, tRowOffsetsPerProcB );

                // send pointers to other procs
                send( mCommunicationList, tAllPointers );

                // send row indices to other procs
                send( mCommunicationList, tAllRows );

                // send column indices to other procs
                send( mCommunicationList, tAllColumns );
            }
            else
            {
                // number of zeros
                receive( 0, mData.mMyNumNnz );

                // allocate memory for swap
                mSwap.set_size( mData.mMyNumNnz  );

                // row offsets
                receive( 0, mData.mMyRowOffset );

                // matrix pointers and indices
                receive( 0, mData.mMyPointers );
                receive( 0, mData.mMyRows );
                receive( 0, mData.mMyColumns );
            }

            mMatrixData = mMyRank == 0 ? mMatrix.data() : mSwap.data() ;
#endif
        }

//------------------------------------------------------------------------------

        PetscErrorCode
        PetscDistributor::initialize_matrix()
        {
#ifdef BELFEM_PETSC
            // create the matrix
            PetscErrorCode aStatus = MatCreate(
                      mData.mComm,
                    & mData.mMat );

            // assert error
            BELFEM_ASSERT( aStatus==0,
              "PETSc has thrown error %i during MatCreate(): %s",
              ( int ) aStatus,
              petsctools_error_message( aStatus ).c_str() );

            // set the size of the matrix
            aStatus = MatSetSizes(
                    mData.mMat,
                    mData.mMyNumRows,
                    mData.mMyNumCols,
                    mData.mNumRows,
                    mData.mNumCols );

            // assert error
            BELFEM_ASSERT( aStatus==0,
              "PETSc has thrown error %i during MatSetSizes(): %s",
              ( int ) aStatus,
              petsctools_error_message( aStatus ).c_str() );

            // set type to parallel mat
            aStatus = MatSetType( mData.mMat, MATMPIAIJ );

            BELFEM_ASSERT( aStatus==0,
              "PETSc has thrown error %i during MatSetType(): %s",
              ( int ) aStatus,
              petsctools_error_message( aStatus ).c_str() );

            // allocate the memory
            aStatus = MatMPIAIJSetPreallocationCSR(
                    mData.mMat,
                    mData.mMyPointers.data(),
                    mData.mMyColumns.data(),
                    NULL );

            BELFEM_ASSERT( aStatus==0,
                          "PETSc has thrown error %i during MatMPIAIJSetPreallocationCSR(): %s",
                          ( int ) aStatus,
                          petsctools_error_message( aStatus ).c_str() );

            return aStatus ;
#else
            return 0 ;
#endif
        }

//------------------------------------------------------------------------------

        void
        PetscDistributor::initialize_vectors()
        {
#ifdef BELFEM_PETSC
            VecCreateMPI( mData.mComm, mData.mMyNumRows, mData.mNumRows, & mData.mRHS );
            VecCreateMPI( mData.mComm, mData.mMyNumRows, mData.mNumRows, & mData.mLHS );


            mLHS.set_size( mData.mMyNumRows );
            mRHS.set_size( mData.mMyNumRows );

            mData.mVectorIndices.set_size( mData.mMyNumRows ) ;
            for( PetscInt k=0; k<mData.mMyNumRows; ++k )
            {
                mData.mVectorIndices( k ) = k + mData.mMyRowOffset ;
            }
#endif
        }


//------------------------------------------------------------------------------

        PetscErrorCode
        PetscDistributor::update_matrix()
        {
#ifdef BELFEM_PETSC

            // enforce zero based indexing for the sparse matrix
            mMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );

            if( mMyRank == 0 )
            {
                // offset in array
                PetscInt tOffset = mNumNnzPerProc( 0 ) ;

                // get raw pointer of vector
                real * tData = mMatrix.data() ;

                // loop over all procs
                for( proc_t p = 1; p < mCommSize; ++p )
                {
                    index_t tNumSamples =  mNumNnzPerProc( p ) ;
                    send( p, tNumSamples, &tData[ tOffset ] );
                    tOffset +=tNumSamples ;
                }
            }
            else
            {
                receive( 0, mSwap );
            }

            // loop over all rows
            PetscErrorCode aStatus ;

            // get  columns
            PetscInt * tColumns = mData.mMyColumns.data() ;

            // get pointers
            PetscInt * tPointers = mData.mMyPointers.data() ;



            PetscInt tNumCols ;
            PetscInt tThisRow ;
            PetscInt tOffset = 0 ;

            // loop over all rows
            for( PetscInt r=0; r<mData.mMyNumRows; ++r )
            {
                // compute num cols
                tNumCols = tPointers[ r+1 ] - tPointers[ r ];

                // get index of this row
                tThisRow = r + mData.mMyRowOffset ;


                // insert columns
                aStatus = MatSetValues(
                        mData.mMat,
                        1,
                        &tThisRow,
                        tNumCols,
                        &tColumns[ tOffset ],
                        &mMatrixData[ tOffset ],
                        INSERT_VALUES );

                BELFEM_ASSERT( aStatus==0,
                    "PETSc has thrown error %i during MatSetValue(): %s",
                    ( int ) aStatus,
                    petsctools_error_message( aStatus ).c_str() );

                // increment offset
                tOffset += tNumCols ;


            }

            aStatus = MatAssemblyBegin( mData.mMat, MAT_FINAL_ASSEMBLY );

            BELFEM_ASSERT( aStatus==0,
              "PETSc has thrown error %i during MatAssemblyBegin(): %s",
              ( int ) aStatus,
              petsctools_error_message( aStatus ).c_str() );


            aStatus = MatAssemblyEnd( mData.mMat, MAT_FINAL_ASSEMBLY );

            BELFEM_ASSERT( aStatus==0,
              "PETSc has thrown error %i during MatAssemblyEnd(): %s",
              ( int ) aStatus,
              petsctools_error_message( aStatus ).c_str() );

            comm_barrier() ;

            return aStatus ;
#else
            return 0 ;
#endif
        }

//------------------------------------------------------------------------------

        void
        PetscDistributor::recover_lhs( Vector <PetscReal> & aLHS )
        {
            petsctools_get_vector( mData.mLHS, mData.mVectorIndices, mLHS );
            this->collect_vector( aLHS, mLHS );
        }


//------------------------------------------------------------------------------

        void
        PetscDistributor::distribute_vector( Vector <PetscReal> & aGlobal, Vector <PetscReal> & aLocal )
        {
#ifdef BELFEM_PETSC
            comm_data_t tDataType   = get_comm_datatype ( ( PetscReal ) 0 );

            if( mMyRank == 0 )
            {
                // copy my own entries
                for( PetscInt k=0; k<mData.mMyNumRows; ++k )
                {
                    aLocal( k ) = aGlobal( k );
                }

                // offset in array
                PetscInt tOffset = mNumRowsPerProc( 0 ) ;

                // get raw pointer of vector
                PetscReal * tData = aGlobal.data() ;

                // loop over all procs
                for( proc_t p = 1; p < gComm.size(); ++p )
                {
                    // create communication tag
                    int tCommTag = comm_tag( mMyRank, p );

                    // get length of vector
                    index_t tLength = mNumRowsPerProc( p );

                    if( tLength > 0 )
                    {
                        // calculate number of individual messages to be sent
                        Vector< int > tLengths = split_message( tLength );
                        index_t tNumMessages = tLengths.length();

                        // Allocate status and request containers
                        MPI_Status*  tStatus  = ( MPI_Status*  ) alloca( tNumMessages * sizeof( MPI_Status  ) );
                        MPI_Request* tRequest = ( MPI_Request* ) alloca( tNumMessages * sizeof( MPI_Request ) );

                        // loop over all messages
                        for( index_t m=0; m<tNumMessages; ++m )
                        {
                            // increment comm tag
                            tCommTag++;

                            // send data
                            MPI_Isend( &tData[ tOffset ],
                                       tLengths( m ),
                                       tDataType,
                                       p,
                                       tCommTag,
                                       gComm.world(),
                                       &tRequest[ m ] );

                            // increment offset
                            tOffset += tLengths( m );

                            // wait until send is complete
                            MPI_Wait( &tRequest[ m ], &tStatus[ m ] );
                        }
                    }
                }
            }
            else
            {
                // get length of vector
                index_t tLength = mData.mMyNumRows ;

                // create communication tag
                int tCommTag = comm_tag( 0, mMyRank );

                if( tLength > 0 )
                {
                    // calculate number of individual messages to be sent
                    Vector< int > tLengths = split_message( tLength );
                    index_t tNumMessages = tLengths.length();

                    // Allocate status and request containers
                    MPI_Status*  tStatus  = ( MPI_Status*  ) alloca( tNumMessages * sizeof( MPI_Status  ) );
                    MPI_Request* tRequest = ( MPI_Request* ) alloca( tNumMessages * sizeof( MPI_Request ) );

                    // offset in array
                    index_t tOffset = 0;

                    // get data
                    PetscReal * tData = aLocal.data() ;

                    // loop over all messages
                    for( index_t m=0; m<tNumMessages; ++m )
                    {
                        // increment comm tag
                        tCommTag++;

                        // receive data
                        MPI_Irecv( & tData[ tOffset ],
                                   tLengths( m ),
                                   tDataType,
                                   0,
                                   tCommTag,
                                   gComm.world(),
                                   &tRequest[ m ] );

                        // increment offset
                        tOffset += tLengths( m );

                        // wait until receive is complete
                        MPI_Wait( &tRequest[ m ], &tStatus[ m ] );
                    }
                }
            }
#endif
        }

        //------------------------------------------------------------------------------

        void
        PetscDistributor::collect_vector( Vector <PetscReal> & aGlobal, Vector <PetscReal> & aLocal )
        {
#ifdef BELFEM_PETSC
            comm_data_t tDataType   = get_comm_datatype ( ( PetscReal ) 0 );

            if( mMyRank == 0 )
            {
                // copy my own entries
                for( PetscInt k=0; k<mData.mMyNumRows; ++k )
                {
                    aGlobal( k ) = aLocal( k );
                }

                // offset in array
                PetscInt tOffset = mData.mMyNumRows ;

                // get raw pointer of vector
                PetscReal * tData = aGlobal.data() ;

                // loop over all procs
                for( proc_t p = 1; p < gComm.size(); ++p )
                {
                    // get length of vector
                    index_t tLength = mNumRowsPerProc( p );

                    if( tLength > 0 )
                    {

                        // get proc ID of source
                        proc_t tSource = mCommunicationList( p );

                        // create communication tag
                        int tCommTag = comm_tag( mMyRank, tSource );

                        // calculate number of individual messages to be sent
                        Vector< int > tLengths = split_message( tLength );
                        index_t tNumMessages = tLengths.length();

                        // Allocate status and request containers
                        MPI_Status*  tStatus  = ( MPI_Status*  ) alloca( tNumMessages * sizeof( MPI_Status  ) );
                        MPI_Request* tRequest = ( MPI_Request* ) alloca( tNumMessages * sizeof( MPI_Request ) );

                        // loop over all messages
                        for( index_t m=0; m<tNumMessages; ++m )
                        {
                            // increment comm tag
                            tCommTag++;

                            // send data
                            MPI_Irecv( &tData[ tOffset ],
                                       tLengths( m ),
                                       tDataType,
                                       tSource,
                                       tCommTag,
                                       gComm.world(),
                                       &tRequest[ m ] );

                            // increment offset
                            tOffset += tLengths( m );

                            // wait until send is complete
                            MPI_Wait( &tRequest[ m ], &tStatus[ m ] );
                        }
                    }
                }
            }
            else
            {
                // get length of vector
                index_t tLength = mData.mMyNumRows ;

                if( tLength > 0 )
                {

                    // create communication tag
                    int tCommTag = comm_tag( 0, mMyRank );

                    // calculate number of individual messages to be sent
                    Vector< int > tLengths = split_message( aLocal.length() );
                    index_t tNumMessages = tLengths.length();

                    // Allocate status and request containers
                    MPI_Status*  tStatus  = ( MPI_Status*  ) alloca( tNumMessages * sizeof( MPI_Status  ) );
                    MPI_Request* tRequest = ( MPI_Request* ) alloca( tNumMessages * sizeof( MPI_Request ) );

                    // offset in array
                    index_t tOffset = 0;

                    // get data
                    PetscReal * tData = aLocal.data() ;

                    // loop over all messages
                    for( index_t m=0; m<tNumMessages; ++m )
                    {
                        // increment comm tag
                        tCommTag++;

                        // receive data
                        MPI_Isend( & tData[ tOffset ],
                                   tLengths( m ),
                                   tDataType,
                                   0,
                                   tCommTag,
                                   gComm.world(),
                                   &tRequest[ m ] );

                        // increment offset
                        tOffset += tLengths( m );

                        // wait until receive is complete
                        MPI_Wait( &tRequest[ m ], &tStatus[ m ] );
                    }
                }
            }
#endif
        }
    }

//------------------------------------------------------------------------------

}