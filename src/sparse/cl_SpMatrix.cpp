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

#include <cstring>    // for std::memcpy
#include <algorithm>  // for std::sort std::fill_n

#include "assert.hpp"
#include "cl_SpMatrix.hpp"
#include "fn_max.hpp"
#include "fn_unique.hpp"
#include "fspblas.hpp"
#include "cl_Logger.hpp"

#include "hdf5_tools.hpp"
#include "cl_HDF5.hpp"
#include "stringtools.hpp"
#include "fn_Graph_sort.hpp"
#include "fspblas.hpp"

#ifdef BELFEM_MPI
#include "commtools.hpp"
#endif

namespace belfem
{


//------------------------------------------------------------------------------

    SpMatrix::SpMatrix( Graph & aGraph,
                        const enum SpMatrixType aType,
                        const index_t aNumRows,
                        const index_t aNumCols ,
                        const bool aSortGraph ) :
            mType( aType )
    {
        if ( aSortGraph ) graph::sort( aGraph );

#ifdef DEBUG
        // guarded by DEBUG, not BELFEM_DEBUG: the latter is defined nowhere,
        // so this check had never run in any configuration ( revived
        // 2026-08-16 ). It is O( V log V ) per construction — debug-build
        // cost only — and it has never been exercised on production graphs,
        // so a first firing is a finding to investigate, not automatically
        // a regression in the caller
        this->check_graph( aGraph );
#endif
        if ( aNumRows == 0 && aNumCols == 0 )
        {
            // set the sizes and tidy indices in graph
            this->set_sizes(
                    aGraph.size(),                         // number of rows
                    this->tidy_graph( aGraph ) // returns number of columns
            );
        }
        else
        {
            // sort indices in graph
            if ( aSortGraph )
            {
                for( graph::Vertex * tVertex : aGraph )
                {
                    tVertex->sort_vertices();
                }
            }
            // take sizes from input
            this->set_sizes( aNumRows, aNumCols );
        }

        // check which kind of matrix this is
        switch ( aType )
        {
            case ( SpMatrixType::CSR ) :
            {
                this->create_csr_indices( aGraph );
                break;
            }
            case( SpMatrixType::CSC ):
            {
                this->create_csc_indices( aGraph );
                break;
            }
            default:
            {
                BELFEM_ERROR( false, "Unknown SpMatrixType" );
                break;
            }
        }

        // allocate the data container
        this->allocate_values();

        // fill container with zeros
        this->fill( 0.0 );

#ifdef BELFEM_NETLIB
        this->allocate_swap();
#endif
        this->set_indexing_base( SpMatrixIndexingBase::Cpp );
    }

//------------------------------------------------------------------------------

    SpMatrix::SpMatrix( SpMatrixType aType,
                  const index_t aNumRows,
                  const index_t aNumCols,
                  const index_t aNumNonZeros,
                  const int_t  * aIndices,
                  const int_t  * aPointers ) :
        mType( aType )
    {
        mNumRows = aNumRows ;
        mNumCols = aNumCols ;
        mNumNonZeros = aNumNonZeros ;

        mPointerSize = ( aType == SpMatrixType::CSC )
                      ? mNumCols + 1
                      : mNumRows + 1;

        // allocate pointer array
        mPointers = ( int_t * ) malloc( ( mPointerSize ) * sizeof( int_t ) );
        std::memcpy( mPointers, aPointers, ( mPointerSize ) * sizeof( int_t ) );


        switch ( aType )
        {
            case ( SpMatrixType::CSR ) :
            {
                // allocate column array
                mColumns = ( int_t * ) malloc( ( mNumNonZeros ) * sizeof( int_t ) );
                std::memcpy( mColumns, aIndices, ( mNumNonZeros ) * sizeof( int_t ) );
                break ;
            }
            case( SpMatrixType::CSC ) :
            {
                // allocate row array
                mRows = ( int_t * ) malloc( ( mNumNonZeros ) * sizeof( int_t ) );
                std::memcpy( mRows, aIndices, ( mNumNonZeros ) * sizeof( int_t ) );
                break ;
            }
            default:
            {
                BELFEM_ERROR( false, "Invalid Matrix type");
            }
        }
        this->allocate_values();
        this->fill( 0.0 );

        // ensure indices are sorted (external arrays may not guarantee this)
        this->sort_entries();

        this->set_indexing_base( SpMatrixIndexingBase::Cpp );
    }

//------------------------------------------------------------------------------

    SpMatrix::SpMatrix( const string & aHDF5Path, const string aLabel )
    {
        this->load( aHDF5Path, aLabel );
    }

//------------------------------------------------------------------------------

    SpMatrix::SpMatrix( const hid_t aParent, const string aLabel )
    {
        herr_t tStatus = 0 ;
        hid_t tGroup = hdf5::open_group( aLabel, aParent );
        this->load( tGroup, tStatus );
        BELFEM_ERROR( hdf5::close_group( tGroup ) == 0, "Something went wrong while trying to load the matrix %s",
                      aLabel.c_str() );
    }

//------------------------------------------------------------------------------

    SpMatrix::SpMatrix( const Matrix< real > & aMatrix, const SpMatrixType aType )
    {
        // set sizes
        this->set_sizes( aMatrix.n_rows(), aMatrix.n_cols() );

        // count nonzeros
        mNumNonZeros = 0 ;
        for( int_t i=0; i<mNumRows; ++i )
        {
            for( int_t j=0; j<mNumCols; ++j )
            {
                if( aMatrix( i, j ) != 0.0 )
                {
                    ++mNumNonZeros;
                }
            }
        }


#ifdef BELFEM_NETLIB
        this->allocate_swap();
#endif

        this->allocate_values() ;

        // set type of matrix
        mType = aType ;

        switch ( aType )
        {
            case ( SpMatrixType::CSR ) :
            {
                mPointerSize = mNumRows + 1;

                // allocate pointer array
                mPointers = ( int_t * ) malloc( ( mPointerSize ) * sizeof( int_t ) );

                // first entry
                mPointers[ 0 ] = 0 ;

                // allocate column array
                mColumns = ( int_t * ) malloc( ( mNumNonZeros ) * sizeof( int_t ) );

                // position in columns
                int_t tStep = 0 ;

                int_t tCount;

                // populate pointers array
                for( int_t i=0; i<mNumRows; ++i )
                {
                    // reset counter
                    tCount = 0 ;

                    // loop over all columns
                    for( int_t j=0; j<mNumCols; ++j )
                    {
                        if( aMatrix( i, j ) != 0.0 )
                        {
                            // write column
                            mColumns[ tStep ] = j ;

                            // write value
                            mValues[ tStep++ ] = aMatrix( i, j );

                            // increment counter
                            ++tCount;
                        }
                    }

                    // count entries
                    mPointers[ i + 1 ] = mPointers[ i ] + tCount ;
                }

                break;
            }
            case( SpMatrixType::CSC ):
            {
                mPointerSize = mNumCols + 1;

                // allocate pointer array
                mPointers = ( int_t * ) malloc( ( mPointerSize ) * sizeof( int_t ) );

                // first entry
                mPointers[ 0 ] = 0 ;

                // allocate row array
                mRows = ( int_t * ) malloc( ( mNumNonZeros ) * sizeof( int_t ) );

                // position in columns
                int_t tStep = 0 ;

                int_t tCount;

                // populate pointers array
                for( int_t j=0; j<mNumCols; ++j )
                {
                    // reset counter
                    tCount = 0 ;

                    // loop over all columns
                    for( int_t i=0; i<mNumRows; ++i )
                    {
                        if( aMatrix( i, j ) != 0.0 )
                        {
                            // write column
                            mRows[ tStep ] = i ;

                            // write value
                            mValues[ tStep++ ] = aMatrix( i, j );

                            // increment counter
                            tCount++;
                        }
                    }

                    // count entries
                    mPointers[ j + 1 ] = mPointers[ j ] + tCount ;
                }

                break;
            }
            default:
            {
                BELFEM_ERROR( false, "Unknown SpMatrixType" );
                break;
            }
        }

        this->set_indexing_base( SpMatrixIndexingBase::Cpp );
    }

    SpMatrix::SpMatrix( SpMatrix * aParent ) :
        mParent( aParent )
    {
        BELFEM_ERROR( aParent != nullptr, "parent matrix must not be null" );
        BELFEM_ERROR( aParent->mParent == nullptr,
                      "parent is itself a child, chains are not supported" );
        BELFEM_ERROR( aParent->mChild == nullptr,
                      "parent already has a child" );
        BELFEM_ERROR( aParent->mPointers != nullptr,
                      "parent matrix has no structure" );

        // copy type, dimensions and structure pointers from the parent
        this->update_from_parent();

        // the child shares the structure but owns its value array
        this->allocate_values();
        this->fill( 0.0 );

#ifdef BELFEM_NETLIB
        this->allocate_swap();
#endif
        this->set_indexing_base( static_cast< SpMatrixIndexingBase >( aParent->indexing_base() ) );

        mParent->set_child( this );
    }

//------------------------------------------------------------------------------

    SpMatrix::~SpMatrix()
    {
        // unlink from the parent so it never touches a dead child
        if ( mParent != nullptr )
        {
            mParent->mChild = nullptr ;
        }

        if ( mChild != nullptr )
        {
            // transfer ownerships
            mChild->update_from_parent() ;
            mChild->mParent = nullptr ;
            mChild = nullptr ;
            mPointerSize = 0;
            mNumRows = 0;
            mNumCols = 0;
            mNumNonZeros = 0;
            if( mValues != nullptr )
            {
                free( mValues );
                mValues = nullptr;
            }
        }
        else
        {
            this->deallocate();
        }
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::set_child( SpMatrix * aChild )
    {
        BELFEM_ERROR( mChild == nullptr , "Child already set for this matrix" );
        mChild = aChild;
    }

    void
    SpMatrix::update_from_parent()
    {
        BELFEM_ERROR( mParent != nullptr , "No parent matrix assigned" );
        mType          =  mParent->type();
        mNumRows        = mParent->n_rows();
        mNumCols        = mParent->n_cols();
        mNumNonZeros    = mParent->number_of_nonzeros();
        mPointers       = mParent->pointers();
        mPointerSize    = mParent->n_pointers() ;
        mRows           = mParent->rows();
        mColumns        = mParent->cols();
        mHaveCooIndices = mParent->have_coo_indices();
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::deallocate()
    {
        // a parent must never free the structure its child still aliases;
        // all public paths ( destructor, load, copy, move ) unlink or refuse first
        BELFEM_ERROR( mChild == nullptr,
                      "deallocate() must not be called on a matrix with an attached child" );

        if ( mParent == nullptr )
        {
            if( mPointers != nullptr )
            {
                free( mPointers );
                mPointers = nullptr;
                mPointerSize = 0;
            }
            if( mRows != nullptr )
            {
                free( mRows );
                mRows = nullptr;
                mNumRows = 0;
            }
            if( mColumns != nullptr )
            {
                free( mColumns );
                mColumns = nullptr;
                mNumCols = 0;
            }
        }
        else
        {
            mPointerSize = 0;
            mNumRows = 0;
            mNumCols = 0;
        }

        if( mValues != nullptr )
        {
            free( mValues );
            mValues = nullptr;
            mNumNonZeros = 0;
        }

        // the flag must not survive its arrays
        mHaveCooIndices = false;
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::check_graph( Graph & aGraph )
    {
        index_t tNumVertices = aGraph.size() ;

        Vector< index_t > tIndices( tNumVertices+1, BELFEM_UINT_MAX );

        index_t tCount = 0 ;

        for( graph::Vertex * tVertex : aGraph )
        {
            tIndices( tCount++ ) = tVertex->index() ;
        }

        unique( tIndices );

        BELFEM_ERROR( tIndices.length() == tNumVertices+1,
                     "Vertex indices in Graph are not unique" );

    }

//------------------------------------------------------------------------------

    index_t
    SpMatrix::tidy_graph( Graph & aGraph )
    {
        // step 1: update indices in graph
        index_t tCount = 0;
        index_t aMaxIndex = 0;

        for( graph::Vertex * tVertex : aGraph )
        {
            tVertex->set_index( tCount++ );
        }

        // step 2: sort connected vertices according to their indices
        for( graph::Vertex * tVertex : aGraph )
        {
            tVertex->sort_vertices();

            // since the connected vertices are sorted, we know that
            // the last connected vertex has the largest index
            if ( tVertex->number_of_vertices() > 0 )
            {
                aMaxIndex = std::max(
                        tVertex->vertex(  tVertex->number_of_vertices() - 1 )->index(),
                        aMaxIndex );
            }
        }

        return aMaxIndex + 1; // << plus 1, because C++ is zero based
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::create_csr_indices( Graph & aGraph )
    {
        // number of vertices
        index_t tSize = aGraph.size();

        mPointerSize = mNumRows + 1;

        // allocate pointer array
        mPointers = ( int_t * ) malloc( ( mPointerSize ) * sizeof( int_t ) );

        // populate pointer array
        std::fill_n( mPointers, mPointerSize, 0 );

        // create pointer array ( step 1 )
        for ( index_t k = 0; k < tSize; ++k )
        {
            mPointers[ aGraph( k )->index() + 1 ] = ( int_t ) aGraph( k )->number_of_vertices();
        }

        // counter to prevent data type overflow
        index_t tCount = 0;

        // create pointer array ( step 2 )
        for ( int_t k = 1; k < mPointerSize; ++k )
        {
            tCount += mPointers[ k ];
            mPointers[ k ] += mPointers[ k - 1 ];
        }

        // set number of nonzeros and check int_t type boundaries
        this->set_nnz( tCount );

        BELFEM_ASSERT( ( index_t ) mPointers[ mNumRows ] == tCount,
            "Something went wrong while creating CSR index" );

        // allocate index vector
        tCount = ( tCount == 0 ) ? 1 : tCount;

        mColumns = ( int_t * ) malloc( tCount * sizeof( int_t ) );

        // reset counter
        tCount = 0;

        // loop over all nodes and create index array
        for( graph::Vertex * tVertex : aGraph )
        {
            int_t tNumVertices = tVertex->number_of_vertices();

            for( int_t k=0; k<tNumVertices; ++k )
            {
                mColumns[ tCount++ ] = tVertex->vertex( k )->index();
            }
        }
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::create_csc_indices(  Graph & aGraph )
    {
        // number of vertices
        index_t tSize = aGraph.size();

        index_t tNumNonzeros = 0;

        // allocate a counting array
        index_t * tCount = ( index_t * ) malloc( mNumCols * sizeof( index_t ) );

        // fill array with zeros
        std::fill_n( tCount, mNumCols, 0 );

        // count vertices per column
        for ( index_t k = 0; k < tSize; ++k )
        {
            int_t tNumberOfVertices =  aGraph( k )->number_of_vertices();

            for( int_t i=0; i<tNumberOfVertices; ++i )
            {
                ++tCount[ aGraph( k )->vertex( i )->index() ];
            }

            tNumNonzeros +=tNumberOfVertices;
        }

        mPointerSize = mNumCols + 1;

        // populate pointer array
        mPointers = ( int_t * ) malloc( ( mPointerSize ) * sizeof( int_t ) );

        mPointers[ 0 ] = 0;
        for( int_t k=0; k<mNumCols; ++k )
        {
            mPointers[ k+1 ] = mPointers[ k ] + tCount[ k ];
        }

        // reset counter
        std::fill_n( tCount, mNumCols, 0 );

        // set number of nonzeros and check int_t type boundaries
        this->set_nnz( tNumNonzeros );

        // populate indices
        tNumNonzeros = ( tNumNonzeros == 0 ) ? 1 : tNumNonzeros;

        mRows = ( int_t * ) malloc( tNumNonzeros * sizeof( int_t ) );

        int_t n = ( int_t ) tSize;

        // count vertices per column
        for ( int_t k = 0; k < n; ++k )
        {
            int_t tNumberOfVertices =  aGraph( k )->number_of_vertices();

            for( int_t i=0; i<tNumberOfVertices; ++i )
            {
                // get column of array
                index_t j = aGraph( k )->vertex( i )->index();

                // write index into array
                mRows[ mPointers[ j ] + tCount[ j ] ] = ( int_t ) aGraph( k )->index();

                // increment counter
                ++tCount[ j ];
            }
        }

        // free counter
        free( tCount );
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::set_sizes(
            const index_t aNumRows,
            const index_t aNumCols )
    {
        // Check for overflow when casting index_t to int_t
        const index_t tMaxSize = static_cast< index_t >( std::numeric_limits< int_t >::max() );

        BELFEM_ERROR( aNumRows <= tMaxSize,
                     "Matrix has too many rows (%lu) for int_t type.\n"
                     "Maximum allowed: %lu\n"
                     "Recommendation: Enable BELFEM_INT64 in CMake configuration to use 64-bit integers.",
                     ( long long unsigned int ) aNumRows,
                     ( long long unsigned int ) tMaxSize );

        BELFEM_ERROR( aNumCols <= tMaxSize,
                     "Matrix has too many columns (%lu) for int_t type.\n"
                     "Maximum allowed: %lu\n"
                     "Recommendation: Enable BELFEM_INT64 in CMake configuration to use 64-bit integers.",
                     ( long long unsigned int ) aNumCols,
                     ( long long unsigned int ) tMaxSize );

        // set data
        mNumRows = ( int_t ) aNumRows;
        mNumCols = ( int_t ) aNumCols;
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::set_nnz( const index_t aNumNonZeros )
    {
#ifdef BELFEM_INT64
        const index_t tMaxNNZ = 9223372036854775807 ;
#else
        const index_t tMaxNNZ = 2147483647 ;
#endif
        // make sure that NNZ is OK
        BELFEM_ERROR( aNumNonZeros < tMaxNNZ,
                     "too many non-zeros in matrix (%lu > %lu )",
                     ( long unsigned int ) aNumNonZeros,
                     ( long unsigned int ) tMaxNNZ );

        // set data
        mNumNonZeros = ( int_t ) aNumNonZeros;
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::allocate_values()
    {
        // One slot beyond the nonzeros. index() / position() return
        // mNumNonZeros for an entry that is not in the sparsity pattern, and
        // the writable operator() and the assembly routines guard that with
        // BELFEM_ASSERT only - so in a release build, assembling into an
        // entry the pattern does not contain used to write one element past
        // the allocation. The extra slot turns that into a harmless dump.
        //
        // The slot is zeroed because assembly accumulates ( += ) into it,
        // which reads it first; fill() only covers the logical length.
        // mNumNonZeros stays the logical length everywhere else - fill(),
        // HDF5 I/O, std::copy and MPI transfers never touch the dump slot.
        mValues = ( real * ) malloc( ( mNumNonZeros + 1 ) * sizeof( real ) );
        BELFEM_ERROR( mValues != nullptr,
                      "failed to allocate value array" );

        mValues[ mNumNonZeros ] = 0.0 ;
    }
//------------------------------------------------------------------------------

#ifdef BELFEM_NETLIB
    /**
     * allocate swap vector
     */
    void
    SpMatrix::allocate_swap()
    {
        mSwapSize = ( mNumRows > mNumCols ) ? mNumRows : mNumCols ;
        mSwap.set_size( mSwapSize );
    }
#endif

//------------------------------------------------------------------------------

    void
    SpMatrix::fill( const real aValue )
    {
        std::fill_n( mValues, mNumNonZeros, aValue );
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::sort_entries()
    {
        // sorting permutes the shared index array together with the own value
        // array — on a linked pair either direction desyncs the sibling
        BELFEM_ERROR( mParent == nullptr, "can't sort a child matrix" );
        BELFEM_ERROR( mChild  == nullptr, "can't sort a parent matrix while a child is attached" );

        // determine which index array to check
        int_t * tIndices = nullptr ;
        int_t   tNumSlices = 0 ;

        if ( mType == SpMatrixType::CSR )
        {
            tIndices   = mColumns ;
            tNumSlices = mNumRows ;
        }
        else if ( mType == SpMatrixType::CSC )
        {
            tIndices   = mRows ;
            tNumSlices = mNumCols ;
        }
        else
        {
            return ;
        }

        // indexing base offset (0 for C++, 1 for Fortran)
        int_t tBase = mPointers[ 0 ] ;

        // reusable work buffer (clear() preserves allocation)
        Cell< std::pair< int_t, real > > tWork ;

        // check each row (CSR) or column (CSC)
        for ( int_t i = 0; i < tNumSlices; ++i )
        {
            int_t tBegin = mPointers[ i ] - tBase ;
            int_t tEnd   = mPointers[ i + 1 ] - tBase ;

            if ( tEnd - tBegin <= 1 ) continue ;

            // check if already sorted
            bool tSorted = true ;
            for ( int_t k = tBegin; k < tEnd - 1; ++k )
            {
                if ( tIndices[ k ] > tIndices[ k + 1 ] )
                {
                    tSorted = false ;
                    break ;
                }
            }

            if ( !tSorted )
            {
                // gather index-value pairs
                tWork.clear() ;
                for ( int_t k = tBegin; k < tEnd; ++k )
                {
                    tWork.push( { tIndices[ k ], mValues[ k ] } ) ;
                }

                // sort by index
                std::sort( tWork.data(),
                           tWork.data() + tWork.size(),
                           []( const std::pair< int_t, real > & a,
                               const std::pair< int_t, real > & b )
                           {
                               return a.first < b.first ;
                           } ) ;

                // scatter back
                int_t k = tBegin ;
                for ( auto & tPair : tWork )
                {
                    tIndices[ k ] = tPair.first ;
                    mValues[ k ]  = tPair.second ;
                    ++k ;
                }
            }
        }
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::set_indexing_base( const enum SpMatrixIndexingBase & aBasis )
    {
        if ( mParent != nullptr )
        {
            if ( mParent->indexing_base() != static_cast< int_t >( aBasis ) )
            {
                mParent->set_indexing_base( aBasis );
            }
        }

        switch( aBasis )
        {
            case( SpMatrixIndexingBase::Cpp ) :
            {
                // test if this is in fortran base
                if( mPointers[ 0 ] == 1 && mParent == nullptr )
                {
                    // decrement all pointers
                    std::for_each( mPointers, mPointers + mPointerSize,
                                   [ ]( int_t & tValue ){ --tValue; } );

                    // decrement all rows
                    if( mRows != nullptr )
                    {
                        std::for_each( mRows, mRows + mNumNonZeros,
                                       [ ]( int_t & tValue ){ --tValue; } );
                    }

                    // decrement all cols
                    if( mColumns != nullptr )
                    {
                        std::for_each( mColumns, mColumns + mNumNonZeros,
                                       [ ]( int_t & tValue ){ --tValue; } );
                    }
                }

                // no index function to select: position() reads the base
                // from mPointers[ 0 ] on every call
                break;
            }
            case( SpMatrixIndexingBase::Fortran ) :
            {
                // test if this is in c++ base
                if( mPointers[ 0 ] == 0 && mParent == nullptr )
                {
                    // decrement all pointers
                    std::for_each( mPointers, mPointers + mPointerSize,
                                   [ ]( int_t & tValue ){ ++tValue; } );

                    // increment all rows
                    if( mRows != nullptr )
                    {
                        std::for_each( mRows, mRows + mNumNonZeros,
                                       [ ]( int_t & tValue ){ ++tValue; } );
                    }

                    // increment all cols
                    if( mColumns != nullptr )
                    {
                        std::for_each( mColumns, mColumns + mNumNonZeros,
                                       [ ]( int_t & tValue ){ ++tValue; } );
                    }

                }
                // no index function to select: position() reads the base
                // from mPointers[ 0 ] on every call
                break;
            }
            default:
            {
                BELFEM_ERROR( false, "Unknown SpMatrixIndexingBase" );
                break;
            }
        }

        if ( mChild != nullptr )
        {
            // no guard possible here: mChild->indexing_base() reads the shared
            // mPointers[ 0 ], which is already converted at this point; the
            // child converts nothing of its own ( mParent != nullptr ) and only
            // forwards the base down the chain
            mChild->set_indexing_base( aBasis );
        }
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::create_coo_indices()
    {
        mHaveCooIndices = true ;
        if ( mParent != nullptr )
        {
            if ( ! mParent->have_coo_indices() )
            {
                mParent->create_coo_indices();
            }

            // always re-alias, regardless of which side created the array:
            // the child never computes or owns coo indices of its own
            switch( mType )
            {
                case( SpMatrixType::CSR ) :
                {
                    mRows = mParent->rows();
                    break ;
                }
                case( SpMatrixType::CSC ) :
                {
                    mColumns = mParent->cols();
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Unknown SpMatrixType" );
                    break;
                }
            }
            return ;
        }

        // this one genuinely rewrites the index arrays, so it needs them in
        // a known base and must put the matrix back as it found it.
        // ( multiply() no longer does this - its kernels take the base. )
        int_t tOldBase = this->indexing_base();
        this->set_indexing_base( SpMatrixIndexingBase::Cpp );

        switch( mType )
        {
            case( SpMatrixType::CSR ) :
            {
                // test if rows already exist
                if( mRows == nullptr )
                {
                    if( mNumNonZeros == 0 )
                    {
                        mRows = ( int_t * ) malloc( 1 * sizeof( int_t ));
                    }
                   else
                    {
                        mRows = ( int_t * ) malloc( mNumNonZeros * sizeof( int_t ));
                    }

                    // populate row indices
                    int_t tCount = 0;
                    for ( int_t k = 0; k < mNumRows; ++k )
                    {
                        int_t tN = mPointers[ k + 1 ] - mPointers[ k ];
                        for ( int_t i = 0; i < tN; ++i )
                        {
                            mRows[ tCount++ ] = k;
                        }
                    }
                }
                break;
            }
            case( SpMatrixType::CSC ) :
            {
                if( mColumns == nullptr )
                {
                    if ( mNumNonZeros == 0 )
                    {
                        mColumns = ( int_t * ) malloc( 1 * sizeof( int_t ));
                    }
                    else
                    {
                        mColumns = ( int_t * ) malloc( mNumNonZeros * sizeof( int_t ));
                    }
                    // populate row indices
                    int_t tCount = 0;
                    for ( int_t k = 0; k < mNumCols; ++k )
                    {
                        int_t tN = mPointers[ k + 1 ] - mPointers[ k ];
                        for ( int_t i = 0; i < tN; ++i )
                        {
                            mColumns[ tCount++ ] = k;
                        }
                    }
                }
                break;
            }
            default:
            {
                BELFEM_ERROR( false, "Unknown SpMatrixType" );
                break;
            }
        }

        // restore indexing base
        if( tOldBase == 1 )
        {
            this->set_indexing_base( SpMatrixIndexingBase::Fortran );
        }

        if ( mChild != nullptr )
        {
            if ( ! mChild->have_coo_indices() )
            {
                mChild->create_coo_indices();
            }
        }
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::free_coo_indices()
    {
        // the flag must fall before any recursion, both here and in the
        // parent branch below — it is the termination condition
        mHaveCooIndices = false ;

        if ( mParent != nullptr )
        {
            // drop the alias; the parent owns and frees the array
            switch( mType )
            {
                case( SpMatrixType::CSR ) :
                {
                    mRows = nullptr ;
                    break ;
                }
                case( SpMatrixType::CSC ) :
                {
                    mColumns = nullptr ;
                    break ;
                }
                default :
                {
                    break ;
                }
            }

            if ( mParent->have_coo_indices() )
            {
                mParent->free_coo_indices();
            }
            return;
        }

        switch( mType )
        {
            case( SpMatrixType::CSR ) :
            {
                if ( mRows != nullptr )
                {
                    free( mRows );
                    mRows = nullptr;
                }
                break;
            }
            case(  SpMatrixType::CSC ) :
            {
                if ( mColumns != nullptr )
                {
                    free( mColumns );
                    mColumns = nullptr;
                }
                break;
            }
            default :
            {
                break;
            }
        }

        if ( mChild != nullptr )
        {
            // unconditional: the child must drop its alias even if its own
            // flag is already false, otherwise it keeps a dangling pointer
            mChild->free_coo_indices();
        }
    }

//------------------------------------------------------------------------------

    /**
     * print_t the matrix to the screen ( for debugging )
     */
    void
    SpMatrix::print( const string aLabel  )
    {
        int_t tCount = 0;

        std::cout<< "SpMatrix " << aLabel << " (" << mNumRows << ", " << mNumCols << ") " << std::endl;

        if ( mType == SpMatrixType::CSR )
        {
            for ( int_t r = 0; r < mNumRows; ++r )
            {
                // get bandwidth
                int_t n = mPointers[ r + 1 ] - mPointers[ r ];

                for ( int_t c = 0; c < n; ++c )
                {
                    std::cout << tCount << " : (" << r << ", " << mColumns[ tCount ] << ")  : " << mValues[ tCount ] << std::endl;
                    ++tCount;
                }
            }
        }
        else if ( mType == SpMatrixType::CSC )
        {
            for ( int_t r = 0; r < mNumCols; ++r )
            {
                // get bandwidth
                int_t n = mPointers[ r + 1 ] - mPointers[ r ];

                for ( int_t c = 0; c < n; ++c )
                {
                    std::cout <<  tCount << " : (" << r << ", " << mRows[ tCount ] << ")  : " << mValues[ tCount ] << std::endl;
                    ++tCount;
                }
            }
        }
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::print2( const string aLabel  )
    {
        std::cout<< "SpMatrix " << aLabel << " (" << mNumRows << ", " << mNumCols << ") " << std::endl;

        for( int_t k=0; k<mPointerSize; ++k )
        {
            std::cout << "P ( " << k << " ) = " <<mPointers[ k ] << std::endl;
        }

        if( mType == SpMatrixType::CSC )
        {
            for ( int_t k = 0; k < mNumNonZeros; ++k )
            {
                std::cout << "R ( " << k << " ) = " << mRows[ k ] << std::endl;
            }
        }
        else if( mType == SpMatrixType::CSR )
        {
            for ( int_t k = 0; k < mNumNonZeros; ++k )
            {
                std::cout << "C ( " << k << " ) = " << mColumns[ k ] << std::endl;
            }
        }

    }

//------------------------------------------------------------------------------

    void
    SpMatrix::save(
            const string & aPath,
            const string   aLabel,
            const enum FileMode aMode)
    {
#ifdef BELFEM_HDF5
        // create a new file
        HDF5 tFile( aPath, aMode, false );

        // get status from file
        herr_t & tStatus = tFile.status();

        // create a group with the specified label
        hid_t tGroup = tFile.create_group( aLabel );

        // save matrix into this group
        this->save( tGroup, tStatus );

        // close HDF5 file
        tFile.close();
#else
        BELFEM_ERROR( false, "Trying to save a sparse matrix to HDF5, but BELFEM is not link against HDF5 libraries." );
#endif
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::save(   hid_t        & aGroup,
            herr_t       & aStatus )
    {
#ifdef BELFEM_HDF5
        // the format label of this file
        string tFormatLabel;

        switch( mType )
        {
            case( SpMatrixType::CSR ) :
            {
                tFormatLabel = "CSR";
                break;
            }
            case( SpMatrixType::CSC ) :
            {
                tFormatLabel = "CSC";
                break;
            }
            default :
            {
                tFormatLabel = "unknown";
                break;
            }
        }

        // save format
        hdf5::save_string_to_file( aGroup, "Format", tFormatLabel, aStatus );

        // save sizes
        hdf5::save_scalar_to_file( aGroup, "NumRows", mNumRows, aStatus );
        hdf5::save_scalar_to_file( aGroup, "NumCols", mNumCols, aStatus );
        // save number of nonzeros
        hdf5::save_scalar_to_file( aGroup, "NumNonZeros", mNumNonZeros, aStatus );

        // save pointers
        hdf5::save_array_to_file( aGroup, "Pointers", mPointers, mPointerSize, aStatus );

        // save indices
        if( mType == SpMatrixType::CSC )
        {
            hdf5::save_array_to_file( aGroup, "Indices", mRows, mNumNonZeros, aStatus );
        }
        else if ( mType == SpMatrixType::CSR )
        {
            hdf5::save_array_to_file( aGroup, "Indices", mColumns, mNumNonZeros, aStatus );
        }

        // save data
        hdf5::save_array_to_file( aGroup, "Values", mValues, mNumNonZeros, aStatus );
#else
        BELFEM_ERROR( false, "Trying to save a sparse matrix to HDF5, but BELFEM is not link against HDF5 libraries." );
#endif
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::load(
            const string & aPath,
            const string aLabel )
    {
#ifdef BELFEM_HDF5
        // create a new file
        HDF5 tFile( aPath, FileMode::OPEN_RDONLY, false );

        hid_t   tGroup = tFile.select_group( aLabel );
        herr_t & tStatus = tFile.status();


        this->load( tGroup, tStatus );

        tFile.close_active_group();
        tFile.close();
#else
        BELFEM_ERROR( false, "Trying to load a sparse matrix from HDF5, but BELFEM is not link against HDF5 libraries." );
#endif
    }

//------------------------------------------------------------------------------

    /**
     * load matrix from a specific group in a hdf5 file
     */
    void
    SpMatrix::load( hid_t & aGroup, herr_t & aStatus )
    {
#ifdef BELFEM_HDF5
        // a linked matrix ( child, or parent with an attached child ) never
        // replaces its structure: the file is verified against the existing
        // arrays and only the values are loaded. Unlinked matrices reload
        // everything as before.
        const bool tLinked = ( mParent != nullptr ) || ( mChild != nullptr );

        if ( ! tLinked )
        {
            this->deallocate();
        }

        string tFormatLabel;

        // load format
        hdf5::load_string_from_file( aGroup, "Format", tFormatLabel, aStatus );

        string tFormat = string_to_upper( tFormatLabel );

        // resolve type
        SpMatrixType tType = SpMatrixType::UNDEFINED ;
        if( tFormat == "CSR" )
        {
            tType = SpMatrixType::CSR;
        }
        else if ( tFormat == "CSC" )
        {
            tType = SpMatrixType::CSC ;
        }
        BELFEM_ERROR( tType != SpMatrixType::UNDEFINED,
                      "unknown sparse matrix format in file: %s",
                      tFormatLabel.c_str() );

        if ( tLinked )
        {
            BELFEM_ERROR( tType == mType,
                          "matrix format of file does not match" );

            // load sizes into locals, never overwriting the linked structure
            int_t tNumRows     = 0 ;
            int_t tNumCols     = 0 ;
            int_t tNumNonZeros = 0 ;
            hdf5::load_scalar_from_file( aGroup, "NumRows", tNumRows, aStatus );
            hdf5::load_scalar_from_file( aGroup, "NumCols", tNumCols, aStatus );
            hdf5::load_scalar_from_file( aGroup, "NumNonZeros", tNumNonZeros, aStatus );

            BELFEM_ERROR( tNumRows == mNumRows
                       && tNumCols == mNumCols
                       && tNumNonZeros == mNumNonZeros,
                          "incompatible matrix adjacency" );

            // compare the pointer array, tolerating an indexing base offset;
            // mismatches are recorded and raised only after the scratch is
            // freed, so the error path does not leak
            int_t * tSwap = ( int_t * ) malloc( mPointerSize * sizeof( int_t ) );
            BELFEM_ERROR( tSwap != nullptr || mPointerSize == 0,
                          "failed to allocate scratch array" );
            hdf5::load_array_from_file( aGroup, "Pointers", tSwap, mPointerSize, aStatus );
            BELFEM_ERROR( aStatus >= 0, "failed to read pointer array from file" );

            int_t tShift = mPointers[ 0 ] - tSwap[ 0 ];

            int_t tMismatch = -1;
            for ( int_t k = 0; k < mPointerSize; ++k )
            {
                if ( tSwap[ k ] + tShift != mPointers[ k ] )
                {
                    tMismatch = k;
                    break;
                }
            }
            free( tSwap );

            BELFEM_ERROR( tShift >= -1 && tShift <= 1,
                          "trying to load incompatible matrix (illegal indexing base)" );
            BELFEM_ERROR( tMismatch < 0,
                          "trying to load incompatible matrix (pointer mismatch at %ld)",
                          ( long ) tMismatch );

            // compare the index array against the own structural array
            tSwap = ( int_t * ) malloc( mNumNonZeros * sizeof( int_t ) );
            BELFEM_ERROR( tSwap != nullptr || mNumNonZeros == 0,
                          "failed to allocate scratch array" );
            hdf5::load_array_from_file( aGroup, "Indices", tSwap, mNumNonZeros, aStatus );
            BELFEM_ERROR( aStatus >= 0, "failed to read index array from file" );

            const int_t * tComp = mType == SpMatrixType::CSC ? mRows : mColumns;

            for ( int_t k = 0; k < mNumNonZeros; ++k )
            {
                if ( tSwap[ k ] + tShift != tComp[ k ] )
                {
                    tMismatch = k;
                    break;
                }
            }
            free( tSwap );

            BELFEM_ERROR( tMismatch < 0,
                          "trying to load incompatible matrix (index mismatch at %ld)",
                          ( long ) tMismatch );

            // structure verified: load the values into the existing buffer
            if ( mValues == nullptr )
            {
                this->allocate_values();
            }
            hdf5::load_array_from_file( aGroup, "Values", mValues, mNumNonZeros, aStatus );

            // the check above proved the indices equal the own, already
            // sorted arrays, so no sort is needed; the indexing base of the
            // shared structure is deliberately left untouched
        }
        else
        {
            mType = tType ;

            // load sizes
            hdf5::load_scalar_from_file( aGroup, "NumRows", mNumRows, aStatus );
            hdf5::load_scalar_from_file( aGroup, "NumCols", mNumCols, aStatus );

            // load number of nonzeros
            hdf5::load_scalar_from_file( aGroup, "NumNonZeros", mNumNonZeros, aStatus );

            // determine pointer size
            if( mType == SpMatrixType::CSC )
            {
                mPointerSize = mNumCols + 1;
            }
            else
            {
                mPointerSize = mNumRows + 1;
            }

            // make sure that array is not used
            if( mPointers != nullptr )
            {
                free( mPointers );
            }
            if( mRows != nullptr )
            {
                free( mRows );
            }
            if( mColumns != nullptr )
            {
                free( mColumns );
            }

            // allocate memory
            mPointers = ( int_t * ) malloc( mPointerSize * sizeof( int_t ) );
            hdf5::load_array_from_file( aGroup, "Pointers", mPointers,  mPointerSize, aStatus );

            // load indices
            if( mType == SpMatrixType::CSC )
            {
                mRows = ( int_t * ) malloc( mNumNonZeros * sizeof( int_t ) );
                hdf5::load_array_from_file( aGroup, "Indices", mRows, mNumNonZeros, aStatus );
            }
            else if ( mType == SpMatrixType::CSR )
            {
                mColumns = ( int_t * ) malloc( mNumNonZeros * sizeof( int_t ) );
                hdf5::load_array_from_file( aGroup, "Indices", mColumns, mNumNonZeros, aStatus );
            }

            if( mValues != nullptr )
            {
                free ( mValues );
            }

            this->allocate_values();

            // load values
            hdf5::load_array_from_file( aGroup, "Values", mValues, mNumNonZeros, aStatus );
        }

        // same-base call: converts nothing ( the arrays are already
        // zero-based ), only forwards the base to a linked child
        if( mPointers[ 0 ] == 0 )
        {
            this->set_indexing_base( SpMatrixIndexingBase::Cpp );
        }
        else
        {
            this->set_indexing_base( SpMatrixIndexingBase::Fortran );
        }

        // ensure indices are sorted (external files may not guarantee this);
        // linked matrices verified their structure above and must not sort
        if ( ! tLinked )
        {
            this->sort_entries();
        }

#ifdef BELFEM_NETLIB
        this->allocate_swap();
#endif
#else
        BELFEM_ERROR( false, "Trying to load a matrix from HDF5, but BELFEM is not link against HDF5 libraries." );
#endif
    }

//------------------------------------------------------------------------------

    // The four index_csr/csc_{zero,one}_based search functions lived here.
    // They are replaced by the single inline SpMatrix::position() in the
    // header, which handles both types and both bases and, unlike a cached
    // member-function pointer, cannot go stale when a parent is rebased.

//------------------------------------------------------------------------------

    void
    SpMatrix::positions_in_slice( const index_t   aSlice,
                                  const int_t   * aCols,
                                  const uint      aNumCols,
                                  int_t         * aPos ) const
    {
        const int_t tBase  = mPointers[ 0 ];
        const bool  tIsCsr = ( mType == SpMatrixType::CSR );

        BELFEM_ASSERT( aSlice < ( index_t ) ( tIsCsr ? mNumRows : mNumCols ),
                      "slice index out of bounds ( %lu >= %lu )",
                      ( long unsigned int ) aSlice,
                      ( long unsigned int ) ( tIsCsr ? mNumRows : mNumCols ) );

#ifndef NDEBUG
        // strictly ascending: within one element the dof indices are unique,
        // so equal neighbors cannot occur and the merge below does not
        // handle them. If this ever fires, dof construction changed and the
        // merge needs a tie case - it must NOT be relaxed to non-strict.
        for ( uint k = 1; k < aNumCols; ++k )
        {
            BELFEM_ASSERT( aCols[ k ] > aCols[ k - 1 ],
                "column list must be strictly ascending ( entry %u: %ld <= %ld )",
                ( unsigned int ) k,
                ( long int ) aCols[ k ],
                ( long int ) aCols[ k - 1 ] );
        }
#endif

        const int_t * tIndices = tIsCsr ? mColumns : mRows ;

        int_t tCursor = mPointers[ aSlice ]     - tBase ;
        const int_t tEnd = mPointers[ aSlice + 1 ] - tBase ;

        // one pass through the slice and one through aCols, advancing
        // whichever is behind: O( slice + cols ) streaming instead of
        // aNumCols independent O( log slice ) probes
        for ( uint k = 0; k < aNumCols; ++k )
        {
            const int_t tTarget = aCols[ k ] + tBase ;

            while ( tCursor < tEnd && tIndices[ tCursor ] < tTarget )
            {
                ++tCursor ;
            }

            aPos[ k ] = ( tCursor < tEnd && tIndices[ tCursor ] == tTarget )
                        ? tCursor
                        : mNumNonZeros ;
        }
    }
//----------------------------------------------------------------------------

#if defined(__clang__) || defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-but-set-variable"
#elif defined(__GNUC__) || defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

#ifdef BELFEM_MKL
    // matrix-vector product y = alpha * op(A) * x + beta * y via the
    // Inspector-Executor sparse BLAS ( the classic mkl_dcsrmm/mkl_dcscmm
    // NIST API was removed in oneMKL 2026 ). The handle wraps the existing
    // arrays without copying them.
    static void
    mkl_sparse_multiply(
            const SpMatrixType   aType,
            const int_t          aNumRows,
            const int_t          aNumCols,
            int_t              * aPointers,
            int_t              * aIndices,
            real               * aValues,
            const real         * aX,
            real               * aY,
            const real           aAlpha,
            const real           aBeta,
            const bool           aTransposedFlag )
    {
        BELFEM_ERROR( aType == SpMatrixType::CSR || aType == SpMatrixType::CSC,
                "Unsupported matrix type." );

        // int_t and MKL_INT have the same width by construction
        // ( see the static_asserts in lapacktools.hpp )
        MKL_INT * tPointers = reinterpret_cast< MKL_INT * >( aPointers );
        MKL_INT * tIndices  = reinterpret_cast< MKL_INT * >( aIndices );

        const sparse_index_base_t tBase = ( aPointers[ 0 ] == 1 ) ?
                SPARSE_INDEX_BASE_ONE : SPARSE_INDEX_BASE_ZERO ;

        sparse_matrix_t tHandle ;
        sparse_status_t tStatus ;

        if ( aType == SpMatrixType::CSR )
        {
            tStatus = mkl_sparse_d_create_csr( &tHandle, tBase,
                    aNumRows, aNumCols,
                    tPointers, tPointers + 1, tIndices, aValues );
        }
        else
        {
            tStatus = mkl_sparse_d_create_csc( &tHandle, tBase,
                    aNumRows, aNumCols,
                    tPointers, tPointers + 1, tIndices, aValues );
        }

        BELFEM_ERROR( tStatus == SPARSE_STATUS_SUCCESS,
                "mkl_sparse_d_create has thrown an error: %i", ( int ) tStatus );

        matrix_descr tDescr ;
        tDescr.type = SPARSE_MATRIX_TYPE_GENERAL ;
        tDescr.mode = SPARSE_FILL_MODE_LOWER ;
        tDescr.diag = SPARSE_DIAG_NON_UNIT ;

        tStatus = mkl_sparse_d_mv(
                aTransposedFlag ? SPARSE_OPERATION_TRANSPOSE
                                : SPARSE_OPERATION_NON_TRANSPOSE,
                aAlpha, tHandle, tDescr, aX, aBeta, aY );

        BELFEM_ERROR( tStatus == SPARSE_STATUS_SUCCESS,
                "mkl_sparse_d_mv has thrown an error: %i", ( int ) tStatus );

        mkl_sparse_destroy( tHandle );
    }

#endif
//------------------------------------------------------------------------------

    void
    SpMatrix::multiply(
            const Vector< real > & aX,
                  Vector< real > & aY,
            const real aAlpha,
            const real aBeta,
            const bool aTransposedFlag )
    {
        int_t tLengthX = mNumCols;
        int_t tLengthY = mNumRows;

        // flip X and Y length if matrix is transposed
        if( aTransposedFlag )
        {
            tLengthX = mNumRows;
            tLengthY = mNumCols;
        }

        BELFEM_ASSERT( (index_t) tLengthX == (index_t) aX.length(),
                      "Number of rows and cols of matrix and vector does not match ( %lu and %lu ).",
                      ( long unsigned int ) tLengthX,
                      ( long unsigned int ) aX.length() );

        BELFEM_ASSERT( (index_t) aY.length() == ( index_t ) tLengthY,
            "Number of cols of vector does not match ( is %lu, but should be %lu )",
                      ( long unsigned int ) aY.length(),
                      ( long unsigned int ) tLengthY );

        // Nothing here converts the indexing base any more. Doing so used to
        // cost two full rewrites of the index arrays per matvec, and - because
        // parent and child share those arrays while caching their own base -
        // it left a sibling reading the wrong base for the duration of the
        // call. MKL detects the base from aPointers[ 0 ] itself; the Fortran
        // kernels take it as an argument ( see the #else branch ).

#ifdef BELFEM_MKL

        mkl_sparse_multiply(
                mType,
                mNumRows,
                mNumCols,
                mPointers,
                mType == SpMatrixType::CSR ? mColumns : mRows,
                mValues,
                aX.data(),
                aY.data(),
                aAlpha,
                aBeta,
                aTransposedFlag );

#else
        // the base is handed to the kernel instead of being converted into
        // the matrix and back
        const int_t tBase = this->indexing_base();

        // the kernels overwrite y, so beta * y must be saved first
        if( aBeta != 0.0 )
        {
            mSwap = aY * aBeta;
        }

        // transposed products use the CSR <-> CSC duality: op(A) keeps the
        // same arrays and calls the sibling kernel with swapped dimensions
        const bool tRunCsr = ( mType == SpMatrixType::CSR ) != aTransposedFlag ;

        const int_t * tNumRows = aTransposedFlag ? &mNumCols : &mNumRows ;
        const int_t * tNumCols = aTransposedFlag ? &mNumRows : &mNumCols ;

        int_t * tIndices = ( mType == SpMatrixType::CSR ) ? mColumns : mRows ;

        switch( mType )
        {
            case( SpMatrixType::CSC ) :
            case( SpMatrixType::CSR ) :
            {
                if ( tRunCsr )
                {
                    matvec_csr(
                        tNumRows,
                        tNumCols,
                        & mNumNonZeros,
                        mValues,
                        tIndices,
                        mPointers,
                        aX.data(),
                        aY.data(),
                        & tBase );
                }
                else
                {
                    matvec_csc(
                        tNumRows,
                        tNumCols,
                        & mNumNonZeros,
                        mValues,
                        tIndices,
                        mPointers,
                        aX.data(),
                        aY.data(),
                        & tBase );
                }
                break;
            }
            default :
            {
                BELFEM_ERROR( false, "Unsupported matrix type.");
            }
        }
        if( aAlpha != 1.0 )
        {
            aY *= aAlpha;
        }
        if( aBeta != 0.0 )
        {
           aY += mSwap;
        }
#endif
        // no base restore: nothing was converted
    }


    void
    SpMatrix::multiply( const Vector< real > & aX,
                        Vector< real > & aY )
    {
        // y = A * x is the extended product at alpha = 1, beta = 0, untransposed,
        // so this forwards rather than carrying a second implementation. It used
        // to switch on mType and call the Fortran kernels directly, with no
        // BELFEM_MKL branch at all - so an MKL build ran BELFEM's own matvec on
        // this path while the overload below ran MKL's, and the two drifted: only
        // the other one validates its dimensions. This path is the ARPACK inverse
        // iteration ( DofMgr_EigenValues ), i.e. once per iteration.
        //
        // The fill is REQUIRED, not defensive, and only on the MKL path. The
        // Fortran kernels zero y themselves ( splinalg.f90 ), so every caller of
        // this overload is entitled to hand us an UNINITIALISED aY, and the eigen
        // driver does exactly that ( a set_size() with no fill value, then straight
        // into here ). mkl_sparse_d_mv computes y := alpha * op( A ) * x + beta * y
        // and the Inspector-Executor sparse BLAS is not documented to short-circuit
        // beta == 0; 0.0 * y is 0 for finite garbage but NaN for a NaN or Inf bit
        // pattern, which uninitialised heap supplies. On the Fortran branch the
        // kernel's own zeroing already covers it, so filling there would be a
        // second pointless O( n ) pass on a hot path.
#ifdef BELFEM_MKL
        aY.fill( 0.0 );
#endif
        this->multiply( aX, aY, 1.0, 0.0, false );
    }

#if defined(__clang__) || defined(__INTEL_LLVM_COMPILER)
#pragma clang diagnostic pop
#elif defined(__GNUC__) || defined(__GNUG__)
#pragma GCC diagnostic pop
#endif

    /**
     * copy operator
     */
    SpMatrix &
    SpMatrix::operator=( const SpMatrix & aMatrix )
    {

        if( this == &aMatrix ) return *this;

        BELFEM_ERROR( mParent == nullptr, "can't copy onto a child matrix" );
        BELFEM_ERROR( mChild  == nullptr, "can't copy onto a parent matrix" );

        // delete current data
        this->deallocate();

        // copy the type
        mType = aMatrix.type();

        // set matrix size
        mNumRows = aMatrix.n_rows();
        mNumCols = aMatrix.n_cols();

        if( mType == SpMatrixType::CSR )
        {
            mPointerSize = mNumRows + 1;
        }
        else if ( mType == SpMatrixType::CSC )
        {
            mPointerSize = mNumCols + 1;
        }

        // set number of nonzeros
        mNumNonZeros = aMatrix.number_of_nonzeros();

        // copy pointers
        mPointers = ( int_t * ) malloc( ( mPointerSize ) * sizeof( int_t ) );
        std::memcpy( mPointers, aMatrix.pointers(), ( mPointerSize ) * sizeof( int_t ) );

        // copy cols
        if( aMatrix.rows() != nullptr )
        {
            mRows = ( int_t * ) malloc( mNumNonZeros * sizeof( int_t ) );
            std::memcpy( mRows, aMatrix.rows(), mNumNonZeros * sizeof( int_t ) );
        }
        else
        {
            mRows = nullptr;
        }

        // copy cols
        if( aMatrix.cols() != nullptr )
        {
            mColumns = ( int_t * ) malloc( mNumNonZeros * sizeof( int_t ) );
            std::memcpy( mColumns, aMatrix.cols(), mNumNonZeros * sizeof( int_t ) );
        }
        else
        {
            mColumns = nullptr;
        }

        // copy data ( allocate_values() adds the sentinel slot )
        this->allocate_values();
        std::memcpy( mValues, aMatrix.data(), mNumNonZeros * sizeof( real ) );

        // the coo array was deep-copied above, so the flag must follow
        mHaveCooIndices = aMatrix.mHaveCooIndices ;

#ifdef BELFEM_NETLIB
        this->allocate_swap();
#endif

        // link
        if( mPointers[ 0 ] == 0 )
        {
            this->set_indexing_base( SpMatrixIndexingBase::Cpp );
        }
        else
        {
            this->set_indexing_base( SpMatrixIndexingBase::Fortran );
        }

        // return ref to this matrix
        return *this;
    }
//----------------------------------------------------------------------------

    SpMatrix &
    SpMatrix::operator=( SpMatrix && aMatrix )
    {
        if ( this == &aMatrix )
        {
            return *this;
        }

        BELFEM_ERROR( mParent == nullptr,         "can't move onto a child matrix" );
        BELFEM_ERROR( mChild  == nullptr,         "can't move onto a parent matrix" );
        BELFEM_ERROR( aMatrix.mParent == nullptr, "can't move from a child matrix" );
        BELFEM_ERROR( aMatrix.mChild  == nullptr, "can't move from a parent matrix" );

        // delete current data
        this->deallocate();

        // move the type
        mType = aMatrix.mType ;

        // move matrix size
        mNumRows = aMatrix.mNumRows ;
        aMatrix.mNumRows = 0 ;

        mNumCols = aMatrix.mNumCols ;
        aMatrix.mNumCols = 0 ;

        // move pointer size
        mPointerSize = aMatrix.mPointerSize ;
        aMatrix.mPointerSize = 0 ;

        // move number of nonzeros
        mNumNonZeros = aMatrix.mNumNonZeros ;
        aMatrix.mNumNonZeros = 0 ;

        // move pointers
        mPointers = aMatrix.mPointers ;
        aMatrix.mPointers = nullptr ;

        // move rows
        mRows = aMatrix.mRows ;
        aMatrix.mRows = nullptr ;

        // move cols
        mColumns = aMatrix.mColumns ;
        aMatrix.mColumns = nullptr ;

        // move data
        mValues = aMatrix.mValues ;
        aMatrix.mValues = nullptr ;

        // move coo flag
        mHaveCooIndices = aMatrix.mHaveCooIndices ;
        aMatrix.mHaveCooIndices = false ;

#ifdef BELFEM_NETLIB
        // move swap
        mSwap = std::move( aMatrix.mSwap ) ;
        mSwapSize = aMatrix.mSwapSize ;
        aMatrix.mSwapSize = 0 ;
#endif

        // move function pointer

        // return ref to this matrix
        return *this ;
    }

    void
    SpMatrix::transpose()
    {
        BELFEM_ERROR( mParent == nullptr, "can't transpose a child matrix" );
        BELFEM_ERROR( mChild  == nullptr, "can't transpose a parent matrix" );

        int_t tNumRows = mNumCols;
        int_t tNumCols = mNumRows;

        if( mType == SpMatrixType::CSC )
        {
            // CSC -> CSR: row indices become column indices
            if( mColumns != nullptr )
            {
                free( mColumns );
            }
            mColumns = mRows;
            mRows = nullptr;
            mType = SpMatrixType::CSR;

            mNumRows = tNumRows;
            mNumCols = tNumCols;
            mPointerSize = mNumRows + 1;
        }
        else if( mType == SpMatrixType::CSR )
        {
            // CSR -> CSC: column indices become row indices
            if( mRows != nullptr )
            {
                free( mRows );
            }
            mRows = mColumns;
            mColumns = nullptr;
            mType = SpMatrixType::CSC;

            mNumRows = tNumRows;
            mNumCols = tNumCols;
            mPointerSize = mNumCols + 1;
        }

        // the branches above free the coo array ( it is the non-structural
        // one of mRows/mColumns ), so the flag must fall with it
        mHaveCooIndices = false;

        // same-base call: converts nothing; position() selects the index
        // array from mType on every call, so nothing needs rebinding
        if ( mPointers != nullptr )
        {
            this->set_indexing_base( this->indexing_base() == 1
                    ? SpMatrixIndexingBase::Fortran
                    : SpMatrixIndexingBase::Cpp );
        }
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::set_type( const SpMatrixType aType )
    {
        BELFEM_ERROR( mParent == nullptr && mChild == nullptr,
                      "can't change the type of a linked matrix" );
        mType = aType;
    }

//------------------------------------------------------------------------------

    size_t
    SpMatrix::memory() const
    {
        // values are always owned; the structure belongs to the parent
        size_t aMemory = mNumNonZeros * sizeof( real ) ;

        if ( mParent == nullptr )
        {
            aMemory += mPointerSize * sizeof( int_t );
            aMemory += mNumNonZeros * sizeof( int_t );

            if ( mHaveCooIndices )
            {
                aMemory += mNumNonZeros * sizeof( int_t );
            }
        }

        return aMemory;
    }

//------------------------------------------------------------------------------
}