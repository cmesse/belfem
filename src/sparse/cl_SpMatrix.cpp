//
// Created by Christian Messe on 2019-01-20.
//

#include <cstring>    // for std::memcpy
#include <algorithm>  // for std::sort std::fill_n

#include "assert.hpp"
#include "cl_SpMatrix.hpp"
#include "fn_max.hpp"
#include "fn_unique.hpp"
#include "cl_Timer.hpp"
#include "cl_Logger.hpp"

#include "fspblas.hpp"
#include "HDF5_Tools.hpp"
#include "cl_HDF5.hpp"
#include "stringtools.hpp"
#include "fn_Graph_sort.hpp"

#ifdef BELFEM_MPI
#include "commtools.hpp"
#endif

namespace belfem
{

//------------------------------------------------------------------------------

    SpMatrix::SpMatrix( Cell<graph::Vertex *> & aGraph,
                        const enum SpMatrixType aType,
                        const index_t aNumRows,
                        const index_t aNumCols) :
            mType( aType )
    {
        belfem::graph::sort( aGraph );

#ifdef BELFEM_DEBUG
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
            for( graph::Vertex * tVertex : aGraph )
            {
                tVertex->sort_vertices();
            }

            // take sizes from input
            this->set_sizes( aNumRows, aNumCols );
        }

        // chech which kind of matrix this is
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

    SpMatrix::SpMatrix( const string & aHDF5Path, const string aLabel )
    {
        this->load( aHDF5Path, aLabel );
    }

//------------------------------------------------------------------------------

    SpMatrix::SpMatrix( const Matrix< real > & aMatrix, const SpMatrixType aType )
    {
        // set sizes
        this->set_sizes( aMatrix.n_rows(), aMatrix.n_cols() );

        // count nonzeros
        mNumNonZeros = 0 ;
        for( int i=0; i<mNumRows; ++i )
        {
            for( int j=0; j<mNumCols; ++j )
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
                mPointers = ( int * ) malloc( ( mPointerSize ) * sizeof( int ) );

                // first entry
                mPointers[ 0 ] = 0 ;

                // allocate row array
                mColumns = ( int * ) malloc( ( mNumNonZeros ) * sizeof( int ) );

                // position in columns
                int tStep = 0 ;

                int tCount;

                // populate pointers array
                for( int i=0; i<mNumRows; ++i )
                {
                    // reset counter
                    tCount = 0 ;

                    // loop over all columns
                    for( int j=0; j<mNumCols; ++j )
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
                mPointers = ( int * ) malloc( ( mPointerSize ) * sizeof( int ) );

                // first entry
                mPointers[ 0 ] = 0 ;

                // allocate row array
                mRows = ( int * ) malloc( ( mNumNonZeros ) * sizeof( int ) );

                // position in columns
                int tStep = 0 ;

                int tCount;

                // populate pointers array
                for( int j=0; j<mNumCols; ++j )
                {
                    // reset counter
                    tCount = 0 ;

                    // loop over all columns
                    for( int i=0; i<mNumRows; ++i )
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

    }

//------------------------------------------------------------------------------

    SpMatrix::~SpMatrix()
    {
        this->deallocate();
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::deallocate()
    {
        if( mPointers != NULL )
        {
            free( mPointers );
            mPointers = nullptr;
            mPointerSize = 0;
        }
        if( mRows != NULL )
        {
            free( mRows );
            mRows = nullptr;
            mNumRows = 0;
        }
        if( mColumns != NULL )
        {
            free( mColumns );
            mColumns = nullptr;
            mNumCols = 0;
        }
        if( mValues != NULL )
        {
            free( mValues );
            mValues = nullptr;
            mNumNonZeros = 0;
        }

#ifdef BELFEM_NETLIB
        if( mSwap != NULL )
        {
            free( mSwap );
            mSwap = nullptr;
            mSwapSize = 0;
        }
#endif

    }

//------------------------------------------------------------------------------

    void
    SpMatrix::check_graph( Cell< graph::Vertex * > & aGraph )
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
    SpMatrix::tidy_graph( Cell< graph::Vertex * > & aGraph )
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
    SpMatrix::create_csr_indices( Cell< graph::Vertex * > & aGraph )
    {
        // number of vertices
        index_t tSize = aGraph.size();

        mPointerSize = mNumRows + 1;

        // allocate pointer array
        mPointers = ( int * ) malloc( ( mPointerSize ) * sizeof( int ) );

        // populate pointer array
        std::fill_n( mPointers, mPointerSize, 0 );

        // create pointer array ( step 1 )
        for ( index_t k = 0; k < tSize; ++k )
        {
            mPointers[ aGraph( k )->index() + 1 ] = ( int ) aGraph( k )->number_of_vertices();
        }

        // counter to prevent data type overflow
        index_t tCount = 0;

        // create pointer array ( step 2 )
        for ( int k = 1; k < mPointerSize; ++k )
        {
            tCount += mPointers[ k ];
            mPointers[ k ] += mPointers[ k - 1 ];
        }

        // set number of nonzeros and check int type boundaries
        this->set_nnz( tCount );

        BELFEM_ASSERT( ( index_t ) mPointers[ mNumRows ] == tCount,
            "Something went wrong while creating CSR index" );

        // allocate index vector
        tCount = ( tCount == 0 ) ? 1 : tCount;

        mColumns = ( int * ) malloc( tCount * sizeof( int ) );

        // reset counter
        tCount = 0;

        // loop over all nodes and create index array
        for( graph::Vertex * tVertex : aGraph )
        {
            index_t tNumVertices = tVertex->number_of_vertices();

            for( uint k=0; k<tNumVertices; ++k )
            {
                mColumns[ tCount++ ] = tVertex->vertex( k )->index();
            }
        }
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::create_csc_indices(  Cell< graph::Vertex * > & aGraph )
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
            uint tNumberOfVertices =  aGraph( k )->number_of_vertices();

            for( uint i=0; i<tNumberOfVertices; ++i )
            {
                ++tCount[ aGraph( k )->vertex( i )->index() ];
            }

            tNumNonzeros +=tNumberOfVertices;
        }

        mPointerSize = mNumCols + 1;

        // populate pointer array
        mPointers = ( int * ) malloc( ( mPointerSize ) * sizeof( int ) );

        mPointers[ 0 ] = 0;
        for( int k=0; k<mNumCols; ++k )
        {
            mPointers[ k+1 ] = mPointers[ k ] + tCount[ k ];
        }

        // reset counter
        std::fill_n( tCount, mNumCols, 0 );

        // set number of nonzeros and check int type boundaries
        this->set_nnz( tNumNonzeros );

        // populate indices
        tNumNonzeros = ( tNumNonzeros == 0 ) ? 1 : tNumNonzeros;

        mRows = ( int * ) malloc( tNumNonzeros * sizeof( int ) );

        int n = ( int ) tSize;

        // count vertices per column
        for ( int k = 0; k < n; ++k )
        {
            uint tNumberOfVertices =  aGraph( k )->number_of_vertices();

            for( uint i=0; i<tNumberOfVertices; ++i )
            {
                // get column of array
                index_t j = aGraph( k )->vertex( i )->index();

                // write index into array
                mRows[ mPointers[ j ] + tCount[ j ] ] = ( int ) aGraph( k )->index();

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
        // make sure that size is OK
        BELFEM_ERROR( aNumRows + 1 < ( index_t ) BELFEM_INT_MAX,
                "value overflow in number of rows of matrix (%lu + 1 > %lu )",
                     ( long unsigned int ) aNumRows,
                     ( long unsigned int ) BELFEM_INT_MAX );

        BELFEM_ERROR( aNumCols + 1 < ( index_t ) BELFEM_INT_MAX,
                     "value overflow in numner cols of matrix (%lu + 1 > %lu )",
                     ( long unsigned int ) aNumCols,
                     ( long unsigned int ) BELFEM_INT_MAX );

        // set data
        mNumRows = ( int ) aNumRows;
        mNumCols = ( int ) aNumCols;
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::set_nnz( const index_t aNumNonZeros )
    {
        // make sure that NNZ is OK
        BELFEM_ERROR( aNumNonZeros < ( index_t )  BELFEM_INT_MAX,
                     "too many non-zeros in matrix (%lu > %lu )",
                     ( long unsigned int ) aNumNonZeros,
                     ( long unsigned int ) BELFEM_INT_MAX );

        // set data
        mNumNonZeros = ( int ) aNumNonZeros;
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::allocate_values()
    {
        mValues = ( real * ) malloc( mNumNonZeros * sizeof( real ) );
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
        mSwap = ( real * ) malloc( mSwapSize * sizeof( real ) );
    }
#endif

//------------------------------------------------------------------------------

    void
    SpMatrix::fill( const real & aValue )
    {
        std::fill_n( mValues, mNumNonZeros, aValue );
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::set_indexing_base( const enum SpMatrixIndexingBase & aBasis )
    {
        switch( aBasis )
        {
            case( SpMatrixIndexingBase::Cpp ) :
            {
                // test if this is in fortran base
                if( mPointers[ 0 ] == 1 )
                {
                    // decrement all pointers
                    std::for_each( mPointers, mPointers + mPointerSize,
                                   [ ]( int & tValue ){ --tValue; } );

                    // decrement all rows
                    if( mRows != nullptr )
                    {
                        std::for_each( mRows, mRows + mNumNonZeros,
                                       [ ]( int & tValue ){ --tValue; } );
                    }

                    // decrement all cols
                    if( mColumns != nullptr )
                    {
                        std::for_each( mColumns, mColumns + mNumNonZeros,
                                       [ ]( int & tValue ){ --tValue; } );
                    }

                }

                // set index function
                switch( mType )
                {
                    case( SpMatrixType::CSR ) :
                    {
                        mIndexFunction = & SpMatrix::index_csr_zero_based;
                        break;
                    }
                    case( SpMatrixType::CSC ) :
                    {
                        mIndexFunction = & SpMatrix::index_csc_zero_based;
                        break;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Unknown matrix type");
                        break;
                    }
                }
                break;
            }
            case( SpMatrixIndexingBase::Fortran ) :
            {
                // test if this is in c++ base
                if( mPointers[ 0 ] == 0 )
                {
                    // decrement all pointers
                    std::for_each( mPointers, mPointers + mPointerSize,
                                   [ ]( int & tValue ){ ++tValue; } );

                    // increment all rows
                    if( mRows != nullptr )
                    {
                        std::for_each( mRows, mRows + mNumNonZeros,
                                       [ ]( int & tValue ){ ++tValue; } );
                    }

                    // increment all cols
                    if( mColumns != nullptr )
                    {
                        std::for_each( mColumns, mColumns + mNumNonZeros,
                                       [ ]( int & tValue ){ ++tValue; } );
                    }

                }
                // set index function
                switch( mType )
                {
                    case( SpMatrixType::CSR ) :
                    {
                        mIndexFunction = & SpMatrix::index_csr_one_based;
                        break;
                    }
                    case( SpMatrixType::CSC ) :
                    {
                        mIndexFunction = & SpMatrix::index_csc_one_based;
                        break;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Unknown matrix type");
                        break;
                    }
                }
                break;
            }
            default:
            {
                BELFEM_ERROR( false, "Unknown SpMatrixIndexingBase" );
                break;
            }
        }
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::create_coo_indices()
    {
        // make sure that we use Cpp indexing
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
                        mRows = ( int * ) malloc( 1 * sizeof( int ));
                    }
                   else
                    {
                        mRows = ( int * ) malloc( mNumNonZeros * sizeof( int ));
                    }

                    // populate row indices
                    int tCount = 0;
                    for ( int k = 0; k < mNumRows; ++k )
                    {
                        int tN = mPointers[ k + 1 ] - mPointers[ k ];
                        for ( int i = 0; i < tN; ++i )
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
                        mColumns = ( int * ) malloc( 1 * sizeof( int ));
                    }
                    else
                    {
                        mColumns = ( int * ) malloc( mNumNonZeros * sizeof( int ));
                    }
                    // populate row indices
                    int tCount = 0;
                    for ( int k = 0; k < mNumCols; ++k )
                    {
                        int tN = mPointers[ k + 1 ] - mPointers[ k ];
                        for ( int i = 0; i < tN; ++i )
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
    }

//------------------------------------------------------------------------------

    void
    SpMatrix::free_coo_indices()
    {
        switch( mType )
        {
            case( SpMatrixType::CSR ) :
            {
                if ( mRows != NULL )
                {
                    free( mRows );
                }
                break;
            }
            case(  SpMatrixType::CSC ) :
            {
                if ( mColumns != NULL )
                {
                    free( mColumns );
                }
                break;
            }
            default :
            {
                break;
            }
        }
    }

//------------------------------------------------------------------------------

    /**
     * print the matrix to the screen ( for debugging )
     */
    void
    SpMatrix::print( const string aLabel  )
    {
        int tCount = 0;

        std::cout<< "SpMatrix " << aLabel << " (" << mNumRows << ", " << mNumCols << ") " << std::endl;

        if ( mType == SpMatrixType::CSR )
        {
            for ( int r = 0; r < mNumRows; ++r )
            {
                // get bandwidth
                int n = mPointers[ r + 1 ] - mPointers[ r ];

                for ( int c = 0; c < n; ++c )
                {
                    std::cout << " (" << r << ", " << mColumns[ tCount ] << ")  : " << mValues[ tCount ] << " [ " << tCount << " ]" << std::endl;
                    ++tCount;
                }
            }
        }
        else if ( mType == SpMatrixType::CSC )
        {
            for ( int r = 0; r < mNumCols; ++r )
            {
                // get bandwidth
                int n = mPointers[ r + 1 ] - mPointers[ r ];

                for ( int c = 0; c < n; ++c )
                {
                    std::cout << " (" << r << ", " << mRows[ tCount ] << ")  : " << mValues[ tCount ] << " [ " << tCount << " ]" << std::endl;
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

        for( int k=0; k<mPointerSize; ++k )
        {
            std::cout << "P ( " << k << " ) = " <<mPointers[ k ] << std::endl;
        }

        if( mType == SpMatrixType::CSC )
        {
            for ( int k = 0; k < mNumNonZeros; ++k )
            {
                std::cout << "R ( " << k << " ) = " << mRows[ k ] << std::endl;
            }
        }
        else if( mType == SpMatrixType::CSR )
        {
            for ( int k = 0; k < mNumNonZeros; ++k )
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
        this->deallocate();

        string tFormatLabel;


        // load format
        hdf5::load_string_from_file( aGroup, "Format", tFormatLabel, aStatus );

        string tFormat = string_to_upper( tFormatLabel );

        // set type and link indices
        if( tFormat == "CSR" )
        {
            mType = SpMatrixType::CSR;
        }
        else if ( tFormat == "CSC" )
        {
            mType = SpMatrixType::CSC ;
        }

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
        if( mValues != nullptr )
        {
            free ( mValues );
        }

        // allocate memory
        mPointers = ( int * ) malloc( mPointerSize * sizeof( int ) );

        hdf5::load_array_from_file( aGroup, "Pointers", mPointers,  mPointerSize, aStatus );


        // load indices
        if( mType == SpMatrixType::CSC )
        {
            mRows = ( int * ) malloc( mNumNonZeros * sizeof( int ) );
            hdf5::load_array_from_file( aGroup, "Indices", mRows, mNumNonZeros, aStatus );

        }
        else if ( mType == SpMatrixType::CSR )
        {
            mColumns = ( int * ) malloc( mNumNonZeros * sizeof( int ) );
            hdf5::load_array_from_file( aGroup, "Indices", mColumns, mNumNonZeros, aStatus );
        }

        mValues = ( real * ) malloc( mNumNonZeros * sizeof( real ) );

        // load values
        hdf5::load_array_from_file( aGroup, "Values", mValues, mNumNonZeros, aStatus );

#ifdef BELFEM_NETLIB
        this->allocate_swap();
#endif
#else
        BELFEM_ERROR( false, "Trying to load a matrix from HDF5, but BELFEM is not link against HDF5 libraries." );
#endif
    }

//------------------------------------------------------------------------------

    index_t
    SpMatrix::index_csr_zero_based(
            const index_t & aRowIndex,
            const index_t & aColIndex ) const
    {
        BELFEM_ASSERT( aRowIndex < ( index_t ) mNumRows, "aRowIndex out of bounds ( %lu >= %lu )",
                      ( long unsigned int ) aRowIndex ,
                      ( long unsigned int ) mNumRows );

        BELFEM_ASSERT( aColIndex < ( index_t ) mNumCols, "aColIndex out of bounds ( %lu >= %lu )",
                      ( long unsigned int ) aColIndex ,
                      ( long unsigned int ) mNumCols );

        auto tBegin = mColumns + mPointers[ aRowIndex ];

        auto tEnd   = mColumns + mPointers[ aRowIndex+1 ];

        // find position in memory
        auto tFound = std::find( tBegin, tEnd, aColIndex );

        if( tFound < tEnd )
        {
            return ( tFound - tBegin ) + mPointers[ aRowIndex ];
        }
        else
        {
            return mNumNonZeros;
        }
    }

//----------------------------------------------------------------------------

    index_t
    SpMatrix::index_csc_zero_based(
            const index_t & aRowIndex,
            const index_t & aColIndex ) const
    {
        BELFEM_ASSERT( aRowIndex < ( index_t ) mNumRows, "aRowIndex out of bounds" );
        BELFEM_ASSERT( aColIndex < ( index_t ) mNumCols, "aColIndex out of bounds" );

        auto tBegin = mRows + mPointers[ aColIndex ];

        auto tEnd   = mRows + mPointers[ aColIndex+1 ];

        // find position in memory
        auto tFound = std::find( tBegin, tEnd, aRowIndex );

        if( tFound < tEnd )
        {
            return ( tFound - tBegin ) + mPointers[ aColIndex ];
        }
        else
        {
            return mNumNonZeros;
        }
    }

//----------------------------------------------------------------------------

    index_t
    SpMatrix::index_csr_one_based(
            const index_t & aRowIndex,
            const index_t & aColIndex ) const
    {
        BELFEM_ASSERT( aRowIndex < ( index_t ) mNumRows, "aRowIndex out of bounds" );
        BELFEM_ASSERT( aColIndex < ( index_t ) mNumCols, "aColIndex out of bounds" );

        auto tBegin = mColumns + mPointers[ aRowIndex ] - 1;

        auto tEnd   = mColumns + mPointers[ aRowIndex+1 ] - 1;

        // find position in memory
        auto tFound = std::find( tBegin, tEnd, aColIndex+1 );

        if( tFound < tEnd )
        {
            return ( tFound - tBegin ) + mPointers[ aRowIndex ] - 1;
        }
        else
        {
            return mNumNonZeros;
        }
    }

//----------------------------------------------------------------------------

    index_t
    SpMatrix::index_csc_one_based(
            const index_t & aRowIndex,
            const index_t & aColIndex ) const
    {
        BELFEM_ASSERT( aRowIndex < ( index_t ) mNumRows, "aRowIndex out of bounds" );
        BELFEM_ASSERT( aColIndex < ( index_t ) mNumCols, "aColIndex out of bounds" );

        auto tBegin = mRows + mPointers[ aColIndex ] - 1;

        auto tEnd   = mRows + mPointers[ aColIndex+1 ] - 1;

        // find position in memory
        auto tFound = std::find( tBegin, tEnd, aRowIndex+1 );


        if( tFound < tEnd )
        {
            return ( tFound-tBegin ) + mPointers[ aColIndex ] - 1;
        }
        else
        {
            return mNumNonZeros;
        }

    }

//----------------------------------------------------------------------------

    void
    SpMatrix::multiply(
            const Vector< real > & aX,
                  Vector< real > & aY,
            const real aAlpha,
            const real aBeta,
            const bool aTransposedFlag )
    {
        int tLengthX = mNumCols;
        int tLengthY = mNumRows;

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

        int  tOne   = 1;

        // Bugfix see Issue #3
        this->set_indexing_base( SpMatrixIndexingBase::Fortran );

#ifdef BELFEM_NETLIB
        int  tTransposedFlag  = ( aTransposedFlag ) ? 1 : 0 ;

        int tParameters[ 5 ] = { 0, 0, 0, 0, 1 };

        // indexing base, zero or one-based
        tParameters[ 3 ] = mPointers[ 0 ] ;

        switch( mType )
        {
            case( SpMatrixType::CSC ) :
            {
                // call fspblas
                dcscmm_(
                        &tTransposedFlag,
                        &mNumRows,
                        &tOne,
                        &mNumCols,
                        &aAlpha,
                        tParameters,
                        mValues,
                        mRows,
                        mPointers,
                        mPointers+1,
                        aX.data(),
                        &tLengthX,
                        &aBeta,
                        aY.data(),
                        &tLengthY,
                        mSwap,
                        &mSwapSize
                );
                break;
            }
            case( SpMatrixType::CSR ) :
            {
                // call fspblas
                dcsrmm_(
                        &tTransposedFlag,
                        &mNumRows,
                        &tOne,
                        &mNumCols,
                        &aAlpha,
                        tParameters,
                        mValues,
                        mColumns,
                        mPointers,
                        mPointers+1,
                        aX.data(),
                        &tLengthX,
                        &aBeta,
                        aY.data(),
                        &tLengthY,
                        mSwap,
                        &mSwapSize
                );
                break;
            }
            default :
            {
                BELFEM_ERROR( false, "Unsupported matrix type.");
                break;
            }
        }
#elif BELFEM_MKL
        char tParameters[] = { 'G', ' ', 'N', 'C' };

        if(  mPointers[ 0 ] == 1 )
        {
            tParameters[ 3 ] = 'F';
        }

        char tTransposedFlag = ( aTransposedFlag ) ? 'Y' : 'N' ;

        switch( mType )
        {
            case( SpMatrixType::CSC ) :
            {
                // call mkl
                mkl_dcscmm_(
                        &tTransposedFlag,
                        &mNumRows,
                        &tOne,
                        &mNumCols,
                        &aAlpha,
                        tParameters,
                        mValues,
                        mRows,
                        mPointers,
                        mPointers+1,
                        aX.data(),
                        &tLengthX,
                        &aBeta,
                        aY.data(),
                        &tLengthY
                );
                break;
            }
            case( SpMatrixType::CSR ) :
            {
                // call mkl
                mkl_dcsrmm_(
                        &tTransposedFlag,
                        &mNumRows,
                        &tOne,
                        &mNumCols,
                        &aAlpha,
                        tParameters,
                        mValues,
                        mColumns,
                        mPointers,
                        mPointers+1,
                        aX.data(),
                        &tLengthX,
                        &aBeta,
                        aY.data(),
                        &tLengthY
                );
                break;
            }
            default :
            {
                BELFEM_ERROR( false, "Unsupported matrix type.");
            }
        }
#endif
    }


    /**
     * copy operator
     */
    SpMatrix &
    SpMatrix::operator=( const SpMatrix & aMatrix )
    {
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
        mPointers = ( int * ) malloc( ( mPointerSize ) * sizeof( int ) );
        std::memcpy( mPointers, aMatrix.pointers(), ( mPointerSize ) * sizeof( int ) );

        // copy cols
        if( aMatrix.rows() != NULL )
        {
            mRows = ( int * ) malloc( mNumNonZeros * sizeof( int ) );
            std::memcpy( mRows, aMatrix.rows(), mNumNonZeros * sizeof( int ) );
        }
        else
        {
            mRows = nullptr;
        }

        // copy cols
        if( aMatrix.cols() != NULL )
        {
            mColumns = ( int * ) malloc( mNumNonZeros * sizeof( int ) );
            std::memcpy( mColumns, aMatrix.cols(), mNumNonZeros * sizeof( int ) );
        }
        else
        {
            mColumns = nullptr;
        }

        // copy data
        mValues = ( real * ) malloc( mNumNonZeros * sizeof( real ) ) ;
        std::memcpy( mValues, aMatrix.data(), mNumNonZeros * sizeof( real ) );

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

    void
    SpMatrix::transpose()
    {
        int tNumRows = mNumCols;
        int tNumCols = mNumRows;

        if( mType == SpMatrixType::CSC )
        {
            if( mRows != nullptr )
            {
                free( mRows );
            }
            mRows = std::move( mColumns );
            mColumns = nullptr;
            mType = SpMatrixType::CSR;

            mNumRows = tNumRows;
            mNumCols = tNumCols;
        }
        else if( mType == SpMatrixType::CSR )
        {
            if( mColumns != nullptr )
            {
                free( mColumns );
            }
            mColumns = std::move( mRows );
            mRows = nullptr;
            mType = SpMatrixType::CSC;

            mNumRows = tNumRows;
            mNumCols = tNumCols;
        }
    }

//------------------------------------------------------------------------------
}