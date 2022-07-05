//
// Created by Christian Messe on 2019-01-20.
//

#ifndef BELFEM_CL_SPMATRIX_HPP
#define BELFEM_CL_SPMATRIX_HPP

#include <cstring>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"
#include "cl_Vector.hpp"
#include "filetools.hpp"
#include "HDF5_Tools.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    enum class SpMatrixType
    {
        CSC        = 0,         // compressed sparse column
        CSR        = 1,         // compressed sparse row
        UNDEFINED  = 2
    };

//------------------------------------------------------------------------------

    enum class SpMatrixIndexingBase
    {
        Cpp        = 0,
        Fortran    = 1
    };

//------------------------------------------------------------------------------

    class SpMatrix
    {
//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

        // type of matrix, CSC or CSR
        SpMatrixType mType;

        // size of matrix
        int mNumRows = 0;
        int mNumCols = 0;

        // number of nonzeros
        int mNumNonZeros = 0;

        // container for pointers
        int * mPointers = nullptr;

        // size of pointer array
        int mPointerSize = 0;

        int * mRows = nullptr;

        int * mColumns = nullptr;

        // values array
        real * mValues = nullptr;

#ifdef BELFEM_NETLIB
        real * mSwap = nullptr;
        int mSwapSize = 0;
#endif

        // pointer with zero value
        real mZero = 0.0;

        // function that finds the position in the array
        index_t
        ( SpMatrix:: * mIndexFunction )(
                const index_t & aRowIndex,
                const index_t & aColIndex ) const;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        // empty constructor
        SpMatrix() = default;

//------------------------------------------------------------------------------

        // standard constructor using a graph
        SpMatrix( Cell<graph::Vertex *> & aGraph,
                  const enum SpMatrixType aType = SpMatrixType::CSC,
                          const index_t   aNumRows = 0,
                          const index_t   aNumCols = 0);

//------------------------------------------------------------------------------

        // file constructor using a HDF5 file
        SpMatrix( const string & aHDF5Path, const string   aLabel="Matrix" );

//------------------------------------------------------------------------------

        // constructor for testing purposes using dense Matrix
        SpMatrix( const Matrix< real > & aMatrix, const SpMatrixType aType = SpMatrixType::CSC );

//------------------------------------------------------------------------------

        ~SpMatrix();

//------------------------------------------------------------------------------
// Access to data containers
//------------------------------------------------------------------------------

        /**
         * return the data type
         */
        const SpMatrixType &
        type() const;

//------------------------------------------------------------------------------

        /**
         * write a specific value into all entries of the value container
         */
         void
         fill( const real & aValue );

//------------------------------------------------------------------------------

        /**
         * number of rows of this matrix
         */
        index_t
        n_rows() const;

//------------------------------------------------------------------------------

        /**
         * number of columns of this matrix
         */
        index_t
        n_cols() const;

//------------------------------------------------------------------------------

        /**
         * number of nonzero values in this matrix
         */
        index_t
        number_of_nonzeros() const;

//------------------------------------------------------------------------------

        /**
         * expose the pointers
         */
        int *
        pointers();

//------------------------------------------------------------------------------

        /**
         * expose the pointers (const version)
         */
        const int *
        pointers() const;

//------------------------------------------------------------------------------

        /**
         * expose the index array
         */
         int *
         indices();

//------------------------------------------------------------------------------

        /**
         * expose the index array ( const version )
         */
        const int *
        indices() const;

//------------------------------------------------------------------------------

        /**
         * expose the row indices
         */
        int *
        rows();

//------------------------------------------------------------------------------

        /**
         * expose the row indices ( const version )
         */
        const int *
        rows() const;

//------------------------------------------------------------------------------

        /**
         * expose the col indices
         */
        int *
        cols();

//------------------------------------------------------------------------------

        /**
         * expose the col indices ( const version )
         */
        const int *
        cols() const;

//------------------------------------------------------------------------------

        /**
         * expose the data container
         */
        real *
        data();

//------------------------------------------------------------------------------

        /**
         * expose the data container ( const version )
         */
        const real *
        data() const;

//------------------------------------------------------------------------------

        /**
         * expose single entry in the data container
         */
        real &
        data( const index_t aIndex );

//------------------------------------------------------------------------------

        /**
         *  expose single entry in the data containe ( const version )
         */
        const real &
        data( const index_t aIndex ) const;


//------------------------------------------------------------------------------

        /*
         * get the index of a specific row and col
         * */
        int
        index( const index_t & aRowIndex,
               const index_t & aColIndex ) const;

//------------------------------------------------------------------------------
// Utilities
//------------------------------------------------------------------------------

        /**
         * change the indexing base
         */
        void
        set_indexing_base( const enum SpMatrixIndexingBase & aBasis );

//------------------------------------------------------------------------------

        /**
         * create addidional indices that are needed by MUMPS
         */
        void
        create_coo_indices();

//------------------------------------------------------------------------------

        /**
        * delete addidional indices that are needed by MUMPS
        */
        void
        free_coo_indices();

//------------------------------------------------------------------------------

        /**
         * print the matrix to the screen ( for debugging )
         */
        void
        print( const string aLabel="SparseMatrix" );

//------------------------------------------------------------------------------

        /**
         * print the container indices on the screen ( for debugging )
         */
        void
        print2( const string aLabel="SparseMatrix" );

//------------------------------------------------------------------------------

        /**
         * returns the basis type of the matrix
         * 0: c++ indexing
         * 1: fortran indexing
         */
        int
        indexing_base() const;

//------------------------------------------------------------------------------

        /**
         * performs a matrix-vector multiplication
         *
         * @param aX
         * @param aY
         */
        void
        multiply( const Vector< real > & aX,
                        Vector< real > & aY,
                  const real aAlpha = 1.0,
                  const real aBeta  = 0.0,
                  const bool aTransposedFlag=false );

//------------------------------------------------------------------------------

        void
        transpose();

//------------------------------------------------------------------------------
// Saving and Loading
//------------------------------------------------------------------------------

        /**
         * save matrix to a hdf5 file
         *
         * @param aPath   path of hdf5 file
         * @param aLabel  title of this matrix
         * @param aMode   file mode
         */
        void
        save(   const string & aPath,
                const string   aLabel="Matrix",
                const enum FileMode aMode=FileMode::NEW );

//------------------------------------------------------------------------------

        /**
         * save matrix to a specific group in a HDF5 file
         */
        void
        save(   hid_t        & aGroup,
                herr_t       & aStatus );

//------------------------------------------------------------------------------

        /**
         * load matrix from a hdf5 file
         *
         * @param aPath   path of hdf5 file
         * @param aLabel  title of this matrix
         */
        void
        load(   const string & aPath,
                const string   aLabel="Matrix"
                );

//------------------------------------------------------------------------------

        /**
         * load matrix from a specific group in a hdf5 file
         */
        void
        load(   hid_t        & aGroup,
                herr_t       & aStatus );


//------------------------------------------------------------------------------
// Operators
//------------------------------------------------------------------------------

        /**
         * access a specific value with write access
         */
        real &
        operator()( const index_t & aRowIndex,
                    const index_t & aColIndex );

//------------------------------------------------------------------------------

        /**
         * access a specific value with read access
         */
        const real &
        operator()( const index_t & aRowIndex,
                    const index_t & aColIndex ) const;
//------------------------------------------------------------------------------

        /**
         * copy operator
         */
        SpMatrix &
        operator=( const SpMatrix & aMatrix );

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------
// Initialization
//------------------------------------------------------------------------------

        /**
         * make sure that all indices are unique ( debug mode only )
         */
         void
         check_graph( Cell<graph::Vertex *> & aGraph );

//------------------------------------------------------------------------------

        /**
         * tidy up the graph to make sure that data are sane
         */
        index_t
        tidy_graph( Cell<graph::Vertex *> & aGraph );

//------------------------------------------------------------------------------

        void
        create_csr_indices( Cell<graph::Vertex *> & aGraph );

//------------------------------------------------------------------------------

        void
        create_csc_indices( Cell<graph::Vertex *> & aGraph );

//------------------------------------------------------------------------------

        /**
         * throw error if matrices are too big
         */
        void
        set_sizes(
                const index_t aNumRows,
                const index_t aNumCols );

        void
        set_nnz( const index_t aNumberOfNonzeros );

//------------------------------------------------------------------------------

        /**
         * free containers
         */
        void
        deallocate();

//------------------------------------------------------------------------------

        /**
         * allocate the memory container
         */
        void
        allocate_values();

//------------------------------------------------------------------------------

#ifdef BELFEM_NETLIB
        /**
         * allocate swap vector
         */
        void
        allocate_swap();
#endif

//------------------------------------------------------------------------------
// Indexing
//------------------------------------------------------------------------------

        index_t
        index_csr_zero_based (
                const index_t & aRowIndex,
                const index_t & aColIndex ) const;

//------------------------------------------------------------------------------

        index_t
        index_csc_zero_based (
                const index_t & aRowIndex,
                const index_t & aColIndex ) const;

//------------------------------------------------------------------------------

        index_t
        index_csr_one_based (
                const index_t & aRowIndex,
                const index_t & aColIndex ) const;

//------------------------------------------------------------------------------

        index_t
        index_csc_one_based (
                const index_t & aRowIndex,
                const index_t & aColIndex ) const;

//------------------------------------------------------------------------------
    };


//------------------------------------------------------------------------------
// external operators
//------------------------------------------------------------------------------

    /**
     * multiply operator
     */
    Vector< real >
    inline operator * ( SpMatrix     & aA,
                        const Vector<real> & aX )
    {
        Vector<real> aY( aA.n_cols(), 0.0 );
        aA.multiply( aX, aY );
        return aY ;
    }

//------------------------------------------------------------------------------

    inline const SpMatrixType &
    SpMatrix::type() const
    {
        return mType;
    }

//------------------------------------------------------------------------------

    inline index_t
    SpMatrix::n_rows() const
    {
        return ( index_t ) mNumRows;
    }

//------------------------------------------------------------------------------

    inline index_t
    SpMatrix::n_cols() const
    {
        return ( index_t ) mNumCols;
    }

//------------------------------------------------------------------------------

    inline index_t
    SpMatrix::number_of_nonzeros() const
    {
        return ( index_t ) mNumNonZeros;
    }

//------------------------------------------------------------------------------

    inline int *
    SpMatrix::indices()
    {
        if( mType == SpMatrixType::CSC )
        {
            return mRows;
        }
        else if ( mType == SpMatrixType::CSR )
        {
            return mColumns;
        }
        else
        {
            return nullptr;
        }
    }

//------------------------------------------------------------------------------

    inline const int *
    SpMatrix::indices() const
    {
        if( mType == SpMatrixType::CSC )
        {
            return mRows;
        }
        else if ( mType == SpMatrixType::CSR )
        {
            return mColumns;
        }
        else
        {
            return nullptr;
        }
    }

//------------------------------------------------------------------------------

    inline int *
    SpMatrix::rows()
    {
        return mRows;
    }

//------------------------------------------------------------------------------

    inline const int *
    SpMatrix::rows() const
    {
        return mRows;
    }

//------------------------------------------------------------------------------

    inline int *
    SpMatrix::cols()
    {
        return mColumns;
    }

//------------------------------------------------------------------------------

    inline const int *
    SpMatrix::cols() const
    {
        return mColumns;
    }

//------------------------------------------------------------------------------

    inline int *
    SpMatrix::pointers()
    {
        return mPointers;
    }

//------------------------------------------------------------------------------

    const inline int *
    SpMatrix::pointers() const
    {
        return mPointers;
    }

//------------------------------------------------------------------------------

    inline real *
    SpMatrix::data()
    {
        return mValues;
    }

//------------------------------------------------------------------------------

    inline const real *
    SpMatrix::data() const
    {
        return mValues;
    }

//------------------------------------------------------------------------------

    inline real &
    SpMatrix::data( const index_t aIndex )
    {
        BELFEM_ASSERT( aIndex < ( index_t ) mNumNonZeros,
            "Index %lu for sparse matix out of bounds ( must be less than %lu )",
                      ( long unsigned int ) aIndex,
                      ( long unsigned int ) mNumNonZeros );

        return mValues[ aIndex ];
    }

//------------------------------------------------------------------------------

    inline const real &
    SpMatrix::data( const index_t aIndex ) const
    {
        BELFEM_ASSERT( aIndex < ( index_t ) mNumNonZeros,
                      "Index %lu for sparse matix out of bounds ( must be less than %lu )",
                      ( long unsigned int ) aIndex,
                      ( long unsigned int ) mNumNonZeros );

        return mValues[ aIndex ];
    }

//------------------------------------------------------------------------------

    inline int
    SpMatrix::indexing_base() const
    {
        return mPointers[ 0 ];
    }

//------------------------------------------------------------------------------

    inline int
    SpMatrix::index( const index_t & aRowIndex, const index_t & aColIndex ) const
    {
        return ( this->*mIndexFunction )( aRowIndex, aColIndex );
    }

//------------------------------------------------------------------------------

    /**
     * access a specific value with write access
     */
    inline real &
    SpMatrix::operator()( const index_t & aRowIndex,
                const index_t & aColIndex )
    {
        int tIndex = ( this->*mIndexFunction )( aRowIndex, aColIndex );

        BELFEM_ASSERT( tIndex < mNumNonZeros, "tried to access zero value in writable mode( %u, %u )",
                      ( unsigned int ) aRowIndex,
                      ( unsigned int ) aColIndex );

        return mValues[ tIndex ];
    }

//------------------------------------------------------------------------------

    inline const real &
    SpMatrix::operator()( const index_t & aRowIndex,
                const index_t & aColIndex ) const
    {
        int tIndex = ( this->*mIndexFunction )( aRowIndex, aColIndex );

        if( tIndex < mNumNonZeros )
        {
            return mValues[ tIndex ];
        }
        else
        {
            return mZero;
        }
    }

//------------------------------------------------------------------------------
}

#endif //BELFEM_CL_SPMATRIX_HPP
