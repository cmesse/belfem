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

#ifndef BELFEM_CL_SPMATRIX_HPP
#define BELFEM_CL_SPMATRIX_HPP

#include <cstring>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"
#include "cl_Vector.hpp"
#include "filetools.hpp"
#include "hdf5_tools.hpp"

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

    /**
     * @brief Sparse matrix in CSR or CSC format.
     *
     * @ingroup grp_sparse
     * @see @ref sparse_sparse_usage_guide
     */
    class SpMatrix
    {
//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

        SpMatrix * mParent = nullptr ;
        SpMatrix * mChild  = nullptr ;

        // type of matrix, CSC or CSR
        SpMatrixType mType;

        // size of matrix
        int_t mNumRows = 0;
        int_t mNumCols = 0;

        // number of nonzeros
        int_t mNumNonZeros = 0;

        // size of pointer array
        int_t mPointerSize = 0;

        // container for pointers
        int_t * mPointers = nullptr;

        int_t * mRows = nullptr;

        int_t * mColumns = nullptr;

        // values array
        real * mValues = nullptr;

        bool mHaveCooIndices = false;

#ifdef BELFEM_NETLIB
        Vector< real > mSwap ;
        int_t mSwapSize = 0;
#endif

        // pointer with zero value
        real mZero = 0.0;

        // Note: the lookup used to dispatch through a member-function
        // pointer over four near-identical search functions, one per
        // ( type, base ) pair. It is now the inline position() below, which
        // reads the base from mPointers[ 0 ] on every call. That is not only
        // faster - it removes a stale-state bug: a parent and its child share
        // the index arrays, and rebasing the parent did not update a child's
        // cached function pointer.

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        // empty constructor
        SpMatrix() = default;

        // raw-owning class: copy and move construction are not supported,
        // use the assignment operators on an empty matrix instead
        SpMatrix( const SpMatrix & ) = delete;
        SpMatrix( SpMatrix && ) = delete;

//------------------------------------------------------------------------------

        // standard constructor using a graph
        SpMatrix( Graph & aGraph,
                  const enum SpMatrixType aType = SpMatrixType::CSC,
                          const index_t   aNumRows = 0,
                          const index_t   aNumCols = 0,
                          const bool aSortGraph = true );

//------------------------------------------------------------------------------

        SpMatrix( SpMatrixType aType,
                  const index_t aNumRows,
                  const index_t aNumCols,
                  const index_t aNumNonZeros,
                  const int_t  * aIndices,
                  const int_t  * aPointers );

//------------------------------------------------------------------------------

        // file constructor using a HDF5 file
        SpMatrix( const string & aHDF5Path, const string   aLabel="Matrix" );

//------------------------------------------------------------------------------

        SpMatrix( const hid_t aParent, const string aLabel="Matrix" );

//------------------------------------------------------------------------------

        // constructor for testing purposes using dense Matrix
        SpMatrix( const Matrix< real > & aMatrix, const SpMatrixType aType = SpMatrixType::CSC );

//------------------------------------------------------------------------------

        // child constructor: shares the sparsity structure of the parent,
        // owns only its value array
        explicit SpMatrix( SpMatrix * aParent );

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
         fill( const real aValue );

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

        index_t
        n_pointers() const ;

//------------------------------------------------------------------------------

        size_t
        memory() const ;

//------------------------------------------------------------------------------

        bool
        have_coo_indices() const ;

//------------------------------------------------------------------------------

        /**
         * sets the type, must be called after set_sizes
         */
        void
        set_type( const SpMatrixType aType );

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
        int_t *
        pointers();

//------------------------------------------------------------------------------

        /**
         * expose the pointers (const version)
         */
        const int_t *
        pointers() const;

//------------------------------------------------------------------------------

        /**
         * expose the index array
         */
         int_t *
         indices();

//------------------------------------------------------------------------------

        /**
         * expose the index array ( const version )
         */
        const int_t *
        indices() const;

//------------------------------------------------------------------------------

        /**
         * expose the row indices
         */
        int_t *
        rows();

//------------------------------------------------------------------------------

        /**
         * expose the row indices ( const version )
         */
        const int_t *
        rows() const;

//------------------------------------------------------------------------------

        /**
         * expose the col indices
         */
        int_t *
        cols();

//------------------------------------------------------------------------------

        /**
         * expose the col indices ( const version )
         */
        const int_t *
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
         *
         * Returns the position into the value array, or number_of_nonzeros()
         * when the entry is not part of the sparsity pattern. Thin wrapper
         * around position().
         * */
        int_t
        index( const index_t aRowIndex,
               const index_t aColIndex ) const;

//------------------------------------------------------------------------------

        /**
         * Position of ( aRowIndex, aColIndex ) in the value array, or
         * mNumNonZeros when the entry is not in the sparsity pattern.
         *
         * Works in either indexing base and for both CSR and CSC; the base is
         * read from mPointers[ 0 ] and both branches are perfectly predicted
         * in any loop that does not alternate between matrices.
         */
        int_t
        position( const index_t aRowIndex,
                  const index_t aColIndex ) const;

//------------------------------------------------------------------------------

        /**
         * Batched lookup: positions of aNumCols entries of one slice
         * ( a row for CSR, a column for CSC ) in a single merge join over the
         * slice, instead of aNumCols independent searches.
         *
         * @param aSlice     row index (CSR) or column index (CSC)
         * @param aCols      the other index of each entry, STRICTLY ascending
         * @param aNumCols   number of entries
         * @param aPos       output, one position per entry; an entry absent
         *                   from the pattern gets mNumNonZeros
         *
         * aCols must be strictly ascending - within one element the dof
         * indices are unique, so equal neighbors cannot occur and the merge
         * join does not handle them. Asserted in debug builds.
         */
        void
        positions_in_slice( const index_t   aSlice,
                            const int_t   * aCols,
                            const uint      aNumCols,
                            int_t         * aPos ) const;

//------------------------------------------------------------------------------
// Utilities
//------------------------------------------------------------------------------

        /**
         * Ensure that the indices within each row (CSR) or column (CSC)
         * are sorted in ascending order. Unsorted entries are reordered
         * in place together with their values. This is required for
         * binary search in the index() methods.
         *
         * Called automatically when loading from HDF5 or constructing
         * from external arrays.
         */
        void
        sort_entries();

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
        * delete additional indices that are needed by MUMPS
        */
        void
        free_coo_indices();

//------------------------------------------------------------------------------

        /**
         * print_t the matrix to the screen ( for debugging )
         */
        void
        print( const string aLabel="SparseMatrix" );

//------------------------------------------------------------------------------

        /**
         * print_t the container indices on the screen ( for debugging )
         */
        void
        print2( const string aLabel="SparseMatrix" );

//------------------------------------------------------------------------------

        /**
         * returns the basis type of the matrix
         * 0: c++ indexing
         * 1: fortran indexing
         */
        int_t
        indexing_base() const;

//------------------------------------------------------------------------------

        /**
         * performs a matrix-vector multiplication
         *
         * c = alpha * A * b + beta * c
         * @param aX
         * @param aY
         * @param aAlpha           scaling factor of the product
         * @param aBeta            scaling factor of the accumulator aY
         * @param aTransposedFlag  multiply with the transpose of the matrix
         */
        void
        multiply( const Vector< real > & aX,
                        Vector< real > & aY,
                  const real aAlpha,
                  const real aBeta,
                  const bool aTransposedFlag=false );

//------------------------------------------------------------------------------

        /**
         * simple matrix-vector multiplication
         */
        void
        multiply( const Vector< real > & aX,
                        Vector< real > & aY );

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

        /**
         * move operator
         */
         SpMatrix &
         operator=( SpMatrix && aMatrix );

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------
// Initialization
//------------------------------------------------------------------------------

        /**
         * make sure that all indices are unique ( debug mode only )
         */
         void
         check_graph( Graph & aGraph );

//------------------------------------------------------------------------------

        /**
         * tidy up the graph to make sure that data are sane
         */
        index_t
        tidy_graph( Graph & aGraph );

//------------------------------------------------------------------------------

        void
        create_csr_indices( Graph & aGraph );

//------------------------------------------------------------------------------

        void
        create_csc_indices( Graph & aGraph );

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
// Parent/Child logic
//------------------------------------------------------------------------------

        void
        set_child( SpMatrix * aChild );

//------------------------------------------------------------------------------

        void
        update_from_parent();

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
        Vector<real> aY( aA.n_rows(), 0.0 );
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
    SpMatrix::n_pointers() const
    {
        return ( index_t ) mPointerSize;
    }

//------------------------------------------------------------------------------

    inline bool
    SpMatrix::have_coo_indices() const
    {
        return mHaveCooIndices;
    }

//------------------------------------------------------------------------------

    inline index_t
    SpMatrix::number_of_nonzeros() const
    {
        return ( index_t ) mNumNonZeros;
    }

//------------------------------------------------------------------------------

    inline int_t *
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

    inline const int_t *
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

    inline int_t *
    SpMatrix::rows()
    {
        return mRows;
    }

//------------------------------------------------------------------------------

    inline const int_t *
    SpMatrix::rows() const
    {
        return mRows;
    }

//------------------------------------------------------------------------------

    inline int_t *
    SpMatrix::cols()
    {
        return mColumns;
    }

//------------------------------------------------------------------------------

    inline const int_t *
    SpMatrix::cols() const
    {
        return mColumns;
    }

//------------------------------------------------------------------------------

    inline int_t *
    SpMatrix::pointers()
    {
        return mPointers;
    }

//------------------------------------------------------------------------------

    const inline int_t *
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

    inline int_t
    SpMatrix::indexing_base() const
    {
        return mPointers != nullptr ? mPointers[ 0 ] : 0;
    }

//------------------------------------------------------------------------------

    inline int_t
    SpMatrix::position( const index_t aRowIndex, const index_t aColIndex ) const
    {
        BELFEM_ASSERT( aRowIndex < ( index_t ) mNumRows, "aRowIndex out of bounds ( %lu >= %lu )",
                      ( long unsigned int ) aRowIndex ,
                      ( long unsigned int ) mNumRows );

        BELFEM_ASSERT( aColIndex < ( index_t ) mNumCols, "aColIndex out of bounds ( %lu >= %lu )",
                      ( long unsigned int ) aColIndex ,
                      ( long unsigned int ) mNumCols );

        // CSR searches a row for a column, CSC a column for a row. Both the
        // type and the base are loop-invariant in every caller, so these
        // branches predict perfectly and the whole function inlines.
        const int_t  tBase  = mPointers[ 0 ];
        const bool   tIsCsr = ( mType == SpMatrixType::CSR );

        const index_t tSlice  = tIsCsr ? aRowIndex : aColIndex ;
        const int_t   tTarget = ( int_t )( tIsCsr ? aColIndex : aRowIndex ) + tBase ;

        const int_t * tIndices = tIsCsr ? mColumns : mRows ;

        const int_t tBegin = mPointers[ tSlice ]     - tBase ;
        const int_t tEnd   = mPointers[ tSlice + 1 ] - tBase ;

        // Binary search, always. A linear scan with an early exit was tried
        // for short slices on the theory that FE rows are only tens of
        // entries long. On the machine it was measured on it lost at every
        // row length ( 8, 16, 24, 32, 50 ): a slice sits in one or two cache
        // lines, so lower_bound costs ~6 well-predicted steps while the scan
        // averages k/2 iterations behind a data-dependent branch.
        //
        // A second reviewer running the same experiment on different hardware
        // did see linear win at 8 entries per row, so the crossover is
        // microarchitecture-dependent and this is not a universal result. The
        // threshold knob was removed rather than guessed at; if it is ever
        // revisited, re-measure on the production nodes rather than trusting
        // either of those numbers. Both agree the difference is small next to
        // what the assembly restructure targets.
        // See devlog/dl20260818_spmatrix_accessor.md.
        //
        // ( indices within a slice are sorted, see sort_entries )
        const int_t * tFound = std::lower_bound( tIndices + tBegin,
                                                 tIndices + tEnd,
                                                 tTarget );

        return ( tFound < tIndices + tEnd && *tFound == tTarget )
               ? ( int_t )( tFound - tIndices )
               : mNumNonZeros ;
    }

//------------------------------------------------------------------------------

    inline int_t
    SpMatrix::index( const index_t aRowIndex, const index_t aColIndex ) const
    {
        return this->position( aRowIndex, aColIndex );
    }

//------------------------------------------------------------------------------

    /**
     * access a specific value with write access
     */
    inline real &
    SpMatrix::operator()( const index_t & aRowIndex,
                const index_t & aColIndex )
    {
        BELFEM_ASSERT( this->indexing_base() == 0,
            "operator() called while matrix is in Fortran indexing mode. "
            "Call set_indexing_base( SpMatrixIndexingBase::Cpp ) first." );

        const int_t tIndex = this->position( aRowIndex, aColIndex );

        BELFEM_ASSERT( tIndex < mNumNonZeros, "tried to access zero value in writable mode( %lu, %lu )",
                      ( long unsigned int ) aRowIndex,
                      ( long unsigned int ) aColIndex );

        return mValues[ tIndex ];
    }

//------------------------------------------------------------------------------

    inline const real &
    SpMatrix::operator()( const index_t & aRowIndex,
                const index_t & aColIndex ) const
    {
        BELFEM_ASSERT( this->indexing_base() == 0,
            "operator() called while matrix is in Fortran indexing mode. "
            "Call set_indexing_base( SpMatrixIndexingBase::Cpp ) first." );

        const int_t tIndex = this->position( aRowIndex, aColIndex );

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
