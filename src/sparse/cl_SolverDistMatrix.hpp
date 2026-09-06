/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_CL_SOLVERDISTMATRIX_HPP
#define BELFEM_CL_SOLVERDISTMATRIX_HPP

#include <set>        // for std::set (garray unique collection)

#include "typedefs.hpp"
#include "commtools.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_SpMatrix.hpp"
#include "cl_OrderedMap.hpp"
#include "cl_SolverParameters.hpp"

namespace belfem
{
    namespace sparse
    {
        class DistMatrix
        {
        protected:

            const proc_t     mCommRank ;
            const proc_t     mCommSize ;
            const bool       mPermutationSwitch ;

            //! root's deck value, broadcast: the parallel orderings are
            //! collective, so every rank must agree on the branch
            const ReorderingMethod mReorderingMethod ;

            SpMatrix * mMatrix = nullptr ;

            Cell< index_t > mForwardPermutation ;
            Cell< index_t > mBackwardPermutation ;
            Cell< index_t > mIndexPermutation ;

            Vector< int_t > mSizes ;
            Vector< int_t > mDist ;
            Vector< int_t > mOffsets ;

            Vector< real > mRhs ;
            Vector< real > mMyRhs ;
            Vector< real > mLhs ;
            Vector< real > mMyLhs ;

            index_t mNumRows ;
            index_t mMyNumRows ;

        public:

            DistMatrix( const SolverParameters * aParams, SpMatrix * aMatrix = nullptr ) ;

            virtual ~DistMatrix();

            virtual void
            distribute_values( SpMatrix * aMatrix ) = 0 ;

            const int_t *
            dist() const
            {
                return mDist.data();
            }

            const real *
            rhs() const
            {
                return mMyRhs.data();
            }

            real *
            lhs()
            {
                return mMyLhs.data();
            }

            index_t
            size() const
            {
                return mNumRows;
            }

            Vector< real > &
            lhs_vector()
            {
                return mMyLhs;
            }

            Vector< real > &
            rhs_vector()
            {
                return mMyRhs;
            }

            void
            distribute_rhs( const Vector< real > & aRhs ) ;

            // as initial guess
            void
            distribute_lhs( const Vector< real > & aLhs ) ;

            void
            collect_lhs( Vector< real > & aLhs ) ;

        protected:

            void
            determine_sizes_and_offsets() ;

            virtual void
            distribute_sparsity_pattern() = 0 ;

        private:

            bool
            set_permutation_switch( const SolverParameters * aParams );

            ReorderingMethod
            set_reordering_method( const SolverParameters * aParams );

            void
            create_matrix( SpMatrix * aMatrix );

            //! collective: every rank enters it. parmetis / ptscotch are
            //! MPI-collective calls, metis runs on root ( a no-op on the
            //! empty non-root graph )
            void
            order_graph( Graph & aGraph );

            void
            compute_index_permutation( const SpMatrix * aMatrix );


        };

        template< typename T >
        class DistMatrixCSR : public DistMatrix
        {
            Cell< Vector< T > > mAllPointers ;
            Cell< Vector< T > > mAllIndices ;

            Vector< T >    mMyPointers ;
            Vector< T >    mMyIndices ;
            Vector< real > mMyValues ;

            // rank 0 only: points into the (permuted) matrix data set by
            // distribute_values(); the solver copies on set/update, so the
            // pointer only has to outlive that call
            const real * mMyValuesData = nullptr ;

        public:

            DistMatrixCSR( const SolverParameters * aParams, SpMatrix * aMatrix ) :
                DistMatrix( aParams, aMatrix )
            {
                this->distribute_sparsity_pattern();
            }

            void
            distribute_values( SpMatrix * aMatrix ) override
            {
                if ( mCommRank == 0 )
                {
                    aMatrix->set_indexing_base( SpMatrixIndexingBase::Cpp );

                    real * tData = mPermutationSwitch ? mMatrix->data() : aMatrix->data();

                    // reorganize data
                    if ( mPermutationSwitch )
                    {
                        real * tValues = aMatrix->data();

                        index_t tCount = 0 ;
                        for ( index_t k : mIndexPermutation )
                        {
                            tData[ tCount++ ] = tValues[ k ] ;
                        }
                    }

                    // per-rank value slices are contiguous ranges
                    // [ ptr[dist(p)], ptr[dist(p+1)] ) = mOffsets ranges,
                    // so we scatter straight from the matrix data
                    mMyValuesData = tData ;

                    comm_barrier();
                    distribute( tData, mOffsets );
                    comm_barrier();
                }
                else
                {
                    comm_barrier();
                    receive( mMyValues );
                    comm_barrier();
                }
            }

            const T *
            pointers() const
            {
                return mMyPointers.data();
            }

            const T *
            indices() const
            {
                return mMyIndices.data();
            }

            const T
            n_rows() const
            {
                return mMyNumRows;
            }

            const real *
            values() const
            {
                if ( mCommRank == 0 )
                {
                    BELFEM_ASSERT( mMyValuesData != nullptr,
                        "values() called before distribute_values()" );

                    // rank 0 owns rows [ 0, dist(1) ), so its slice starts
                    // at the beginning of the matrix data
                    return mMyValuesData ;
                }
                else
                {
                    return mMyValues.data();
                }
            }

        protected:

            void
            distribute_sparsity_pattern() override
            {
                if ( mCommRank == 0 )
                {
                    const int_t * tPtrs = mMatrix->pointers();
                    const int_t * tIdx  = mMatrix->indices();

                    mAllPointers.set_size( mCommSize, {} );
                    mAllIndices.set_size( mCommSize, {} );

                    for ( proc_t p = 0; p < mCommSize; ++p )
                    {
                        Vector< T >    & tPointers = mAllPointers( p );
                        Vector< T >    & tIndices  = mAllIndices( p );

                        int_t a = mDist( p );
                        int_t b = mDist( p+1 );
                        int_t n = b - a;

                        // Count total nnz
                        int_t nnz = tPtrs[b] - tPtrs[a];

                        tPointers.set_size( n + 1, 0 );
                        tIndices.set_size( nnz );

                        // local index
                        int_t k = 0 ;
                        for ( int_t i = a; i < b; ++i )
                        {
                            // local row
                            int_t r = i - a ;

                            // nnz per row
                            nnz = tPtrs[i+1] - tPtrs[i];

                            tPointers( r+1 ) = tPointers( r ) + nnz;

                            for ( int_t j = tPtrs[i]; j < tPtrs[i+1]; ++j )
                            {
                                tIndices( k++ ) = tIdx[j];
                            }
                        }
                    }

                    comm_barrier();

                    share( mDist );

                    mMyPointers = std::move( mAllPointers( mCommRank ) );
                    mMyIndices  = std::move( mAllIndices( mCommRank ) );
                    mMyNumRows = mMyPointers.length() - 1;

                    distribute( mAllPointers );
                    distribute( mAllIndices );
                    comm_barrier() ;

                    mAllPointers.clear();
                    mAllIndices.clear();

                }
                else
                {
                    comm_barrier();
                    receive( mDist );
                    receive( mMyPointers );
                    receive( mMyIndices );
                    comm_barrier() ;

                    mMyNumRows = mMyPointers.length() - 1;
                    mMyValues.set_size( mMyIndices.length(), 0 );
                    mMyLhs.set_size( mMyNumRows, 0.0 );
                    mMyRhs.set_size( mMyNumRows, 0.0 );
                }
            }

        };

        template< typename T >
        class DistMatrixAIJ : public DistMatrix
        {
            Cell< Vector< T > > mAllDiagonalPointers ;
            Cell< Vector< T > > mAllOffDiagonalPointers ;
            Cell< Vector< T > > mAllDiagonalIndices ;
            Cell< Vector< T > > mAllOffDiagonalIndices ;
            Cell< Vector< T > > mAllGarrays ;

            Cell< Vector< real > > mAllDiagonalValues ;
            Cell< Vector< real > > mAllOffDiagonalValues ;

            Vector< T > mMyDiagonalPointers ;
            Vector< T > mMyOffDiagonalPointers ;
            Vector< T > mMyDiagonalIndices ;
            Vector< T > mMyOffDiagonalIndices ;
            Vector< T > mMyGarray ;

            Vector< real > mMyDiagonalValues ;
            Vector< real > mMyOffDiagonalValues ;

            const bool       mUseLocalIndices ;
        public:

            DistMatrixAIJ( const SolverParameters * aParams, SpMatrix * aMatrix ) :
                DistMatrix( aParams, aMatrix ),
                mUseLocalIndices( aParams->type() == SolverType::STRUMPACK )
            {
                this->distribute_sparsity_pattern();
            }

            void
            distribute_values( SpMatrix * aMatrix ) override
            {
                if ( mCommRank == 0 )
                {
                    const int_t * tPtrs = mPermutationSwitch ? mMatrix->pointers() : aMatrix->pointers();
                    const int_t * tIdx  = mPermutationSwitch ? mMatrix->indices() : aMatrix->indices();
                    real  * tData = mPermutationSwitch ? mMatrix->data() : aMatrix->data();

                    // reorganize data
                    if ( mPermutationSwitch )
                    {
                        real * tValues = aMatrix->data();

                        index_t tCount = 0 ;
                        for ( index_t k : mIndexPermutation )
                        {
                            tData[ tCount++ ] = tValues[ k ] ;
                        }
                    }

                    index_t tCount = 0 ;
                    for ( proc_t p=0; p<mCommSize; ++p )
                    {
                        Vector< real > & tDiagValues = mAllDiagonalValues( p );
                        Vector< real > & tOffDiagValues = mAllOffDiagonalValues( p );
                        int_t a = mDist( p );
                        int_t b = mDist( p+1 );

                        index_t d = 0 ;
                        index_t o = 0 ;

                        for ( int_t i=a; i<b; ++i )
                        {
                            int_t n = tPtrs[i+1] - tPtrs[i];

                            for ( int_t j=0; j<n; ++j )
                            {
                                int_t k = tIdx[ tPtrs[i] + j ];

                                if ( k < a || k >= b )
                                {
                                    tOffDiagValues( o++ ) = tData[ tCount++ ];
                                }
                                else
                                {
                                    tDiagValues( d++ ) = tData[ tCount++ ];
                                }
                            }
                        }
                    }

                    comm_barrier() ;
                    distribute( mAllDiagonalValues );
                    distribute( mAllOffDiagonalValues );
                    comm_barrier() ;
                }
                else
                {
                    comm_barrier();
                    receive( mMyDiagonalValues );
                    receive( mMyOffDiagonalValues );
                    comm_barrier() ;
                }
            }

            T *
            diagonal_pointers()
            {
                return mMyDiagonalPointers.data();
            }

            T *
            offdiagonal_pointers()
            {
                return mMyOffDiagonalPointers.data();
            }

            T *
            diagonal_indices()
            {
                return mMyDiagonalIndices.data();
            }

            T *
            offdiagonal_indices()
            {
                return mMyOffDiagonalIndices.data();
            }


            T *
            garray()
            {
                return mMyGarray.data();
            }

            T
            garray_size() const
            {
                return mMyGarray.length();
            }

            T
            n_rows() const
            {
                return mMyNumRows;
            }

            real *
            diagonal_values()
            {
                if ( mCommRank == 0 )
                {
                    return mAllDiagonalValues( 0 ).data();
                }
                else
                {
                    return mMyDiagonalValues.data();
                }
            }

            real *
            offdiagonal_values()
            {
                if ( mCommRank == 0 )
                {
                    return mAllOffDiagonalValues( 0 ).data();
                }
                else
                {
                    return mMyOffDiagonalValues.data();
                }
            }

        protected:

            void
            distribute_sparsity_pattern() override
            {
                if ( mCommRank == 0 )
                {
                    mAllDiagonalPointers.set_size( mCommSize, {} );
                    mAllOffDiagonalPointers.set_size( mCommSize, {} );
                    mAllDiagonalIndices.set_size( mCommSize, {} );
                    mAllOffDiagonalIndices.set_size( mCommSize, {} );
                    mAllDiagonalValues.set_size( mCommSize, {} );
                    mAllOffDiagonalValues.set_size( mCommSize, {} );
                    mAllGarrays.set_size( mCommSize, {} );

                    const int_t * tPtrs = mMatrix->pointers() ;

                    const int_t * tIdx  = mMatrix->indices() ;

                    for ( proc_t p = 0; p < mCommSize; ++p )
                    {
                        int_t a = mDist( p );
                        int_t b = mDist( p+1 );
                        int_t n = b - a;

                        // STEP 1: Collect unique off-diagonal columns for this rank
                        // Use std::set for automatic sorting and uniqueness
                        std::set< int_t > tUniqueOffDiagCols;

                        // Count total nnz and diagonal/off-diagonal entries
                        //int_t nnz = tPtrs[b] - tPtrs[a];
                        int_t nd = 0 ;
                        int_t no = 0 ;

                        // Count total nnz and diagonal/off-diagonal entries
                        for ( int_t i = a; i < b; ++i )
                        {
                            for ( int_t j=tPtrs[i]; j<tPtrs[i+1]; ++j )
                            {
                                int_t k = tIdx[j];
                                if ( a <=k && k < b )
                                {
                                    ++nd;
                                }
                                else
                                {
                                    ++no;
                                    tUniqueOffDiagCols.insert( k );
                                }
                            }
                        }

                        // STEP 2: Build garray from unique off-diagonal columns
                        // std::set is already sorted, so we can copy directly
                        Vector< int_t > & tGarray = mAllGarrays( p );
                        tGarray.set_size( tUniqueOffDiagCols.size() );
                        int_t garrayIdx = 0;
                        OrderedMap< int_t, int_t > tGlobalToLocal;  // Map global col -> local garray index
                        for ( int_t globalCol : tUniqueOffDiagCols )
                        {
                            tGarray( garrayIdx ) = globalCol;
                            tGlobalToLocal[ globalCol ] = garrayIdx++ ;
                        }

                        // STEP 3: Allocate MPIAIJ format (for set_MPIAIJ_matrix)

                        Vector< T > & tDiagPointers    = mAllDiagonalPointers( p );
                        Vector< T > & tOffDiagPointers = mAllOffDiagonalPointers( p );
                        Vector< T > & tDiagIndices     = mAllDiagonalIndices( p );
                        Vector< T > & tOffDiagIndices  = mAllOffDiagonalIndices( p );
                        Vector< real > & tDiagValues      = mAllDiagonalValues( p );
                        Vector< real > & tOffDiagValues   = mAllOffDiagonalValues( p );

                        tDiagPointers.set_size( n + 1, 0 );
                        tOffDiagPointers.set_size( n + 1, 0 );
                        tDiagIndices.set_size( nd );
                        tOffDiagIndices.set_size( no );
                        tDiagValues.set_size( nd, 0.0 );
                        tOffDiagValues.set_size( no, 0.0 );

                        // STEP4: Build the format
                        int_t d = 0 ;
                        int_t o = 0 ;
                        for ( int_t i=a; i<b; ++i )
                        {
                            // local row index
                            int_t r = i - a ;

                            nd = 0 ;
                            no = 0 ;

                            // loop over the columns
                            for ( int_t j=tPtrs[i]; j<tPtrs[i+1]; ++j )
                            {
                                int_t k = tIdx[j];
                                if ( a <=k && k < b )
                                {
                                    tDiagIndices( d++ ) = k - a ;
                                    ++nd ;
                                }
                                else
                                {
                                    // Off-diagonal block: convention depends on solver
                                    if ( mUseLocalIndices )
                                    {
                                        // STRUMPACK: use local garray index
                                        auto it = tGlobalToLocal.find( k );
                                        BELFEM_ASSERT( it != tGlobalToLocal.end(),
                                                     "Failed to find global column %lu in garray map for rank %lu",
                                                     (long unsigned int) k, (long unsigned int) p );
                                        tOffDiagIndices( o++ ) = it->second;
                                    }
                                    else
                                    {
                                        // PETSc: use global column index
                                        tOffDiagIndices( o++ ) = k;
                                    }
                                    ++no ;
                                }
                            }

                            tDiagPointers( r+1 ) = tDiagPointers( r ) + nd;
                            tOffDiagPointers( r+1 ) = tOffDiagPointers( r ) + no;
                        }
                    }

                    comm_barrier();
                    share( mDist );

                    mMyDiagonalPointers = std::move( mAllDiagonalPointers( mCommRank ) );
                    mMyOffDiagonalPointers = std::move( mAllOffDiagonalPointers( mCommRank ) );
                    mMyDiagonalIndices = std::move( mAllDiagonalIndices( mCommRank ) );
                    mMyOffDiagonalIndices = std::move( mAllOffDiagonalIndices( mCommRank ) );
                    mMyGarray = std::move( mAllGarrays( mCommRank ) );

                    distribute( mAllDiagonalPointers );
                    distribute( mAllOffDiagonalPointers );
                    distribute( mAllDiagonalIndices );
                    distribute( mAllOffDiagonalIndices );
                    distribute( mAllGarrays );

                    comm_barrier();

                    mAllDiagonalPointers.clear();
                    mAllOffDiagonalPointers.clear();
                    mAllDiagonalIndices.clear();
                    mAllOffDiagonalIndices.clear();
                    mAllGarrays.clear();
                }
                else
                {
                    comm_barrier();
                    receive( mDist );
                    receive( mMyDiagonalPointers );
                    receive( mMyOffDiagonalPointers );
                    receive( mMyDiagonalIndices );
                    receive( mMyOffDiagonalIndices );
                    receive( mMyGarray );
                    comm_barrier();

                    mMyDiagonalValues.set_size( mMyDiagonalIndices.length(), 0 );
                    mMyOffDiagonalValues.set_size( mMyOffDiagonalIndices.length(), 0 );
                    mMyNumRows = mMyDiagonalPointers.length() - 1;
                    mMyLhs.set_size( mMyNumRows, 0.0 );
                    mMyRhs.set_size( mMyNumRows, 0.0 );

                }
            }
        };

    }
}
#endif //BELFEM_CL_SOLVERDISTMATRIX_HPP