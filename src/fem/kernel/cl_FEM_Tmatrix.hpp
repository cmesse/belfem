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

#ifndef BELFEM_CL_FEM_TMATRIX_HPP
#define BELFEM_CL_FEM_TMATRIX_HPP
#include "typedefs.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace fem
    {
        class Tmatrix
        {
            const uint mNumRows ;
            const uint mNumCols  ;
            const uint mNumNonZeros  ;

            // data are stored in CSR format
            uint * mPointers = nullptr ; // for rows
            uint * mIndices  = nullptr ; // for columns
            real * mValues   = nullptr ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Tmatrix( const Matrix< real > & aMatrix );

//------------------------------------------------------------------------------

            ~Tmatrix() ;

//------------------------------------------------------------------------------

            /**
             * performs the operation B_j = T_ij * A_i  ( B = T^T A )
             */
            void
            project( const Vector< real > & aA, Vector< real > & aB ) const ;

//------------------------------------------------------------------------------

            /**
             * performs the operation B_il = T_ji * A_jk * T_kl
             */
            void
            project( const Matrix< real > & aA, Matrix< real > & aB ) const ;

//------------------------------------------------------------------------------

            /**
             * for debugging
             */
             void
             get_matrix( Matrix< real > & aMatrix ) const ;

//------------------------------------------------------------------------------

            /**
             * for debugging
             */
            void
            print() const ;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            uint
            count_nnz( const Matrix< real > & aMatrix );

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_FEM_TMATRIX_HPP
