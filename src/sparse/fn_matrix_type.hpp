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

#ifndef FN_MATRIX_TYPE_HPP
#define FN_MATRIX_TYPE_HPP
#include "en_SolverEnums.hpp"
#include "cl_SpMatrix.hpp"
#include "assert.hpp"

namespace belfem
{
    inline
    SpMatrixType matrix_type( const SolverType aType )
    {
        switch( aType )
        {
            case SolverType::UMFPACK :
            {
                return SpMatrixType::CSC ;
            }
            case SolverType::SUPERLU :
            {
                return SpMatrixType::CSC ;
            }
            case SolverType::MUMPS :
            {
                return SpMatrixType::CSC ;
            }
            case( SolverType::PARDISO ) :
            {
                return SpMatrixType::CSR ;
            }
            case( SolverType::PETSc ) :
            {
                return SpMatrixType::CSR ;
            }
            case( SolverType::STRUMPACK ) :
            {
                return SpMatrixType::CSR ;
            }
            default:
            {
                BELFEM_ERROR( false, "Unknown solver type");
                return SpMatrixType::UNDEFINED ;
            }
        }
    }
}
#endif //FN_MATRIX_TYPE_HPP
