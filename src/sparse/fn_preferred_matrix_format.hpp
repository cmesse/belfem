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

#ifndef FN_PREFERRED_MATRIX_FORMAT_HPP
#define FN_PREFERRED_MATRIX_FORMAT_HPP
#include "cl_SpMatrix.hpp"
#include "en_SolverEnums.hpp"
namespace belfem
{
    inline SpMatrixType
    preferred_matrix_format( SolverType aSolver )
    {
        switch( aSolver )
        {
            case( SolverType::UMFPACK ) :
            case( SolverType::SUPERLU ) :
            {
                return SpMatrixType::CSC ;
            }
            case( SolverType::STRUMPACK ) :
            case( SolverType::PARDISO ) :
            case( SolverType::PETSc ) :
            case( SolverType::MUMPS ) :
            {
                return SpMatrixType::CSR ;
            }
            default:
            {
                return SpMatrixType::CSR ;
            }
        }
    }
}
#endif //FN_PREFERRED_MATRIX_FORMAT_HPP
