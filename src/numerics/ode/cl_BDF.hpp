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

#ifndef BELFEM_CL_BDF_HPP
#define BELFEM_CL_BDF_HPP

#include "cl_Matrix.hpp"
#include "cl_ShiftRegister.hpp"
#include "cl_Vector.hpp"

#include "typedefs.hpp"

namespace belfem
{
    namespace ode
    {
        class BDF
        {
            ShiftRegister< real > & mH ;

            Vector< real > mCoefficients ;

            // normalized time values
            Vector< real > mS ;

            Matrix< real > mVandermonde ;

            // pivot vector
            Vector< int_t > mPivot ;

            // for vector integration
            Vector< real > mF ;

        public:

            BDF( ShiftRegister< real > & aH ) ;

            void
            compute_coefficients();

            const Vector< real > &
            coefficients() const ;

            real
            eval( ShiftRegister< real > & aY , const real aF, const bool aUpdateCoefficients = true );

            real
            deval( const ShiftRegister< real > & aY, const bool aUpdateCoefficients = true );

            const Vector< real > &
            eval( ShiftRegister< Vector< real > > & aY , const Vector< real > & aF, const bool aUpdateCoefficients = true );

            const Vector< real > &
            deval( const ShiftRegister< Vector< real > > & aY, const bool aUpdateCoefficients = true );


        };


        inline const Vector< real > &
        BDF::coefficients() const
        {
            return mCoefficients ;
        }

    }
}
#endif // BELFEM_CL_BDF_HPP
