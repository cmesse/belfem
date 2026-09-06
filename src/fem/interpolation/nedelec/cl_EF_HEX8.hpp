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

#ifndef BELFEM_CL_EF_HEX8_HPP
#define BELFEM_CL_EF_HEX8_HPP

#include "nedelec/cl_EF_EdgeFunction.hpp"
#include "cl_FEM_Element.hpp"

namespace belfem
{
    namespace fem
    {
        class Element;

        class EF_HEX8 : public EdgeFunction
        {
            Matrix< real > mX ; // node coordinates
            InterpolationFunction * mHex8  = nullptr;

            Cell< Vector< real > > mF ;
            Cell< Matrix< real > > mFxi ;

            Cell< Matrix< real > > mHex8Nxi ;

            //! per-point Voigt second derivatives of the trilinear Lagrange
            //! basis ( full 6 x 8 from d2NdXi2; rows 0-2 are identically
            //! zero for trilinear N, G() reads the mixed rows 5 = xi-eta,
            //! 4 = xi-zeta, 3 = eta-zeta — the order pinned by
            //! test_LagrangeInterpolation ). Filled in precompute() only.
            Cell< Matrix< real > > mHex8Nxi2 ;

            Matrix< real > mEx ;
            Matrix< real > mEy ;
            Matrix< real > mEz ;

            Matrix< real > mNx ;

            // help vectors for volume estimation
            Vector< real > mU ;
            Vector< real > mV ;

            // edge directions
            real mS[ 12 ];

            uint mLastIndex = BELFEM_UINT_MAX ;

            real mApproxVolume = BELFEM_QUIET_NAN ;
            bool mSmallVolume = false ;

            const Vector< uint > mXLengthA = { 0, 3, 4, 7 };
            const Vector< uint > mXLengthB = { 1, 2, 5, 6 };

            const Vector< uint > mYLengthA = { 0, 1, 4, 5 };
            const Vector< uint > mYLengthB = { 2, 3, 6, 7 };

            const Vector< uint > mZLengthA = { 0, 1, 2, 3 };
            const Vector< uint > mZLengthB = { 4, 5, 6, 7 };

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * constructor
             */
            EF_HEX8();

            /**
             * destructor
             */
            ~EF_HEX8() override;

//------------------------------------------------------------------------------

            void
            precompute( const Matrix< real > & aXi ) override ;

//------------------------------------------------------------------------------

            /**
             * connect function with element
             */
            void
            link( Element *aElement ) override;

//------------------------------------------------------------------------------

            const Matrix <real> &
            E( const uint aIndex ) override ;

//------------------------------------------------------------------------------

            const Matrix <real> &
            C( const uint aIndex ) override ;

//------------------------------------------------------------------------------

            /**
              * gradient operator for h-field
              * ( point-dependent; carries the inverse-map Hessian term the
              *   curl cancels — layout contract in EdgeFunction::mGrad )
              */
            const Matrix< real > &
            G( const uint aIndex ) override ;

//------------------------------------------------------------------------------

            void
            update_nabla( const uint aIndex ) override ;

//------------------------------------------------------------------------------
        private:

            void
            estimate_volume();

        };
    }
}
#endif //BELFEM_CL_EF_HEX8_HPP
