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

#ifndef BELFEM_CL_EF_HEX8TS_HPP
#define BELFEM_CL_EF_HEX8TS_HPP


#include "nedelec/cl_EF_EdgeFunction.hpp"
#include "cl_FEM_Element.hpp"

namespace belfem
{
    namespace fem
    {
        class Element;

        class EF_HEX8TS : public EdgeFunction
        {
            Matrix< real > mX ; // node coordinates
            InterpolationFunction * mQuad4  = nullptr;

            Cell< Vector< real > > mF ;
            Cell< Matrix< real > > mFxi ;

            Cell< Matrix< real > > mHex8Nxi ;
            Cell< Matrix< real > > mQuad4Nxi ;

            Matrix< real > mEx ;
            Matrix< real > mEy ;
            Matrix< real > mEz ;

            Matrix< real > mNx ;
            Matrix< real > mG ;
            Matrix< real > mInvG ;
            Matrix< real > mJ2 ;

            // edge directions
            real mS[ 8 ];

            // help vectors
            Vector< real > mU ;
            Vector< real > mV ;
            Vector< real > mW ;

            real mSurface = BELFEM_QUIET_NAN ;
            real mThickness  = BELFEM_QUIET_NAN ;
            real mVolume  = BELFEM_QUIET_NAN ;

            uint mLastIndex = BELFEM_UINT_MAX ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * constructor
             */
            EF_HEX8TS();

            /**
             * destructor
             */
            ~EF_HEX8TS() override;

//------------------------------------------------------------------------------

            /**
             * connect function with element
             */
            void
            link( Element *aElement ) override;

//------------------------------------------------------------------------------

            void
            precompute( const Matrix< real > & aXi ) override ;

//------------------------------------------------------------------------------

            const Matrix <real> &
            E( const uint aIndex ) override ;

//------------------------------------------------------------------------------

            const Matrix <real> &
            C( const uint aIndex ) override ;

//------------------------------------------------------------------------------

            /**
              * gradient operator for h-field,
              * layout contract in EdgeFunction::mGrad.
              *
              * G chains the reference derivatives through the class's own
              * per-point nablas held fixed — exactly the ingredients C()
              * uses, so antisym( G ) reproduces the implemented C on every
              * geometry. Scope by mid-surface shape:
              *  - affine ( parallelogram ) quads: the map is affine and G
              *    is the exact gradient; on rectangles additionally
              *    div( e_k ) == 0.
              *  - planar non-parallelogram quads: the nablas are exact
              *    gradients and C is the exact curl, but G omits the
              *    symmetric inverse-map Hessian term F_k Hess( xi_a ) —
              *    G is NOT the full gradient there.
              *  - warped quads: nablas follow the thin-shell update_nabla
              *    convention ( mid-surface Jacobian only ); E, C and G all
              *    share that approximation.
              */
            const Matrix< real > &
            G( const uint aIndex ) override ;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            update_nabla( const uint aIndex ) override ;

//------------------------------------------------------------------------------
        };
    }
}

#endif //BELFEM_CL_EF_HEX8TS_HPP
