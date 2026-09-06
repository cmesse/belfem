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

#ifndef BELFEM_CL_EF_HEX8TB_HPP
#define BELFEM_CL_EF_HEX8TB_HPP

#include "nedelec/cl_EF_EdgeFunction.hpp"
#include "cl_FEM_Element.hpp"

namespace belfem
{
    namespace fem
    {
        class Element;

        /**
         * Edge function for the HEX8TB side-connector wall element.
         *
         * Four longitudinal edge dofs on edges (0,1), (3,2), (4,5), (7,6),
         * all pointing in +xi (the side-curve direction). The element is a
         * degenerate cuboid: the metric uses the exact width (eta) and
         * thickness (zeta) stored on the mesh blocks, plus the
         * midpoint-to-midpoint length (xi). The node positions only orient
         * the frame; they never set the thin dimensions, which are physical
         * collapsed-layer data. The Jacobian is therefore constant per
         * element and is computed once in link().
         */
        class EF_HEX8TB : public EdgeFunction
        {
            Cell< Vector< real > > mF ;
            Cell< Matrix< real > > mFxi ;

            // edge directions
            real mS[ 4 ];

            // help vectors
            Vector< real > mU ; // midpoint of face {0,3,4,7} at xi = -1, then tangent t
            Vector< real > mV ; // midpoint of face {1,2,6,5} at xi = +1, then normal n
            Vector< real > mW ; // scratch, then binormal b

            //! constant curl factor nabla_eta x nabla_xi = -4/(w*L) * n
            Vector< real > mP ;

            //! constant curl factor nabla_zeta x nabla_xi = 4/(d*L) * b
            Vector< real > mQ ;

            real mThickness  = BELFEM_QUIET_NAN ;
            real mWidth      = BELFEM_QUIET_NAN ;
            real mLength     = BELFEM_QUIET_NAN ;
            real mCrossSection = BELFEM_QUIET_NAN ;

            real mVolume  = BELFEM_QUIET_NAN ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * constructor
             */
            EF_HEX8TB();

            /**
             * destructor
             */
            ~EF_HEX8TB() override = default ;

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
              * Exact on the class's own convention: the wall metric is the
              * imposed exact cuboid frame ( t, b, n ) from link(), so the
              * map is affine and grad( e_k ) = s_k ( grad F_k ) (x) nabla xi
              * with all frame factors element constants. The frame is
              * orthonormal, hence div( e_k ) == 0 identically.
              */
            const Matrix< real > &
            G( const uint aIndex ) override ;

//------------------------------------------------------------------------------

            /**
             * the Jacobian is constant per element and set in link();
             * this override only exists because the Calculator calls it
             * before reading det_J()
             */
            void
            update_nabla( const uint aIndex ) override ;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline void
        EF_HEX8TB::update_nabla( const uint /* aIndex */ )
        {
            // nothing to do: mJ, mInvJ, mDetJ are element constants,
            // computed once in link()
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_EF_HEX8TB_HPP
