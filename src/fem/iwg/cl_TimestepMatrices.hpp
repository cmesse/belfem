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

#ifndef BELFEM_CL_FEM_TIMESTEPMATRICES_HPP
#define BELFEM_CL_FEM_TIMESTEPMATRICES_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Bitset.hpp"

namespace belfem
{
    namespace fem
    {
        /**
         * @brief Flags indicating which matrices/vectors have been set by the IWG
         *
         * The IWG declares in its constructor which matrices it will populate;
         * initialize() sizes exactly those, reset() zeroes them, and the
         * time-stepping scheme reads the flags to know what to include in assembly.
         */
        enum class MatrixFlag : index_t
        {
            M            = 0,   ///< Mass matrix
            D            = 1,   ///< Damping matrix (reserved for second-order)
            K            = 2,   ///< Stiffness matrix
            F            = 3,   ///< Force vector
            dMdX_times_x = 4,   ///< Contracted mass derivative: (dM_ik/dx_j) * x_k
            dMdX_times_h = 5,   ///< Contracted mass derivative with history: (dM_ik/dx_j) * h_k
            dKdX_times_x = 6,   ///< Contracted stiffness derivative
            dFdX         = 7    ///< Force vector Jacobian df_i/dx_j
        };

        /**
         * @class TimestepMatrices
         * @brief Container for element-level DENSE matrices in transient nonlinear FEM
         *
         * This class serves as interface between the physics-aware IWG and
         * the time-stepping scheme (BDF). All matrices are DENSE and sized
         * n_e × n_e where n_e is the number of element DOFs.
         *
         * The global SPARSE Jacobian is assembled from these element contributions
         * in a separate assembler class.
         *
         * Workflow:
         *   1. IWG constructor sets the flags for the matrices it produces
         *   2. link_to_group() calls initialize( n_e ), which sizes the flagged matrices
         *   3. per element: reset(), then the IWG fills M, K, f and the derivative blocks
         *   4. assemble_dJdx() combines the flagged blocks with the BDF coefficients
         *   5. Global assembler accumulates element J^e into sparse global Jacobian
         *
         * Key design decisions:
         *   - All matrices are DENSE (element-level, typically small: 4-27 DOFs)
         *   - Derivatives are stored as contracted n_e × n_e matrices (dX/dx · x),
         *     not as n_e × n_e × n_e tensors. The IWG computes this contraction directly.
         *   - Uses Bitset<8> to track which matrices are populated.
         *
         * @ingroup grp_fem_iwg
         * @see @ref fem_iwg_iwg_usage_guide
         */
        class TimestepMatrices
        {
            //--------------------------------------------------------------
            // Physics matrices (populated by IWG) - all DENSE, n_e × n_e
            //--------------------------------------------------------------

            Matrix< real > mM;    ///< Mass matrix M(x), dense
            Matrix< real > mD;    ///< Damping matrix D(x), dense (reserved for second-order)
            Matrix< real > mK;    ///< Stiffness matrix K(x), dense
            Vector< real > mF;    ///< Force vector f(x), dense

            // Contracted derivatives: (dX_ik/dx_j) * v_k for some vector v
            // These are dense n_e × n_e matrices, computed by IWG
            Matrix< real > mdMdX_times_x;   ///< (dM_ik/dx_j) * x_k
            Matrix< real > mdMdX_times_h;   ///< (dM_ik/dx_j) * h_k (history correction)
            Matrix< real > mdKdX_times_x;   ///< (dK_ik/dx_j) * x_k
            Matrix< real > mdFdX;           ///< df_i/dx_j

            //--------------------------------------------------------------
            // Assembled element system (computed by time-stepping scheme)
            //--------------------------------------------------------------

            Matrix< real > mJ;    ///< Jacobian
            Matrix< real > mdJdx;    ///< Newton correction term to the jacobian

            //--------------------------------------------------------------
            // Flags
            //--------------------------------------------------------------

            Bitset< 8 > mFlags;

        public:

            //--------------------------------------------------------------
            // Constructor / Destructor
            //--------------------------------------------------------------

            TimestepMatrices() = default;
            ~TimestepMatrices() = default;

            /**
             * @brief Initialize all matrices to given size
             * @param aNumDofs Number of element DOFs (typically 4-27 for common elements)
             */
            void
            initialize( const index_t aNumDofs );

            void
            reset();

            //--------------------------------------------------------------
            // Flag management
            //--------------------------------------------------------------

            /**
             * @brief Clear all flags ( unused; the IWG sets its flags once in the constructor )
             */
            void
            reset_flags()
            {
                mFlags.reset();
            }

            /**
             * @brief Set a flag after populating corresponding matrix
             */
            void
            set_flag( const MatrixFlag aFlag )
            {
                mFlags.set( static_cast< index_t >( aFlag ) );
            }

            /**
             * @brief Check if a matrix has been populated
             */
            bool
            has_flag( const MatrixFlag aFlag ) const
            {
                return mFlags.test( static_cast< index_t >( aFlag ) );
            }

            /**
             * @brief Direct access to flags bitset
             */
            const Bitset< 8 > &
            flags() const
            {
                return mFlags;
            }

            //--------------------------------------------------------------
            // Const accessors
            //--------------------------------------------------------------

            const Matrix< real > & M() const { return mM; }
            const Matrix< real > & D() const { return mD; }
            const Matrix< real > & K() const { return mK; }
            const Vector< real > & f() const { return mF; }

            const Matrix< real > & dMdx_times_x() const { return mdMdX_times_x; }
            const Matrix< real > & dMdx_times_h() const { return mdMdX_times_h; }
            const Matrix< real > & dKdx_times_x() const { return mdKdX_times_x; }
            const Matrix< real > & dfdx()         const { return mdFdX; }

            const Matrix< real > & J() const { return mJ; }
            const Matrix< real > & dJdx() const { return mdJdx; }

            //--------------------------------------------------------------
            // Non-const accessors (for IWG to fill)
            //--------------------------------------------------------------

            Matrix< real > & M() { return mM; }
            Matrix< real > & D() { return mD; }
            Matrix< real > & K() { return mK; }
            Vector< real > & f() { return mF; }

            Matrix< real > & dMdx_times_x() { return mdMdX_times_x; }
            Matrix< real > & dMdx_times_h() { return mdMdX_times_h; }
            Matrix< real > & dKdx_times_x() { return mdKdX_times_x; }
            Matrix< real > & dfdx()         { return mdFdX; }

            //--------------------------------------------------------------
            // Assembly methods
            //--------------------------------------------------------------

            void
            assemble_J( const real adt );

            /**
             * @brief Assemble the Newton correction dJ/dx
             * @param adt    length of the current timestep
             * @param aAlpha BDF coefficient of the current time level; scales
             *        the (dM/dx)·x term (Hairer & Wanner 1996, II.4). One for
             *        all schemes that scale M by one (BDF1, Crank-Nicolson, ...)
             */
            void
            assemble_dJdx( const real adt, const real aAlpha = 1.0 );


        };

    } // namespace fem
} // namespace belfem

#endif // BELFEM_CL_FEM_TIMESTEPMATRICES_HPP