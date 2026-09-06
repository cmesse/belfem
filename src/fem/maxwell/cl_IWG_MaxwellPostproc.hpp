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

#ifndef CL_IWG_MAXWELLPOSTPROC_HPP
#define CL_IWG_MAXWELLPOSTPROC_HPP

#include "cl_IWG.hpp"
#include "cl_Maxwell_FieldList.hpp"
#include "cl_FEM_Dof.hpp"
#include "en_Maxwell_Formulations.hpp"

namespace belfem
{
    namespace fem
    {
        class IWG_MaxwellPostproc : public IWG
        {
            // defines which formulation we have
            const maxwell::Formulation mFormulation ;

            bool mUseEdges = false ;

            void
            ( *mFunKF )
                      ( Calculator     * aCalc,
                        Matrix< real > & aK,
                        Vector< real > & aF );
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_MaxwellPostproc(
              const maxwell::Formulation aFormulation,
              const ModelDimensionality aDimensionality,
                const bool aHigherOrder = false );

            void
            create_custom_vectors_and_matrices( Calculator * aCalc ) override ;

            void
            compute_jacobian_and_rhs(
                Element *aElement,
                Matrix<real> &aJacobian,
                Vector<real> &aRHS) override ;
        };
    }
}

#endif //CL_IWG_MAXWELLPOSTPROC_HPP
