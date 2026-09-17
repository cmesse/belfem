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

#include "mt_maxwell_background.hpp"

#include "cl_FEM_Block.hpp"
#include "cl_FEM_DofManagerBase.hpp"
#include "cl_IWG.hpp"
#include "cl_IWG_Maxwell.hpp"
#include "cl_FEM_Kernel.hpp"
#include "fn_inv2.hpp"
#include "fn_inv3.hpp"
#include "fn_inv.hpp"
#include "fn_det.hpp"
#include "fn_trans.hpp"

namespace belfem
{
    namespace fem
    {
        namespace maxwell
        {
            void
            background_phi( Calculator * aCalc, TimestepMatrices * aMatrices )
            {

                //Get the boundary condition
                real h = aCalc->group()->parent()->parent()->boundary_condition(aCalc->group()->id())->value() ;

                Vector< real > & H = aCalc->vector("h");
                const Vector< real > & r = aCalc->group()->parent()->parent()->boundary_condition(aCalc->group()->id())->direction();

                H = h*r ;

                const Vector< real > & w =  aCalc->integration()->weights();

                for ( uint k=0; k<aCalc->num_intpoints(); ++k )
                {
                    const Matrix< real > & B = aCalc->Bm( k );
                    aMatrices->K() += w( k ) * trans( B ) * B * aCalc->dS();
                    aMatrices->f()   -= w( k ) * trans( B ) * H * aCalc->dS();
                }

            }

        }
    }
}