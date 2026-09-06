//
// Created by christian on 10/23/24.
//

#ifndef MT_MAXWELL_H_HPP
#define MT_MAXWELL_H_HPP

#include "typedefs.hpp"
#include "constants.hpp"
#include "globals.hpp"

#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

#include "cl_FEM_Element.hpp"
#include "cl_FEM_Calculator.hpp"
#include "cl_TimestepMatrices.hpp"


namespace belfem
{
    namespace fem
    {
        namespace maxwell
        {
            inline void
            save_resistivity( Calculator * aCalc , const real aRhoM )
            {
                aCalc->mesh()->field_data( "element_rho" )( aCalc->element()->element()->index() ) = aRhoM  ;
            }

            // ohmic dissipation of one element, added into the mesh-wide
            // "dotQ" global. This is a POWER in watts ( rho |j|^2 dV ), not an
            // energy: an AC loss per cycle is the time integral of it, which
            // nothing in the tree takes. The global is zeroed before every
            // assembly by Controller::reset_dotQ(), so what reaches the Exodus
            // file is the value from that timestep's LAST assembly
            inline void
            save_dotQ( Calculator * aCalc , const real aDotQ )
            {
                aCalc->mesh()->global_variable( "dotQ" )->value() += aDotQ ;
            }

            inline real
            get_resistivity( Calculator * aCalc, const index_t aIndex=BELFEM_UINT_MAX )
            {
                return std::clamp(  aCalc->mesh()->field_data( "element_rho" )
                    ( aIndex == BELFEM_UINT_MAX ? aCalc->element()->element()->index() : aIndex ), gRhoMin, gRhoMax );
            }

            /**
             * Collapsed h-kernels: all per-point material math is delegated
             * to aCalc->maxwell(), whose constructor-time dispatch covers the
             * bulk/thin-shell, metal/alloy/HTS, defect, and piecewise axes.
             * Dispatch in IWG_Maxwell::link_to_group() is by solver algorithm
             * and mu constancy only:
             *
             *   Picard                       -> h_picard
             *   Newton, mu constant          -> h_newton_mu0
             *   Newton, mu field-dependent   -> h_newton_mu
             */
            void
            h_picard( Calculator * aCalc, TimestepMatrices * aMatrices );

            void
            h_newton_mu0( Calculator * aCalc, TimestepMatrices * aMatrices );

            void
            h_newton_mu( Calculator * aCalc, TimestepMatrices * aMatrices );

            /**
             * facet-based thin-shell interface / stabilization kernel
             */
            void
            h_ghost( Calculator * aCalc, TimestepMatrices * aMatrices );

            void
            h_side_connector( Calculator * aCalc, TimestepMatrices * aMatrices );

            void
            h_side_connector_newton( Calculator * aCalc, TimestepMatrices * aMatrices );
        }
    }
}

#endif //MT_MAXWELL_H_HPP
