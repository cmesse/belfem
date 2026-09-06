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

#ifndef BELFEM_CL_FEMTWOTERMINALS_HPP
#define BELFEM_CL_FEMTWOTERMINALS_HPP

#include "cl_TwoTerminals.hpp"
#include "Component_Enums.hpp"
#include "cl_Vector.hpp"
#include "cl_BDF.hpp"

namespace belfem
{
    namespace electronics
    {
//-----------------------------------------------------------------------------

        /**
         * @brief Terminal pair coupling the circuit to the finite-element problem.
         *
         * @ingroup grp_circuit
         * @see @ref circuit_circuit_usage_guide
         */
        class FEMTwoTerminals : public TwoTerminals
        {

            //! Current value from previous FEM problem
            real mIn = 0.0 ;

            //! Voltage value from previous FEM problem
            real mVn = 0.0 ;

            //! Time discretized resistance
            //real mRFEM = 1e-10 ;

            //real miFEMh = 0.0 ;

            //! BDF object
            ode::BDF * mBDF  = nullptr;

            //! Time Step Shift Register
            ShiftRegister< real > mH ;

            //! Current Shift Register
            ShiftRegister< real > mI ;

            //! Voltage Shift Register
            ShiftRegister< real > mV ;

            //! Order
            //uint mOrder = 1 ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            FEMTwoTerminals(Cell < ElectricNode* > & aTerminals, string aLabel = "") ;

//------------------------------------------------------------------------------

            ComponentType
            component_type() const override ;

//-----------------------------------------------------------------------------

            void
            compute_current() override ;

//------------------------------------------------------------------------------

            void
            shift( const real aTime, const real aDeltaTime ) override ;

//------------------------------------------------------------------------------

            void
            shift_back() override ;

//------------------------------------------------------------------------------

            void
            set_IV( const real I, const real V ) ;

//------------------------------------------------------------------------------

            void
            save_state( hid_t aGroup, const string & aPrefix ) override ;

//------------------------------------------------------------------------------

            void
            load_state( hid_t aGroup, const string & aPrefix ) override ;

//------------------------------------------------------------------------------

            /*real
            dIdV() const ;*/

//------------------------------------------------------------------------------

        };

        inline ComponentType
        FEMTwoTerminals::component_type() const
        {
            return ComponentType::TERMINALPAIR ;
        }

//-----------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_FEMTWOTERMINALS_HPP
