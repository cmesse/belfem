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

#ifndef BELFEM_CL_ELECTRICNODE_HPP
#define BELFEM_CL_ELECTRICNODE_HPP

#include "typedefs.hpp"
#include "cl_Graph_Vertex.hpp"

namespace belfem
{
    namespace electronics
    {
//-----------------------------------------------------------------------------

        /**
         * @brief Circuit node; the unknowns of the modified nodal analysis.
         *
         * @ingroup grp_circuit
         * @see @ref circuit_circuit_usage_guide
         */
        class ElectricNode : public graph::Vertex
        {

            //!Voltage of this node
            real mVoltage;

//-----------------------------------------------------------------------------

        public:
            ElectricNode() ;

//-----------------------------------------------------------------------------

            void
            set_voltage(const real aVoltage) ;

//-----------------------------------------------------------------------------

            real
            get_voltage() const ;

//-----------------------------------------------------------------------------

        };
    }

}

#endif //BELFEM_CL_ELECTRICNODE_HPP
