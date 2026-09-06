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

#ifndef BELFEM_CL_ELECTRICALCIRCUITFACTORY_HPP
#define BELFEM_CL_ELECTRICALCIRCUITFACTORY_HPP

#include <memory>

#include "typedefs.hpp"

#include "cl_InputFile.hpp"
#include "cl_FEM_PhysicalBoundaryCondition.hpp"
#include "cl_ElectricalCircuit.hpp"

namespace belfem
{
    namespace electronics
    {
//-----------------------------------------------------------------------------
        /**
         * @brief Builds the circuit from the input deck's circuit{} section.
         *
         * @ingroup grp_circuit
         * @see @ref circuit_circuit_usage_guide
         */
        class ElectricalCircuitFactory
        {

            const proc_t mCommRank ;

            //! unique_ptr members so a refusal that throws mid-construction
            //! ( the dtor of a half-built object never runs ) cannot leak
            std::unique_ptr< const InputFile > mInputFile ;

            Cell < fem::PhysicalBoundaryCondition * > & mPhysicalBoundaryConditions ;

            std::unique_ptr< ElectricalCircuit > mCircuit ;

//-----------------------------------------------------------------------------

        public:

            ElectricalCircuitFactory( const string & aInputFile, Cell< fem::PhysicalBoundaryCondition * > & aPhysicalBoundaryConditions ) ;

            ~ElectricalCircuitFactory() ;

//-----------------------------------------------------------------------------

            void
            read_circuit( const input::Section * aSection ) ;

//-----------------------------------------------------------------------------

            Circuit *
            circuit() ;

//-----------------------------------------------------------------------------

        private:

//-----------------------------------------------------------------------------

            /**
             * the hybrid path ( circuit { file : ... } ): the netlist holds
             * the lumped topology, the deck holds terminal pairs and output,
             * with node references resolved as netlist node NAMES
             */
            void
            read_circuit_from_netlist( const input::Section * aSection ) ;

            /**
             * one terminal pair: the lumped component plus its boundary
             * condition -- shared by the classic and the netlist path
             */
            void
            read_terminal_pair( const input::Section * aComponentSection,
                                const id_t aSectionIndex,
                                const index_t aNodePlus,
                                const index_t aNodeMinus,
                                const string & aLabel,
                                uint & aCount ) ;

            /**
             * the output block; aNetlist resolves voltage node names on the
             * hybrid path, nullptr keeps the classic integer grammar
             */
            void
            read_output( const input::Section * aSection,
                         const class NgspiceCircuitFactory * aNetlist ) ;

//-----------------------------------------------------------------------------

        };
    }
}

#endif //BELFEM_CL_ELECTRICALCIRCUITFACTORY_HPP
