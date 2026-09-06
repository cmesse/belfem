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

#ifndef BELFEM_CL_NGSPICECIRCUITFACTORY_HPP
#define BELFEM_CL_NGSPICECIRCUITFACTORY_HPP

#include <memory>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"
#include "cl_NetlistParser.hpp"

namespace belfem
{
    class SourceFunction;

    namespace electronics
    {
//-----------------------------------------------------------------------------

        class ElectricalCircuit;

//-----------------------------------------------------------------------------

        /**
         * Builds an ElectricalCircuit from a lexed netlist
         * ( todo/ngspice_parser_plan.md §9.2, decisions O5/O6/O8/O9 ).
         *
         * Node map ( O9, frozen for restart ): the ground aliases "0" and
         * "gnd" are one node and land on the LAST BELFEM index; every other
         * node packs into 0 .. N-2 in order of first appearance in the
         * netlist ( element cards and directive n+/n- keys alike ). A deck
         * that never references ground is invalid ( manual, "Ground node" ).
         * The map is logged at construction and exposed via node_index()
         * so the deck-side terminal-pair reader resolves node NAMES through
         * the same map ( O5 ) -- never raw integers.
         *
         * Components are created in netlist line order ( elements and
         * directives merged by line ), which freezes the unknown-current
         * layout for restart ( O9 ). Labels are the case-folded instance
         * names ( O5 ).
         *
         * v1 element support, refusals are hard errors naming the card:
         *  - R/C/L: value from the positional field or the r=/c=/l= keyword
         *    ( exactly one ); nonzero; `ic=` and every other keyword refused
         *    ( value-changing, plan §12.2 ). BDF order from the
         *    `* belfem: order <instance> <n>` directive ( O8 ), default 1,
         *    valid 1..6 ( the BDF kernel's own bound ).
         *  - V/I: bare value or `DC <value>` -> constant source;
         *    `SIN(vo va freq [td theta phase])` -> sine, where vo/td/theta
         *    must be zero or absent ( BELFEM sine cannot represent them );
         *    phase is in DEGREES per SPICE and converted to radians.
         *    PULSE/PWL/EXP/SFFM/AM are refused until Phase 6.
         *  - D: `<model>` referencing a `.model <name> d(...)` card;
         *    Is from `is=` ( default 1e-14 A per SPICE ), Vt = n * 0.026 V
         *    with `n` from the model ( default 1 ); any other model
         *    parameter, an area factor, or a non-d model type is refused.
         *  - directives: `superconductor <name> n+= n-= Ic= n= Ec= length=`,
         *    `switch <name> n+= n-= state=open|closed t_switch=`,
         *    `order <instance> <n>`; unknown kinds or keys are refused.
         *
         * Ownership: the factory produces and hands over -- circuit()
         * transfers the ElectricalCircuit to the caller ( it can be called
         * once ); an unclaimed circuit dies with the factory.
         */
        class NgspiceCircuitFactory
        {
            const string            mSource;      // for error messages

            //! unique_ptr for exception safety ( a refusal mid-build must
            //! not leak the half-built circuit ) and so the class is not
            //! copyable -- a raw owning pointer would double-delete on the
            //! implicit copy. One-off setup allocation, allowed smart use
            std::unique_ptr< ElectricalCircuit > mCircuit;

            //! canonical node name ( "gnd" folded onto "0" ) -> BELFEM index
            Map< string, index_t >  mNodeMap;
            uint                    mNumberOfNodes = 0;

//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            /**
             * build from a lexed netlist
             */
            NgspiceCircuitFactory( const NetlistParser & aNetlist );

            /**
             * convenience: lex the file, then build
             */
            NgspiceCircuitFactory( const string & aPath );

            ~NgspiceCircuitFactory();

//-----------------------------------------------------------------------------

            /**
             * hand the circuit over to the caller, who owns it from here;
             * hard-errors on a second call
             */
            ElectricalCircuit *
            circuit();

            /**
             * resolve a netlist node name ( case-folded; ground aliases
             * "0"/"gnd" allowed ) to its BELFEM index -- the O5 lookup for
             * the deck-side terminal-pair and output readers. Hard-errors
             * on a name the netlist does not contain.
             */
            index_t
            node_index( const string & aName ) const;

            uint
            number_of_nodes() const
            {
                return mNumberOfNodes;
            }

//-----------------------------------------------------------------------------
        private:
//-----------------------------------------------------------------------------

            void
            build( const NetlistParser & aNetlist );

            void
            build_node_map( const NetlistParser & aNetlist );

            //! order directives, checked against the netlist's L/C cards
            void
            collect_orders( const NetlistParser & aNetlist,
                            Map< string, uint > & aOrders ) const;

            void
            create_element( const NetlistElement & aElement,
                            const NetlistParser & aNetlist,
                            const Map< string, uint > & aOrders );

            void
            create_from_directive( const NetlistDirective & aDirective );

            //! the single R/C/L value: positional or keyword, exactly one
            real
            rcl_value( const NetlistElement & aElement,
                       const char * aKeyword ) const;

            //! build the SourceFunction for a V/I card ( the created
            //! source component takes ownership )
            SourceFunction *
            source_function( const NetlistElement & aElement ) const;

            index_t
            node_index_checked( const string & aName,
                                const string & aCard,
                                const index_t aLine ) const;

            //! true if aName is a ground alias ( "0" or "gnd" )
            static bool
            is_ground( const string & aName );

//-----------------------------------------------------------------------------
        };

//-----------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_NGSPICECIRCUITFACTORY_HPP
