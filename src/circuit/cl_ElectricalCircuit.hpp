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

#ifndef CL_ELECTRICALCIRCUIT_HPP
#define CL_ELECTRICALCIRCUIT_HPP

#include "typedefs.hpp"
#include "cl_ElectricNode.hpp"
#include "cl_TwoTerminals.hpp"
#include "cl_FEMTwoTerminals.hpp"
#include "cl_Resistor.hpp"
#include "cl_SourceFunction.hpp"
#include "cl_Cell.hpp"
#include "cl_SpMatrix.hpp"
#include "cl_Vector.hpp"
#include "cl_Solver.hpp"
#include "cl_Circuit.hpp"

namespace belfem
{
    namespace electronics
    {
//-----------------------------------------------------------------------------
        /**
         * @brief The circuit itself: MNA assembly, Jacobian, solve and timestep.
         *
         * @ingroup grp_circuit
         * @see @ref circuit_circuit_usage_guide
         */
        class ElectricalCircuit : public Circuit
        {

            //! Nodes in the circuit
            Cell < ElectricNode * > mNodes ;

            //! Components in the circuit
            Cell < Component * > mComponents ;

            //! Components with unknown currents in the circuit
            Cell < Component * > mComponentsUnknown ;

            //! Components that are solved by a FEM problem
            Cell < FEMTwoTerminals* > mTerminalPairs ;

            //! Vertices in the circuit for matrix
            Cell < graph::Vertex* > mVertices ;

            //! Number of nodes
            uint mNumberOfNodes;

            //! Number of components
            uint mNumberOfComponents = 0;

            //! Number of components with currents as unknown
            uint mNumberOfUnknownCurrents = 0;

            //! Jacobian matrix
            SpMatrix * mJ = nullptr ;

            //! Matrix of the circuit
            SpMatrix * mMNA = nullptr ;

            //! Right-hand side / Newton residual vector
            Vector< real > mRHS;

            real mRHSnorm ;

            //! Solution vector
            Vector< real > mX;

            //! Previous solution
            Vector< real > mPrevX ;

            //! Solver
            Solver mSolver { SolverType::SUPERLU } ;

            //! List of component labels for current output
            Cell< string > mOutputCurrents ;

            //! List of the node voltages for output
            Vector< id_t > mOutputVoltages ;

            //! Output file
            string mOutputFile = "" ;

            //! Time step
            real mDeltaTime ;

            //! Current simulation time
            real mTime = 0.0 ;

            //! Relaxation factor
            real mOmega = 1.0;

//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            ElectricalCircuit( const uint aNumberOfNodes) ;

            ~ElectricalCircuit() override ;

            // NON-COPYABLE, NON-MOVABLE. This class owns every entry of
            // mNodes and mComponents plus the two matrices mJ and mMNA,
            // and deletes all of them in its destructor, so the implicit
            // copy would be a shallow pointer copy and the second
            // destructor a double free. Copy assignment is already
            // implicitly deleted through the by-value mSolver member,
            // whose own assignment is deleted -- these declarations state
            // the invariant instead of leaving it to that chain
            ElectricalCircuit( const ElectricalCircuit & ) = delete ;
            ElectricalCircuit( ElectricalCircuit && ) = delete ;
            ElectricalCircuit & operator=( const ElectricalCircuit & ) = delete ;
            ElectricalCircuit & operator=( ElectricalCircuit && ) = delete ;

//-----------------------------------------------------------------------------

            real
            get_voltage_on_node( const index_t aIndex) const ;

//-----------------------------------------------------------------------------

            real
            get_current_on_component( const index_t aIndex) const ;

//-----------------------------------------------------------------------------

            Component *
            component( const index_t aIndex) const ;

//-----------------------------------------------------------------------------

            Cell < FEMTwoTerminals * > &
            terminal_pairs( ) ;

//-----------------------------------------------------------------------------

            uint
            number_of_nodes() const ;

//-----------------------------------------------------------------------------

            uint
            number_of_components() const ;

//-----------------------------------------------------------------------------

            void
            set_timestep( const real aDeltaTime ) override ;

//-----------------------------------------------------------------------------

            void
            set_file( const string aFile ) ;

//-----------------------------------------------------------------------------

            void
            set_output_currents( const Cell < string > aCurrents ) ;

//-----------------------------------------------------------------------------

            void
            set_output_voltages( const Vector < id_t > aVoltages ) ;

//-----------------------------------------------------------------------------


            void
            set_omega( const real aOmega ) override ;

//-----------------------------------------------------------------------------

            void
            shift() override ;

//-----------------------------------------------------------------------------

            void
            shift_back() override ;

//-----------------------------------------------------------------------------


            void
            create_resistor( const real aValue, const index_t aNIndex1, const index_t aNIndex2, const string aLabel = "") ;

//-----------------------------------------------------------------------------

            void
            create_inductor( const real aValue, const uint aOrder, const index_t aNIndex1, const index_t aNIndex2, const string aLabel = "") ;

//-----------------------------------------------------------------------------

            void
            create_capacitor( const real aValue, const uint aOrder, const index_t aNIndex1, const index_t aNIndex2, const string aLabel = "") ;

//-----------------------------------------------------------------------------

            void
            create_voltage_source( SourceFunction * aFunction, const index_t aNIndex1, const index_t aNIndex2, const string aLabel = "" ) ;

//-----------------------------------------------------------------------------

            void
            create_current_source( SourceFunction * aFunction, const index_t aNIndex1, const index_t aNIndex2, const string aLabel = "" ) ;

//-----------------------------------------------------------------------------

            void
            create_switch( const bool aIsClosed , const real aTimeSwitch, const index_t aNIndex1, const index_t aNIndex2, const string aLabel = "" ) ;

//-----------------------------------------------------------------------------

            void
            create_diode( const real aIs, const real aVt, const index_t aNIndex1, const index_t aNIndex2, const string aLabel = "" ) ;

//-----------------------------------------------------------------------------

            void
            create_superconductor( const real aIc, const real aN, const real aEc, const real aLength, const index_t aNIndex1, const index_t aNIndex2, const string aLabel = "" ) ;

//-----------------------------------------------------------------------------

            void
            create_terminal_pair( const index_t aNIndex1, const index_t aNIndex2, const string aLabel = "" ) ;

//-----------------------------------------------------------------------------

            void
            compute_adjacency() ;

//-----------------------------------------------------------------------------

            /**
             * number of vertices in the adjacency graph handed to the SpMatrix:
             * the non-ground nodes, then one branch vertex per unknown-current
             * component. NOT number_of_nodes(), which counts ground and no
             * branches. Valid as soon as the circuit is built
             */
            uint
            number_of_graph_vertices() const ;

//-----------------------------------------------------------------------------

            /**
             * adjacency-list length of one vertex of that graph, i.e. how many
             * entries it contributes to the sparsity pattern. Zero until
             * compute_adjacency() has filled the graph
             */
            uint
            degree_of_graph_vertex( const index_t aIndex ) const ;

//-----------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs() override ;

//-----------------------------------------------------------------------------

            void
            compute_MNA_matrix() override ;

//-----------------------------------------------------------------------------

            void
            solve() override ;

//-----------------------------------------------------------------------------

            real
            residual() const override ;

//-----------------------------------------------------------------------------

            void
            save_timestep() override ;

//-----------------------------------------------------------------------------

            void
            init_output_file() ;

//-----------------------------------------------------------------------------

            real
            current( const index_t aIndex ) const override ;

//-----------------------------------------------------------------------------

            real
            voltage( const index_t aIndex ) const override ;

            void
            set_current_and_voltage( const index_t aIndex, const real aI, const real aV ) override ;

//-----------------------------------------------------------------------------

            void
            save_state( hid_t aFile ) override ;

            void
            load_state( hid_t aFile ) override ;

        private:

            void
            update_components();

//-----------------------------------------------------------------------------
        };

        inline real
        ElectricalCircuit::current( const index_t aIndex ) const
        {
            return mTerminalPairs( aIndex )->get_current() ;
        }

//-----------------------------------------------------------------------------

        inline real
        ElectricalCircuit::voltage( const index_t aIndex ) const
        {
            return mTerminalPairs( aIndex )->get_node_plus()->get_voltage()-
                   mTerminalPairs( aIndex )->get_node_minus()->get_voltage() ;
        }

//-----------------------------------------------------------------------------

        inline void
        ElectricalCircuit::set_current_and_voltage( const index_t aIndex, const real aI, const real aV )
        {
            mTerminalPairs( aIndex )->set_IV( aI, aV ) ;
        }
    }
}

#endif //CL_ELECTRICALCIRCUIT_HPP
