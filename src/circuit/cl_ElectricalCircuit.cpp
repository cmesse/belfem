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

#include "cl_ElectricalCircuit.hpp"
#include "cl_Component.hpp"
#include "cl_Resistor.hpp"
#include "cl_Inductor.hpp"
#include "cl_Capacitor.hpp"
#include "cl_VoltageSource.hpp"
#include "cl_CurrentSource.hpp"
#include "cl_Switch.hpp"
#include "cl_Diode.hpp"
#include "cl_Superconductor.hpp"
#include "fn_max.hpp"
#include "fn_norm.hpp"
#include "fn_sprint.hpp"

namespace belfem
{
    namespace electronics
    {
//----------------------------------------------------------------------------

        ElectricalCircuit::ElectricalCircuit( const uint aNumberOfNodes ) :
                mNumberOfNodes(aNumberOfNodes)
        {
            //Initialize the nodes
            mNodes = Cell< ElectricNode *>(mNumberOfNodes, nullptr) ;

            //Set the indices of the nodes
            for (uint j = 0; j < mNumberOfNodes ; ++j)
            {
                mNodes(j) = new ElectricNode() ;
                mNodes(j)->set_index(j) ;
                if (j != mNumberOfNodes-1)
                {
                    mVertices.push(mNodes(j)) ;
                }
            }
        }

        ElectricalCircuit::~ElectricalCircuit()
        {
            for ( ElectricNode * tNode : mNodes)
            {
                delete tNode;
            }
            mNodes.clear() ;
            mNumberOfNodes = 0 ;

            for ( Component * tComponent : mComponents )
            {
                delete tComponent;
            }
            mComponents.clear() ;
            mNumberOfComponents = 0 ;

            if(mJ != nullptr)
            {
                delete mJ ;
            }

            if(mMNA != nullptr)
            {
                delete mMNA ;
            }
        }

//----------------------------------------------------------------------------

        real
        ElectricalCircuit::get_voltage_on_node( const index_t aIndex) const
        {
            return mNodes(aIndex)->get_voltage() ;
        }

//-----------------------------------------------------------------------------

        real
        ElectricalCircuit::get_current_on_component( const index_t aIndex) const
        {
            return mComponents(aIndex)->get_current() ;
        }

//-----------------------------------------------------------------------------

        Component *
        ElectricalCircuit::component( const index_t aIndex) const
        {
            return mComponents(aIndex)  ;
        }

//-----------------------------------------------------------------------------

        Cell < FEMTwoTerminals * > &
        ElectricalCircuit::terminal_pairs( )
        {
            return mTerminalPairs ;
        }

//-----------------------------------------------------------------------------

        uint
        ElectricalCircuit::number_of_nodes() const
        {
            return mNumberOfNodes ;
        }

//-----------------------------------------------------------------------------

        uint
        ElectricalCircuit::number_of_components() const
        {
            return mNumberOfComponents ;
        }

//-----------------------------------------------------------------------------

        uint
        ElectricalCircuit::number_of_graph_vertices() const
        {
            return mVertices.size() ;
        }

//-----------------------------------------------------------------------------

        uint
        ElectricalCircuit::degree_of_graph_vertex( const index_t aIndex ) const
        {
            return mVertices( aIndex )->number_of_vertices() ;
        }

//-----------------------------------------------------------------------------

        void
        ElectricalCircuit::set_timestep( const real aDeltaTime )
        {
            mDeltaTime = aDeltaTime ;

            //Also update the timestep for the required components ;
            for(Component * tComponent : mComponents)
            {
                if(tComponent->component_type() == ComponentType::INDUCTOR ||
                   tComponent->component_type() == ComponentType::CAPACITOR )
                {
                    tComponent->set_timestep(aDeltaTime) ;
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        ElectricalCircuit::set_file( const string aFile )
        {
            mOutputFile = aFile ;
        }

//-----------------------------------------------------------------------------

        void
        ElectricalCircuit::set_output_currents( const Cell < string > aCurrents )
        {
            mOutputCurrents = aCurrents ;
        }

//-----------------------------------------------------------------------------

        void
        ElectricalCircuit::set_output_voltages( const Vector < id_t > aVoltages )
        {
            mOutputVoltages = aVoltages ;
        }

//-----------------------------------------------------------------------------

        void
        ElectricalCircuit::set_omega( const real aOmega )
        {
            mOmega = aOmega ;
        }

//-----------------------------------------------------------------------------


        void
        ElectricalCircuit::shift()
        {

            //Set the previous solution vector
            mPrevX = mX ;

            //update the time
            mTime += mDeltaTime ;

            //shift the required components
            for(Component * tComponent : mComponents)
            {
                if(tComponent->component_type() == ComponentType::INDUCTOR ||
                   tComponent->component_type() == ComponentType::CAPACITOR ||
                   tComponent->component_type() == ComponentType::TERMINALPAIR ||
                   tComponent->component_type() == ComponentType::CURRENTSOURCE ||
                   tComponent->component_type() == ComponentType::VOLTAGESOURCE ||
                   tComponent->component_type() == ComponentType::SWITCH)
                {
                    tComponent->shift(mTime, mDeltaTime) ;
                }
                /*else
                {
                    tComponent->set_current(0) ;
                }*/
            }
        }

//-----------------------------------------------------------------------------


        void
        ElectricalCircuit::shift_back()
        {

            //Return to the previous solution
            mX = mPrevX ;

            //Return to the prevous time step
            mTime -= mDeltaTime ;

            this->update_components();

            // Reset the currents to the other components
            for(Component * tComponent : mComponents)
            {
                switch ( tComponent->component_type() )
                {
                    case(ComponentType::RESISTOR) :
                    case(ComponentType::DIODE):
                    case(ComponentType::SUPERCONDUCTOR):
                    {
                        tComponent->compute_current() ;
                        break ;
                    }
                    case(ComponentType::INDUCTOR) :
                    case(ComponentType::CAPACITOR) :
                    case(ComponentType::TERMINALPAIR) :
                    {
                        tComponent->shift_back() ;
                        tComponent->compute_current() ;
                        break ;
                    }
                    case(ComponentType::SWITCH) :
                    {
                        // undo a latch that fired during the rejected attempt;
                        // the switch current is an unknown-current dof already
                        // restored by update_components()
                        tComponent->shift_back() ;
                        break ;
                    }
                    case(ComponentType::UNDEFINED) :
                    {
                        BELFEM_ERROR(false, "Undefined component in the circuit") ;
                        break ;
                    }

                    default:
                        break ;
                }
            }
        }

        void
        ElectricalCircuit::update_components()
        {
            // Reset the voltages at each node
            for (uint j = 0 ; j < mNumberOfNodes-1 ; ++j)
            {
                mNodes(j)->set_voltage(mX(j)) ;
            }

            // Reset the unknown currents to the corresponding components
            for (uint j = mNumberOfNodes-1 ; j < mNumberOfNodes+mNumberOfUnknownCurrents-1; ++j)
            {
                mComponentsUnknown(j-(mNumberOfNodes-1))->set_current(mX(j)) ;
            }

        }

//-----------------------------------------------------------------------------

        void
        ElectricalCircuit::create_resistor( const real aValue, const index_t aNIndex1, const index_t aNIndex2, string aLabel)
        {
            Cell < ElectricNode* > tNodes = {mNodes(aNIndex1), mNodes(aNIndex2)} ;
            mComponents.push(new Resistor(aValue,tNodes, aLabel)) ;
            ++mNumberOfComponents ;
        }

//-----------------------------------------------------------------------------

        void
        ElectricalCircuit::create_inductor( const real aValue, const uint aOrder, const index_t aNIndex1, const index_t aNIndex2, string aLabel)
        {
            Cell < ElectricNode* > tNodes = {mNodes(aNIndex1), mNodes(aNIndex2)} ;
            mComponents.push(new Inductor(aValue, aOrder,tNodes, aLabel)) ;
            ++mNumberOfComponents ;
        }

//-----------------------------------------------------------------------------

        void
        ElectricalCircuit::create_capacitor( const real aValue, const uint aOrder, const index_t aNIndex1, const index_t aNIndex2, string aLabel)
        {
            Cell < ElectricNode* > tNodes = {mNodes(aNIndex1), mNodes(aNIndex2)} ;
            mComponents.push(new Capacitor(aValue, aOrder,tNodes, aLabel)) ;
            ++mNumberOfComponents ;
        }

//-----------------------------------------------------------------------------

        void
        ElectricalCircuit::create_voltage_source( SourceFunction * aFunction, const index_t aNIndex1, const index_t aNIndex2, string aLabel )
        {
            Cell < ElectricNode* > tNodes = {mNodes(aNIndex1), mNodes(aNIndex2)} ;
            VoltageSource * tVSource = new VoltageSource(aFunction,tNodes, aLabel) ;
            mComponents.push(tVSource) ;

            //The voltage source is also an unknown current value
            tVSource->set_index(mNumberOfNodes+mNumberOfUnknownCurrents-1) ;
            mVertices.push(tVSource) ;
            mComponentsUnknown.push(tVSource) ;

            //Increment the number of components and unknown currents
            ++mNumberOfComponents ;
            ++mNumberOfUnknownCurrents ;
        }

//-----------------------------------------------------------------------------

        void
        ElectricalCircuit::create_current_source( SourceFunction * aFunction, const index_t aNIndex1, const index_t aNIndex2, string aLabel )
        {
            Cell < ElectricNode* > tNodes = {mNodes(aNIndex1), mNodes(aNIndex2)} ;
            mComponents.push(new CurrentSource(aFunction,tNodes, aLabel)) ;
            ++mNumberOfComponents ;
        }

//-----------------------------------------------------------------------------

        void
        ElectricalCircuit::create_switch( const bool aIsClosed, const real aTimeSwitch, const index_t aNIndex1, const index_t aNIndex2, string aLabel )
        {
            Cell < ElectricNode* > tNodes = {mNodes(aNIndex1), mNodes(aNIndex2)} ;
            Switch * tSwitch = new Switch(aIsClosed, aTimeSwitch, tNodes, aLabel) ;
            mComponents.push(tSwitch) ;

            //The switch is also an unknown current value
            tSwitch->set_index(mNumberOfNodes+mNumberOfUnknownCurrents-1) ;
            mVertices.push(tSwitch) ;
            mComponentsUnknown.push(tSwitch) ;

            //Increment the number of components and unknown currents
            ++mNumberOfComponents ;
            ++mNumberOfUnknownCurrents ;
        }

//-----------------------------------------------------------------------------

        void
        ElectricalCircuit::create_diode( const real aIs, const real aVt, const index_t aNIndex1, const index_t aNIndex2, string aLabel )
        {
            Cell < ElectricNode* > tNodes = {mNodes(aNIndex1), mNodes(aNIndex2)} ;
            mComponents.push(new Diode(aIs, aVt,tNodes, aLabel)) ;
            ++mNumberOfComponents ;
        }

//-----------------------------------------------------------------------------

        void
        ElectricalCircuit::create_superconductor( const real aIc, const real aN, const real aEc, const real aLength, const index_t aNIndex1, const index_t aNIndex2, string aLabel )
        {
            Cell < ElectricNode* > tNodes = {mNodes(aNIndex1), mNodes(aNIndex2)} ;
            mComponents.push(new Superconductor(aIc, aN, aEc, aLength,tNodes, aLabel)) ;
            ++mNumberOfComponents ;
        }

//-----------------------------------------------------------------------------

        void
        ElectricalCircuit::create_terminal_pair( const index_t aNIndex1, const index_t aNIndex2, const string aLabel )
        {
            Cell < ElectricNode* > tNodes = {mNodes(aNIndex1), mNodes(aNIndex2)} ;
            FEMTwoTerminals * tFEMTwoTerminals = new FEMTwoTerminals( tNodes, aLabel ) ;
            mComponents.push(tFEMTwoTerminals) ;

            //The voltage source is also an unknown current value
            //tFEMTwoTerminals->set_index(mNumberOfNodes+mNumberOfUnknownCurrents-1) ;
            //mVertices.push(tFEMTwoTerminals) ;
            //mComponentsUnknown.push(tFEMTwoTerminals) ;

            //Increment the number of components and unknown currents
            ++mNumberOfComponents ;
            //++mNumberOfUnknownCurrents ;

            mTerminalPairs.push(tFEMTwoTerminals) ;
        }

//-----------------------------------------------------------------------------

        void
        ElectricalCircuit::compute_adjacency()
        {
            Vector < uint > tNumComponentsPerNode = Vector < uint >(mNumberOfNodes-1,0) ;
            // Loop over all components to determine how many components are connected to each node (not counting the ground)
            for(uint c = 0; c < mComponents.size(); ++c)
            {
                for(uint n = 0; n < mComponents(c)->number_of_terminals(); ++n)
                {
                    //Update the number of components per node (except the ground node)
                    if (mComponents(c)->node(n)->index() < mNumberOfNodes-1)
                    {
                        tNumComponentsPerNode(mComponents(c)->node(n)->index()) += 1 ;
                    }
                }
            }

            //Allocate memory for storing which elements are connected to each node
            Matrix< uint > tElementsPerNode = Matrix< uint >(mNumberOfNodes-1, belfem::max(tNumComponentsPerNode));
            tNumComponentsPerNode.fill(0) ;

            // Populate the ElementsPerNode matrix by mapping nodes to their connected elements
            for(uint c = 0; c < mComponents.size(); ++c)
            {
                for(uint n = 0; n < mComponents(c)->number_of_terminals(); ++n)
                {
                    //Update the number of components per node (except the ground node)
                    if (mComponents(c)->node(n)->index() < mNumberOfNodes-1)
                    {
                        uint j = mComponents(c)->node(n)->index() ;
                        tNumComponentsPerNode(j) += 1 ;
                        tElementsPerNode(j,tNumComponentsPerNode(j)-1) = c;
                    }
                }
            }

            // Build the exact neighbor set each node receives. This is the
            // single source for both the allocation below and the population
            // that follows: the two used to filter differently ( a denylist
            // here, an allowlist there ), and neither counted the node-side
            // insert of an unknown-current vertex marked (*) below. A source
            // or switch tied to ground contributed nothing to the old size
            // yet still inserted its branch vertex, which is the overflow.
            Cell< Cell< index_t > > tNodeNeighbors( mNumberOfNodes-1, Cell< index_t >() ) ;

            for(uint n = 0 ; n < mNumberOfNodes-1; ++n)
            {
                Cell< index_t > & tNeighbors = tNodeNeighbors( n ) ;

                for(uint c = 0 ; c < tNumComponentsPerNode(n); ++c)
                {
                    uint d = tElementsPerNode(n, c) ;

                    for(uint j = 0 ; j < mComponents(d)->number_of_terminals(); ++j)
                    {
                        if (mComponents(d)->node(j)->index() != mNumberOfNodes-1 &&
                            mComponents(d)->node(j)->index() != n &&
                            (mComponents(d)->component_type()==ComponentType::RESISTOR ||
                             mComponents(d)->component_type()==ComponentType::CAPACITOR||
                             mComponents(d)->component_type()==ComponentType::INDUCTOR ||
                             mComponents(d)->component_type()==ComponentType::TERMINALPAIR ||
                             mComponents(d)->component_type()==ComponentType::DIODE ||
                             mComponents(d)->component_type()==ComponentType::SUPERCONDUCTOR))
                        {
                            tNeighbors.push(mComponents(d)->node(j)->index()) ;
                        }
                    }
                }

                //Remove duplicate nodes (when components are in parallel
                unique( tNeighbors ) ;
            }

            // Count the unknown-current vertices each node receives at (*).
            // These are NODE COUNTS ONLY -- never fold the branch vertices into
            // tNodeNeighbors: the first voltage source carries index
            // mNumberOfNodes-1, the same number as ground, so mNodes( index )
            // would silently resolve to the ground node. Iterating
            // mComponentsUnknown rather than testing component types keeps this
            // matched to (*) if that list ever gains a type.
            Cell< uint > tNumUnknownPerNode( mNumberOfNodes-1, 0 ) ;

            for (uint u = 0 ; u < mNumberOfUnknownCurrents; ++u)
            {
                Component * tComponent = mComponentsUnknown(u) ;

                for (uint j = 0; j < tComponent->number_of_terminals(); ++j)
                {
                    if(tComponent->node(j)->index() < mNumberOfNodes-1)
                    {
                        // not unique()d on purpose: two sources on one node are
                        // two distinct branch dofs and each gets its own slot
                        tNumUnknownPerNode(tComponent->node(j)->index()) += 1 ;
                    }
                }
            }

            // Initialize the vertex container on each nodes, sized to exactly
            // what gets inserted: the node itself, its unique neighbors, and
            // one slot per unknown-current vertex
            for(uint n = 0 ; n < mNumberOfNodes-1; ++n)
            {
                mNodes(n)->init_vertex_container(
                        1 + tNodeNeighbors(n).size() + tNumUnknownPerNode(n) ) ;
            }

            // Populate the vertices for each node from the same sets that sized
            // the containers above
            for(uint n = 0 ; n < mNumberOfNodes-1; ++n)
            {
                mNodes(n)->insert_vertex(mNodes(n)) ;

                for (uint c = 0 ; c < tNodeNeighbors(n).size(); ++c)
                {
                    mNodes(n)->insert_vertex(mNodes(tNodeNeighbors(n)(c))) ;
                }
            }

            // Also add the adjacency for the unknown current components

            //Initialize the Number of Nodes per Unknown current components
            Vector < uint > tNumNodesPerUnknown = Vector < uint >(mNumberOfUnknownCurrents,0) ;
            for (uint n = mNumberOfNodes-1; n < mNumberOfNodes+mNumberOfUnknownCurrents-1 ; ++n)
            {
                Component * tComponent = mComponentsUnknown(n-(mNumberOfNodes-1));
                for (uint j = 0; j < tComponent->number_of_terminals(); ++j)
                {
                    if(tComponent->node(j)->index() < mNumberOfNodes-1)
                    {
                        tNumNodesPerUnknown(n-(mNumberOfNodes-1)) += 1 ;
                    }
                }

                // a switch also inserts itself as the diagonal term below, and
                // the graph is built once, so that slot must always be reserved
                if (tComponent->component_type() == ComponentType::SWITCH)
                {
                    tNumNodesPerUnknown(n-(mNumberOfNodes-1)) += 1 ;
                }
            }

            //Initialize the vertex container for the components with unknown currents
            for (uint n = mNumberOfNodes-1; n < mNumberOfNodes+mNumberOfUnknownCurrents-1 ; ++n)
            {
                mVertices(n)->init_vertex_container(tNumNodesPerUnknown(n-(mNumberOfNodes-1))) ;
            }

            //Populate the vertices with components with unknown currents
            for (uint n = mNumberOfNodes-1; n < mNumberOfNodes+mNumberOfUnknownCurrents-1 ; ++n)
            {
                Component * tComponent = mComponentsUnknown(n-(mNumberOfNodes-1));
                for (uint j = 0; j < tComponent->number_of_terminals(); ++j)
                {
                    if(tComponent->node(j)->index() < mNumberOfNodes-1)
                    {
                        mVertices(n)->insert_vertex(tComponent->node(j)) ;
                        mNodes(tComponent->node(j)->index())->insert_vertex(mVertices(n)) ; // (*)
                    }
                }

                // Add also the switch to its own vertices, as the diagonal term
                if (tComponent->component_type() == ComponentType::SWITCH)
                {
                    mVertices(n)->insert_vertex(mVertices(n)) ;
                }
            }

            // Vertex holds no capacity member, so insert_vertex cannot check
            // itself ( its assertion guards the counter WIDTH, not the buffer ).
            // Assert the drift here instead, where both numbers are in scope.
            for(uint n = 0 ; n < mNumberOfNodes-1; ++n)
            {
                BELFEM_ASSERT( mNodes(n)->number_of_vertices()
                                   == 1 + tNodeNeighbors(n).size() + tNumUnknownPerNode(n),
                               "Node %lu received %u vertices but was sized for %lu",
                               ( long unsigned int ) n,
                               ( unsigned int ) mNodes(n)->number_of_vertices(),
                               ( long unsigned int ) ( 1 + tNodeNeighbors(n).size()
                                                         + tNumUnknownPerNode(n) ) );
            }
            for (uint n = mNumberOfNodes-1; n < mNumberOfNodes+mNumberOfUnknownCurrents-1 ; ++n)
            {
                BELFEM_ASSERT( mVertices(n)->number_of_vertices()
                                   == tNumNodesPerUnknown(n-(mNumberOfNodes-1)),
                               "Branch vertex %lu received %u vertices but was sized for %u",
                               ( long unsigned int ) n,
                               ( unsigned int ) mVertices(n)->number_of_vertices(),
                               ( unsigned int ) tNumNodesPerUnknown(n-(mNumberOfNodes-1)) );
            }
        }

//-----------------------------------------------------------------------------

        void
        ElectricalCircuit::compute_MNA_matrix()
        {
            // Compute the adjacency to initialize the sparse matrix
            if (mMNA == nullptr)
            {
                this->compute_adjacency() ;
                // the delete keeps a retry from leaking a prior mJ; the nulling
                // keeps that same delete from leaving a dangling pointer if the
                // allocation below throws, which the destructor would then free
                // a second time
                if ( mJ != nullptr ) { delete mJ ; mJ = nullptr ; }
                mJ = new SpMatrix(mVertices) ;

                // we already know that mMNA is nullptr, so we don't need to check again
                mMNA = new SpMatrix(mVertices) ;
                mRHS = Vector< real > (mNumberOfNodes+mNumberOfUnknownCurrents-1, 0.0) ;

                // load_state() may already have restored mX/mPrevX from a memdump
                // before the first stamp — only a cold start zero-initializes them
                uint tNumDofs = mNumberOfNodes+mNumberOfUnknownCurrents-1 ;
                if ( mX.length() != tNumDofs || mPrevX.length() != tNumDofs )
                {
                    mX = Vector< real > (tNumDofs, 0.0) ;
                    mPrevX = Vector< real > (tNumDofs, 0.0) ;
                }
            }
            mMNA->fill(0.0) ;

            uint tCount = 0; //Count for number of unknown current components.

            // Populate the matrix and RHS depending on the components and their connectivities
            for(Component * tComponent : mComponents)
            {
                switch ( tComponent->component_type() )
                {
                    case(ComponentType::RESISTOR) :
                    {
                        if(tComponent->get_node_plus()->index() != mNumberOfNodes-1 &&
                           tComponent->get_node_minus()->index() == mNumberOfNodes-1)
                        {
                            mMNA->operator()(tComponent->get_node_plus()->index(),
                                             tComponent->get_node_plus()->index()) += 1.0/tComponent->get_value();
                        }
                        else if(tComponent->get_node_minus()->index() != mNumberOfNodes-1 &&
                                tComponent->get_node_plus()->index() == mNumberOfNodes-1)
                        {
                            mMNA->operator()(tComponent->get_node_minus()->index(),
                                             tComponent->get_node_minus()->index()) += 1.0/tComponent->get_value();
                        }
                        else
                        {
                            mMNA->operator()( tComponent->get_node_plus()->index(),
                                              tComponent->get_node_plus()->index()) += 1.0 / tComponent->get_value();
                            mMNA->operator()( tComponent->get_node_minus()->index(),
                                              tComponent->get_node_minus()->index()) += 1.0 / tComponent->get_value();
                            mMNA->operator()( tComponent->get_node_plus()->index(),
                                              tComponent->get_node_minus()->index()) -= 1.0 / tComponent->get_value();
                            mMNA->operator()( tComponent->get_node_minus()->index(),
                                              tComponent->get_node_plus()->index()) -= 1.0 / tComponent->get_value();
                        }
                        break ;
                    }

                    case(ComponentType::INDUCTOR) :
                    case(ComponentType::CAPACITOR) :
                    {
                        if(tComponent->get_node_plus()->index() != mNumberOfNodes-1 &&
                           tComponent->get_node_minus()->index() == mNumberOfNodes-1)
                        {
                            mMNA->operator()(tComponent->get_node_plus()->index(),
                                             tComponent->get_node_plus()->index()) += 1.0/tComponent->get_discretized_resistance();
                        }
                        else if(tComponent->get_node_minus()->index() != mNumberOfNodes-1 &&
                                tComponent->get_node_plus()->index() == mNumberOfNodes-1)
                        {
                            mMNA->operator()(tComponent->get_node_minus()->index(),
                                             tComponent->get_node_minus()->index()) += 1.0/tComponent->get_discretized_resistance();
                        }
                        else
                        {
                            mMNA->operator()( tComponent->get_node_plus()->index(),
                                              tComponent->get_node_plus()->index()) += 1.0 / tComponent->get_discretized_resistance();
                            mMNA->operator()( tComponent->get_node_minus()->index(),
                                              tComponent->get_node_minus()->index()) += 1.0 / tComponent->get_discretized_resistance();
                            mMNA->operator()( tComponent->get_node_plus()->index(),
                                              tComponent->get_node_minus()->index()) -= 1.0 / tComponent->get_discretized_resistance();
                            mMNA->operator()( tComponent->get_node_minus()->index(),
                                              tComponent->get_node_plus()->index()) -= 1.0 / tComponent->get_discretized_resistance();
                        }
                        break ;
                    }

                    case(ComponentType::VOLTAGESOURCE) :
                    {
                        if(tComponent->get_node_plus()->index() != mNumberOfNodes-1)
                        {
                            mMNA->operator()(tComponent->get_node_plus()->index(),
                                             tCount+mNumberOfNodes-1) += 1.0;
                            mMNA->operator()(tCount+mNumberOfNodes-1,
                                             tComponent->get_node_plus()->index()) += 1.0 ;
                        }
                        if(tComponent->get_node_minus()->index() != mNumberOfNodes-1)
                        {
                            mMNA->operator()(tComponent->get_node_minus()->index(),
                                             tCount+mNumberOfNodes-1) -= 1.0;
                            mMNA->operator()(tCount+mNumberOfNodes-1,
                                             tComponent->get_node_minus()->index()) -= 1.0 ;
                        }
                        ++tCount ;
                        break ;
                    }

                    case(ComponentType::SWITCH) :
                    {
                        // the switch occupies an unknown-current dof; its stamps are
                        // voltage-state dependent and live in compute_jacobian_and_rhs,
                        // but the shared dof counter must advance here too
                        ++tCount ;
                        break ;
                    }

                    case(ComponentType::UNDEFINED) :
                    {
                        BELFEM_ERROR(false, "Undefined component in the circuit") ;
                        break ;
                    }

                    default :
                    {
                        break ;
                    }
                }
            }

        }

//-----------------------------------------------------------------------------

        void
        ElectricalCircuit::compute_jacobian_and_rhs()
        {
            mRHS.fill(0.0) ;
            mJ->fill(0.0) ;

            uint tCount = 0; //Count for number of unknown current components.

            // Populate the matrix and RHS depending on the components and their connectivities
            for(Component * tComponent : mComponents)
            {
                switch ( tComponent->component_type() )
                {
                    case(ComponentType::RESISTOR) :
                    {
                        if(tComponent->get_node_plus()->index() != mNumberOfNodes-1 &&
                           tComponent->get_node_minus()->index() == mNumberOfNodes-1)
                        {
                            mRHS(tComponent->get_node_plus()->index()) += tComponent->get_node_plus()->get_voltage()/tComponent->get_value() ;
                        }
                        else if(tComponent->get_node_minus()->index() != mNumberOfNodes-1 &&
                                tComponent->get_node_plus()->index() == mNumberOfNodes-1)
                        {
                            mRHS(tComponent->get_node_minus()->index()) += tComponent->get_node_minus()->get_voltage()/tComponent->get_value() ;
                        }
                        else
                        {
                            mRHS(tComponent->get_node_plus()->index()) += tComponent->get_node_plus()->get_voltage()/tComponent->get_value() ;
                            mRHS(tComponent->get_node_minus()->index()) += tComponent->get_node_minus()->get_voltage()/tComponent->get_value() ;
                            mRHS(tComponent->get_node_plus()->index()) -= tComponent->get_node_minus()->get_voltage()/tComponent->get_value() ;
                            mRHS(tComponent->get_node_minus()->index()) -= tComponent->get_node_plus()->get_voltage()/tComponent->get_value() ;
                        }
                        break ;
                    }

                    case(ComponentType::INDUCTOR) :
                    case(ComponentType::CAPACITOR) :
                    {
                        if(tComponent->get_node_plus()->index() != mNumberOfNodes-1 &&
                           tComponent->get_node_minus()->index() == mNumberOfNodes-1)
                        {
                            mRHS(tComponent->get_node_plus()->index()) += tComponent->get_node_plus()->get_voltage()/tComponent->get_discretized_resistance() ;
                        }
                        else if(tComponent->get_node_minus()->index() != mNumberOfNodes-1 &&
                                tComponent->get_node_plus()->index() == mNumberOfNodes-1)
                        {
                            mRHS(tComponent->get_node_minus()->index()) += tComponent->get_node_minus()->get_voltage()/tComponent->get_discretized_resistance() ;
                        }
                        else
                        {
                            mRHS(tComponent->get_node_plus()->index()) += tComponent->get_node_plus()->get_voltage()/tComponent->get_discretized_resistance() ;
                            mRHS(tComponent->get_node_minus()->index()) += tComponent->get_node_minus()->get_voltage()/tComponent->get_discretized_resistance() ;
                            mRHS(tComponent->get_node_plus()->index()) -= tComponent->get_node_minus()->get_voltage()/tComponent->get_discretized_resistance() ;
                            mRHS(tComponent->get_node_minus()->index()) -= tComponent->get_node_plus()->get_voltage()/tComponent->get_discretized_resistance() ;
                        }

                        if(tComponent->get_node_plus()->index() != mNumberOfNodes-1)
                        {
                            mRHS(tComponent->get_node_plus()->index()) += tComponent->get_discretized_source();
                        }

                        if(tComponent->get_node_minus()->index() != mNumberOfNodes-1)
                        {
                            mRHS(tComponent->get_node_minus()->index()) -= tComponent->get_discretized_source();
                        }
                        break ;
                    }

                    case(ComponentType::VOLTAGESOURCE) :
                    {
                        if(tComponent->get_node_plus()->index() != mNumberOfNodes-1)
                        {
                            mRHS(tCount+mNumberOfNodes-1) += tComponent->get_node_plus()->get_voltage() ;
                            mRHS(tComponent->get_node_plus()->index()) += tComponent->get_current() ;
                        }
                        if(tComponent->get_node_minus()->index() != mNumberOfNodes-1)
                        {
                            mRHS(tCount+mNumberOfNodes-1) -= tComponent->get_node_minus()->get_voltage() ;
                            mRHS(tComponent->get_node_minus()->index()) -= tComponent->get_current() ;
                        }
                        mRHS(tCount+mNumberOfNodes-1) -= tComponent->get_value() ;
                        ++tCount ;
                        break ;
                    }

                    case(ComponentType::CURRENTSOURCE) :
                    case(ComponentType::TERMINALPAIR) :
                    {
                        if(tComponent->get_node_plus()->index() != mNumberOfNodes-1)
                        {
                            mRHS(tComponent->get_node_plus()->index()) += tComponent->get_current();
                        }

                        if(tComponent->get_node_minus()->index() != mNumberOfNodes-1)
                        {
                            mRHS(tComponent->get_node_minus()->index()) -= tComponent->get_current();
                        }
                        break ;
                    }

                    case(ComponentType::SWITCH) :
                    {

                        if(tComponent->get_node_plus()->index() != mNumberOfNodes-1)
                        {
                            mJ->operator()( tComponent->get_node_plus()->index(),
                                            tCount + mNumberOfNodes - 1 ) += 1.0;
                            mRHS(tComponent->get_node_plus()->index()) += tComponent->get_current() ;

                            // Switch is closed, meaning that we set the voltage difference at 0
                            if (tComponent->is_closed())
                            {
                                mJ->operator()(tCount+mNumberOfNodes-1,
                                               tComponent->get_node_plus()->index()) += 1.0 ;
                                mRHS(tCount+mNumberOfNodes-1) += tComponent->get_node_plus()->get_voltage() ;
                            }
                        }
                        if (tComponent->get_node_minus()->index() != mNumberOfNodes-1)
                        {
                            mJ->operator()( tComponent->get_node_minus()->index(),
                                            tCount + mNumberOfNodes - 1 ) -= 1.0;
                            mRHS(tComponent->get_node_minus()->index()) -= tComponent->get_current() ;

                            // Switch is closed, meaning that we set the voltage difference at 0
                            if (tComponent->is_closed())
                            {
                                mJ->operator()(tCount+mNumberOfNodes-1,
                                               tComponent->get_node_minus()->index()) -= 1.0 ;
                                mRHS(tCount+mNumberOfNodes-1) -= tComponent->get_node_minus()->get_voltage() ;
                            }
                        }

                        // Switch is open, meaning that the current is set to 0
                        if ( !tComponent->is_closed() )
                        {
                            mJ->operator()(tCount+mNumberOfNodes-1,
                                           tCount+mNumberOfNodes-1) += 1.0 ;
                            mRHS(tCount+mNumberOfNodes-1) += tComponent->get_current() ;
                        }
                        ++tCount ;
                        break ;
                    }
                    case(ComponentType::DIODE) :
                    case(ComponentType::SUPERCONDUCTOR) :
                    {
                        if(tComponent->get_node_plus()->index() != mNumberOfNodes-1 &&
                           tComponent->get_node_minus()->index() == mNumberOfNodes-1)
                        {
                            mJ->operator()(tComponent->get_node_plus()->index(),
                                           tComponent->get_node_plus()->index()) += tComponent->dIdV();
                        }
                        else if(tComponent->get_node_minus()->index() != mNumberOfNodes-1 &&
                                tComponent->get_node_plus()->index() == mNumberOfNodes-1)
                        {
                            mJ->operator()(tComponent->get_node_minus()->index(),
                                           tComponent->get_node_minus()->index()) += tComponent->dIdV();
                        }
                        else
                        {
                            mJ->operator()( tComponent->get_node_plus()->index(),
                                            tComponent->get_node_plus()->index()) += tComponent->dIdV();
                            mJ->operator()( tComponent->get_node_minus()->index(),
                                            tComponent->get_node_minus()->index()) += tComponent->dIdV();
                            mJ->operator()( tComponent->get_node_plus()->index(),
                                            tComponent->get_node_minus()->index()) -= tComponent->dIdV();
                            mJ->operator()( tComponent->get_node_minus()->index(),
                                            tComponent->get_node_plus()->index()) -= tComponent->dIdV();

                        }

                        if(tComponent->get_node_plus()->index() != mNumberOfNodes-1)
                        {
                            mRHS(tComponent->get_node_plus()->index()) += tComponent->get_current();
                        }

                        if(tComponent->get_node_minus()->index() != mNumberOfNodes-1)
                        {
                            mRHS(tComponent->get_node_minus()->index()) -= tComponent->get_current();
                        }
                        break ;
                    }

                    case(ComponentType::UNDEFINED) :
                    {
                        BELFEM_ERROR(false, "Undefined component in the circuit") ;
                        break ;
                    }
                }
            }

            // Add the MNA matrix to the Jacobian
            for (uint i = 0; i < mJ->number_of_nonzeros(); ++i)
            {
                mJ->data(i)+=mMNA->data(i) ;
            }
        }

//-----------------------------------------------------------------------------

        void
        ElectricalCircuit::solve()
        {
            //Solve the matrix system and update the voltages and currents of the components
            Vector< real > tdX = Vector< real > (mNumberOfNodes+mNumberOfUnknownCurrents-1, 0.0) ;

            // compute the residual as r = A * x - b and write it into RHS vector
            //mJ->multiply( mX, mRHS, 1.0, -1.0 );

            mSolver.solve(*mJ,tdX,mRHS) ;

            mX-=(mOmega*tdX) ;

            this->update_components();

            // Set the currents to the other components
            for(Component * tComponent : mComponents)
            {
                switch ( tComponent->component_type() )
                {
                    case(ComponentType::RESISTOR) :
                    case(ComponentType::INDUCTOR) :
                    case(ComponentType::CAPACITOR) :
                    case(ComponentType::TERMINALPAIR) :
                    case(ComponentType::DIODE):
                    case(ComponentType::SUPERCONDUCTOR):
                    {
                        tComponent->compute_current() ;
                        break ;
                    }
                    case(ComponentType::UNDEFINED) :
                    {
                        BELFEM_ERROR(false, "Undefined component in the circuit") ;
                        break ;
                    }

                    default:
                        break ;
                }
            }

            //Compute the norm of the actual RHS (not the residual)
            Vector< real > tRHS = Vector< real > (mNumberOfNodes+mNumberOfUnknownCurrents-1, 0.0) ;
            mJ->multiply( mX, tRHS, 1.0, -1.0 );
            mRHSnorm = norm(tRHS) ;

        }

//-----------------------------------------------------------------------------

        real
        ElectricalCircuit::residual() const
        {
            return norm(mRHS)/mRHSnorm ;
        }

//-----------------------------------------------------------------------------

        void
        ElectricalCircuit::save_timestep()
        {
            std::ofstream tFile;
            tFile.open (mOutputFile, std::ios_base::app);

            tFile << mTime << " " ;

            for(string tString : mOutputCurrents)
            {
                uint tCount = 0 ;
                for(Component * tComp : mComponents)
                {
                    if (tComp->get_label() == tString)
                    {
                        tFile << tComp->get_current() << " " ;
                        ++tCount ;
                    }
                }
                BELFEM_ERROR(tCount > 0, "%s not defined in the circuit for output current", tString.c_str()) ;
            }

            for(index_t tInd : mOutputVoltages)
            {
                BELFEM_ERROR(tInd < mNodes.size(), "Node %d not defined in the circuit for output voltage", tInd) ;
                tFile << mNodes(tInd)->get_voltage() << " " ;
            }
            tFile << " \n";

            tFile.close() ;
        }

//-----------------------------------------------------------------------------

        void
        ElectricalCircuit::init_output_file()
        {
            std::ofstream tFile;
            tFile.open (mOutputFile);

            tFile << "t " ;
            for(string tString : mOutputCurrents)
            {
                tFile << "I" << tString << " " ;
            }

            for(index_t tInd : mOutputVoltages)
            {
                tFile << "V" << tInd << " " ;
            }
            tFile << " \n";

            tFile.close() ;
        }

        void
        ElectricalCircuit::save_state( hid_t aFile )
        {
#ifdef BELFEM_HDF5
            herr_t tStatus = 0 ;
            hdf5::save_scalar_to_file( aFile, "time", mTime, tStatus );
            hdf5::save_scalar_to_file( aFile, "delta_time", mDeltaTime, tStatus );
            hdf5::save_vector_to_file( aFile, "x", mX, tStatus );
            hdf5::save_vector_to_file( aFile, "prev_x", mPrevX, tStatus );

            // v2: per-component state ( histories, currents, latches ),
            // datasets prefixed by creation index
            uint tNumComponents = mComponents.size() ;
            hdf5::save_scalar_to_file( aFile, "n_components", tNumComponents, tStatus );
            for ( index_t c = 0; c < mComponents.size(); ++c )
            {
                string tPrefix = sprint( "c%03u_", ( unsigned int ) c );
                uint tType = ( uint ) mComponents(c)->component_type() ;
                hdf5::save_scalar_to_file( aFile, tPrefix + "type", tType, tStatus );
                mComponents(c)->save_state( aFile, tPrefix );
            }
#endif
        }

        void
        ElectricalCircuit::load_state( hid_t aFile )
        {
#ifdef BELFEM_HDF5
            herr_t tStatus = 0 ;

            // validate the file BEFORE touching any state, so a refused
            // restart leaves the circuit as constructed.
            // An old dump without n_components would silently cold-start the
            // L/C histories — refuse it instead
            BELFEM_ERROR( hdf5::dataset_exists( aFile, "n_components" ),
                "Circuit state in restart file has no per-component data ( pre-v2 format ). Delete the memdump to restart with a cold circuit." ) ;

            uint tNumComponents = 0 ;
            hdf5::load_scalar_from_file( aFile, "n_components", tNumComponents, tStatus );
            BELFEM_ERROR( tNumComponents == ( uint ) mComponents.size(),
                "Circuit state in restart file has %u components, the circuit has %lu. Delete the memdump to restart with a cold circuit.",
                tNumComponents,
                ( long unsigned int ) mComponents.size() ) ;

            for ( index_t c = 0; c < mComponents.size(); ++c )
            {
                uint tType = 0 ;
                hdf5::load_scalar_from_file( aFile,
                        sprint( "c%03u_type", ( unsigned int ) c ), tType, tStatus );
                BELFEM_ERROR( tType == ( uint ) mComponents(c)->component_type(),
                    "Component %lu in the restart file has a different type than the circuit ( component order changed? ). Delete the memdump to restart with a cold circuit.",
                    ( long unsigned int ) c ) ;
            }

            hdf5::load_scalar_from_file( aFile, "time", mTime, tStatus );
            hdf5::load_scalar_from_file( aFile, "delta_time", mDeltaTime, tStatus );
            hdf5::load_vector_from_file( aFile, "x", mX, tStatus );
            hdf5::load_vector_from_file( aFile, "prev_x", mPrevX, tStatus );

            // guard against a state vector from an edited circuit
            uint tNumDofs = mNumberOfNodes - 1 + mNumberOfUnknownCurrents ;
            BELFEM_ERROR( mX.length() == tNumDofs && mPrevX.length() == tNumDofs,
                "Circuit state in restart file does not match the circuit (expect %u dofs, file has x: %lu, prev_x: %lu). Delete the memdump to restart with a cold circuit.",
                tNumDofs,
                ( long unsigned int ) mX.length(),
                ( long unsigned int ) mPrevX.length() ) ;

            this->update_components();

            // seed the component timesteps and BDF1 companions from the
            // restored delta_time, so that a component whose dumped history
            // is EMPTY ( dump before the first shift ) is still coherent for
            // a stamp without a prior shift; non-empty loaders overwrite the
            // seed with the exact dumped-state companions
            this->set_timestep( mDeltaTime );

            for ( index_t c = 0; c < mComponents.size(); ++c )
            {
                mComponents(c)->load_state( aFile, sprint( "c%03u_", ( unsigned int ) c ) );
            }
#endif
        }

//-----------------------------------------------------------------------------
    }
}
