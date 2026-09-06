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

#ifndef BELFEM_CL_COMPONENT_HPP
#define BELFEM_CL_COMPONENT_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"
#include "cl_ElectricNode.hpp"
#include "Component_Enums.hpp"
#include "hdf5_types.hpp"

namespace belfem
{
    namespace electronics
    {
//-----------------------------------------------------------------------------
        /**
         * @brief Base class for every circuit component.
         *
         * @ingroup grp_circuit
         * @see @ref circuit_circuit_usage_guide
         */
        class Component : public graph::Vertex
        {
        protected:

            //! Number of terminals
            const uint mNumTerminals ;

            //! Label
            string mLabel ;

            //! Nodes at each terminal ( non-owning )
            ElectricNode ** mTerminals = nullptr ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Component( const uint aNumTerminals, Cell < ElectricNode* > & aTerminals, string aLabel = "") ;

            ~Component() override;

//------------------------------------------------------------------------------

            virtual ComponentType
            component_type() const ;

//-----------------------------------------------------------------------------

            uint
            number_of_terminals() const ;

//------------------------------------------------------------------------------

            ElectricNode *
            node( index_t tIndex ) const ;

//------------------------------------------------------------------------------

            virtual real
            get_value() const ;

//------------------------------------------------------------------------------

            virtual real
            dIdV() const ;

//------------------------------------------------------------------------------

            virtual real
            get_discretized_resistance() const ;

//------------------------------------------------------------------------------

            virtual real
            get_discretized_source() const ;

//------------------------------------------------------------------------------

            virtual real
            get_current() const ;

//------------------------------------------------------------------------------

            string
            get_label() const ;

//------------------------------------------------------------------------------

            virtual void
            set_current( const real aCurrent) ;

//------------------------------------------------------------------------------

            virtual void
            compute_current() ;

//------------------------------------------------------------------------------

            virtual void
            set_timestep( const real aDeltaTime ) ;

//------------------------------------------------------------------------------

            virtual void
            shift( const real aTime, const real aDeltaTime ) ;

//------------------------------------------------------------------------------

            virtual void
            shift_back() ;

//------------------------------------------------------------------------------

            /**
             * write this component's time-stepping state into the restart
             * group, dataset names prefixed with aPrefix. Default: no state
             * ( resistors, diodes and superconductors recompute from voltages )
             */
            virtual void
            save_state( hid_t aGroup, const string & aPrefix ) ;

//------------------------------------------------------------------------------

            /**
             * restore the state written by save_state(). The loaded object is
             * coherent for both legal continuations: the controller triad
             * ( set_timestep - shift - stamp ) and a stamp without shift.
             * Calling shift_back() before the first post-restore shift() is
             * forbidden
             */
            virtual void
            load_state( hid_t aGroup, const string & aPrefix ) ;

//------------------------------------------------------------------------------

            virtual bool
            is_closed() const ;

//------------------------------------------------------------------------------

            virtual void
            switch_state() ;

//------------------------------------------------------------------------------

            virtual ElectricNode *
            get_node_plus() const ;

//------------------------------------------------------------------------------

            virtual ElectricNode *
            get_node_minus() const ;

//------------------------------------------------------------------------------

        };
    }
}

#endif //BELFEM_CL_COMPONENT_HPP
