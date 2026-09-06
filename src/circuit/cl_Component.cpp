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

#include "cl_Component.hpp"
#include "commtools.hpp"

namespace belfem
{
    namespace electronics
    {
//----------------------------------------------------------------------------

        Component::Component( const uint aNumTerminals, Cell < ElectricNode * > & aTerminals, string aLabel) :
                graph::Vertex(), mNumTerminals(aNumTerminals), mLabel(aLabel)
        {
            BELFEM_ERROR(aTerminals.size() == aNumTerminals,"Unmatching number of terminals") ;
            mTerminals = ( ElectricNode** ) std::malloc( mNumTerminals * sizeof( ElectricNode * ) );
            for (uint j = 0 ; j < mNumTerminals; ++j)
            {
                mTerminals[j] = aTerminals(j) ;
            }
        }

//----------------------------------------------------------------------------

        Component::~Component()
        {
            if (mTerminals != nullptr)
            {
                std::free( mTerminals) ;
            }
        }

//----------------------------------------------------------------------------

        ComponentType
        Component::component_type() const
        {
            BELFEM_ERROR( false, "Invalid call to abstract Component class : component_type()");
            return ComponentType::UNDEFINED ;
        }

//-----------------------------------------------------------------------------

        uint
        Component::number_of_terminals() const
        {
            return mNumTerminals ;
        }

//------------------------------------------------------------------------------

        ElectricNode *
        Component::node( index_t tIndex ) const
        {
            return mTerminals[tIndex] ;
        }

//------------------------------------------------------------------------------

        real
        Component::get_value() const
        {
            BELFEM_ERROR( false, "Invalid call to abstract Component class : get_value()");
            return 0.0 ;
        }

//------------------------------------------------------------------------------

        real
        Component::dIdV() const
        {
            BELFEM_ERROR( false, "Invalid call to abstract Component class : dIdV()");
            return 0.0 ;
        }

//------------------------------------------------------------------------------

        real
        Component::get_discretized_resistance() const
        {
            BELFEM_ERROR( false, "Invalid call to abstract Component class : get_discretized_resistance()");
            return 0.0 ;
        }

//------------------------------------------------------------------------------

        real
        Component::get_discretized_source() const
        {
            BELFEM_ERROR( false, "Invalid call to abstract Component class : get_discretized_source()");
            return 0.0 ;
        }

//------------------------------------------------------------------------------

        real
        Component::get_current() const
        {
            BELFEM_ERROR( false, "Invalid call to abstract Component class : get_current()");
            return 0.0 ;
        }

//------------------------------------------------------------------------------

        string
        Component::get_label() const
        {
            return mLabel ;
        }

//------------------------------------------------------------------------------

        void
        Component::set_current( const belfem::real aCurrent )
        {
            BELFEM_ERROR( false, "Invalid call to abstract Component class : set_current()");
            return ;
        }

//------------------------------------------------------------------------------

        void
        Component::compute_current()
        {
            BELFEM_ERROR( false, "Invalid call to abstract Component class : compute_current()");
            return ;
        }

//------------------------------------------------------------------------------

        void
        Component::set_timestep( const belfem::real aDeltaTime )
        {
            BELFEM_ERROR( false, "Invalid call to abstract Component class : set_timestep()");
            return ;
        }

//------------------------------------------------------------------------------

        void
        Component::shift( const real aTime, const real aDeltaTime )
        {
            BELFEM_ERROR( false, "Invalid call to abstract Component class : shift()");
            return ;
        }

//------------------------------------------------------------------------------

        void
        Component::shift_back()
        {
            BELFEM_ERROR( false, "Invalid call to abstract Component class : shift_back()");
            return ;
        }

//------------------------------------------------------------------------------

        void
        Component::save_state( hid_t aGroup, const string & aPrefix )
        {
            // default: stateless in the restart file
        }

//------------------------------------------------------------------------------

        void
        Component::load_state( hid_t aGroup, const string & aPrefix )
        {
            // default: stateless in the restart file
        }

//------------------------------------------------------------------------------

        bool
        Component::is_closed() const
        {
            BELFEM_ERROR( false, "Invalid call to abstract Component class : is_closed()");
            return 0 ;
        }

//------------------------------------------------------------------------------

        void
        Component::switch_state()
        {
            BELFEM_ERROR( false, "Invalid call to abstract Component class : switch_state()");
            return ;
        }

//------------------------------------------------------------------------------

        ElectricNode *
        Component::get_node_plus() const
        {
            BELFEM_ERROR( false, "Invalid call to abstract Component class : get_node_plus()");
            return 0 ;
        }

//------------------------------------------------------------------------------

        ElectricNode *
        Component::get_node_minus() const
        {
            BELFEM_ERROR( false, "Invalid call to abstract Component class : get_node_minus()");
            return 0 ;
        }

//------------------------------------------------------------------------------
    }
}
