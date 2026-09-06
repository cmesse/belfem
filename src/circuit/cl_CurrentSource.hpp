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

#ifndef BELFEM_CL_CURRENTSOURCE_HPP
#define BELFEM_CL_CURRENTSOURCE_HPP

#include "cl_TwoTerminals.hpp"
#include "cl_SourceFunction.hpp"

namespace belfem
{
    namespace electronics
    {
//-----------------------------------------------------------------------------

        class CurrentSource : public TwoTerminals
        {

            //! Function
            SourceFunction * mFunction ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            CurrentSource( SourceFunction * aFunction, Cell < ElectricNode* > & aTerminals, string aLabel = "") ;

//------------------------------------------------------------------------------

            ~CurrentSource() override ;

//------------------------------------------------------------------------------

            ComponentType
            component_type() const override ;

//-----------------------------------------------------------------------------

            real
            get_value() const override ;

//-----------------------------------------------------------------------------

            void
            shift( const real aTime, const real aDeltaTime ) override ;

//------------------------------------------------------------------------------

            void
            save_state( hid_t aGroup, const string & aPrefix ) override ;

//------------------------------------------------------------------------------

            void
            load_state( hid_t aGroup, const string & aPrefix ) override ;

//------------------------------------------------------------------------------


        };

        inline ComponentType
        CurrentSource::component_type() const
        {
            return ComponentType::CURRENTSOURCE ;
        }

//-----------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_CURRENTSOURCE_HPP
