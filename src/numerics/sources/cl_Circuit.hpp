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

#ifndef BELFEM_CL_CIRCUIT_HPP
#define BELFEM_CL_CIRCUIT_HPP
#include "typedefs.hpp"
#include "hdf5_types.hpp"

namespace belfem
{
    class Circuit
    {
        public:

        Circuit() = default ;
        virtual ~Circuit() = default ;

        // NON-COPYABLE, NON-MOVABLE. The base owns nothing itself and is
        // abstract, so it cannot be sliced into a value. This is not what
        // repairs the ElectricalCircuit double free -- the deletes on the
        // derived class do that. It buys one narrower thing: a future
        // subclass that owns allocations cannot stay IMPLICITLY copyable
        // by forgetting to declare its own. A hand-written copy
        // constructor that default-initializes this base is still
        // possible, and is still that subclass's problem
        Circuit( const Circuit & ) = delete ;
        Circuit( Circuit && ) = delete ;
        Circuit & operator=( const Circuit & ) = delete ;
        Circuit & operator=( Circuit && ) = delete ;

        virtual void
        set_timestep( const real aTimestep ) = 0 ;

        virtual void
        set_omega( const real aOmega ) = 0 ;

        virtual void
        compute_jacobian_and_rhs() = 0 ;

        virtual void
        compute_MNA_matrix() = 0 ;

        virtual void
        shift() = 0 ;

        virtual void
        shift_back() = 0 ;

        virtual void
        solve() = 0 ;

        virtual real
        residual() const = 0 ;

        virtual real
        current( const index_t aIndex ) const = 0 ;

        virtual real
        voltage( const index_t aIndex ) const = 0 ;

        virtual void
        set_current_and_voltage( const index_t aIndex, const real aI, const real aV ) = 0 ;

        virtual void
        save_timestep() = 0 ;

        virtual void
        save_state( hid_t aFile ) = 0 ;

        virtual void
        load_state( hid_t aFile ) = 0 ;

    };
}
#endif // BELFEM_CL_CIRCUIT_HPP
