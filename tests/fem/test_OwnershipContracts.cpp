/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

/*
 * Ownership contracts that are checked by the compiler.
 *
 * Each class below owns raw resources and frees them in its destructor, so
 * a copy would free them twice. All four special members must stay deleted.
 * The move traits matter as much as the copy traits: a class with a
 * user-declared destructor has no implicit move, so a move request binds to
 * the copy constructor, and a later "= default" move would copy the owning
 * flags and pointers just the same. The gate is the build: a regression
 * turns this file into a compile error. ElectricalCircuit carries the same
 * four assertions in the circuit suite, which has its include paths.
 */

#include <gtest/gtest.h>
#include <type_traits>

#include "typedefs.hpp"
#include "cl_Solver.hpp"
#include "cl_IF_ElementMapper.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_Controller.hpp"

namespace
{
    template< typename T >
    constexpr bool
    is_pinned_in_place()
    {
        return ! std::is_copy_constructible_v< T >
            && ! std::is_copy_assignable_v< T >
            && ! std::is_move_constructible_v< T >
            && ! std::is_move_assignable_v< T > ;
    }
}

static_assert( is_pinned_in_place< belfem::fem::Kernel >(),
               "Kernel owns dof managers, IWGs, materials and boundary conditions; it must not be copied or moved" );

static_assert( is_pinned_in_place< belfem::Solver >(),
               "Solver owns its wrapper; it must not be copied or moved" );

static_assert( is_pinned_in_place< belfem::fem::ElementMapper >(),
               "ElementMapper owns its interpolation functions; it must not be copied or moved" );

// the clock accessors return references, not copies; which member each one
// refers to is not something a type check can see
static_assert( std::is_same_v< decltype( std::declval< const belfem::fem::Controller & >().time() ),
                               const belfem::real & >,
               "Controller::time() must return a reference" );

static_assert( std::is_same_v< decltype( std::declval< const belfem::fem::Controller & >().time_thermal() ),
                               const belfem::real & >,
               "Controller::time_thermal() must return a reference" );

TEST( Ownership, ContractsHoldAtCompileTime )
{
    // the assertions above are the test; this case records in the suite's
    // output that the file was compiled into it
    EXPECT_TRUE( is_pinned_in_place< belfem::fem::Kernel >() );
    EXPECT_TRUE( is_pinned_in_place< belfem::Solver >() );
    EXPECT_TRUE( is_pinned_in_place< belfem::fem::ElementMapper >() );
}
