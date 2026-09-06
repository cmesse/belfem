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

#ifndef BELFEM_FN_FEM_GHOST_SWITCH_HPP
#define BELFEM_FN_FEM_GHOST_SWITCH_HPP

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_InputFile.hpp"
#include "cl_Input_Section.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        /**
         * The nonlinear section the penalty blocks are read from: the same
         * alias rule Controller::set_params applies to every key in it --
         * `nonlinear magnetic` wins over `nonlinear`, and a block under the
         * losing section is silently ignored. nullptr when neither exists
         * or when aSolver itself is nullptr
         */
        inline const input::Section *
        winning_nonlinear_section( const input::Section * aSolver )
        {
            if ( aSolver == nullptr )
            {
                return nullptr ;
            }
            if ( aSolver->section_exists( "nonlinear magnetic" ) )
            {
                return aSolver->section( "nonlinear magnetic" );
            }
            if ( aSolver->section_exists( "nonlinear" ) )
            {
                return aSolver->section( "nonlinear" );
            }
            return nullptr ;
        }

//------------------------------------------------------------------------------

        /**
         * The Nitsche ghost switch as the deck states it, read from the
         * `solver` section:
         *
         *   -1  no `nitsche ghost penalty` block  -> ghost OFF ( opt-in )
         *    0  `eta : 0`                          -> ghost OFF
         *   >0  `eta : <value>`                    -> ghost ON with that eta
         *
         * A present block without `eta` is a setup error: `k_reg` alone
         * cannot say whether the ghost is wanted, and an empty block that
         * silently meant "off" would be indistinguishable from a typo.
         *
         * ONE reader for three consumers -- the thin-shell factory ( which
         * decides whether duplicate interface dofs and ghost facets are
         * created at all ), the controller ( which writes eta into the IWG )
         * and the mesh cache tag ( which must miss the cache when the switch
         * flips ) -- so the three cannot disagree about what "off" means.
         * Setup-path code, called a handful of times per run
         */
        inline real
        read_ghost_eta( const input::Section * aSolver )
        {
            const input::Section * tNonLinear = winning_nonlinear_section( aSolver );

            if (    tNonLinear == nullptr
                 || ! tNonLinear->section_exists( "nitsche ghost penalty" ) )
            {
                return -1.0 ;
            }

            const input::Section * tGhost = tNonLinear->section( "nitsche ghost penalty" );

            BELFEM_ERROR( tGhost->key_exists( "eta" ),
                "a 'nitsche ghost penalty' block requires the key 'eta': a positive value switches "
                "the ghost coupling of the thin-shell layers on, 'eta : 0' switches it off "
                "( so does omitting the block ). 'k_reg' alone does not say which." );

            // dimensionless: a dimensioned value ( "eta : 4 mOhm ;" ) is
            // rejected instead of silently SI-scaled
            const value tEta = tGhost->get_value( "eta", "-" );

            BELFEM_ERROR( tEta.first >= 0.0, "eta must not be negative, but is %g",
                ( double ) tEta.first );

            return tEta.first ;
        }

//------------------------------------------------------------------------------

        //! true if the deck asks for duplicate interface dofs and ghost facets
        inline bool
        ghost_facets_requested( const input::Section * aSolver )
        {
            return read_ghost_eta( aSolver ) > 0.0 ;
        }

//------------------------------------------------------------------------------

        //! the same, from the whole deck; a deck without a solver section
        //! has no ghost
        inline bool
        ghost_facets_requested( const InputFile & aInputFile )
        {
            return aInputFile.section_exists( "solver" ) ?
                ghost_facets_requested( aInputFile.section( "solver" ) ) : false ;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_FN_FEM_GHOST_SWITCH_HPP
