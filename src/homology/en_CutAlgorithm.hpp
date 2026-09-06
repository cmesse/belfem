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

#ifndef BELFEM_EN_CUTALGORITHM_HPP
#define BELFEM_EN_CUTALGORITHM_HPP
#include "typedefs.hpp"
#include "stringtools.hpp"

namespace belfem
{
    namespace mesh
    {
        enum class CutAlgorithm
        {
            Pellikka = 0,
            CCR = 1,
            BeltedTree = 2,
            PellikkaGeneralized = 3,
            UNDEFINED = 4
        };

        inline string
        to_string( CutAlgorithm aCutAlgorithm )
        {
            switch ( aCutAlgorithm )
            {
                case CutAlgorithm::Pellikka:
                    return "pellikka";
                case CutAlgorithm::CCR:
                    return "ccr";
                case CutAlgorithm::BeltedTree:
                    return "belted tree";
                case CutAlgorithm::PellikkaGeneralized:
                    return "generalized pellikka";
                default:
                    return "undefined";
            }
        }

        inline CutAlgorithm
        to_cut_algorithm( const string & aCutAlgorithm )
        {
            uint tNumAlgorithms = static_cast< uint >( CutAlgorithm::UNDEFINED );
            string tCutAlgorithm = string_to_lower( aCutAlgorithm );

            for ( uint a=0; a<tNumAlgorithms; ++a )
            {
                if ( to_string( static_cast< CutAlgorithm >( a ) ) == tCutAlgorithm )
                {
                    return static_cast< CutAlgorithm >( a );
                }
            }
            return CutAlgorithm::UNDEFINED;
        }

    }
}

#endif //BELFEM_EN_CUTALGORITHM_HPP
