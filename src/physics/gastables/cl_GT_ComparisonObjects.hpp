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

#ifndef BELFEM_CL_GT_COMPARISONOBJECTS_HPP
#define BELFEM_CL_GT_COMPARISONOBJECTS_HPP

#include "cl_GT_HeatPoly.hpp"
#include "cl_GT_TransportPoly.hpp"
namespace belfem
{
    namespace gastables
    {
        // comparision object for heat polynomial
        static struct OpHeatPoly
        {
            bool
            operator()( const HeatPoly * aA, const HeatPoly * aB )
            {
                return aA->T_min() < aB->T_min();
            }
        } fHeatPoly;

        // comparision object for transport polynomial
        static struct OpTransportPoly
        {
            bool
            operator()( const TransportPoly * aA, const TransportPoly * aB )
            {
                return aA->T_min() < aB->T_min();
            }
        } fTransportPoly;

    }
}
#endif //BELFEM_CL_GT_COMPARISONOBJECTS_HPP
