//
// Created by Christian Messe on 25.08.19.
//

#ifndef BELFEM_CL_GT_COMPARISONOBJECTS_HPP
#define BELFEM_CL_GT_COMPARISONOBJECTS_HPP

#include "cl_GT_HeatPoly.hpp"
#include "cl_GT_TransportPoly.hpp"
namespace belfem
{
    namespace gastables
    {
        // comparision object for heat polynomial
        struct
        {
            bool
            operator()( const HeatPoly * aA, const HeatPoly * aB )
            {
                return aA->T_min() < aB->T_min();
            }
        } fHeatPoly;

        // comparision object for transport polynomial
        struct
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
