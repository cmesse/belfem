//
// Created by Christian Messe on 26.08.19.
//

#ifndef BELFEM_CL_GT_HEATPOLYEMPTY_HPP
#define BELFEM_CL_GT_HEATPOLYEMPTY_HPP

#include "GT_globals.hpp"
#include "cl_GT_HeatPoly.hpp"

namespace belfem
{
    namespace gastables
    {
        /**
          * An empty polynomial that returns only zeros.
          * Creted if no data has been found in thermo.inp
          */
        class HeatPolyEmpty : public HeatPoly
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            HeatPolyEmpty() :
                HeatPoly( 0.0, gTmax, 0.0, 0.0, Vector<real> ( 1, 0.0 ) )
            {

            }

            ~HeatPolyEmpty() = default;

//------------------------------------------------------------------------------

            real
            Cp( const real aT ) const
            {
                return 0.0;
            }

//------------------------------------------------------------------------------

            real
            H( const real aT ) const
            {
                return 0.0;
            }

//------------------------------------------------------------------------------

            real
            S( const real aT ) const
            {
                return 0.0;
            }

//------------------------------------------------------------------------------

            real
            dCpdT( const real aT ) const
            {
                return 0.0;
            }

//------------------------------------------------------------------------------

            real
            d2CpdT2( const real aT ) const
            {
                return 0.0;
            }

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_GT_HEATPOLYEMPTY_HPP
