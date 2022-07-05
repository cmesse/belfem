//
// Created by Christian Messe on 2019-08-22.
//

#ifndef BELFEM_CL_GT_HEATPOLYGLUE_HPP
#define BELFEM_CL_GT_HEATPOLYGLUE_HPP


#include "cl_Vector.hpp"
#include "cl_GT_HeatPoly.hpp"

namespace belfem
{
    namespace gastables
    {
        /**
          * A special polynomial that smoothly connects two other polynomials
          */
        class HeatPolyGlue : public HeatPoly
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            HeatPolyGlue(
                    const real             aTmin,
                    const real             aTmax,
                    const real             aEnthalpyConstant,
                    const real             aEntropyConstant,
                    const Vector< real > & aCoefficients );

//------------------------------------------------------------------------------

            ~HeatPolyGlue() = default;

//------------------------------------------------------------------------------

            real
            Cp( const real & aT ) const;

//------------------------------------------------------------------------------

            real
            H( const real & aT ) const;

//------------------------------------------------------------------------------

            real
            S( const real & aT ) const;

//------------------------------------------------------------------------------

            real
            dSdT( const real & aT ) const;

//------------------------------------------------------------------------------

            real
            dCpdT( const real & aT ) const;

//------------------------------------------------------------------------------

            real
            d2CpdT2( const real & aT ) const;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */
#endif //BELFEM_CL_GT_HEATPOLYGLUE_HPP
