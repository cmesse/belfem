//
// Created by Christian Messe on 2019-08-19.
//

#ifndef BELFEM_CL_GT_HEATPOLYCUSTOM_HPP
#define BELFEM_CL_GT_HEATPOLYCUSTOM_HPP

#include "cl_Vector.hpp"
#include "cl_GT_HeatPoly.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        /**
         * The polynomial with arbitrary coefficients.
         */
        class HeatPolyCustom : public HeatPoly
        {
            // coefficients
            const Vector<real> mExponents;

            const uint mNumberOfCoeffs;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            HeatPolyCustom(
                    const real           aTmin,
                    const real           aTmax,
                    const real           aEnthalpyConstant,
                    const real           aEntropyConstant,
                    const Vector<real> & aCoefficients,
                    const Vector<real> & aExponents );

//------------------------------------------------------------------------------

            ~HeatPolyCustom() = default;

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
    } /* namespace gastables */
} /* namespace belfem */
#endif //BELFEM_CL_GT_HEATPOLYCUSTOM_HPP
