//
// Created by Christian Messe on 2019-08-18.
//

#ifndef BELFEM_CL_GT_THERMOPOLY_HPP
#define BELFEM_CL_GT_THERMOPOLY_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace gastables
    {
        class HeatPoly
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------
            // minumum temperature of this poly
            real mTmin;

            // maximum temperature
            real mTmax;

            // offset for enthalpy
            real mEnthalpyConstant;

            // offset for entropy
            real mEntropyConstant;

            // coefficients
            const Vector <real> mCoefficients;


//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            HeatPoly(
                    const real aTmin,
                    const real aTmax,
                    const real aEnthalpyConstant,
                    const real aEntropyConstant,
                    const Vector <real> & aCoefficients );

//------------------------------------------------------------------------------

            virtual ~HeatPoly() = default;

//------------------------------------------------------------------------------

            virtual real
            Cp( const real & aT ) const;

//------------------------------------------------------------------------------

            virtual real
            H( const real & aT ) const;

//------------------------------------------------------------------------------

            virtual real
            S( const real & aT ) const;

//------------------------------------------------------------------------------

            virtual real
            dSdT( const real & aT ) const;

//------------------------------------------------------------------------------

            virtual real
            dCpdT( const real & aT ) const;

//------------------------------------------------------------------------------

            virtual real
            d2CpdT2( const real & aT ) const;

//------------------------------------------------------------------------------

            inline const real &
            T_min() const;

//------------------------------------------------------------------------------

            inline const real &
            T_max() const;

//------------------------------------------------------------------------------

            void
            set_T_min( const real & aTmin );

//------------------------------------------------------------------------------

            void
            set_T_max( const real & aTmax );

//------------------------------------------------------------------------------

            void
            set_enthalpy_constant( const real & aEnthalpyConstant );

//------------------------------------------------------------------------------

            void
            set_entropy_constant( const real & aEntropyConstant );

//------------------------------------------------------------------------------

            real
            enthalpy_constant() const;
//------------------------------------------------------------------------------

            real
            entropy_constant() const;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        const real &
        HeatPoly::T_min() const
        {
            return mTmin;
        }

//------------------------------------------------------------------------------

        const real &
        HeatPoly::T_max() const
        {
            return mTmax;
        }

//------------------------------------------------------------------------------


    }
}
#endif //BELFEM_CL_GT_THERMOPOLY_HPP
