//
// Created by Christian Messe on 25.08.19.
//

#ifndef BELFEM_CL_GT_TRANSPORTPOLY_HPP
#define BELFEM_CL_GT_TRANSPORTPOLY_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------

        enum class TransportPolyType
        {
            VISCOSITY,
            CONDUCTIVITY,
            UNDEFINED
        };

//------------------------------------------------------------------------------

        class TransportPoly
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------
            const enum TransportPolyType mType;

            // unit scale
            const real mScale;

            // miniumum temperature
            real mTmin;

            // maximum temperature
            real mTmax;

            // coefficients
            const Vector <real> mCoefficients;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            TransportPoly(
                    const enum TransportPolyType aType,
                    const real aTmin,
                    const real aTmax,
                    const Vector <real> & aCoefficients );

            virtual ~TransportPoly() = default;

//------------------------------------------------------------------------------

            virtual real
            rawpoly( const real & aT ) const;

//------------------------------------------------------------------------------

            virtual real
            drawpoly( const real & aT ) const;

//------------------------------------------------------------------------------

            virtual real
            ddrawpoly( const real & aT ) const;

//------------------------------------------------------------------------------

            virtual real
            eval( const real & aT ) const;

//------------------------------------------------------------------------------

            virtual real
            deval( const real & aT ) const;

//------------------------------------------------------------------------------

            virtual real
            ddeval( const real & aT ) const;

//------------------------------------------------------------------------------

            inline const real &
            T_min() const;

//------------------------------------------------------------------------------

            inline const real &
            T_max() const;

//------------------------------------------------------------------------------

            inline const TransportPolyType &
            type();

//------------------------------------------------------------------------------

            void
            set_T_min( const real & aTmin );

//------------------------------------------------------------------------------

            void
            set_T_max( const real & aTmax );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        const real &
        TransportPoly::T_min() const
        {
            return mTmin;
        }

//------------------------------------------------------------------------------

        const real &
        TransportPoly::T_max() const
        {
            return mTmax;
        }
//------------------------------------------------------------------------------

        inline const TransportPolyType &
        TransportPoly::type()
        {
            return mType;
        }

//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */
#endif //BELFEM_CL_GT_TRANSPORTPOLY_HPP
