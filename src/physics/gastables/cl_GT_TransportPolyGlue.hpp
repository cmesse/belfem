//
// Created by Christian Messe on 26.08.19.
//

#ifndef BELFEM_CL_GT_TRANSPORTPOLYGLUE_HPP
#define BELFEM_CL_GT_TRANSPORTPOLYGLUE_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_GT_TransportPoly.hpp"
namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------
        class TransportPolyGlue: public TransportPoly
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            TransportPolyGlue(
                    const enum TransportPolyType aType,
                    const real aTmin,
                    const real aTmax,
                    const Vector<real> & aCoefficients );

//------------------------------------------------------------------------------

            ~TransportPolyGlue() = default;

//------------------------------------------------------------------------------

            real
            rawpoly( const real aT ) const;

//------------------------------------------------------------------------------

            real
            drawpoly( const real aT ) const;

//------------------------------------------------------------------------------

            real
            ddrawpoly( const real aT ) const;

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */


#endif //BELFEM_CL_GT_TRANSPORTPOLYGLUE_HPP
