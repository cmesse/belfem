//
// Created by Christian Messe on 26.08.19.
//

#ifndef BELFEM_CL_GT_TRANSPORTPOLYEMPTY_HPP
#define BELFEM_CL_GT_TRANSPORTPOLYEMPTY_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_GT_TransportPoly.hpp"
namespace belfem
{
    namespace gastables
    {
//------------------------------------------------------------------------------
        class TransportPolyEmpty: public TransportPoly
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            TransportPolyEmpty( const enum TransportPolyType aType  );

//------------------------------------------------------------------------------

            ~TransportPolyEmpty() = default;

//------------------------------------------------------------------------------

            real
            rawpoly( const real & aT ) const;

//------------------------------------------------------------------------------

            real
            drawpoly( const real & aT ) const;

//------------------------------------------------------------------------------

            real
            ddrawpoly( const real & aT ) const;

//------------------------------------------------------------------------------

            real
            eval( const real &  aT ) const;

//------------------------------------------------------------------------------

            real
            deval( const real &  aT ) const;

//------------------------------------------------------------------------------

            real
            ddeval( const real &  aT ) const;

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace gastables */
} /* namespace belfem */

#endif //BELFEM_CL_GT_TRANSPORTPOLYEMPTY_HPP
