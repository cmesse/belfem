//
// Created by christian on 2/9/23.
//

#ifndef BELFEM_CL_IWG_MAXWELL_THERMAL2D_HPP
#define BELFEM_CL_IWG_MAXWELL_THERMAL2D_HPP

#include "cl_IWG.hpp"

namespace belfem
{
    namespace fem
    {
        class IWG_Maxwell_Thermal2D : public IWG
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_Maxwell_Thermal2D( const uint aNumberOfDimensions=1 );

//------------------------------------------------------------------------------

            ~IWG_Maxwell_Thermal2D() = default;

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

        };
    }
}

#endif //BELFEM_CL_IWG_MAXWELL_THERMAL2D_HPP
