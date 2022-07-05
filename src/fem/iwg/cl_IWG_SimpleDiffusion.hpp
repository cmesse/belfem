//
// Created by Christian Messe on 08.11.19.
//

#ifndef BELFEM_CL_IWG_SIMPLEDIFFUSION_HPP
#define BELFEM_CL_IWG_SIMPLEDIFFUSION_HPP

#include "cl_IWG.hpp"

namespace belfem
{
    namespace fem
    {
        class IWG_SimpleDiffusion : public IWG
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_SimpleDiffusion( const uint aNumberOfDimensions=2 );

//------------------------------------------------------------------------------

            ~IWG_SimpleDiffusion() = default;

//------------------------------------------------------------------------------

            void
            compute_jacobian(
                    Element        * aElement,
                    Matrix< real > & aJacobian );

//------------------------------------------------------------------------------

        };
    }
}
#endif //BELFEM_CL_IWG_SIMPLEDIFFUSION_HPP
