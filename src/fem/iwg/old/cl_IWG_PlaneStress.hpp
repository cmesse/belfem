//
// Created by Christian Messe on 15.11.19.
//

#ifndef BELFEM_CL_IWG_PlaneStress_HPP
#define BELFEM_CL_IWG_PlaneStress_HPP

#include "cl_IWG.hpp"

namespace belfem
{
    namespace fem
    {
        class IWG_PlaneStress : public IWG
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_PlaneStress();

//------------------------------------------------------------------------------

            ~IWG_PlaneStress() = default;

//------------------------------------------------------------------------------

            void
            compute_jacobian(
                    Element        * aElement,
                    Matrix< real > & aJacobian );

//------------------------------------------------------------------------------

        };
    }
}

#endif //BELFEM_CL_IWG_PlaneStress_HPP
