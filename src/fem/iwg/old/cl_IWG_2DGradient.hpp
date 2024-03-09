 //
// Created by Christian Messe on 16.11.19.
//

#ifndef BELFEM_CL_IWG_2DGRADIENT_HPP
#define BELFEM_CL_IWG_2DGRADIENT_HPP

#include "cl_IWG.hpp"


namespace belfem
{
    namespace fem
    {
        class IWG_2DGradient : public IWG
        {

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_2DGradient();

//------------------------------------------------------------------------------

            ~IWG_2DGradient() = default;

//------------------------------------------------------------------------------

            void
            compute_jacobian(
                    Element        * aElement,
                    Matrix< real > & aJacobian );

//------------------------------------------------------------------------------

            void
            compute_rhs(
                    Element        * aElement,
                    Matrix< real > & aRHS );


//------------------------------------------------------------------------------
        };
    }
}


#endif //BELFEM_CL_IWG_2DGRADIENT_HPP
