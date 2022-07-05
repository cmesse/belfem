//
// Created by Christian Messe on 11.11.19.
//

#ifndef BELFEM_CL_IWG_SURFACEGRADIENT_HPP
#define BELFEM_CL_IWG_SURFACEGRADIENT_HPP


#include "cl_IWG.hpp"

namespace belfem
{
    namespace fem
    {
        class IWG_SurfaceGradient : public IWG
        {

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_SurfaceGradient();

//------------------------------------------------------------------------------

            ~IWG_SurfaceGradient() = default;

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

            void
            allocate_work_matrices( Group * aGroup ) ;

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* end namespace fem */
} /* end namespace belfem */

#endif //BELFEM_CL_IWG_GRADIENT_HPP
