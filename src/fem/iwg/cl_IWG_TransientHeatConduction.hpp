//
// Created by Christian Messe on 03.08.20.
//

#ifndef BELFEM_CL_IWG_TRANSIENTHEATCONDUCTION_HPP
#define BELFEM_CL_IWG_TRANSIENTHEATCONDUCTION_HPP


#include "cl_IWG_StationaryHeatConduction.hpp"

#include "fn_dot.hpp"

namespace belfem
{
    namespace fem
    {
        class IWG_TransientHeatConduction :
                public IWG_StationaryHeatConduction
        {


//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_TransientHeatConduction (
                    const uint aNumberOfDimensions,
                    const IwgType aType=IwgType::TransientHeatConduction );

//------------------------------------------------------------------------------

            virtual ~IWG_TransientHeatConduction () = default;

//------------------------------------------------------------------------------
// Functions called by Field during assembly
//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            shift_fields();

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            void
            allocate_work_matrices( Group    * aGroup );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    /* end class TransientHeatConcuction */
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
    } /* end namespace fem */
} /* end namespace belfem */
#endif //BELFEM_CL_IWG_TRANSIENTHEATCONDUCTION_HPP
