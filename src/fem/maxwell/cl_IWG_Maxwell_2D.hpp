//
// Created by christian on 2/16/23.
//

#ifndef BELFEM_CL_IWG_MAXWELL_2D_HPP
#define BELFEM_CL_IWG_MAXWELL_2D_HPP

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
//#include "fn_norm.hpp"
//#include "fn_dot.hpp"
//#include "fn_trans.hpp"
#include "cl_IWG_Timestep.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Group.hpp"
//#include "en_SolverEnums.hpp"
//#include "en_FEM_DomainType.hpp"
//#include "cl_FEM_DomainGroup.hpp"
//#include "cl_FEM_BoundaryCondition.hpp"
//#include "en_FEM_MagfieldBcType.hpp"
#include "cl_EF_EdgeFunction.hpp"
//#include "cl_Maxwell_FieldList.hpp"

namespace belfem
{
    namespace fem
    {
        class IWG_Maxwell_2D //: public IWG_Timestep
        {
            // delete me
            Group * mGroup = nullptr ;

            //! the edge shape function
            EdgeFunction * mEdgeFunction = nullptr ;

            void
            ( IWG_Maxwell_2D::*mFunComputeMatrices )
                    ( Element * aElement,
                      Matrix< real > & aM,
                      Matrix< real > & aK,
                      Vector< real > & aF );

            void
            ( IWG_Maxwell_2D::*mFunComputeTimestep )
                    ( Element * aElement,
                      Matrix< real > & aM, // in: Mass matrix out : jacobian matrix
                      Matrix< real > & aK, // in: Stiffness matrix
                      Vector< real > & aF  // in: Load vector, out: RHS
                      );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_Maxwell_2D() = default;

            ~IWG_Maxwell_2D() = default;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            /*
             * help function, computes the inverse of the jacobian matrix
             */
            const Matrix< real > &
            inv_J_tri3( Element * aElement );

//------------------------------------------------------------------------------

            /*
             * help function, computes the B-operator matrix
             */
            const Matrix< real > &
            B_tri3( Element * aElement ) ;

//------------------------------------------------------------------------------

            void
            vol_phi_faraday_tri3( Element * aElement,
                                  Matrix< real > & aM,
                                  Matrix< real > & aK,
                                  Vector< real > & aF );


//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* end namespace fem */
} /* end namespace belfem */

#endif //BELFEM_CL_IWG_MAXWELL_2D_HPP
