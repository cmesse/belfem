//
// Created by christian on 12/3/21.
//

#ifndef BELFEM_CL_EF_TRI3_HPP
#define BELFEM_CL_EF_TRI3_HPP

#include "nedelec/cl_EF_EdgeFunction.hpp"
#include "cl_FEM_Element.hpp"

namespace belfem
{
    namespace fem
    {
        class Element;

        class EF_TRI3 : public EdgeFunction
        {
            // edge directions
            real mS[ 3 ];

            // derivatives
            real mNablaXi[ 2 ];
            real mNablaEta[ 2 ];
            real mNablaZeta[ 2 ];

            // coefficients for shape function
            Matrix< real > mG ;
            Matrix< real > mH ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * constructor
             */
            EF_TRI3();

//------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            ~EF_TRI3() override = default;

//------------------------------------------------------------------------------

            void
            precompute( const Matrix< real > & aXi ) override;

//------------------------------------------------------------------------------

            /**
             * connect function with element
             */
            void
            link( Element * aElement ) override;

//------------------------------------------------------------------------------

            /**
             * interpolation operator for h-field
             */
            const Matrix< real > &
            E( const uint aIndex ) override;

//------------------------------------------------------------------------------

            /**
              * curl operator for h-field
              */
            const Matrix< real > &
            C( const uint aIndex = 0 ) override;

//------------------------------------------------------------------------------

            /**
              * gradient operator for h-field
              * ( constant for this element, filled in link();
              *   layout contract in EdgeFunction::mGrad )
              */
            const Matrix< real > &
            G( const uint aIndex = 0 ) override ;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline const Matrix <real> &
        EF_TRI3::C( const uint aIndex )
        {
            return mC ;
        }

//------------------------------------------------------------------------------

        inline const Matrix <real> &
        EF_TRI3::G( const uint aIndex )
        {
            return mGrad ;
        }

//------------------------------------------------------------------------------
    } /* end namespace fem */
}  /* end namespace belfem */

#endif //BELFEM_CL_EF_TRI3_HPP