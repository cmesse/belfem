//
// Created by claude on 1/11/25.
//

#ifndef BELFEM_CL_EF_QUAD4TS_HPP
#define BELFEM_CL_EF_QUAD4TS_HPP

#include "nedelec/cl_EF_EdgeFunction.hpp"
#include "cl_FEM_Element.hpp"

namespace belfem
{
    namespace fem
    {
        class Element;

        class EF_QUAD4TS : public EdgeFunction
        {
            // edge directions
            real mS[ 2 ];

            // derivatives for line element
            real mNablaXi[ 2 ];
            real mNablaEta[ 2 ];

            // normal to the line (in 2D space)
            real mN[ 2 ];
            real mDetGram ;
            real mThickness ;
            real mLength ;

            // coefficients for shape function
            Matrix< real > mF ;
            Matrix< real > mG ;

            // pseudo inverse (1D jacobian for line)
            Matrix< real > mPseudoInvJ ;

            // gram matrix (1x1 for line)
            Matrix< real > mGram ;

            // coefficient matrix
            Matrix< real > mCoeffs ;

            // parameter matrix for curl
            Cell< Matrix< real > > mCurlPars ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * constructor
             */
            EF_QUAD4TS();

//------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            ~EF_QUAD4TS() override = default;

//------------------------------------------------------------------------------

            void
            precompute( const Matrix< real > & aXi ) override ;

//------------------------------------------------------------------------------

            /**
             * connect function with element
             */
            void
            link( Element * aElement ) override ;

//------------------------------------------------------------------------------

            /**
             * interpolation operator for h-field
             */
            const Matrix< real > &
            E( const uint aIndex ) override ;

//------------------------------------------------------------------------------

            /**
              * curl operator for h-field
              */
            const Matrix< real > &
            C( const uint aIndex = 0 ) override ;

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

        inline const Matrix< real > &
        EF_QUAD4TS::G( const uint aIndex )
        {
            return mGrad ;
        }

//------------------------------------------------------------------------------
    } /* end namespace fem */
}  /* end namespace belfem */

#endif //BELFEM_CL_EF_QUAD4TS_HPP
