//
// Created by christian on 12/3/21.
//

#ifndef BELFEM_CL_EF_TRI3_HPP
#define BELFEM_CL_EF_TRI3_HPP

#include "cl_EF_EdgeFunction.hpp"
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
            ~EF_TRI3() = default;

//------------------------------------------------------------------------------

            void
            precompute( const Matrix< real > & aXi );

//------------------------------------------------------------------------------

            /**
             * connect function with element
             */
            void
            link( Element * aElement,
                  const bool aOperatorCurlH,
                  const bool aOperatorGradPhi,
                  const bool aOperatorCurlA );

//------------------------------------------------------------------------------

            /**
             * interpolation operator for h-field
             */
            const Matrix< real > &
            E( const uint aIndex );

//------------------------------------------------------------------------------

            /**
              * curl operator for h-field
              */
            const Matrix< real > &
            C( const uint aIndex = 0 );

//------------------------------------------------------------------------------

            /**
             * gradient operator for phi-field
             */
            const Matrix< real > &
            B( const uint aIndex = 0 ) ;

//------------------------------------------------------------------------------

            /**
              * curl operator for a-field
              */
            const Matrix< real > &
            CA( const uint aIndex = 0 ) ;

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
        EF_TRI3::B( const uint aIndex )
        {
            return mB ;
        }

//------------------------------------------------------------------------------

        inline const Matrix <real> &
        EF_TRI3::CA( const uint aIndex )
        {
            return mCurlA ;
        }

//------------------------------------------------------------------------------
    } /* end namespace fem */
}  /* end namespace belfem */

#endif //BELFEM_CL_EF_TRI3_HPP