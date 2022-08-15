//
// Created by christian on 8/11/22.
//

#ifndef BELFEM_CL_EF_LINE3_HPP
#define BELFEM_CL_EF_LINE3_HPP

#include "cl_EF_EdgeFunction.hpp"
#include "cl_FEM_Element.hpp"
#include "fn_polyval.hpp"

namespace belfem
{
    namespace fem
    {
        class Element;

        /**
         * edge function for thin shell
         */
        class EF_LINE3 : public EdgeFunction
        {
            // coordinates of edge
            Vector< real > mX ;
            Vector< real > mY ;

            // scaling parameters
            Cell< Vector< real > > mNxi ;

            Cell< Matrix< real > > mE ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * constructor
             */
            EF_LINE3();

//------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            ~EF_LINE3();

//------------------------------------------------------------------------------

            /**
             * @param aXi
             */
            void
            precompute( const Matrix< real > & aXi ) ;

//------------------------------------------------------------------------------

            /**
             * connect function with element
             * ( here, we don't care about the switches )
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
        private:
//------------------------------------------------------------------------------

            real
            compute_det_J( const uint aIndex );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline const Matrix< real > &
        EF_LINE3::E( const uint aIndex )
        {
            mAbsDetJ = this->compute_det_J( aIndex );
            return mE( aIndex );
        }

//--------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_EF_LINE3_HPP
