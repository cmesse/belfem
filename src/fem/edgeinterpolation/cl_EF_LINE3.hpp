//
// Created by christian on 8/11/22.
//

#ifndef BELFEM_CL_EF_LINE3_HPP
#define BELFEM_CL_EF_LINE3_HPP

#include "cl_EF_EdgeFunction.hpp"
#include "cl_FEM_Element.hpp"
#include "fn_polyval.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"

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

            // coordinates of edge dofs
            Vector< real > mXi;

            // scaling parameters
            Cell< Vector< real > > mNxi ;
            Cell< Matrix< real > > mE ;

            // coefficients for edge function on dof 0
            Vector< real > mC0 ;

            // coefficients for edge function on dof 1
            Vector< real > mC1 ;

            // shape function of volume element
            InterpolationFunction * mVolumeLagrange = nullptr ;
            Cell< Matrix< real > >  mVolumeNxi ;

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

//-----------------------------------------------------------------------------

            const Matrix< real > &
            vol_N_xi( const uint aIndex = 0 ) ;

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
//-----------------------------------------------------------------------------

        inline const Matrix< real > &
        EF_LINE3::vol_N_xi( const uint aIndex )
        {
            return mVolumeNxi( aIndex );
        }

//--------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_EF_LINE3_HPP
