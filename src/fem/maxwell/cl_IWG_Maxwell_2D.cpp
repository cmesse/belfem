//
// Created by christian on 2/16/23.
//

#include "cl_IWG_Maxwell_2D.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        const Matrix< real > &
        IWG_Maxwell_2D::inv_J_tri3( Element * aElement )
        {
            Matrix< real > & tJ    = mGroup->work_J() ;
            Matrix< real > & aInvJ = mGroup->work_invJ() ;

            // Jacobian
            tJ( 0, 0 ) =
                      aElement->element()->node( 0 )->x()
                    - aElement->element()->node( 2 )->x() ;

            tJ( 1, 0 ) =
                      aElement->element()->node( 1 )->x()
                    - aElement->element()->node( 2 )->x() ;

            tJ( 0, 1 ) =
                      aElement->element()->node( 0 )->y()
                    - aElement->element()->node( 2 )->y() ;

            tJ( 1, 1 ) =
                      aElement->element()->node( 1 )->y()
                    - aElement->element()->node( 2 )->y() ;

            mGroup->work_det_J() =
                      tJ( 0, 0 ) * tJ( 1, 1 )
                    - tJ( 0, 1 ) * tJ( 1, 0 );

            aInvJ( 0, 0 ) =   tJ( 1, 1 );
            aInvJ( 1, 0 ) = - tJ( 1, 0 );
            aInvJ( 0, 1 ) = - tJ( 0, 1 );
            aInvJ( 1, 1 ) =   tJ( 0, 0 );
            aInvJ /= mGroup->work_det_J() ;

            return aInvJ ;
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        IWG_Maxwell_2D::B_tri3( Element * aElement )
        {
            Matrix< real > & aB = mGroup->work_B() ;
            const Matrix< real > & tInvJ = this->inv_J_tri3( aElement );

            aB( 0, 0 ) =  tInvJ( 0, 0 );
            aB( 1, 0 ) =  tInvJ( 1, 0 );
            aB( 0, 1 ) =  tInvJ( 0, 1 );
            aB( 1, 1 ) =  tInvJ( 1, 1 );
            aB( 0, 2 ) = -tInvJ( 0, 0 )
                                            -tInvJ( 0, 1 );
            aB( 1, 2 ) = -tInvJ( 1, 0 )
                                            -tInvJ( 1, 1 );
            return aB ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_2D::vol_phi_faraday_tri3(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK,
                Vector< real > & aF )
        {
        }

//------------------------------------------------------------------------------
} /* end namespace fem */
} /* end namespace belfem */