//
// Created by christian on 2/17/23.
//

#include "cl_IWG_Maxwell_3D.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        const Matrix< real > &
        IWG_Maxwell_3D::inv_J_tri3( Element * aElement )
        {
            Matrix< real > & tJ    = mGroup->work_J() ;
            Matrix< real > & aInvJ = mGroup->work_invJ() ;

            // Jacobian
            tJ( 0, 0 ) =
                      aElement->element()->node( 0 )->x()
                    - aElement->element()->node( 3 )->x() ;

            tJ( 1, 0 ) =
                      aElement->element()->node( 1 )->x()
                    - aElement->element()->node( 3 )->x() ;

            tJ( 2, 0 ) =
                      aElement->element()->node( 2 )->x()
                    - aElement->element()->node( 3 )->x() ;

            tJ( 0, 1 ) =
                      aElement->element()->node( 0 )->y()
                    - aElement->element()->node( 3 )->y() ;

            tJ( 1, 1 ) =
                      aElement->element()->node( 1 )->y()
                    - aElement->element()->node( 3 )->y() ;

            tJ( 2, 1 ) =
                      aElement->element()->node( 2 )->y()
                    - aElement->element()->node( 3 )->y() ;

            tJ( 0, 2 ) =
                      aElement->element()->node( 0 )->z()
                    - aElement->element()->node( 3 )->z() ;

            tJ( 1, 2 ) =
                      aElement->element()->node( 1 )->z()
                    - aElement->element()->node( 3 )->z() ;

            tJ( 2, 2 ) =
                      aElement->element()->node( 2 )->z()
                    - aElement->element()->node( 3 )->z() ;

            mGroup->work_det_J() =
                    tJ( 0, 0 )*tJ( 1, 1 )*tJ( 2, 2 )
                  - tJ( 0, 0 )*tJ( 1, 2 )*tJ( 2, 1 )
                  - tJ( 0, 1 )*tJ( 1, 0 )*tJ( 2, 2 )
                  + tJ( 0, 1 )*tJ( 1, 2 )*tJ( 2, 0 )
                  + tJ( 0, 2 )*tJ( 1, 0 )*tJ( 2, 1 )
                  - tJ( 0, 2 )*tJ( 1, 1 )*tJ( 2, 0 ) ;

            aInvJ( 0, 0 ) =
                      tJ( 1, 1 )*tJ( 2, 2 )
                    - tJ( 1, 2 )*tJ( 2, 1 ) ;

            aInvJ( 1, 0 ) =
                    - tJ( 1, 0 )*tJ( 2, 2 )
                    + tJ( 1, 2 )*tJ( 2, 0 ) ;

            aInvJ( 2, 0 ) =
                      tJ( 1, 0 )*tJ( 2, 1 )
                    - tJ( 1, 1 )*tJ( 2, 0 ) ;

            aInvJ( 0, 1 ) =
                    - tJ( 0, 1 )*tJ( 2, 2 )
                    + tJ( 0, 2 )*tJ( 2, 1 ) ;

            aInvJ( 1, 1 ) =
                      tJ( 0, 0 )*tJ( 2, 2 )
                    - tJ( 0, 2 )*tJ( 2, 0 ) ;

            aInvJ( 2, 1 ) =
                    - tJ( 0, 0 )*tJ( 2, 1 )
                    + tJ( 0, 1 )*tJ( 2, 0 ) ;

            aInvJ( 0, 2 ) =
                      tJ( 0, 1 )*tJ( 1, 2 )
                    - tJ( 0, 2 )*tJ( 1, 1 ) ;

            aInvJ( 1, 2 ) =
                    - tJ( 0, 0 )*tJ( 1, 2 )
                    + tJ( 0, 2 )*tJ( 1, 0 ) ;

            aInvJ( 2, 2 ) =
                      tJ( 0, 0 )*tJ( 1, 1 )
                    - tJ( 0, 1 )*tJ( 1, 0 ) ;

            aInvJ /= mGroup->work_det_J() ;

            return aInvJ ;
        }
    }
}