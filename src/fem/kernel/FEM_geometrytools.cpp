//
// Created by christian on 2/28/23.
//

#include "FEM_geometrytools.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        const Vector< real > &
        normal_tri( Group * aGroup, const uint aMasterIndex )
        {
            Vector< real > & aN = aGroup->work_normal() ;
            const Matrix< real > & tX = aGroup->work_Xm() ;

            switch( aMasterIndex )
            {
                case( 0 ) :
                {
                    aN( 0 ) = tX( 1, 1 )- tX( 0, 1 );
                    aN( 1 ) = tX( 0, 0 ) - tX( 1, 0 );

                    break;
                }
                case( 1 ) :
                {
                    aN( 0 ) = tX( 2, 1 )-tX( 1, 1 );
                    aN( 1 ) = tX( 1, 0 )-tX( 2, 0 );
                    break;
                }
                case( 2 ) :
                {
                    aN( 0 ) = tX( 0, 1 )-tX( 2, 1 );
                    aN( 1 ) = tX( 2, 0 )-tX( 0, 0 );
                    break;
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid index" );
                }
            }

            // this value will contain the length of the side
            aGroup->work_det_J() = std::sqrt(
                    aN( 0 )* aN( 0 ) + aN( 1 )* aN( 1 ) );

            // now let's norm the vector
            aN /= aGroup->work_det_J() ;

            // finally, we must adapt this value, since along the edge
            // we integrate from -1 to 1 rather than 0 to 1
            aGroup->work_det_J() *= 0.5 ;

            // now we can return the normal
            return aN ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        normal_tri( Group * aGroup, const uint aMasterIndex, const uint aIndex )
        {
            Vector< real > & aN = aGroup->work_normal() ;

            Matrix< real > & tJ = aGroup->work_J();

            // compute the Jacobian matrix
            tJ = aGroup->master_integration( aMasterIndex )
                         ->dNdXi( aIndex )* aGroup->work_Xm() ;

            switch( aMasterIndex )
            {
                case( 0 ) :
                {
                    aN( 0 ) = tJ( 1, 1 ) - tJ( 0, 1 );
                    aN( 1 ) = tJ( 0, 0 ) - tJ( 1, 0 );
                    break;
                }
                case( 1 ) :
                {
                    aN( 0 ) = - tJ( 1, 1 ) ;
                    aN( 1 ) =   tJ( 1, 0 ) ;
                    break;
                }
                case( 2 ) :
                {
                    aN( 0 ) =  tJ( 0, 1 ) ;
                    aN( 1 ) = -tJ( 0, 0 ) ;
                    break;
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid index for facet" );
                }
            }

            // if this edge was straight, this would be the length of this side
            aGroup->work_det_J() = std::sqrt(
                    aN( 0 )* aN( 0 ) + aN( 1 )* aN( 1 ) );

            // now let's norm the vector
            aN /= aGroup->work_det_J() ;

            // finally, we must adapt this value, since along the edge
            // we integrate from -1 to 1 rather than 0 to 1
            aGroup->work_det_J() *= 0.5 ;

            // now we can return the normal
            return aN ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        normal_tet( Group * aGroup, const uint aMasterIndex )
        {
            Vector< real >       & aN = aGroup->work_normal() ;
            const Matrix< real > & tX = aGroup->work_Xm() ;

            switch( aMasterIndex )
            {
                case( 0 ) :
                {
                    aN( 0 ) = (tX(0,1)-tX(1,1))*(tX(0,2)-tX(3,2))-(tX(0,1)-tX(3,1))*(tX(0,2)-tX(1,2));
                    aN( 1 ) = (tX(0,0)-tX(3,0))*(tX(0,2)-tX(1,2))-(tX(0,0)-tX(1,0))*(tX(0,2)-tX(3,2));
                    aN( 2 ) = (tX(0,0)-tX(1,0))*(tX(0,1)-tX(3,1))-(tX(0,0)-tX(3,0))*(tX(0,1)-tX(1,1));
                    break ;
                }
                case( 1 ) :
                {
                    aN( 0 ) = (tX(1,1)-tX(2,1))*(tX(1,2)-tX(3,2))-(tX(1,1)-tX(3,1))*(tX(1,2)-tX(2,2));
                    aN( 1 ) = (tX(1,0)-tX(3,0))*(tX(1,2)-tX(2,2))-(tX(1,0)-tX(2,0))*(tX(1,2)-tX(3,2));
                    aN( 2 ) = (tX(1,0)-tX(2,0))*(tX(1,1)-tX(3,1))-(tX(1,0)-tX(3,0))*(tX(1,1)-tX(2,1));
                    break ;
                }
                case( 2 ) :
                {
                    aN( 0 ) = (tX(1,1)-tX(2,1))*(tX(1,2)-tX(3,2))-(tX(1,1)-tX(3,1))*(tX(1,2)-tX(2,2));
                    aN( 1 ) = (tX(1,0)-tX(3,0))*(tX(1,2)-tX(2,2))-(tX(1,0)-tX(2,0))*(tX(1,2)-tX(3,2));
                    aN( 2 ) = (tX(1,0)-tX(2,0))*(tX(1,1)-tX(3,1))-(tX(1,0)-tX(3,0))*(tX(1,1)-tX(2,1));
                    break ;
                }
                case( 3 ) :
                {
                    aN( 0 ) = (tX(0,1)-tX(2,1))*(tX(0,2)-tX(1,2))-(tX(0,1)-tX(1,1))*(tX(0,2)-tX(2,2));
                    aN( 1 ) = (tX(0,0)-tX(1,0))*(tX(0,2)-tX(2,2))-(tX(0,0)-tX(2,0))*(tX(0,2)-tX(1,2));
                    aN( 2 ) = (tX(0,0)-tX(2,0))*(tX(0,1)-tX(1,1))-(tX(0,0)-tX(1,0))*(tX(0,1)-tX(2,1));
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid index for facet" );
                }
            }

            // if this edge was straight, this would be the length of this side
            aGroup->work_det_J() = std::sqrt(
                    aN( 0 )* aN( 0 ) + aN( 1 )* aN( 1 ) + aN( 2 ) * aN( 2 ) );

            // now let's norm the vector
            aN /= aGroup->work_det_J() ;

            // finally, we must adapt this value, since along the edge
            // we integrate from -1 to 1 rather than 0 to 1
            aGroup->work_det_J() *= 0.5 ;

            // now we can return the normal
            return aN ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        normal_tet( Group * aGroup, const uint aMasterIndex, const uint aIndex )
        {
            Vector< real > & aN = aGroup->work_normal() ;

            Matrix< real > & tJ = aGroup->work_J();

            // compute the Jacobian matrix
            tJ = aGroup->master_integration( aMasterIndex )
                         ->dNdXi( aIndex )* aGroup->work_Xm() ;

            switch( aMasterIndex )
            {
                case( 0 ) :
                {
                    aN( 0 ) = tJ(0,2)*(tJ(0,1)-tJ(1,1))-tJ(0,1)*(tJ(0,2)-tJ(1,2));
                    aN( 1 ) = tJ(0,0)*(tJ(0,2)-tJ(1,2))-tJ(0,2)*(tJ(0,0)-tJ(1,0));
                    aN( 2 ) = tJ(0,1)*(tJ(0,0)-tJ(1,0))-tJ(0,0)*(tJ(0,1)-tJ(1,1));
                    break ;
                }
                case( 1 ) :
                {
                    aN( 0 ) = tJ(1,2)*(tJ(1,1)-tJ(2,1))-tJ(1,1)*(tJ(1,2)-tJ(2,2));
                    aN( 1 ) = tJ(1,0)*(tJ(1,2)-tJ(2,2))-tJ(1,2)*(tJ(1,0)-tJ(2,0));
                    aN( 2 ) = tJ(1,1)*(tJ(1,0)-tJ(2,0))-tJ(1,0)*(tJ(1,1)-tJ(2,1));
                    break ;
                }
                case( 2 ) :
                {
                    aN( 0 ) = tJ(1,2)*(tJ(1,1)-tJ(2,1))-tJ(1,1)*(tJ(1,2)-tJ(2,2));
                    aN( 1 ) = tJ(1,0)*(tJ(1,2)-tJ(2,2))-tJ(1,2)*(tJ(1,0)-tJ(2,0));
                    aN( 2 ) = tJ(1,1)*(tJ(1,0)-tJ(2,0))-tJ(1,0)*(tJ(1,1)-tJ(2,1));
                    break ;
                }
                case( 3 ) :
                {
                    aN( 0 ) = tJ(1,2)*(tJ(1,1)-tJ(2,1))-tJ(1,1)*(tJ(1,2)-tJ(2,2));
                    aN( 1 ) = tJ(1,0)*(tJ(1,2)-tJ(2,2))-tJ(1,2)*(tJ(1,0)-tJ(2,0));
                    aN( 2 ) = tJ(1,1)*(tJ(1,0)-tJ(2,0))-tJ(1,0)*(tJ(1,1)-tJ(2,1));
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid index for facet" );
                }
            }
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        inv_J_tri3( Group * aGroup, Element * aElement )
        {
            Matrix< real > & aInvJ = aGroup->work_invJ() ;

            const Matrix< real > & tX = aGroup->work_X() ;

            aInvJ( 0, 0 ) = tX( 1, 1 ) - tX( 2, 1 );

            aInvJ( 0, 1 ) = tX( 2, 1 ) - tX( 0, 1 );

            aInvJ( 1, 0 ) = tX( 2, 0 ) - tX( 1, 0 );

            aInvJ( 1, 1 ) = tX( 0, 0 ) - tX( 2, 0 );

            aGroup->work_det_J() =
                    aInvJ( 0, 0 ) * aInvJ( 1, 1 )
                  - aInvJ( 0, 1 ) * aInvJ( 1, 0 );

            aInvJ /= aGroup->work_det_J() ;

            return aInvJ ;
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        inv_J_tet4( Group * aGroup )
        {
            Matrix< real > & aInvJ = aGroup->work_invJ() ;

            // grab node coordinates from element
            const Matrix< real > & tX = aGroup->work_X() ;

            // work vector
            real * tW = aGroup->work_geo().data() ;

            tW[  0 ] = tX( 1, 1 ) - tX( 0, 1 ); // y1-y0
            tW[  1 ] = tX( 2, 1 ) - tX( 0, 1 ); // y2-y0
            tW[  2 ] = tX( 3, 1 ) - tX( 0, 1 ); // y3-y0
            tW[  3 ] = tX( 2, 1 ) - tX( 1, 1 ); // y2-y1
            tW[  4 ] = tX( 3, 1 ) - tX( 1, 1 ); // y3-y1
            tW[  5 ] = tX( 3, 1 ) - tX( 2, 1 ); // y3-y2

            tW[  6 ] = tX( 1, 2 ) - tX( 0, 2 ); // z1-z0
            tW[  7 ] = tX( 2, 2 ) - tX( 0, 2 ); // z2-z0
            tW[  8 ] = tX( 3, 2 ) - tX( 0, 2 ); // z3-z0
            tW[  9 ] = tX( 2, 2 ) - tX( 1, 2 ); // z2-z1
            tW[ 10 ] = tX( 3, 2 ) - tX( 1, 2 ); // z3-z1
            tW[ 11 ] = tX( 3, 2 ) - tX( 2, 2 ); // z3-z2

            // determinant
            aGroup->work_det_J() =
				tX(0,0)*(tX(2,1)*tW[10]-tX(1,1)*tW[11]-tX(3,1)*tW[9])
               +tX(1,0)*(tX(0,1)*tW[11]-tX(2,1)*tW[8]+tX(3,1)*tW[7])
               +tX(2,0)*(tX(1,1)*tW[8]-tX(0,1)*tW[10]-tX(3,1)*tW[6])
               +tX(3,0)*(tX(0,1)*tW[9]-tX(1,1)*tW[7]+tX(2,1)*tW[6]);

            // inverse jacobian
            aInvJ(0,0)=tX(2,1)*tW[10]-tX(1,1)*tW[11]-tX(3,1)*tW[9];
            aInvJ(1,0)=tX(1,0)*tW[11]-tX(2,0)*tW[10]+tX(3,0)*tW[9];
            aInvJ(2,0)=tX(2,0)*tW[4]-tX(1,0)*tW[5]-tX(3,0)*tW[3];
            aInvJ(0,1)=tX(0,1)*tW[11]-tX(2,1)*tW[8]+tX(3,1)*tW[7];
            aInvJ(1,1)=tX(2,0)*tW[8]-tX(0,0)*tW[11]-tX(3,0)*tW[7];
            aInvJ(2,1)=tX(0,0)*tW[5]-tX(2,0)*tW[2]+tX(3,0)*tW[1];
            aInvJ(0,2)=tX(1,1)*tW[8]-tX(0,1)*tW[10]-tX(3,1)*tW[6];
            aInvJ(1,2)=tX(0,0)*tW[10]-tX(1,0)*tW[8]+tX(3,0)*tW[6];
            aInvJ(2,2)=tX(1,0)*tW[2]-tX(0,0)*tW[4]-tX(3,0)*tW[0];

            aInvJ /= aGroup->work_det_J() ;

            return aInvJ ;
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        inv_J_2d( Group * aGroup )
        {
            const Matrix< real > & tJ = aGroup->work_J() ;
            Matrix< real > & aInvJ = aGroup->work_invJ() ;

            aInvJ( 0, 0 ) =  tJ( 1, 1 );
            aInvJ( 1, 0 ) = -tJ( 1, 0 );
            aInvJ( 0, 1 ) = -tJ( 0, 1 );
            aInvJ( 1, 1 ) =  tJ( 0, 0 );

            aGroup->work_det_J() =  tJ(0,0)*tJ(1,1) - tJ(0,1)*tJ(1,0) ;

            aInvJ /= aGroup->work_det_J() ;

            return aInvJ ;
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        inv_J_3d( Group * aGroup )
        {
            const Matrix< real > & tJ = aGroup->work_J() ;
            Matrix< real > & aInvJ = aGroup->work_invJ() ;

            aInvJ(0,0) =  tJ(1,1)*tJ(2,2)-tJ(1,2)*tJ(2,1);
            aInvJ(1,0) = -tJ(1,0)*tJ(2,2)+tJ(1,2)*tJ(2,0);
            aInvJ(2,0) =  tJ(1,0)*tJ(2,1)-tJ(1,1)*tJ(2,0);

            aInvJ(0,1) = -tJ(0,1)*tJ(2,2)+tJ(0,2)*tJ(2,1);
            aInvJ(1,1) =  tJ(0,0)*tJ(2,2)-tJ(0,2)*tJ(2,0);
            aInvJ(2,1) = -tJ(0,0)*tJ(2,1)+tJ(0,1)*tJ(2,0);

            aInvJ(0,2) =  tJ(0,1)*tJ(1,2)-tJ(0,2)*tJ(1,1);
            aInvJ(1,2) = -tJ(0,0)*tJ(1,2)+tJ(0,2)*tJ(1,0);
            aInvJ(2,2) =  tJ(0,0)*tJ(1,1)-tJ(0,1)*tJ(1,0);

            aGroup->work_det_J() = tJ(0,0)*tJ(1,1)*tJ(2,2)
                                 - tJ(0,0)*tJ(1,2)*tJ(2,1)
                                 - tJ(0,1)*tJ(1,0)*tJ(2,2)
                                 + tJ(0,1)*tJ(1,2)*tJ(2,0)
                                 + tJ(0,2)*tJ(1,0)*tJ(2,1)
                                 - tJ(0,2)*tJ(1,1)*tJ(2,0) ;

            aInvJ /= aGroup->work_det_J() ;

            return aInvJ ;
        }

//------------------------------------------------------------------------------
    }
}
