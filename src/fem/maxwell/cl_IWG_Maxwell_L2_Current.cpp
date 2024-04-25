//
// Created by christian on 3/18/22.
//

#include "cl_IWG_Maxwell_L2_Current.hpp"
#include "fn_det.hpp"
#include "fn_trans.hpp"
#include "fn_inv.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Maxwell_L2_Current::IWG_Maxwell_L2_Current(
                const ElementType aElementType,
                const real            aAlpha ) :
                IWG_Maxwell( aElementType,
                        IwgType::MAXWELL_L2,
                        IwgMode::Direct,
                        SymmetryMode::PositiveDefiniteSymmetric,
                        SideSetDofLinkMode::FacetOnly,
                        true ),
                mAlpha( aAlpha )
        {
            switch( mElementType )
            {
                case( ElementType::TRI3 ) :
                {
                    mFields.Superconductor = { "jz" };

                    mFunLink = & IWG_Maxwell_L2_Current::link_jacobian_2d ;

                    mEdgeDofMultiplicity = 1 ;
                    mFaceDofMultiplicity = 0 ;
                    mNumberOfRhsDofsPerEdge = 1 ;
                    mNumberOfRhsDofsPerFace = 0 ;
                    mNumberOfRhsEdgeDofsPerElement = 3 ;

                    break ;
                }
                case( ElementType::TRI6 ) :
                {
                    mFields.Superconductor = { "jz" };

                    mFunLink = & IWG_Maxwell_L2_Current::link_jacobian_2d ;

                    mEdgeDofMultiplicity = 2 ;
                    mFaceDofMultiplicity = 2 ;
                    mNumberOfRhsDofsPerEdge = 2 ;
                    mNumberOfRhsDofsPerFace = 2 ;
                    mNumberOfRhsEdgeDofsPerElement = 8 ;
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid Element type") ;
                }
            }
        }

//------------------------------------------------------------------------------

        IWG_Maxwell_L2_Current::~IWG_Maxwell_L2_Current()
        {

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2_Current::allocate_work_matrices( Group * aGroup )
        {
            IWG::allocate_work_matrices( aGroup );

            aGroup->work_J().set_size( mNumberOfSpatialDimensions,
                                       mNumberOfSpatialDimensions, 0.0 );

            if( mAlpha != 0.0 )
            {
                mGroup->work_B().set_size( mNumberOfSpatialDimensions, mNumberOfDofsPerElement, 0.0 );
            }

            mGroup->work_nedelec().set_size( mNumberOfRhsEdgeDofsPerElement, 0.0 );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2_Current::link_jacobian_2d( Group * aGroup )
        {
            mEdgeFunction->link( aGroup );

            if( mAlpha == 0.0 )
            {
                mFunJacobian = & IWG_Maxwell_L2_Current::jacobian_2d_alpha0 ;
            }
            else
            {
                mFunJacobian = & IWG_Maxwell_L2_Current::jacobian_2d_alpha1 ;
            }
            mFunRHS = & IWG_Maxwell_L2_Current::rhs_2d ;

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2_Current::jacobian_2d_alpha0( Element * aElement, Matrix< real > & aK )
        {
            aElement->get_node_coors( mGroup->node_coords() );

            // reset the matrix
            aK.fill( 0.0 );

            // integration weights
            const Vector< real > & tW = mGroup->integration()->weights() ;

            // from the geometry jacobian
            Matrix< real > & tJ = mGroup->work_J() ;

            for( uint k=0; k<mGroup->integration()->number_of_integration_points(); ++k )
            {
                tJ = mGroup->integration()->dNdXi( k ) * mGroup->node_coords() ;

                const Matrix< real > & tN = mGroup->integration()->N( k );

                aK += ( tW( k ) *  std::abs( det( tJ ) ) )
                      * ( trans( tN ) * tN  );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2_Current::jacobian_2d_alpha1( Element * aElement, Matrix< real > & aK )
        {
            aElement->get_node_coors( mGroup->node_coords() );

            // reset the matrix
            aK.fill( 0.0 );

            // integration weights
            const Vector< real > & tW = mGroup->integration()->weights() ;

            // from the geometry jacobian
            Matrix< real > & tJ = mGroup->work_J() ;
            Matrix< real > & tB = mGroup->work_B() ;

            real tDetJ ;

            for( uint k=0; k<mGroup->integration()->number_of_integration_points(); ++k )
            {
                tJ = mGroup->integration()->dNdXi( k ) * mGroup->node_coords() ;
                tDetJ = std::abs( det( tJ ) );

                const Matrix< real > & tN = mGroup->integration()->N( k );

                tB = inv( tJ ) * mGroup->integration()->dNdXi( k );

                aK += ( tW( k ) * tDetJ )
                      * ( trans( tN ) * tN + mAlpha * tDetJ * trans( tB ) * tB );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2_Current::rhs_2d(Element * aElement, Vector< real > & aRHS )
        {
            // link edge function with element
            mEdgeFunction->link( aElement, true, false, false );

            // integration weights
            const Vector< real > & tW = mGroup->integration()->weights() ;

            Vector< real > & tH = mGroup->work_nedelec();
            this->collect_nedelec_data( aElement, NedelecField::H, tH );

            aRHS.fill( 0.0 );
            for( uint k=0; k<mGroup->integration()->number_of_integration_points(); ++k )
            {
                const Matrix< real > & tC = mEdgeFunction->C( k );

                aRHS += ( tW( k ) * mEdgeFunction->abs_det_J() )
                        * trans( mGroup->integration()->N( k ) ) * tC * tH ;
            }
        }

//------------------------------------------------------------------------------
    }
}