//
// Created by christian on 2/16/23.
//

#include "cl_IWG_Maxwell_2D.hpp"
#include "fn_trans.hpp"
#include "FEM_geometry.hpp"
#include "fn_inv.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"
namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Maxwell_2D::IWG_Maxwell_2D(
                const IwgType aType,
                const IwgMode aMode,
                const SymmetryMode aSymmetryMode,
                const SideSetDofLinkMode aSideSetDofLinkMode,
                const bool         aUseEdges ) :
                IWG_TimestepOld(
                        aType,
                        aMode,
                        aSymmetryMode,
                        DofMode::BlockSpecific,
                        aSideSetDofLinkMode ),
                mNumberOfDimensions( 2 ),
                mUseEdges( aUseEdges ),
                mFields( mDofMap, mDofFields, mOtherFields, mHiddenFields, mAllFields )
        {
            mNumberOfSpatialDimensions = mNumberOfDimensions ;
            mNumberOfDerivativeDimensions =  mNumberOfDimensions ;

            // make sure that this is a backwards implicit scheme
            // BELFEM_ERROR( this->theta() == 1.0,
            //               "Maxwell IWGs must be Backward Implicit (theta=1)" );

            // flag telling that we run the computation on the blocks
            mComputeJacobianOnBlock   = true ;

            // flag telling that we run the computations on the sidesets
            mComputeJacobianOnSideset = true ;
        }


        void
        IWG_Maxwell_2D::allocate_work_matrices( Group * aGroup )
        {
            IWG::allocate_work_matrices( aGroup );

            aGroup->work_invJ().set_size( mNumberOfDimensions, mNumberOfDimensions );
            aGroup->work_J().set_size( mNumberOfDimensions, mNumberOfDimensions );

            // check matrix function for superconductors
            switch ( aGroup->domain_type() )
            {
                case( DomainType::Air ) :
                {
                    aGroup->work_phi().set_size( mNumberOfNodesPerElement );
                    aGroup->work_psi().set_size( mNumberOfNodesPerElement );
                    aGroup->work_B().set_size( mNumberOfDimensions,
                                               mNumberOfNodesPerElement );

                }
                default :
                {
                    BELFEM_ERROR( false,
                                  "Unknown Domain type for %s %lu",
                                  aGroup->type() == GroupType::BLOCK ? "block" : "sideset",
                                  ( long unsigned int ) aGroup->id() );
                }
            }
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        IWG_Maxwell_2D::B_tri3( Element * aElement )
        {
            Matrix< real > & aB = mGroup->work_B() ;
            const Matrix< real > & tInvJ = inv_J_tri3( mGroup );

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
        IWG_Maxwell_2D::vol_phi_gauss_tri3_mu0(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK,
                Vector< real > & aF )
        {
            // reset matrices
            aM.fill( 0.0 );
            aF.fill( 0.0 );

            // node coordinates of element
            this->collect_node_coords( aElement, mGroup->work_X() );

            // gradient matrix
            const Matrix< real > & tB = this->B_tri3( aElement );

            aK = ( 0.5 * constant::mu0 * this->abs_det_J() ) *
                 trans( tB ) * tB ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_2D::vol_phi_gauss_tri3_mut(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK,
                Vector< real > & aF )
        {
            // reset matrices
            aM.fill( 0.0 );
            aF.fill( 0.0 );

            // node coordinates of element
            this->collect_node_coords( aElement, mGroup->work_X() );

            // gradient matrix
            const Matrix< real > & tB = this->B_tri3( aElement );

            aK = ( 0.5 * mGroup->material()->mu_s( norm( tB* mGroup->work_phi()) )
                    * this->abs_det_J() ) *
                 trans( tB ) * tB ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_2D::vol_phi_faraday_tri3_mu0(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK,
                Vector< real > & aF )
        {
            // reset matrices
            aK.fill( 0.0 );
            aF.fill( 0.0 );

            // node coordinates of element
            this->collect_node_coords( aElement, mGroup->work_X() );

            // gradient matrix
            const Matrix< real > & tB = this->B_tri3( aElement );

            aM = ( 0.5 * constant::mu0 * this->abs_det_J() ) *
                 trans( tB ) * tB ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_2D::vol_phi_faraday_tri3_mut(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK,
                Vector< real > & aF )
        {
            // reset matrices
            aF.fill( 0.0 );

            // node coordinates of element
            this->collect_node_coords( aElement, mGroup->work_X() );

            // gradient matrix
            const Matrix< real > & tB = this->B_tri3( aElement );

            // mu value from last timestep
            real tMu0 = mGroup->material()->mu_s( norm( tB* mGroup->work_psi()) );
            real tMu1 = mGroup->material()->mu_s( norm( tB* mGroup->work_phi()) );

            aM = 0.5 * this->abs_det_J() * trans( tB ) * tB ;
            aK = ( tMu1 - tMu0 ) / mDeltaTime * aM ;
            aM *= tMu1 ;
        }



//------------------------------------------------------------------------------

        void
        IWG_Maxwell_2D::vol_phi_gauss_tri6_mut(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK,
                Vector< real > & aF )
        {
            // reset matrices
            aM.fill( 0.0 );
            aK.fill( 0.0 );
            aF.fill( 0.0 );

            // node coordinates of element
            this->collect_node_coords( aElement, mGroup->work_X() );

            Matrix< real > & tB = mGroup->work_B();

            // integration weights
            const Vector< real > & tW = mGroup->integration_weights() ;

            if( aElement->element()->is_curved() )
            {
                Matrix< real > & tJ = mGroup->work_J() ;

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // Jacobian Matrix
                    tJ = mGroup->integration()->dNdXi( k ) * mGroup->work_X() ;

                    // gradient operator
                    tB = inv_J_2d( mGroup ) * mGroup->integration()->dNdXi( k );

                    // contribution to mass matrix
                    aK += ( mGroup->material()->mu_s( dot( tB, mGroup->work_phi() ))
                            * tW( k ) * this->abs_det_J() ) * trans( tB ) * tB  ;
                }
            }
            else
            {
                // matrix iJ is constant
                const Matrix< real > & tinvJ = inv_J_tri3( mGroup );

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    tB = tinvJ * mGroup->integration()->dNdXi( k );

                    // contribution to mass matrix
                    aK += ( mGroup->material()->mu_s( norm( tB* mGroup->work_phi()) )
                            * tW( k ) ) * trans( tB ) * tB  ;
                }

                aK *= this->abs_det_J() ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_2D::vol_phi_gauss_tri6_mu0(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK,
                Vector< real > & aF )
        {
            // reset matrices
            aM.fill( 0.0 );
            aK.fill( 0.0 );
            aF.fill( 0.0 );

            // node coordinates of element
            this->collect_node_coords( aElement, mGroup->work_X() );

            Matrix< real > & tB = mGroup->work_B();

            // integration weights
            const Vector< real > & tW = mGroup->integration_weights() ;

            if( aElement->element()->is_curved() )
            {
                Matrix< real > & tJ = mGroup->work_J() ;

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // Jacobian Matrix
                    tJ = mGroup->integration()->dNdXi( k ) * mGroup->work_X() ;

                    // gradient operator
                    tB = inv_J_2d( mGroup ) * mGroup->integration()->dNdXi( k );
                    // contribution to mass matrix
                    aK += ( tW( k ) * this->abs_det_J() ) * trans( tB ) * tB  ;
                }
                aK *= constant::mu0 ;
            }
            else
            {
                // matrix iJ is constant
                const Matrix< real > & tinvJ = inv_J_tri3( mGroup );

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    tB = tinvJ * mGroup->integration()->dNdXi( k );

                    // contribution to stiffness matrix
                    aK += tW( k ) * trans( tB ) * tB  ;
                }

                aK *= constant::mu0 * this->abs_det_J() ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_2D::vol_phi_faraday_tri6_mu0(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK,
                Vector< real > & aF )
        {
            // reset matrices
            aM.fill( 0.0 );
            aK.fill( 0.0 );
            aF.fill( 0.0 );

            // node coordinates of element
            this->collect_node_coords( aElement, mGroup->work_X() );

            Matrix< real > & tB = mGroup->work_B();

            // integration weights
            const Vector< real > & tW = mGroup->integration_weights() ;

            if( aElement->element()->is_curved() )
            {
                Matrix< real > & tJ = mGroup->work_J() ;

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // Jacobian Matrix
                    tJ = mGroup->integration()->dNdXi( k ) * mGroup->work_X() ;

                    // gradient operator
                    tB = inv_J_2d( mGroup ) * mGroup->integration()->dNdXi( k );

                    // contribution to mass matrix
                    aM += tW( k ) * trans( tB ) * tB * this->abs_det_J();
                }
                aM *= constant::mu0 ;
            }
            else
            {
                // matrix iJ is constant
                const Matrix< real > & tinvJ = inv_J_tri3( mGroup );

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    tB = tinvJ * mGroup->integration()->dNdXi( k );

                    aM += tW( k ) * trans( tB ) * tB ;
                }
                aM *= constant::mu0 * this->abs_det_J() ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_2D::vol_phi_faraday_tri6_mut(
                Element * aElement,
                Matrix< real > & aM,
                Matrix< real > & aK,
                Vector< real > & aF )
        {
            // reset matrices
            aM.fill( 0.0 );
            aK.fill( 0.0 );
            aF.fill( 0.0 );

            // node coordinates of element
            this->collect_node_coords( aElement, mGroup->work_X() );

            Matrix< real > & tB = mGroup->work_B();

            // integration weights
            const Vector< real > & tW = mGroup->integration_weights() ;

            if( aElement->element()->is_curved() )
            {
                Matrix< real > & tJ = mGroup->work_J() ;

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // Jacobian Matrix
                    tJ = mGroup->integration()->dNdXi( k ) * mGroup->work_X() ;

                    // gradient operator
                    tB = inv_J_2d( mGroup ) * mGroup->integration()->dNdXi( k );

                    real tMu0 = mGroup->material()->mu_s( std::abs( dot( tB, mGroup->work_psi() ) ) );
                    real tMu1 = mGroup->material()->mu_s( std::abs( dot( tB, mGroup->work_phi() ) ) );
                    real tdMudt = ( tMu1 - tMu0 ) / mDeltaTime ;

                    // contribution to mass matrix
                    aM += ( tW( k ) * tMu1 * this->abs_det_J() ) * trans( tB ) * tB  ;

                    // contribution to stiffness matrix
                    aK += ( tW( k ) * tdMudt * this->abs_det_J() ) * trans( tB ) * tB  ;
                }
            }
            else
            {
                // matrix iJ is constant
                const Matrix< real > & tinvJ = inv_J_tri3( mGroup );

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    tB = tinvJ * mGroup->integration()->dNdXi( k );

                    real tMu0 = mGroup->material()->mu_s( std::abs( dot( tB, mGroup->work_psi() ) ) );
                    real tMu1 = mGroup->material()->mu_s( std::abs( dot( tB, mGroup->work_phi() ) ) );
                    real tdMudt = ( tMu1 - tMu0 ) / mDeltaTime ;

                    // contribution to mass matrix
                    aM += ( tW( k ) * tMu1 ) * trans( tB ) * tB  ;
                }

                aM *= this->abs_det_J() ;
                aK *= this->abs_det_J() ;
            }
        }


//------------------------------------------------------------------------------


} /* end namespace fem */
} /* end namespace belfem */