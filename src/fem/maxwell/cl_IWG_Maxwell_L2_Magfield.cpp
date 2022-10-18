//
// Created by christian on 3/23/22.
//

#include "cl_IWG_Maxwell_L2_Magfield.hpp"
#include "fn_det.hpp"
#include "fn_trans.hpp"
#include "fn_inv.hpp"

namespace belfem
{
    namespace fem
    {

//------------------------------------------------------------------------------

        IWG_Maxwell_L2_Magfield::IWG_Maxwell_L2_Magfield(
                const ElementType aElementType,
                const Magfield_L2_Mode aMode,
                const real            aAlpha  ) :
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
                    mEdgeDofMultiplicity = 1 ;
                    mFaceDofMultiplicity = 0 ;
                    mNumberOfRhsDofsPerEdge = 1 ;
                    mNumberOfRhsDofsPerFace = 0 ;
                    mNumberOfRhsEdgeDofsPerElement = 3 ;
                    break ;
                }
                case( ElementType::TRI6 ) :
                {
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

            if( mNumberOfDimensions == 2 )
            {
                mFields.Superconductor = { "bx", "by" };
                mFields.Air = { "bx", "by" };
                mFields.Ferro = { "bx", "by" };
                switch( aMode )
                {
                    case( Magfield_L2_Mode::HA ) :
                    {
                        mFunLink = & IWG_Maxwell_L2_Magfield::link_jacobian_2d_ha ;
                        break ;
                    }
                    case( Magfield_L2_Mode::HPHI ) :
                    {
                        mFields.Cut = { "lambda_x", "lambda_y" };
                        mLambdaDofMultiplicity = 1 ;

                        mFunLink = & IWG_Maxwell_L2_Magfield::link_jacobian_2d_hphi ;

                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Invalid projection mode");
                    }
                }
            }
            else
            {
                BELFEM_ERROR( false, "Invalid dimension");
            }
        }

//------------------------------------------------------------------------------

        IWG_Maxwell_L2_Magfield::~IWG_Maxwell_L2_Magfield()
        {

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2_Magfield::link_jacobian_2d_hphi( Group * aGroup )
        {
            mEdgeFunction->link( aGroup );

            switch ( aGroup->domain_type() )
            {
                case( DomainType::Conductor ) :
                {
                    mFunJacobian = mAlpha == 0 ?
                                   & IWG_Maxwell_L2_Magfield::compute_jacobian_2d_alpha0 :
                                   & IWG_Maxwell_L2_Magfield::compute_jacobian_2d_alpha1 ;

                    mFunRHS      = & IWG_Maxwell_L2_Magfield::compute_rhs_2d_h ;

                    break ;
                }
                case( DomainType::Air ) :
                {
                    mFunJacobian = mAlpha == 0 ?
                                   & IWG_Maxwell_L2_Magfield::compute_jacobian_2d_alpha0 :
                                   & IWG_Maxwell_L2_Magfield::compute_jacobian_2d_alpha1 ;

                    mFunRHS      = & IWG_Maxwell_L2_Magfield::compute_rhs_2d_phi ;

                    break ;
                }
                case( DomainType::Ferro ) :
                {
                    mFunJacobian = mAlpha == 0 ?
                                   &IWG_Maxwell_L2_Magfield::compute_jacobian_2d_alpha0 :
                                   &IWG_Maxwell_L2_Magfield::compute_jacobian_2d_alpha1;
#ifdef BELFEM_FERRO_HPHIA
                    mFunRHS = &IWG_Maxwell_L2_Magfield::compute_rhs_2d_a;
#else
                    mFunRHS      = & IWG_Maxwell_L2_Magfield::compute_rhs_2d_phi_ferro ;
#endif
                    break;
                }
                case( DomainType::Cut ) :
                {
                    mFunJacobian = & IWG_Maxwell_L2_Magfield::compute_jacobian_cut ;

                    mFunRHS = & IWG_Maxwell_L2_Magfield::compute_rhs_cut ;
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "invalid domain type");
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2_Magfield::link_jacobian_2d_ha( Group * aGroup )
        {
            mEdgeFunction->link( aGroup );

            switch ( aGroup->domain_type() )
            {
                case( DomainType::Conductor ) :
                {
                    mFunJacobian = mAlpha == 0 ?
                                   & IWG_Maxwell_L2_Magfield::compute_jacobian_2d_alpha0 :
                                   & IWG_Maxwell_L2_Magfield::compute_jacobian_2d_alpha1 ;

                    mFunRHS      = & IWG_Maxwell_L2_Magfield::compute_rhs_2d_h ;

                    break ;
                }
                case( DomainType::Ferro ) :
                case( DomainType::Air ) :
                case( DomainType::Coil ) :
                {
                    mFunJacobian = mAlpha == 0 ?
                                   & IWG_Maxwell_L2_Magfield::compute_jacobian_2d_alpha0 :
                                   & IWG_Maxwell_L2_Magfield::compute_jacobian_2d_alpha1;

                    mFunRHS = & IWG_Maxwell_L2_Magfield::compute_rhs_2d_a;

                    break;
                }
                default :
                {
                    BELFEM_ERROR( false, "invalid domain type");
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2_Magfield::allocate_work_matrices( Group * aGroup )
        {
            IWG::allocate_work_matrices( aGroup );

            aGroup->work_J().set_size( mNumberOfSpatialDimensions,
                                       mNumberOfSpatialDimensions, 0.0 );
            // mGroup->work_N().set_size( mNumberOfSpatialDimensions, mNumberOfDofsPerElement, 0.0 );

            if( mAlpha != 0.0 )
            {
                mGroup->work_dNdX().set_size( mNumberOfSpatialDimensions, mNumberOfNodesPerElement, 0.0 );
                mGroup->work_B().set_size( mNumberOfSpatialDimensions, mNumberOfDofsPerElement, 0.0 );
            }

            // for nodal scalar field
            mGroup->work_phi().set_size( mNumberOfNodesPerElement, 0.0 );

            // for nedelec field
            mGroup->work_nedelec().set_size( mNumberOfRhsEdgeDofsPerElement, 0.0 );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2_Magfield::compute_jacobian_2d_alpha0(
                Element * aElement, Matrix< real > & aJacobian )
        {
            aElement->get_node_coors( mGroup->node_coords() );

            // reset the matrix
            aJacobian.fill( 0.0 );

            // integration weights
            const Vector< real > & tW = mGroup->integration_weights() ;

            // from the geometry jacobian
            Matrix< real > & tJ = mGroup->work_J() ;

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                tJ = mGroup->dNdXi( k ) * mGroup->node_coords() ;
                const Matrix< real > & tN = mGroup->Nvector( k );

                aJacobian += ( tW( k ) * std::abs( det( tJ ) ) )
                      * trans( tN ) * tN ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2_Magfield::compute_jacobian_2d_alpha1(
                Element * aElement, Matrix< real > & aJacobian )
        {
            aElement->get_node_coors( mGroup->node_coords() );

            // reset the matrix
            aJacobian.fill( 0.0 );

            // integration weights
            const Vector< real > & tW = mGroup->integration_weights() ;

            // from the geometry jacobian
            Matrix< real > & tJ    = mGroup->work_J() ;
            Matrix< real > & tdNdx = mGroup->work_dNdX();
            Matrix< real > & tB    = mGroup->work_B() ;

            real tDetJ ;
            uint tCount ;

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                tJ = mGroup->dNdXi( k ) * mGroup->node_coords() ;
                tDetJ = std::abs( det( tJ ) );

                const Matrix< real > & tN = mGroup->Nvector( k );

                tdNdx = inv( tJ ) * mGroup->dNdXi( k );

                // build stiffness matrix
                tCount = 0 ;
                for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                {
                    tB( 0, tCount++ ) = tdNdx( 0, i );
                    tB( 1, tCount++ ) = tdNdx( 1, i );
                }


                aJacobian += ( tW( k ) * tDetJ )
                      * ( trans( tN ) * tN + mAlpha * tDetJ * trans( tB ) * tB );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2_Magfield::compute_jacobian_cut(
                Element        * aElement,
                Matrix< real > & aJacobian )
        {
            aJacobian.fill( 0.0 );
            aJacobian( 0, 4 ) =  1.0 ;
            aJacobian( 1, 4 ) = -1.0 ;
            aJacobian( 2, 5 ) =  1.0 ;
            aJacobian( 3, 5 ) = -1.0 ;
            aJacobian( 4, 0 ) =  1.0 ;
            aJacobian( 4, 1 ) = -1.0 ;
            aJacobian( 5, 2 ) =  1.0 ;
            aJacobian( 5, 3 ) = -1.0 ;
        }


//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2_Magfield::compute_rhs_2d_h(
                Element * aElement, Vector< real > & aRHS )
        {
            // link edge function with element
            mEdgeFunction->link( aElement, true, false, false );

            // integration weights
            const Vector< real > & tW = mGroup->integration_weights() ;

            Vector< real > & tH = mGroup->work_nedelec();
            this->collect_nedelec_data( aElement, NedelecField::H, tH );

            aRHS.fill( 0.0 );
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                aRHS += ( tW( k ) * mEdgeFunction->abs_det_J() )
                        * trans( mGroup->Nvector( k ) )
                        * mEdgeFunction->E( k ) * tH ;
            }

            aRHS *= mGroup->material()->mu_r() * constant::mu0 ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2_Magfield::compute_rhs_2d_phi(
                Element * aElement, Vector< real > & aRHS )
        {
            // link edge function with element
            mEdgeFunction->link( aElement, false, true, false );

            // integration weights
            const Vector< real > & tW = mGroup->integration_weights() ;

            Vector< real > & tPhi = mGroup->work_phi();
            this->collect_node_data( aElement, "phi", tPhi );

            aRHS.fill( 0.0 );
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                aRHS += ( tW( k ) * mEdgeFunction->abs_det_J() )
                        * trans( mGroup->Nvector( k ) )
                        * mEdgeFunction->B( k ) * tPhi ;
            }

            aRHS *= -constant::mu0 ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2_Magfield::compute_rhs_2d_phi_ferro(
                Element * aElement, Vector< real > & aRHS )
        {
            // link edge function with element
            mEdgeFunction->link( aElement, false, true, false );

            // integration weights
            const Vector< real > & tW = mGroup->integration_weights() ;

            Vector< real > & tPhi = mGroup->work_phi();
            this->collect_node_data( aElement, "phi", tPhi );

            aRHS.fill( 0.0 );
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                mH =mEdgeFunction->B( k ) * tPhi  ;

                real tMu = mMaterial->mu_s( norm( mH ) );
                aRHS -= ( tW( k ) * mEdgeFunction->abs_det_J() * tMu )
                        * trans( mGroup->Nvector( k ) )
                        * mH ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2_Magfield::compute_rhs_2d_a(
                Element * aElement, Vector< real > & aRHS )
        {
            // link edge function with element
            mEdgeFunction->link( aElement, false, false, true );

            // integration weights
            const Vector< real > & tW = mGroup->integration_weights() ;

            Vector< real > & tAz = mGroup->work_phi();
            this->collect_node_data( aElement, "az", tAz );

            aRHS.fill( 0.0 );
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                aRHS += ( tW( k ) * mEdgeFunction->abs_det_J() )
                        * trans( mGroup->Nvector( k ) )
                        * mEdgeFunction->CA( k ) * tAz ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_L2_Magfield::compute_rhs_cut(
                Element        * aElement,
                Vector< real > & aRHS )
        {
            aRHS.fill( 0.0 );
        }

//------------------------------------------------------------------------------

    }
}