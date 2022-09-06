//
// Created by Christian Messe on 27.01.22.
//

#include "cl_IWG_Maxwell_HPhi_Tet4.hpp"

#include "fn_trans.hpp"
#include "cl_Timer.hpp"
#include "fn_det.hpp"
#include "fn_inv.hpp"
#include "fn_crossmat.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Maxwell_HPhi_Tet4::IWG_Maxwell_HPhi_Tet4() :
        IWG_Maxwell( ElementType::TET4,
                     IwgType::MAXWELL_HPHI_TET4,
                     IwgMode::Iterative,
                     SymmetryMode::Unsymmetric,
                     SideSetDofLinkMode::MasterAndSlave,
                     true )
        {
            // dof fields
            mFields.Superconductor = { "edge_h" };
            mFields.Ferro = { "ax", "ay", "az" };
            mFields.Air = { "phi" };
            mFields.InterfaceScAir = { "lambda_x", "lambda_y", "lambda_z" };

            mFields.Cut = { "lambda" };

            // non-dof fields
            mFields.MagneticFieldDensity     = { "bx", "by" };
            mFields.CurrentDensity = {  "jz" };
            mFields.CurrentBC = { "lambda_I" };

            mEdgeDofMultiplicity = 1 ;
            mFaceDofMultiplicity = 0 ;
            mLambdaDofMultiplicity = 1 ;

            mNumberOfRhsDofsPerEdge = 1 ;
            mNumberOfRhsDofsPerFace = 0 ;

            mThinShellDofLinkMode = SideSetDofLinkMode::MasterAndSlave ;
        }

//------------------------------------------------------------------------------

       IWG_Maxwell_HPhi_Tet4::~IWG_Maxwell_HPhi_Tet4()
       {

       }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tet4::link_jacobian_function( Group * aGroup )
        {
            // link edge function with group
            mEdgeFunction->link( aGroup );

            switch( aGroup->domain_type() )
            {
                case( DomainType::SuperConductor ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tet4::compute_jacobian_and_rhs_superconductor ;
                    break ;
                }
                case( DomainType::FerroMagnetic ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tet4::compute_jacobian_and_rhs_ferro ;
                    break ;
                }
                case( DomainType::Air ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tet4::compute_jacobian_and_rhs_air ;
                    break ;
                }
                case( DomainType::InterfaceScAir ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tet4::compute_jacobian_and_rhs_scair ;
                    break ;
                }
                case( DomainType::InterfaceScFm ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tet4::compute_jacobian_and_rhs_scfm ;
                    break ;
                }
                case( DomainType::InterfaceFmAir ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tet4::compute_jacobian_and_rhs_fmair ;
                    break ;
                }
                case( DomainType::Cut ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tet4::compute_jacobian_and_rhs_cut ;
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid group type" );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tet4::compute_jacobian_and_rhs_scair(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // collect the node coordinates for the surface
            aElement->get_node_coors( mGroup->work_X() );

            const Vector< real > & tW = mGroup->integration_weights();
            const IntegrationData * tNodeFunction = this->interface_data( aElement );

            // face normal, also computes scaling expression
            const Vector< real > & tn = this->normal_straight_3d( aElement );

            // hack to make matrix less ill-conditioned
            mGroup->work_det_J() *= constant::mu0 ;

            // expression for gradient operatpor B
            Matrix< real >  & tB   = mGroup->work_B() ;

            mGroup->work_J() = tNodeFunction->dNdXi( 0 ) * mGroup->node_coords() ;
            mGroup->work_invJ() = inv( mGroup->work_J() );
            tB = mGroup->work_invJ() * tNodeFunction->dNdXi( 0 );
            tB *= mEdgeFunction->sum_w() * mGroup->work_det_J() ;

            // expression ∫ (n x B) dS
            Matrix< real  > & tnxB = mGroup->work_Sigma() ;
            crossmat( tn, tB, tnxB );

            // expression ∫ (n x E) dS
            Matrix< real >  & tnxE = mGroup->work_Tau() ;
            tnxE.fill( 0.0 );

            // expression ∫ ( N' · n ' · E) dS
            Matrix< real >  & tNnE = mGroup->work_Chi() ;
            tNnE.fill( 0.0 );

            real tWork ;

            // perform integration
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // get the node function
                const Vector< real > & tN = tNodeFunction->phi( k );

                // get the edge function
                const Matrix < real > & tE = mEdgeFunction->E( k );

                // for expression ∫ (n x E) dS
                crossmat( tn, tE, tW( k ), tnxE );

                for( uint j=0; j<6; ++j )
                {
                    tWork = tn( 0 ) * tE( 0, j )
                             + tn( 1 ) * tE( 1, j )
                             + tn( 2 ) * tE( 2, j );
                    for( uint i=0; i<4; ++i )
                    {
                        tNnE( i, j ) += tW( k ) * tN( i ) * tWork ;
                    }
                }
            }
            tnxE *=  mGroup->work_det_J() ;
            tNnE *= -mGroup->work_det_J() ; // <-- don't forget the times minus one!

            // manual assembly of matrix
            aRHS.fill( 0.0 );
            aJacobian.fill( 0.0 );
            Matrix< real > & tK = aJacobian ;

            uint J ;
            for( uint j=0; j<6; ++j )
            {
                tK(  6,  j ) = tNnE( 0, j );
                tK(  7,  j ) = tNnE( 1, j );
                tK(  8,  j ) = tNnE( 2, j );
                tK(  9,  j ) = tNnE( 3, j );
                tK( 10,  j ) = tnxE( 0, j );
                tK( 11,  j ) = tnxE( 1, j );
                tK( 12,  j ) = tnxE( 2, j );
            }
            for( uint j=0; j<4; ++j )
            {
                J = j + 6 ;
                tK( 10, J ) = tnxB( 0, j );
                tK( 11, J ) = tnxB( 1, j );
                tK( 12, J ) = tnxB( 2, j );
            }
            for( uint j=0; j<3; ++j )
            {
                J = j + 10 ;
                tK( 0, J ) = tnxE( j, 0 );
                tK( 1, J ) = tnxE( j, 1 );
                tK( 2, J ) = tnxE( j, 2 );
                tK( 3, J ) = tnxE( j, 3 );
                tK( 4, J ) = tnxE( j, 4 );
                tK( 5, J ) = tnxE( j, 5 );
                tK( 6, J ) = tnxB( j, 0 );
                tK( 7, J ) = tnxB( j, 1 );
                tK( 8, J ) = tnxB( j, 2 );
                tK( 9, J ) = tnxB( j, 3 );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tet4::compute_jacobian_and_rhs_fmair(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            const IntegrationData * tSlave =  mGroup->slave_integration(
                    aElement->facet()->slave_index() ) ;

            // compute Jacobian and B-function for air
            aElement->slave()->get_node_coors( mGroup->node_coords() );

            mGroup->work_J() = tSlave->dNdXi( 0 ) * mGroup->node_coords() ;

            mGroup->work_B() = inv( mGroup->work_J() ) * tSlave->dNdXi( 0 );

            Matrix< real > & tM    = aJacobian ;
            Matrix< real > & tK    = mGroup->work_K();
            const Vector< real > & tn    = this->normal_straight_3d( aElement );
            const Matrix< real > & tB    = mGroup->work_B() ;
            Matrix< real > & tnxB  = mGroup->work_Sigma() ;

            // again, we use a mathematical trick since n and B are constant
            crossmat( tn, tB, tnxB );
            tnxB *= mGroup->work_det_J() * constant::mu0 / 6.0 ;

            tM.fill( 0.0 );

            switch( aElement->facet()->master_index() )
            {
                case( 0 ) :
                {
                    tM( 12, 0 ) = tnxB( 0, 0 ) ;
                    tM( 13, 0 ) = tnxB( 0, 1 ) ;
                    tM( 14, 0 ) = tnxB( 0, 2 ) ;
                    tM( 15, 0 ) = tnxB( 0, 3 ) ;
                    tM( 12, 1 ) = tnxB( 0, 0 ) ;
                    tM( 13, 1 ) = tnxB( 0, 1 ) ;
                    tM( 14, 1 ) = tnxB( 0, 2 ) ;
                    tM( 15, 1 ) = tnxB( 0, 3 ) ;
                    tM( 12, 3 ) = tnxB( 0, 0 ) ;
                    tM( 13, 3 ) = tnxB( 0, 1 ) ;
                    tM( 14, 3 ) = tnxB( 0, 2 ) ;
                    tM( 15, 3 ) = tnxB( 0, 3 ) ;
                    tM( 12, 4 ) = tnxB( 1, 0 ) ;
                    tM( 13, 4 ) = tnxB( 1, 1 ) ;
                    tM( 14, 4 ) = tnxB( 1, 2 ) ;
                    tM( 15, 4 ) = tnxB( 1, 3 ) ;
                    tM( 12, 5 ) = tnxB( 1, 0 ) ;
                    tM( 13, 5 ) = tnxB( 1, 1 ) ;
                    tM( 14, 5 ) = tnxB( 1, 2 ) ;
                    tM( 15, 5 ) = tnxB( 1, 3 ) ;
                    tM( 12, 7 ) = tnxB( 1, 0 ) ;
                    tM( 13, 7 ) = tnxB( 1, 1 ) ;
                    tM( 14, 7 ) = tnxB( 1, 2 ) ;
                    tM( 15, 7 ) = tnxB( 1, 3 ) ;
                    tM( 12, 8 ) = tnxB( 2, 0 ) ;
                    tM( 13, 8 ) = tnxB( 2, 1 ) ;
                    tM( 14, 8 ) = tnxB( 2, 2 ) ;
                    tM( 15, 8 ) = tnxB( 2, 3 ) ;
                    tM( 12, 9 ) = tnxB( 2, 0 ) ;
                    tM( 13, 9 ) = tnxB( 2, 1 ) ;
                    tM( 14, 9 ) = tnxB( 2, 2 ) ;
                    tM( 15, 9 ) = tnxB( 2, 3 ) ;
                    tM( 12,11 ) = tnxB( 2, 0 ) ;
                    tM( 13,11 ) = tnxB( 2, 1 ) ;
                    tM( 14,11 ) = tnxB( 2, 2 ) ;
                    tM( 15,11 ) = tnxB( 2, 3 ) ;
                    break ;
                }
                case( 1 ) :
                {
                    tM( 12, 1 ) = tnxB( 0, 0 ) ;
                    tM( 13, 1 ) = tnxB( 0, 1 ) ;
                    tM( 14, 1 ) = tnxB( 0, 2 ) ;
                    tM( 15, 1 ) = tnxB( 0, 3 ) ;
                    tM( 12, 2 ) = tnxB( 0, 0 ) ;
                    tM( 13, 2 ) = tnxB( 0, 1 ) ;
                    tM( 14, 2 ) = tnxB( 0, 2 ) ;
                    tM( 15, 2 ) = tnxB( 0, 3 ) ;
                    tM( 12, 3 ) = tnxB( 0, 0 ) ;
                    tM( 13, 3 ) = tnxB( 0, 1 ) ;
                    tM( 14, 3 ) = tnxB( 0, 2 ) ;
                    tM( 15, 3 ) = tnxB( 0, 3 ) ;
                    tM( 12, 5 ) = tnxB( 1, 0 ) ;
                    tM( 13, 5 ) = tnxB( 1, 1 ) ;
                    tM( 14, 5 ) = tnxB( 1, 2 ) ;
                    tM( 15, 5 ) = tnxB( 1, 3 ) ;
                    tM( 12, 6 ) = tnxB( 1, 0 ) ;
                    tM( 13, 6 ) = tnxB( 1, 1 ) ;
                    tM( 14, 6 ) = tnxB( 1, 2 ) ;
                    tM( 15, 6 ) = tnxB( 1, 3 ) ;
                    tM( 12, 7 ) = tnxB( 1, 0 ) ;
                    tM( 13, 7 ) = tnxB( 1, 1 ) ;
                    tM( 14, 7 ) = tnxB( 1, 2 ) ;
                    tM( 15, 7 ) = tnxB( 1, 3 ) ;
                    tM( 12, 9 ) = tnxB( 2, 0 ) ;
                    tM( 13, 9 ) = tnxB( 2, 1 ) ;
                    tM( 14, 9 ) = tnxB( 2, 2 ) ;
                    tM( 15, 9 ) = tnxB( 2, 3 ) ;
                    tM( 12,10 ) = tnxB( 2, 0 ) ;
                    tM( 13,10 ) = tnxB( 2, 1 ) ;
                    tM( 14,10 ) = tnxB( 2, 2 ) ;
                    tM( 15,10 ) = tnxB( 2, 3 ) ;
                    tM( 12,11 ) = tnxB( 2, 0 ) ;
                    tM( 13,11 ) = tnxB( 2, 1 ) ;
                    tM( 14,11 ) = tnxB( 2, 2 ) ;
                    tM( 15,11 ) = tnxB( 2, 3 ) ;
                    break ;
                }
                case( 2 ) :
                {
                    tM( 12, 0 ) = tnxB( 0, 0 ) ;
                    tM( 13, 0 ) = tnxB( 0, 1 ) ;
                    tM( 14, 0 ) = tnxB( 0, 2 ) ;
                    tM( 15, 0 ) = tnxB( 0, 3 ) ;
                    tM( 12, 2 ) = tnxB( 0, 0 ) ;
                    tM( 13, 2 ) = tnxB( 0, 1 ) ;
                    tM( 14, 2 ) = tnxB( 0, 2 ) ;
                    tM( 15, 2 ) = tnxB( 0, 3 ) ;
                    tM( 12, 3 ) = tnxB( 0, 0 ) ;
                    tM( 13, 3 ) = tnxB( 0, 1 ) ;
                    tM( 14, 3 ) = tnxB( 0, 2 ) ;
                    tM( 15, 3 ) = tnxB( 0, 3 ) ;
                    tM( 12, 4 ) = tnxB( 1, 0 ) ;
                    tM( 13, 4 ) = tnxB( 1, 1 ) ;
                    tM( 14, 4 ) = tnxB( 1, 2 ) ;
                    tM( 15, 4 ) = tnxB( 1, 3 ) ;
                    tM( 12, 6 ) = tnxB( 1, 0 ) ;
                    tM( 13, 6 ) = tnxB( 1, 1 ) ;
                    tM( 14, 6 ) = tnxB( 1, 2 ) ;
                    tM( 15, 6 ) = tnxB( 1, 3 ) ;
                    tM( 12, 7 ) = tnxB( 1, 0 ) ;
                    tM( 13, 7 ) = tnxB( 1, 1 ) ;
                    tM( 14, 7 ) = tnxB( 1, 2 ) ;
                    tM( 15, 7 ) = tnxB( 1, 3 ) ;
                    tM( 12, 8 ) = tnxB( 2, 0 ) ;
                    tM( 13, 8 ) = tnxB( 2, 1 ) ;
                    tM( 14, 8 ) = tnxB( 2, 2 ) ;
                    tM( 15, 8 ) = tnxB( 2, 3 ) ;
                    tM( 12,10 ) = tnxB( 2, 0 ) ;
                    tM( 13,10 ) = tnxB( 2, 1 ) ;
                    tM( 14,10 ) = tnxB( 2, 2 ) ;
                    tM( 15,10 ) = tnxB( 2, 3 ) ;
                    tM( 12,11 ) = tnxB( 2, 0 ) ;
                    tM( 13,11 ) = tnxB( 2, 1 ) ;
                    tM( 14,11 ) = tnxB( 2, 2 ) ;
                    tM( 15,11 ) = tnxB( 2, 3 ) ;
                    break ;
                }
                case( 3 ) :
                {
                    tM( 12, 0 ) = tnxB( 0, 0 ) ;
                    tM( 13, 0 ) = tnxB( 0, 1 ) ;
                    tM( 14, 0 ) = tnxB( 0, 2 ) ;
                    tM( 15, 0 ) = tnxB( 0, 3 ) ;
                    tM( 12, 1 ) = tnxB( 0, 0 ) ;
                    tM( 13, 1 ) = tnxB( 0, 1 ) ;
                    tM( 14, 1 ) = tnxB( 0, 2 ) ;
                    tM( 15, 1 ) = tnxB( 0, 3 ) ;
                    tM( 12, 2 ) = tnxB( 0, 0 ) ;
                    tM( 13, 2 ) = tnxB( 0, 1 ) ;
                    tM( 14, 2 ) = tnxB( 0, 2 ) ;
                    tM( 15, 2 ) = tnxB( 0, 3 ) ;
                    tM( 12, 4 ) = tnxB( 1, 0 ) ;
                    tM( 13, 4 ) = tnxB( 1, 1 ) ;
                    tM( 14, 4 ) = tnxB( 1, 2 ) ;
                    tM( 15, 4 ) = tnxB( 1, 3 ) ;
                    tM( 12, 5 ) = tnxB( 1, 0 ) ;
                    tM( 13, 5 ) = tnxB( 1, 1 ) ;
                    tM( 14, 5 ) = tnxB( 1, 2 ) ;
                    tM( 15, 5 ) = tnxB( 1, 3 ) ;
                    tM( 12, 6 ) = tnxB( 1, 0 ) ;
                    tM( 13, 6 ) = tnxB( 1, 1 ) ;
                    tM( 14, 6 ) = tnxB( 1, 2 ) ;
                    tM( 15, 6 ) = tnxB( 1, 3 ) ;
                    tM( 12, 8 ) = tnxB( 2, 0 ) ;
                    tM( 13, 8 ) = tnxB( 2, 1 ) ;
                    tM( 14, 8 ) = tnxB( 2, 2 ) ;
                    tM( 15, 8 ) = tnxB( 2, 3 ) ;
                    tM( 12, 9 ) = tnxB( 2, 0 ) ;
                    tM( 13, 9 ) = tnxB( 2, 1 ) ;
                    tM( 14, 9 ) = tnxB( 2, 2 ) ;
                    tM( 15, 9 ) = tnxB( 2, 3 ) ;
                    tM( 12,10 ) = tnxB( 2, 0 ) ;
                    tM( 13,10 ) = tnxB( 2, 1 ) ;
                    tM( 14,10 ) = tnxB( 2, 2 ) ;
                    tM( 15,10 ) = tnxB( 2, 3 ) ;
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid Facet index");
                };
            }

            tK = -trans( tM );
            aRHS = tM * this->collect_q0_aphi_3d( aElement );
            aJacobian += tK * this->timestep() ;
        }

//------------------------------------------------------------------------------
    }
}
