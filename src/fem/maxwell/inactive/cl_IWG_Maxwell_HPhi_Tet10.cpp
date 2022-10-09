//
// Created by christian on 12/17/21.
//

#include "cl_IWG_Maxwell_HPhi_Tet10.hpp"
#include "fn_trans.hpp"
#include "fn_inv.hpp"
#include "fn_crossmat.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Maxwell_HPhi_Tet10::IWG_Maxwell_HPhi_Tet10() :
                IWG_Maxwell( ElementType::TET10,
                                 IwgType::MAXWELL_HPHI_TET10 ,
                                 IwgMode::Iterative,
                                 SymmetryMode::Unsymmetric,
                                 SideSetDofLinkMode::MasterAndSlave,
                                 true )
        {
            // dof fields
            mFields.Superconductor = { "edge_h", "face_h" };
            mFields.Air = { "phi" };
            mFields.Ferro = { "ax", "ay", "az" };
            mFields.Cut = { "lambda" };

            // non-dof fields
            mFields.MagneticFieldDensity = { "bx", "by", "bz" };
            mFields.CurrentDensity = { "jx", "jy", "jz" };
            mFields.CurrentBC = { "lambda_I" };

            mFields.InterfaceScAir = { "lambda_x", "lambda_y", "lambda_z" };

            mNumberOfRhsDofsPerEdge = 2 ;
            mNumberOfRhsDofsPerFace = 2 ;

            mEdgeDofMultiplicity    = 2 ;
            mFaceDofMultiplicity    = 2 ;
            mLambdaDofMultiplicity = 1 ;

            // mNumberOfEdgeDofsPerElement = 20 ;
            // mNumberOfRhsEdgeDofsPerElement = 20 ;ve ;
        }

//------------------------------------------------------------------------------

        IWG_Maxwell_HPhi_Tet10::~IWG_Maxwell_HPhi_Tet10()
        {

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tet10::link_jacobian_function( Group * aGroup )
        {
            // link edge function with group
            mEdgeFunction->link( aGroup );

            switch( aGroup->domain_type() )
            {
                case( DomainType::Conductor ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tet10::compute_jacobian_and_rhs_superconductor ;
                    break ;
                }
                case( DomainType::Air ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tet10::compute_jacobian_and_rhs_air ;
                    break ;
                }
                case( DomainType::Ferro ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tet10::compute_jacobian_and_rhs_ferro ;
                    break ;
                }
                case( DomainType::InterfaceScAir  ):
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tet10::compute_jacobian_and_rhs_scair ;
                    break ;
                }
                case( DomainType::InterfaceScFm ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tet10::compute_jacobian_and_rhs_scfm ;
                    break ;
                }
                case( DomainType::InterfaceFmAir ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tet10::compute_jacobian_and_rhs_fmair ;
                    break ;
                }
                case( DomainType::Cut ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tet10::compute_jacobian_and_rhs_cut ;
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
        IWG_Maxwell_HPhi_Tet10::compute_jacobian_and_rhs_scair(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // collect node coordinates of surface element
            aElement->get_node_coors( mGroup->work_X() );

            const Vector< real > & tW = mGroup->integration_weights();
            const IntegrationData * tNodeFunction = this->interface_data( aElement );

            // reset matrices
            aJacobian.fill( 0.0 );
            aRHS.fill( 0.0 );

            Matrix< real > & tK = aJacobian ;

            // - - - - - - - - - - grab work matrices

            // geometry Jacobian (transposed)
            Matrix< real > & tJ = mGroup->work_J() ;

            // inverse of geometry Jacobian (transposed)
            Matrix< real > & tInvJ = mGroup->work_invJ() ;

            // the gradient operator matrix
            Matrix< real > & tB = mGroup->work_B() ;

            // expression for n x E
            Matrix< real >  & tnxE = mGroup->work_Tau() ;

            // expression for n x B
            Matrix< real > & tnxB = mGroup->work_Sigma() ;

            // expression for ( N' · n ' · E)
            // Gradient operator matrix

            tnxE.fill( 0.0 );
            tnxB.fill( 0.0 );

            // scaling parameter
            real tOmega ;
            real tBeta ;

            if( aElement->element()->is_curved() )
            {

                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // get the normal vector and compute integration element
                    // ( also computes dS )
                    const Vector < real > & tn = this->normal_curved_3d( aElement, k );

                    // compute scaling parameter
                    tOmega = tW( k ) * constant::mu0 * mGroup->work_det_J() ;

                    // get edge function
                    const Matrix < real > & tE = mEdgeFunction->E( k );

                    // add component to container
                    crossmat( tn, tE, tOmega, tnxE );

                    // compute the geometry Jacobian (transposed)
                    tJ = tNodeFunction->dNdXi( k ) * mGroup->node_coords() ;

                    // compute the inverse
                    tInvJ = inv( tJ );

                    // compute the gradient operator
                    tB = tInvJ * tNodeFunction->dNdXi( k );

                    // add to component
                    crossmat( tn, tB, tOmega , tnxB );

                    // get node function
                    const Vector< real > & tN = tNodeFunction->phi( k );

                    // expression for N'n'E
                    for( uint j=0; j<20; ++j )
                    {
                        tBeta =   tn( 0 ) * tE( 0, j )
                                + tn( 1 ) * tE( 1, j )
                                + tn( 2 ) * tE( 2, j );
                        tBeta *= tOmega ;

                        for( uint i=0; i<10; ++i )
                        {
                            tK( i+20, j ) -= tBeta * tN( i );
                        }
                    }
                }
            }
            else
            {
                // get the normal vector and compute integration element
                // ( also computes dS )
                const Vector < real > & tn = this->normal_straight_3d( aElement );

                // compute the geometry Jacobian (transposed)
                tJ = tNodeFunction->dNdXi( 0 ) * mGroup->node_coords() ;

                // compute the inverse
                tInvJ = inv( tJ );

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // compute scaling parameter
                    tOmega = tW( k ) * constant::mu0 * mGroup->work_det_J() ;

                    // get edge function
                    const Matrix < real > & tE = mEdgeFunction->E( k );


                    // add component to container
                    crossmat( tn, tE, tOmega, tnxE );

                    // compute the gradient operator
                    tB = tInvJ * tNodeFunction->dNdXi( k );

                    // add component to container
                    crossmat( tn, tB, tOmega , tnxB );

                    // get node function
                    const Vector< real > & tN = tNodeFunction->phi( k );

                    // expression for N'n'E
                    for( uint j=0; j<20; ++j )
                    {
                        tBeta =     tn( 0 ) * tE( 0, j )
                                  + tn( 1 ) * tE( 1, j )
                                  + tn( 2 ) * tE( 2, j );
                        tBeta *= tOmega ;

                        for( uint i=0; i<10; ++i )
                        {
                            tK( i+20, j ) -= tBeta * tN( i );
                        }
                    }
                }
            }

            // finalize matrix
            uint J ;
            for( uint j=0; j<20; ++j )
            {
                tK( 30, j ) = tnxE( 0, j );
                tK( 31, j ) = tnxE( 1, j );
                tK( 32, j ) = tnxE( 2, j );
            }
            for( uint j=0; j<10; ++j )
            {
                J = j + 20 ;
                tK( 30, J ) = tnxB( 0, j );
                tK( 31, J ) = tnxB( 1, j );
                tK( 32, J ) = tnxB( 2, j );
            }
            for( uint j = 0; j<3; ++j )
            {
                J = j + 30 ;
                for( uint i=0; i<20; ++i )
                {
                    tK( i, J ) = tnxE( j, i );
                }
                for( uint i=0; i<10; ++i )
                {
                    tK( i+20, J ) = tnxB( j, i );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tet10::compute_jacobian_and_rhs_fmair(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // shape function for ferro element
            const IntegrationData * tMaster =  mGroup->master_integration(
                    aElement->facet()->master_index() ) ;

            // shape function for air element
            const IntegrationData * tSlave =  mGroup->slave_integration(
                    aElement->facet()->slave_index() * 3 + aElement->facet()->orientation_on_slave() ) ;

            // get node coords, needed to compute geometry Jacobian
            aElement->slave()->get_node_coors( mGroup->node_coords() );

            // get the node coords for the surface
            aElement->get_node_coors( mGroup->work_X() );

            // the B-Matrix
            Matrix< real > & tK = mGroup->work_K() ;
            Matrix< real > & tM = aJacobian ;
            Matrix< real > & tJ = mGroup->work_J() ;
            Matrix< real > & tB = mGroup->work_B() ;
            Matrix< real > & tnxB = mGroup->work_Sigma() ;
            const Vector< real > & tW = mGroup->integration_weights() ;

            tK.fill( 0.0 );

            uint I ;  // row of element matrix
            uint J ;  // column of element matrix
            real tOmega ; // integration weight
            if( aElement->element()->is_curved() )
            {
                // grab node coords for normal vector
                aElement->get_node_coors( mGroup->work_X() );

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // compute the normal
                    const Vector< real > & tn = this->normal_curved_3d( aElement, k );

                    // compute the integration constant
                    tOmega = -tW( k ) * mGroup->work_det_J() ;

                    // get the shape function
                    const Vector< real > & tN = tMaster->phi( k );

                    tJ = tSlave->dNdXi( k ) * mGroup->node_coords() ;
                    tB = inv( tJ ) * tSlave->dNdXi( k );
                    crossmat( tn, tB, tnxB );

                    for( uint j=0; j<10; ++j )
                    {
                        J = j + 30 ;
                        I = 0 ;
                        for ( uint d=0; d<3; ++d )
                        {
                            for ( uint i = 0; i < 10; ++i )
                            {
                                tK( I++, J ) += tOmega * tN( i ) * tnxB( d, j );
                            }
                        }
                    }
                }
            }
            else
            {
                // get the normal
                const Vector< real > & tn = this->normal_straight_3d( aElement );

                // compute the Jacobian
                tJ = tSlave->dNdXi( 0 ) * mGroup->node_coords() ;
                Matrix< real > & tInvJ = mGroup->work_invJ();
                tInvJ = inv( tJ );

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    const Vector< real > & tN = tMaster->phi( k );
                    tB = tInvJ * tSlave->dNdXi( k );
                    crossmat( tn, tB, tnxB );
                    tOmega = -tW( k ) * mGroup->work_det_J() ;

                    for( uint j=0; j<10; ++j )
                    {
                        J = j + 30 ;
                        I = 0 ;
                        for ( uint d=0; d<3; ++d )
                        {
                            for ( uint i = 0; i < 10; ++i )
                            {
                                tK( I++, J ) += tOmega * tN( i ) * tnxB( d, j );
                            }
                        }
                    }
                }
            }
            tM = - trans( tK );
            aRHS = tM * this->collect_q0_aphi_3d( aElement );
            aJacobian += tK * this->timestep() ;
        }

//------------------------------------------------------------------------------
    }
}