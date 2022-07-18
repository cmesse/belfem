//
// Created by christian on 3/23/22.
//

#include "cl_IWG_Maxwell_HA_Tri6.hpp"
#include "fn_trans.hpp"
#include "cl_FEM_MaxwellBoundaryConditionMagfield.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Maxwell_HA_Tri6::IWG_Maxwell_HA_Tri6() :
                IWG_Maxwell( ElementType::TRI6,
                             IwgType::MAXWELL_HA_TRI6,
                             IwgMode::Iterative,
                             SymmetryMode::Unsymmetric,
                             SideSetDofLinkMode::MasterAndSlave,
                             true )
        {
            // dof fields
            mFields.Superconductor = { "edge_h", "face_h" };
            mFields.Ferro = { "az" };
            mFields.Air = { "az" };
            mFields.Coil = { "az" };
            mFields.BoundaryAir = { "lambda" };

            // non-dof fields
            mFields.MagneticFieldDensity =  { "bx", "by", "bz" };
            mFields.CurrentDensity = { "jz" };
            mFields.CurrentBC = { "lambda_I" };

            mEdgeDofMultiplicity = 2 ;
            mFaceDofMultiplicity = 2 ;
            mLambdaDofMultiplicity = 1 ;

            mNumberOfRhsDofsPerEdge = 2 ;
            mNumberOfRhsDofsPerFace = 2 ;

        }
//------------------------------------------------------------------------------

        IWG_Maxwell_HA_Tri6::~IWG_Maxwell_HA_Tri6()
        {

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HA_Tri6::link_jacobian_function( Group * aGroup )
        {
            // link edge function with group
            mEdgeFunction->link( aGroup );

            switch( aGroup->domain_type() )
            {
                case( DomainType::SuperConductor ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HA_Tri6::compute_jacobian_and_rhs_superconductor ;
                    break ;
                }
                case( DomainType::Air ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HA_Tri6::compute_jacobian_and_rhs_air ;
                    break ;
                }
                case( DomainType::Coil ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HA_Tri6::compute_jacobian_and_rhs_coil ;
                    break ;
                }
                case( DomainType::FerroMagnetic ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HA_Tri6::compute_jacobian_and_rhs_ferro ;
                    break ;
                }
                case( DomainType::InterfaceScAir ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HA_Tri6::compute_jacobian_and_rhs_interface ;
                    break ;
                }
                case( DomainType::InterfaceScFm ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HA_Tri6::compute_jacobian_and_rhs_interface ;
                    break ;
                }
                case( DomainType::InterfaceFmAir ) :
                {
                    BELFEM_ERROR( false, "Ferro-Air-Interface must be blacklisted for H-A formulation" );
                    mFunJacobian = nullptr ;
                    break ;
                }
                case( DomainType::Boundary ) :
                {
                    BELFEM_ERROR( false, "Boundary conditions are strong for H-A formulation" );
                    mFunJacobian = nullptr ;
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid group type" );
                }
            }
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HA_Tri6::compute_jacobian_and_rhs_air(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {

            // link edge function with element
            mEdgeFunction->link( aElement, false, false, true );


            aRHS.fill( 0.0 );
            aJacobian.fill( 0.0 );

            // get the integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // get curl matrix ( constant for this element )
                const Matrix< real > & tC = mEdgeFunction->CA( k );

                aJacobian += ( tW( k ) * mEdgeFunction->abs_det_J() ) *
                        trans( tC ) * tC ;
            }
            aJacobian *= mDeltaTime * constant::nu0 ;
        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HA_Tri6::compute_jacobian_and_rhs_coil(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {

            // link edge function with element
            mEdgeFunction->link( aElement, false, false, true );

            // get the integration weights
            const Vector< real > & tW = mGroup->integration_weights() ;

            // get currents
            Vector< real > & tJz = mGroup->work_phi() ;


            this->collect_node_data( aElement, "jz", tJz );

            aJacobian.fill( 0.0 );
            aRHS.fill( 0.0 );

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // get curl matrix ( constant for this element )
                const Matrix< real > & tC = mEdgeFunction->CA( k );

                aJacobian += ( tW( k ) * mEdgeFunction->abs_det_J() ) *
                             trans( tC ) * tC ;
            }
            aJacobian *= mDeltaTime * constant::nu0 ;

            // compute the right hand side
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                aRHS += tW( k ) * mGroup->n( k ) * dot( mGroup->n( k ), tJz ) ;
            }
            aRHS *= mEdgeFunction->abs_det_J() * mDeltaTime ;
        }

//------------------------------------------------------------------------------
    }
}