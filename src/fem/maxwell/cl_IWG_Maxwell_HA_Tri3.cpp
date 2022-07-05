//
// Created by Christian Messe on 03.03.22.
//

#include "cl_IWG_Maxwell_HA_Tri3.hpp"
#include "fn_trans.hpp"
#include "cl_FEM_MaxwellBoundaryConditionMagfield.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Maxwell_HA_Tri3::IWG_Maxwell_HA_Tri3() :
                IWG_Maxwell( ElementType::TRI3,
                             IwgType::MAXWELL_HA_TRI3,
                             IwgMode::Iterative,
                             SymmetryMode::Unsymmetric,
                             SideSetDofLinkMode::MasterAndSlave,
                             true )
        {
            // dof fields
            mFields.Superconductor = { "edge_h" };
            mFields.Ferro = { "az" };
            mFields.Air = { "az" };
            mFields.Coil = { "az" };
            mFields.BoundaryAir = { "lambda" };

            // non-dof fields
            mFields.MagneticFieldDensity     = { "bx", "by" };
            mFields.CurrentDensity = { "jz" };
            mFields.CurrentBC = { "lambda_I" };

            mEdgeDofMultiplicity = 1 ;
            mFaceDofMultiplicity = 0 ;
            mLambdaDofMultiplicity = 1 ;

            mNumberOfRhsDofsPerEdge = 1 ;
            mNumberOfRhsDofsPerFace = 0 ;

            mNumberOfLayersPerShell = 3 ;

        }
//------------------------------------------------------------------------------

        IWG_Maxwell_HA_Tri3::~IWG_Maxwell_HA_Tri3()
        {

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HA_Tri3::link_jacobian_function( Group * aGroup )
        {
            // link edge function with group
            mEdgeFunction->link( aGroup );

            switch( aGroup->domain_type() )
            {
                case( DomainType::SuperConductor ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HA_Tri3::compute_jacobian_and_rhs_superconductor ;
                    break ;
                }
                case( DomainType::Air ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HA_Tri3::compute_jacobian_and_rhs_air ;
                    break ;
                }
                case( DomainType::Coil ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HA_Tri3::compute_jacobian_and_rhs_coil ;
                    break ;
                }
                case( DomainType::FerroMagnetic ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HA_Tri3::compute_jacobian_and_rhs_ferro ;
                    break ;
                }
                case( DomainType::InterfaceScAir ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HA_Tri3::compute_jacobian_and_rhs_interface ;
                    break ;
                }
                case( DomainType::InterfaceScFm ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HA_Tri3::compute_jacobian_and_rhs_interface ;
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
        IWG_Maxwell_HA_Tri3::compute_jacobian_and_rhs_air(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {

            // link edge function with element
            mEdgeFunction->link( aElement, false, false, true );

            // get curl matrix ( constant for this element )
            const Matrix< real > & tC = mEdgeFunction->CA();

            aJacobian = ( mEdgeFunction->sum_w() * mEdgeFunction->abs_det_J() * constant::nu0 * mDeltaTime )
                    * trans( tC ) * tC ;

            aRHS.fill( 0.0 );

        }

//------------------------------------------------------------------------------

        inline void
        IWG_Maxwell_HA_Tri3::compute_jacobian_and_rhs_coil(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {

            // link edge function with element
            mEdgeFunction->link( aElement, false, false, true );

            // get the integration weights
            const Vector< real > & tW     = mGroup->integration_weights() ;

            // get currents
            Vector< real > & tJz     = mGroup->work_phi() ;


            this->collect_node_data( aElement, "jz", tJz );

            // get curl matrix ( constant for this element )
            const Matrix< real > & tC = mEdgeFunction->CA();

            // the matrix is the same as for air
            aJacobian = ( mEdgeFunction->sum_w() * mEdgeFunction->abs_det_J() * constant::nu0 * mDeltaTime  )
                        * trans( tC ) * tC ;

            // compute the right hand side
            aRHS.fill( 0.0 );

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                aRHS += tW( k ) * mGroup->n( k ) * dot( mGroup->n( k ), tJz ) ;
            }
            aRHS *= mEdgeFunction->abs_det_J() * mDeltaTime ;
        }

//-----------------------------------------------------------------------------
    }
}