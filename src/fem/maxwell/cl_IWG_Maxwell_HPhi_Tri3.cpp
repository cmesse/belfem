//
// Created by Christian Messe on 12.01.22.
//

#include "cl_IWG_Maxwell_HPhi_Tri3.hpp"
#include "cl_FEM_MaxwellBoundaryConditionMagfield.hpp"
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

        IWG_Maxwell_HPhi_Tri3::IWG_Maxwell_HPhi_Tri3() :
            IWG_Maxwell( ElementType::TRI3,
                         IwgType::MAXWELL_HPHI_TRI3,
                         IwgMode::Iterative,
                         SymmetryMode::Unsymmetric,
                         SideSetDofLinkMode::MasterAndSlave,
                         true )
        {
            // dof fields
            mFields.Superconductor = { "edge_h" };
            mFields.Ferro = { "az" };
            mFields.Air = { "phi" };
            mFields.InterfaceScAir = { "lambda" };
            mFields.Cut = { "lambda" };
            mFields.ThinShell = { "lambda_m", "lambda_s" };

            // non-dof fields
            mFields.MagneticFieldDensity     = { "bx", "by" };
            mFields.CurrentDensity = {  "jz" };
            mFields.CurrentBC = { "lambda_I" };

            mEdgeDofMultiplicity = 1 ;
            mFaceDofMultiplicity = 0 ;
            mLambdaDofMultiplicity = 1 ;

            mNumberOfRhsDofsPerEdge = 1 ;
            mNumberOfRhsDofsPerFace = 0 ;

        }

//------------------------------------------------------------------------------

        IWG_Maxwell_HPhi_Tri3::~IWG_Maxwell_HPhi_Tri3()
        {

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tri3::link_jacobian_function( Group * aGroup )
        {
            // link edge function with group
            mEdgeFunction->link( aGroup );

            switch( aGroup->domain_type() )
            {
                case( DomainType::SuperConductor ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_superconductor ;
                    break ;
                }
                case( DomainType::Air ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_air ;
                    break ;
                }
                case( DomainType::FerroMagnetic ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_ferro ;
                    break ;
                }
                case( DomainType::InterfaceScAir ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_scair ;
                    break ;
                }
                case( DomainType::InterfaceScFm ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_scfm ;
                    break ;
                }
                case( DomainType::InterfaceFmAir ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_fmair ;
                    break ;
                }
                case( DomainType::Boundary ) :
                {
                    switch( mGroup->boundary_condition()->physics() )
                    {
                        case ( BoundaryConditionPhysics::Magfield ) :
                        {
                            // get subtype for magfield bc
                            switch ( reinterpret_cast< const MaxwellBoundaryConditionMagfield * >(
                                    mGroup->boundary_condition())->subtype())
                            {
                                case ( MagfieldBcType::Wave ) :
                                {
                                    mFunJacobian =
                                            &IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_wave;
                                    break;
                                }
                                default :
                                {
                                    BELFEM_ERROR( false, "Invalid boundary condition subtype" );
                                }
                            }
                            break;
                        }
                        default :
                        {
                            BELFEM_ERROR( false, "Invalid boundary condition type" );
                        }
                    }
                    break ;

                }
                case( DomainType::Cut ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_cut ;
                    break ;
                }
                case( DomainType::ThinShell ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_thinshell ;
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
         IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_scair(
                 Element        * aElement,
                 Matrix< real > & aJacobian,
                 Vector< real > & aRHS )
         {

            const IntegrationData * tNodeFunction = this->interface_data( aElement );

            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights() ;

            // crossed expressions
            Vector< real > & tNxB  = mGroup->work_sigma() ;
            Vector< real > & tNxE  = mGroup->work_tau() ;
            Matrix< real > & tIntE = mGroup->work_Tau() ;

            // compute the normal
            const Vector< real > & tn = this->normal_straight_2d( aElement );

            // reset matrix
            aJacobian.fill( 0.0 );

            // we can use a little mathematical hack here
            // by integrating E first
            tIntE.fill( 0.0 );

            Matrix< real > & tK = aJacobian ;
            tK.fill( 0.0 );

            // scaling parameter
            real tScale = mGroup->work_det_J() * constant::mu0 * this->timestep() ;

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                const Matrix< real > & tE = mEdgeFunction->E( k );
                const Vector< real > & tN = tNodeFunction->phi( k );

                tIntE += tW( k ) * mEdgeFunction->E( k );

                // perpendicular coupling
                for( uint i=0; i<3; ++i )
                {
                    for( uint j=0; j<3; ++j )
                    {
                        tK( i+3, j ) -= tW( k ) * tN( i ) *
                                ( tn( 0 ) * tE( 0, j )
                                + tn( 1 ) * tE( 1, j ) );
                    }
                }
            }

            // scale perpendicular coupling
            for( uint i=0; i<3; ++i )
            {
                for( uint j=0; j<3; ++j )
                {
                    tK( i+3, j ) *= tScale ;
                }
            }

            // compute Jacobian
            mGroup->work_J() = tNodeFunction->dNdXi( 0 ) * mGroup->node_coords() ;

            // compute the B-Matrix ( is constant for linear elements )
            mGroup->work_B() = inv( mGroup->work_J() ) * tNodeFunction->dNdXi( 0 );

            // compute the cross product
            crossmat( tn, mGroup->work_B(), tNxB );
            crossmat( tn, tIntE, tNxE );

            // now we scale the parameters
            tNxB *= tScale * mEdgeFunction->sum_w() ;
            tNxE *= tScale ;

            // and feed them into the matrices
            aRHS.fill( 0.0 );

            aJacobian( 6, 0 ) = tNxE( 0 );
            aJacobian( 6, 1 ) = tNxE( 1 );
            aJacobian( 6, 2 ) = tNxE( 2 );
            aJacobian( 6, 3 ) = tNxB( 0 );
            aJacobian( 6, 4 ) = tNxB( 1 );
            aJacobian( 6, 5 ) = tNxB( 2 );
            aJacobian( 0, 6 ) = tNxE( 0 );
            aJacobian( 1, 6 ) = tNxE( 1 );
            aJacobian( 2, 6 ) = tNxE( 2 );
            aJacobian( 3, 6 ) = tNxB( 0 );
            aJacobian( 4, 6 ) = tNxB( 1 );
            aJacobian( 5, 6 ) = tNxB( 2 );
         }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_fmair(
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
            const Vector< real > & tn    = this->normal_straight_2d( aElement );
            const Matrix< real > & tB    = mGroup->work_B() ;
                  Vector< real > & tnxB  = mGroup->work_sigma() ;

            // again, we use a mathematical trick since n and B are constant

            crossmat( tn, tB, tnxB );
            tnxB *= mGroup->work_det_J() * constant::mu0 ; ;


            tM.fill( 0.0 );
            switch( aElement->facet()->master_index() )
            {
                case( 0 ) :
                {
                    tM( 3, 0 ) = tnxB( 0 );
                    tM( 4, 0 ) = tnxB( 1 );
                    tM( 5, 0 ) = tnxB( 2 );
                    tM( 3, 1 ) = tnxB( 0 );
                    tM( 4, 1 ) = tnxB( 1 );
                    tM( 5, 1 ) = tnxB( 2 );
                    break ;
                }
                case( 1 ) :
                {
                    tM( 3, 1 ) = tnxB( 0 );
                    tM( 4, 1 ) = tnxB( 1 );
                    tM( 5, 1 ) = tnxB( 2 );
                    tM( 3, 2 ) = tnxB( 0 );
                    tM( 4, 2 ) = tnxB( 1 );
                    tM( 5, 2 ) = tnxB( 2 );
                    break ;
                }
                case( 2 ) :
                {
                    tM( 3, 0 ) = tnxB( 0 );
                    tM( 4, 0 ) = tnxB( 1 );
                    tM( 5, 0 ) = tnxB( 2 );
                    tM( 3, 2 ) = tnxB( 0 );
                    tM( 4, 2 ) = tnxB( 1 );
                    tM( 5, 2 ) = tnxB( 2 );
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid Facet index");
                };
            }
            tK = -trans( tM );
            aRHS = tM * this->collect_q0_aphi_2d( aElement );
            aJacobian += tK * this->timestep() ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_wave(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {

            // get length
            real tHx = aElement->element()->node( 1 )->x() - aElement->element()->node( 0 )->x();
            real tHy = aElement->element()->node( 1 )->y() - aElement->element()->node( 0 )->y();
            real tL = std::sqrt( tHx * tHx + tHy * tHy );

            // the right hand side of the matrix
            aJacobian.fill( 0.0 );

            Vector< real > & tPhi = mGroup->work_phi() ;

            switch ( aElement->facet()->master_index() )
            {
                case ( 0 ) :
                {
                    aJacobian( 0, 0 ) = 2.0;
                    aJacobian( 1, 0 ) = 1.0;
                    aJacobian( 0, 1 ) = 1.0;
                    aJacobian( 1, 1 ) = 2.0;
                    break;
                }
                case ( 1 ) :
                {
                    aJacobian( 1, 1 ) = 2.0;
                    aJacobian( 2, 1 ) = 1.0;
                    aJacobian( 1, 2 ) = 1.0;
                    aJacobian( 2, 2 ) = 2.0;
                    break;
                }
                case ( 2 ) :
                {
                    aJacobian( 0, 0 ) = 2.0;
                    aJacobian( 2, 0 ) = 1.0;
                    aJacobian( 0, 2 ) = 1.0;
                    aJacobian( 2, 2 ) = 2.0;
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Invalid side index" );
                }
            }
            aJacobian *= tL * mGroup->boundary_condition()->penalty() / 6.;
            //aJacobian *= constant::mu0;

            // get the field data (actually negative h )
            tHx = mGroup->boundary_condition()->data( 0 );
            tHy = mGroup->boundary_condition()->data( 1 );

            // compute potential for the individual nodes

            tPhi( 0 ) = tHx * aElement->facet()->master()->node( 0 )->x()
                            + tHy * aElement->facet()->master()->node( 0 )->y();
            tPhi( 1 ) = tHx * aElement->facet()->master()->node( 1 )->x()
                            + tHy * aElement->facet()->master()->node( 1 )->y();
            tPhi( 2 ) = tHx * aElement->facet()->master()->node( 2 )->x()
                            + tHy * aElement->facet()->master()->node( 2  )->y();
            aRHS = aJacobian * tPhi ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_thinshell(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            const IntegrationData * tNodeFunction = this->interface_data( aElement );

            // mass matrix
            Matrix< real > & tM = aJacobian ;

            // stiffness matrix
            Matrix< real > & tK = mGroup->work_K() ;

            // geometry jacobian for master and slave
            Matrix< real > & tJ =  mGroup->work_J() ;

            // node coordinates for master and slave
            Matrix< real > & tX = mGroup->work_X() ;

            // gradient operator for master and slave
            Matrix< real > &tB = mGroup->work_B() ;

            // container for expression cross( n, B )
            Vector< real > & tnxB = mGroup->work_sigma() ;

            // sub mass matrix
            Matrix< real > & tMlayer = mGroup->work_M() ;

            // sub stiffness matrix
            Matrix< real > & tKlayer = mGroup->work_M() ;

            // the sign of the edge
            const real tSign = aElement->edge_direction( 0 ) ? 1.0 : -1.0 ;

            // reset values
            tM.fill( 0.0 );
            tK.fill( 0.0 );
            aRHS.fill( 0.0 );

            this->print_dofs( aElement );
            // std::cout << "check " << aRHS.length() << " " << tM.n_cols() << " " << tK.n_cols() << std::endl ;


            // compute the normal
            const Vector< real > & tn = this->normal_straight_2d( aElement );
            const real tElementLength = mGroup->work_det_J() + mGroup->work_det_J();

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // lambda-condition for master element, air side
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // compute the B-Matrix ( is constant for linear elements )
            aElement->master()->get_node_coors( tX );
            tJ = tNodeFunction->dNdXi( 0 ) * tX ;
            tB = inv( tJ ) * tNodeFunction->dNdXi( 0 );


            crossmat( tn, mGroup->work_B(), tnxB );

            // number of nodes
            const uint n = tX.n_rows() ;

            // get the column
            uint i = 0 ;
            uint j = aRHS.length() - 2 ;

            for( uint k=0; k<n; ++k)
            {
                tM( i, j ) = tnxB( k ) ;
                tM( j, i ) = tnxB( k ) ;
                ++i ;
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // lambda-condition for slave element, air side
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // next column
            j++;
            aElement->slave()->get_node_coors( tX );
            tJ = tNodeFunction->dNdXi( 0 ) * tX ;
            tB = inv( tJ ) * tNodeFunction->dNdXi( 0 );
            for( uint k=0; k<n; ++k)
            {
                tM( i, j ) = tnxB( k ) ;
                tM( j, i ) = tnxB( k ) ;
                ++i ;
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // lambda-condition for master element, sheet side
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            j--;
            tM( i, j ) = tSign ;
            tM( j, i ) = tSign ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // lambda-condition for slave element, sheet side
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            i = j-1 ;
            j++;
            tM( i, j ) = tSign ;
            tM( j, i ) = tSign ;

            // scale mass matrix with element length
            // ( thus simplifying an integral over constant)
            tM *= tElementLength ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // contribution for mass matrix
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            tMlayer( 0, 0 ) = 2.0 ;
            tMlayer( 1, 0 ) = 1.0 ;
            tMlayer( 0, 1 ) = 1.0 ;
            tMlayer( 1, 1 ) = 2.0 ;
            tMlayer *= tElementLength * constant::nu0 / 6.0 ;

            // need thicknesses here
            i = 2 * mNumberOfNodesPerElement ;

            uint tNumLayers = mGroup->number_of_thin_shell_layers() ;

            // loop over all layers
            for( uint l=0; l<tNumLayers; ++l )
            {

                // get thickness of thin shell
                real t = mGroup->thin_shell_thickness( l );

                tM(i,i) += tMlayer(0,0) * t ;
                tM(i+1,i) += tMlayer(1,0) * t ;
                tM(i,i+1) += tMlayer(0,1) * t ;
                tM(i+1,i+1) += tMlayer(1,1) * t ;

                ++i ;

            }

            // not doing  the laplace flux thing for now (Alves doesn't do it either)

            // now we need the edges
            //tM.print("M");
            //std::cout << " l " << tElementLength << " " << mNumberOfLayersPerShell << std::endl ;


            this->compute_layer_stiffness( aElement, 2, tKlayer );

            //mGroup->master_integration( )
            exit( 0 );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tri3::compute_layer_stiffness(
                Element * aElement,
                const uint aLayer,
                Matrix< real > & aK )
        {
            // reset matrix
            aK.fill( 0.0 );

            // grab ghost element
            mesh::Facet * tFacet = mMesh->ghost_facet( aElement->id(), aLayer );

            std::cout << "edge " << tFacet->edge( 0 )->id() << std::endl ;

            std::cout << "nodes" << tFacet->edge( 0 )->node( 0 ) ->id() << " " <<tFacet->edge( 0 )->node( 1 )->id() << std::endl ;
            aK.print( "K_layer" );
        }

//------------------------------------------------------------------------------
    }
}