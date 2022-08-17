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
#include "fn_cross.hpp"

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
            mFields.ThinShell = { "lambda_m", "lambda_s", "lambda_z" };

            // non-dof fields
            mFields.MagneticFieldDensity     = { "bx", "by", "bz" };
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
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // containers and symbols
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // mass matrix
            Matrix< real > & tM = aJacobian ;

            // stiffness matrix
            Matrix< real > & tK = mGroup->work_K() ;

            // the sign of the edge, in this case, we use the physical tag
            const real tSign = aElement->element()->physical_tag() == 1 ? 1.0 : -1.0 ;

            BELFEM_ASSERT( tSign == ( aElement->edge_direction( 0 ) ? 1.0 : -1.0 ) ,
                           "Invalid Edge Direction for Element %lu",
                           ( long unsigned int ) aElement->id() );

            // container for geometry jacobian for master and slave
            Matrix< real > & tJ =  mGroup->work_J() ;

            // container for gradient operator
            Matrix< real > & tBm = mGroup->work_B() ;
            Matrix< real > & tBs = mGroup->work_D() ;

            // coordinates for master
            Matrix< real > & tXm = mGroup->work_Xi() ;

            // coordinates for slave
            Matrix< real > & tXs = mGroup->work_Eta() ;

            // expression cross( n, B )
            Vector< real > tPhi = mGroup->work_phi() ;

            // expression cross( n, B )
            Vector< real > tnxB = mGroup->work_psi() ;

            // mass matrix for individual layer
            Matrix< real > & tMlayer = mGroup->work_M() ;

            // stiffness matrix for individual layer
            Matrix< real > & tKlayer = mGroup->work_L() ;

            // data for magnetic field (H-data)
            Vector< real > & tHt = mGroup->work_sigma() ;

            // additional indices
            uint p ;
            uint q ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // reset matrices
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            tM.fill( 0.0 );
            tK.fill( 0.0 );
            aRHS.fill( 0.0 );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // geometry considerations and integration data
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // compute the normal
            const Vector< real > & tn = this->normal_straight_2d( aElement );


            // compute the element length
            const real tElementLength = mGroup->work_det_J() + mGroup->work_det_J();

            // collect coordinates from master
            this->collect_node_coords( aElement->master(), tXm );

            // collect coordinates from slave
            this->collect_node_coords( aElement->slave(), tXs );

            // get integration data from master
            const IntegrationData * tMaster =
                    mGroup->master_integration( aElement->facet()->master_index() );

            // get integration data from slave
            const IntegrationData * tSlave =
                    mGroup->slave_integration( aElement->facet()->slave_index() );

            // integration weights on master
            const Vector< real > & tWm = tMaster->weights() ;

            // integration weights on slave
            const Vector< real > & tWs = tSlave->weights() ;

            BELFEM_ASSERT( tWm.length() == mNumberOfIntegrationPoints, "number of integraiton points on master does not match" );
            BELFEM_ASSERT( tWs.length() == mNumberOfIntegrationPoints, "number of integraiton points on slave does not match" );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // gradient operator master
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // compute jacobian ( it is constant in first order, so we can set k=0
            tJ = tMaster->dNdXi( 0 ) * tXm ;

            // gradient operator
            tBm = inv( tJ ) * tMaster->dNdXi( 0 ) ;

            // get phi field
            this->collect_node_data( aElement->master(), "phi", tPhi );

            // delete me
            Vector< real > tH0 ( tBm * tPhi );
            tH0 *= -1 ;

            real tHn  = dot( tn, tH0 );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // lambda-condition for master element, air side
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // compute the B-Matrix ( is constant for linear elements )
            crossmat( tn, tBm, tnxB );

            // get the column
            p = 0 ;
            q = aRHS.length() - 3 ; // #lambda_z

            for( uint i=0; i<mNumberOfNodesPerElement; ++i )
            {
                tM( p, q ) = tnxB( i ) ;
                tM( q, p ) = tnxB( i ) ;
                ++p ;
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // gradient operator slave
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // compute jacobian ( it is constant in first order, so we can set k=0
            tJ = tSlave->dNdXi( 0 ) * tXs ;

            // gradient operator
            tBs = inv( tJ ) * tSlave->dNdXi( 0 ) ;

            // get phi field
            this->collect_node_data( aElement->slave(), "phi", tPhi );

            // real tHns = -dot( tn.vector_data(), tBs * tPhi );
            //Vector< real > tH1 = tBs * tPhi ;
            //tH1 *= -1 ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // lambda-condition for slave element, air side
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // next column
            p = mNumberOfNodesPerElement ;
            q = aRHS.length() - 2 ; // #lambda_z

            crossmat( tn, tBs, tnxB );

            // loop over all nodes
            for( uint i=0; i<mNumberOfNodesPerElement; ++i )
            {
                tM( p, q ) = tnxB( i ) ;
                tM( q, p ) = tnxB( i ) ;
                ++p ;
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // lambda-condition for master element, sheet side
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            q--;
            tM( p, q ) = tSign ;
            tM( q, p ) = tSign ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // lambda-condition for slave element, sheet side
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            p = q-1 ;
            q++;
            tM( p, q ) = tSign ;
            tM( q, p ) = tSign ;


            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // contribution for mass and stiffness matrix
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            tMlayer( 0, 0 ) = 2.0 ;
            tMlayer( 1, 0 ) = 1.0 ;
            tMlayer( 0, 1 ) = 1.0 ;
            tMlayer( 1, 1 ) = 2.0 ;
            tMlayer *= constant::mu0 / 6.0 ;

            // need thicknesses here
            p = 2 * mNumberOfNodesPerElement ;

            uint tNumLayers = mGroup->number_of_thin_shell_layers() ;

            // loop over all layers
            for( uint l=0; l<tNumLayers; ++l )
            {

                this->collect_edge_data_from_layer( aElement, "edge_h", l, tHt );

                // get thickness of thin shell
                real t = mGroup->thin_shell_thickness( l );

                tM(p,p)     += tMlayer(0,0) * t ;
                tM(p+1,p)   += tMlayer(1,0) * t ;
                tM(p,p+1)   += tMlayer(0,1) * t ;
                tM(p+1,p+1) += tMlayer(1,1) * t ;

                this->compute_layer_stiffness( aElement, l, tHt, tHn,tKlayer );
                tK(p,p)     += tKlayer(0,0) ;
                tK(p+1,p)   += tKlayer(1,0) ;
                tK(p,p+1)   += tKlayer(0,1) ;
                tK(p+1,p+1) += tKlayer(1,1) ;

                ++p ;

            }

            // scale mass matrix with element length
            // ( thus simplifying an integral over constant)
            tM *= tElementLength ;
            tK *= tElementLength ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // additional BC
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            p = aRHS.length() - 1 ;
            q = 0 ;

            // normal component master
            for( uint i=0; i<mNumberOfNodesPerElement; ++i )
            {
                tM( p, q )  =  (
                           tn( 0 ) * tBm( 0, i)
                          + tn( 1 ) * tBm( 1, i ) ) * tElementLength ;
                tM( q, p ) = tM( p, q ) ;
                ++q ;
            }

            // normal component slave
            for( uint i=0; i<mNumberOfNodesPerElement; ++i )
            {
                tM( p, q )  =  -(
                        tn( 0 ) * tBs( 0, i)
                                       + tn( 1 ) * tBs( 1, i ) ) * tElementLength ;
                tM( q, p ) = tM( p, q ) ;
                ++q ;
            }


            /*const IntegrationData * tIntMaster = mGroup->master_integration( aElement->facet()->master_index() );
            const IntegrationData * tIntSlave = mGroup->slave_integration( aElement->facet()->slave_index() );

            const Vector< real > & tW = mGroup->integration_weights();

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                {
                    q = mNumberOfNodesPerElement ;
                    for ( uint j = 0; j < mNumberOfNodesPerElement; ++j )
                    {
                        tM( i, q ) += tW( k ) * tIntMaster->N( k )( 0, i )
                                      * ( tn( 0 ) * tBs( 0, j )
                                          + tn( 1 ) * tBs( 1, j )) * 0.5 * tElementLength ;
                        ++q ;
                    }
                    aRHS( i ) += tW( k ) * 0.5 * tElementLength * tIntMaster->N( k )( 0, i ) * tHns ;

                }
                p = mNumberOfNodesPerElement ;
                for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                {
                    for ( uint j = 0; j < mNumberOfNodesPerElement; ++j )
                    {
                        tM( p, j ) -= tW( k ) * tIntSlave->N( k )( 0, i )
                                      * ( tn( 0 ) * tBm( 0, j )
                                          + tn( 1 ) * tBm( 1, j )) * 0.5 * tElementLength ;
                    }
                    ++p ;
                }
            } */


            // compute the right hand side
            aRHS += tM *  this->collect_q0_hphi_thinshell( aElement ) ;

            // finalize the Jacobian
            aJacobian += mDeltaTime * mGroup->work_K() ;




        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tri3::compute_layer_stiffness(
                Element * aElement,
                const uint aLayer,
                const Vector < real > & aHt,
                const real aHn,
                Matrix< real > & aK )
        {
            // reset matrix
            aK.fill( 0.0 );

            // grab material
            const Material * tMaterial = mGroup->thin_shell_material( aLayer );

            // grab layer thickness
            const real tThickness = mGroup->thin_shell_thickness( aLayer );

            // for the curl-Operator
            Matrix< real > & tC = mGroup->work_C();

            // integration data
            const IntegrationData * tInteg = mGroup->thinshell_integration();

            // get the ingeration weights
            const Vector< real > & tW = tInteg->weights() ;
            uint tNumIntpoints = tW.length() ;

            // geometry Jacobian
            real tJ = 0.5 * tThickness ;

            // magnetic field
            real tB ;

            // todo: compute temperature
            real tT = 4.0 ;

            // current density is constant for linear element
            real tJz = (  aHt( 1 ) - aHt( 0 ) )  / tThickness ;

            // curl operator is constant for linear element
            tC = tInteg->dNdXi( 0 ) ;
            tC /= tJ ;
            real tHt ;

            for( uint k=0; k<tNumIntpoints; ++k )
            {
                // get shape function

                const Matrix< real > & tN = tInteg->N( k );


                // magnetic field
                tHt = tN( 0, 0 ) * aHt( 0 ) + tN( 0, 1 ) * aHt( 1 );

                tB = constant::mu0 * sqrt( tHt * tHt + aHn * aHn );

                //real tRho = std::min( tMaterial->rho_el( tJz, tT, tB ) + 1e-12 , 1e-3 );
                real tRho = tMaterial->rho_el( tJz, tT, tB ) ;

                aK += ( tW( k ) * tRho ) * trans( tC ) * tC ;

                // std::cout << k << " " << tMaterial->label() << " B=" << tB << " J=" << tJz << " rho " << tRho << std::endl ;

            }
            aK *= tJ ;
        }

//------------------------------------------------------------------------------
    }
}