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


//#define HPHI_TRI3_LAMBDAN

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
#ifdef HPHI_TRI3_LAMBDAN
            mFields.ThinShell = { "lambda_m", "lambda_s", "lambda_n" };
#else
            mFields.ThinShell = { "lambda_m", "lambda_s" };
#endif
            mFields.SymmetryAir = { "lambda_n" };
            mFields.SymmetryFerro = { "lambda_n" };
            mFields.AntiSymmetryAir = { "lambda_n" };
            mFields.AntiSymmetryFerro = { "lambda_n" };

            //mFields.ThinShell = { "lambda_m", "lambda_s" };
            // non-dof fields
            mFields.MagneticFieldDensity     = { "bx", "by", "bz" };
            mFields.CurrentDensity = {  "jz" };
            mFields.CurrentBC = { "lambda_I" };
            mFields.Ghost = { "elementEJ", "elementJ", "elementJJc", "elementRho" };

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
                case( DomainType::Conductor ) :
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
                case( DomainType::Ferro ) :
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
                case( DomainType::InterfaceFmFm ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_fmfm ;
                    break ;
                }
                case( DomainType::InterfaceFmAir ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_fmair ;
                    break ;
                }
                case ( DomainType::SymmetryAir ) :
                {
                    mFunJacobian =
                            &IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_symmetry_air ;
                    break;
                }
                case ( DomainType::AntiSymmetryAir ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_antisymmetry_air ;
                    break;
                }
                case( DomainType::SymmetryFerro ) :
                {
                    mFunJacobian =
                            &IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_symmetry_ferro ;
                    break;
                }
                case( DomainType::AntiSymmetryFerro ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_antisymmetry_ferro ;
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
                                            & IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_wave ;
                                    break;
                                }
                                case ( MagfieldBcType::Symmetry ) :
                                {
                                    mFunJacobian =
                                            & IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_symmetry_air ;
                                    break;
                                }
                                case ( MagfieldBcType::AntiSymmetry ) :
                                {
                                    mFunJacobian =
                                            & IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_antisymmetry_air ;
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
            Matrix< real > & tB = mGroup->work_B() ;
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

            // compute the B-Matrix ( is constant for linear elements )
            tB =  inv( tNodeFunction->dNdXi( 0 ) * mGroup->work_Xs()  )
                    * tNodeFunction->dNdXi( 0 );

            // compute the cross product
            crossmat( tn, tB, tNxB );
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
            // make sure that the facet is correctly oriented
            BELFEM_ASSERT( static_cast< DomainType >( aElement->slave()->element()->physical_tag() ) == DomainType::Ferro,
                           "Master element of facet %lu must be of DomainType::Ferro", ( long unsigned int ) aElement->id() );

            BELFEM_ASSERT( static_cast< DomainType >( aElement->master()->element()->physical_tag() ) == DomainType::Air,
                           "Master element of facet %lu must be of DomainType::Air", ( long unsigned int ) aElement->id() );

            const IntegrationData * tMaster =  mGroup->master_integration(
                    aElement->facet()->master_index() ) ;

            // compute Jacobian and B-function for air
            this->collect_node_coords( aElement->master(), mGroup->work_Xm() );
            this->collect_node_coords( aElement->slave(), mGroup->work_Xs() );

            mGroup->work_B() = inv( tMaster->dNdXi( 0 ) * mGroup->work_Xm() )
                    * tMaster->dNdXi( 0 );

                  Matrix< real > & tM    = aJacobian ;
                  Matrix< real > & tK    = mGroup->work_K();
            const Vector< real > & tn    = this->normal_straight_2d( aElement );
            const Matrix< real > & tB    = mGroup->work_B() ;
                  Vector< real > & tnxB  = mGroup->work_sigma() ;

            // again, we use a mathematical trick since n and B are constant

            crossmat( tn, tB, tnxB );
            tnxB *= -mGroup->work_det_J()  ; // flip sign, because n should point out of ferro

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
        IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_symmetry_air(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // that this works is a coincidence that we shamelessly exploit
            this->compute_jacobian_and_rhs_antisymmetry_ferro( aElement,
                                                               aJacobian,
                                                               aRHS ) ;
            aJacobian *= constant::mu0 ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_antisymmetry_air(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // that this works is a coincidence that we shamelessly exploit
            this->compute_jacobian_and_rhs_symmetry_ferro( aElement, aJacobian, aRHS );

            aJacobian *= constant::mu0 ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_symmetry_ferro(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            this->compute_jacobian_and_rhs_symmmetry_a_tri3(
                    aElement,
                    aJacobian,
                    aRHS );

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tri3::compute_jacobian_and_rhs_antisymmetry_ferro(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            this->compute_jacobian_and_rhs_antisymmmetry_a_tri3(
                    aElement,
                    aJacobian, aRHS );
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
            //const real tSign = aElement->element()->physical_tag() == 1 ? 1.0 : -1.0 ;

            BELFEM_ASSERT( aElement->element()->physical_tag() == 1 ? 1.0 : -1.0 == ( aElement->edge_direction( 0 ) ? 1.0 : -1.0 ) ,
                           "Invalid Edge Direction for Element %lu",
                           ( long unsigned int ) aElement->id() );

            // container for gradient operator
            Matrix< real > & tBm = mGroup->work_B() ;
            Matrix< real > & tBs = mGroup->work_D() ;

            // coordinates for master
            Matrix< real > & tXm = mGroup->work_Xm() ;

            // coordinates for slave
            Matrix< real > & tXs = mGroup->work_Xs() ;

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

            // stiffness matrix
            Matrix< real > & tA = tK;

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

            BELFEM_ASSERT( tMaster->weights().length() == mNumberOfIntegrationPoints, "number of integraiton points on master does not match" );
            BELFEM_ASSERT( tSlave->weights().length() == mNumberOfIntegrationPoints, "number of integraiton points on slave does not match" );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // gradient operator master
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // gradient operators
            tBm = inv( tMaster->dNdXi( 0 ) * tXm ) * tMaster->dNdXi( 0 ) ;
            tBs = inv( tSlave->dNdXi( 0 ) * tXs ) * tSlave->dNdXi( 0 ) ;

            this->collect_node_data( aElement->slave(), "phi", tPhi );
            real tHn  = -0.5 * dot( tn.vector_data(), tBs * tPhi );

            // get phi field
            this->collect_node_data( aElement->master(), "phi", tPhi );

            // compute normal field (needed for material properties in stiffness matrix)
            tHn  -= 0.5 * dot( tn.vector_data(), tBm * tPhi );

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

            // loop over all layers
            for( uint l=0; l<mNumberOfThinShellLayers; ++l )
            {

                this->collect_edge_data_from_layer( aElement, "edge_h", l, tHt );

                //Vector< real > tH( 2 );
                //tH = tBm * tPhi ;

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


            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // lambda-condition for master element, air side
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



#ifdef HPHI_TRI3_LAMBDAN
            q = aRHS.length() - 3 ;
#else
            q = aRHS.length() - 2 ;
#endif
            // compute the B-Matrix ( is constant for linear elements )
            crossmat( tn, tBm, tnxB );


            // get the column
            p = 0 ;

            for( uint i=0; i<mNumberOfNodesPerElement; ++i )
            {
                tA( p, q ) = tnxB( i ) ;
                tA( q, p ) = tnxB( i ) ;
                ++p ;
            }


            // next column
            q += 1 ;

            crossmat( tn, tBs, tnxB );
            // loop over all nodes
            for( uint i=0; i<mNumberOfNodesPerElement; ++i )
            {
                tA( p, q ) = tnxB( i ) ;
                tA( q, p ) = tnxB( i ) ;
                ++p ;
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // lambda-condition for master element, sheet side
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            q-= 1 ;

            tA( p, q ) = 1.0 ;
            tA( q, p ) = 1.0 ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // lambda-condition for slave element, sheet side
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef HPHI_TRI3_LAMBDAN
            p = aRHS.length() - 4 ;
#else
            p = aRHS.length() - 3 ;
#endif
            q+= 1 ;
            tA( p, q ) = 1.0 ;
            tA( q, p ) = 1.0 ;

#ifdef HPHI_TRI3_LAMBDAN
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // additional BC
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            p = aRHS.length() - 1 ;
            q = 0 ;

            // normal component master
            for( uint i=0; i<mNumberOfNodesPerElement; ++i )
            {
                tA( p, q )  =  (
                            tn( 0 ) * tBm( 0, i )
                          + tn( 1 ) * tBm( 1, i ) );
                tA( q, p ) = tA( p, q ) ;
                ++q ;
            }

            // normal component slave
            for( uint i=0; i<mNumberOfNodesPerElement; ++i )
            {
                tA( p, q )  =  -(
                        tn( 0 ) * tBs( 0, i )
                      + tn( 1 ) * tBs( 1, i ) ) ;
                tA( q, p ) = tA( p, q ) ;
                ++q ;
            }

#endif
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
            }*/

            // interface matrix
            Matrix< real > & tL = mGroup->work_H() ;
            this->collect_edge_data_from_layer( aElement, "edge_h", 0, tHt );
            this->compute_thinshell_interface( true, aElement, 0, tHt, tHn,tL );
            p = 6 ;
            for( uint j=0; j<2; ++j )
            {
                for( uint i=0; i<3; ++i )
                {
                    tK( i, j+p) = tL( i, j );
                }
            }
            this->collect_edge_data_from_layer( aElement, "edge_h", mNumberOfThinShellLayers-1, tHt );
            this->compute_thinshell_interface( false, aElement, mNumberOfThinShellLayers-1, tHt, tHn,tL );
            p = 5 + mNumberOfThinShellLayers ;
            for( uint j=0; j<2; ++j )
            {
                for( uint i=0; i<3; ++i )
                {
                    tK( i+3, j+p ) = tL( i, j );
                }
            }

            // scale mass matrix with element length
            // ( thus simplifying an integral over constant)
            tM *= tElementLength ;
            tK *= tElementLength ;

            // compute the right hand side
            aRHS += tM *  this->collect_q0_thinshell( aElement ) ;

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
            real tH ;

            // todo: compute #temperature
            real tT = 10.0 ;

            // current density is constant for linear element
            real tJz = (  aHt( 1 ) - aHt( 0 ) )  / tThickness ;

            // curl operator is constant for linear element
            tC = tInteg->dNdXi( 0 ) ;
            tC /= tJ ;
            real tHt ;

            // value for ac losses
            real tEJ = 0.0 ;
            real tJJc = 0.0 ;
            for( uint k=0; k<tNumIntpoints; ++k )
            {
                // get shape function

                const Matrix< real > & tN = tInteg->N( k );

                // in-plane component
                tHt = tN( 0, 0 ) * aHt( 0 ) + tN( 0, 1 ) * aHt( 1 );

                // field angle and B-Field
                tH = std::sqrt( tHt*tHt + aHn*aHn ) ;
                real tAlpha = tH < BELFEM_EPSILON ? 0 : std::abs( std::asin( tHt / tH ) );

                //real tRho = std::min( tMaterial->rho_el( tJz, tT, tB ) + 1e-12 , 1e-3 );
                real tRho = std::max( tMaterial->rho_el( tJz, tT, constant::mu0*tH, tAlpha ), BELFEM_EPS );


                aK += ( tW( k ) * tRho ) * trans( tC ) * tC ;

                // std::cout << k << " " << tMaterial->label() << " B=" << tB << " J=" << tJz << " rho " << tRho << std::endl ;

                tEJ += tW( k ) * tRho ;
                tJJc += tW( k ) * tJz / tMaterial->j_crit( tT, constant::mu0*tH, tAlpha );
            }

            // scale the value of the stiffness matrix
            aK  *= tJ ;

            // compute element length and scale the value of EJ to the cross section
            tEJ *= tJz * tJz * 0.5 ;
            tJJc *= 0.5 ;

            // compute the field index
            index_t tIndex = mGhostElementMap( aElement->id()
                    *  mGroup->number_of_thin_shell_layers() + aLayer )->index() ;

            mMesh->field_data( "elementEJ")( tIndex ) = tEJ ;
            mMesh->field_data( "elementJ")( tIndex )  = tJz ;
            mMesh->field_data( "elementJJc")( tIndex )  = tJJc ;
        }


//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tri3::compute_thinshell_interface(
                const bool aMaster,
                Element * aElement,
                const uint aLayer,
                const Vector < real > & aHt,
                const real              aHn,
                Matrix< real > & aL )
        {
            // grab material
            const Material * tMaterial = mGroup->thin_shell_material( aLayer );

            // grab layer thickness
            const real tThickness = mGroup->thin_shell_thickness( aLayer );

            // for the curl-Operator
            Matrix< real > & tC = mGroup->work_C();

            // integration data
            const IntegrationData * tInteg = mGroup->thinshell_integration();

            // normal vector
            const Vector< real > & tN = mNormal2D ;

            // geometry Jacobian
            real tJ = 0.5 * tThickness ;

            // todo: compute #temperature
            real tT = 10.0 ;

            // current density is constant for linear element
            real tJz = (  aHt( 1 ) - aHt( 0 ) )  / tThickness ;

            // curl operator is constant for linear element
            tC = tInteg->dNdXi( 0 ) ;
            tC /= tJ ;

            // container for gradient operator, this has already been computed
            const Matrix< real > & tB = aMaster ? mGroup->work_B() : mGroup->work_D() ;

            // in-plane component
            const real tHt = aMaster ? aHt( 0 ) : aHt( 1 );

            // field angle and B-Field
            real tH = std::sqrt( tHt*tHt + aHn*aHn ) ;
            real tAlpha = tH < BELFEM_EPSILON ? 0 : std::abs( std::asin( tHt / tH ) );

            // resistivity
            real tRho = std::max( tMaterial->rho_el( tJz, tT, constant::mu0*tH, tAlpha ), BELFEM_EPS );

            for( uint j=0; j<2; ++j )
            {
                for( uint i=0; i<3; ++i )
                {
                    // note, the we multiply this with the element length later on
                    aL( i, j ) = tRho *
                            ( tN( 1 ) * tB( 0, i )
                            - tN( 0 ) * tB( 1, i ) ) * tC( 0, j ) ;
                }
            }

            if( aMaster )
            {
                aL *= -1. ;
            }
        }

//------------------------------------------------------------------------------

        // specifically for thin shells
        const Vector< real > &
        IWG_Maxwell_HPhi_Tri3::collect_q0_thinshell( Element * aElement )
        {
            // grab the output vector
            Vector< real > & aQ0 = mGroup->work_nedelec() ;

            BELFEM_ASSERT( mGroup->domain_type() == DomainType::ThinShell,
                           "function IWG_Maxwell_HPhi_Tri3::collect_q0_thinshell can only be applied to a thin shell" );

            // grab field data from mesh
            const Vector< real > & tPhi  = mMesh->field_data( "phi0" );
            const Vector< real > & tHe   = mMesh->field_data("edge_h0" );


            uint tLambdaCount = 0 ;

            // loop over all dofs
            for( uint k=0; k<mNumberOfDofsPerElement; ++k )
            {
                // grab dof
                Dof * tDof = aElement->dof( k );

                switch( tDof->mesh_basis()->entity_type() )
                {
                    case( EntityType::NODE ) :
                    {
                        aQ0( k ) = tPhi( tDof->mesh_basis()->index() );
                        break ;
                    }
                    case( EntityType::EDGE ) :
                    {
                        aQ0( k ) = tHe( tDof->mesh_basis()->index() );
                        break ;
                    }
                    case( EntityType::FACET ) :
                    {
                        // get lambda field
                        const Vector< real > & tL = mMesh->field_data( mFields.ThinShellLast( tLambdaCount++ ) );
                        aQ0( k ) = tL( tDof->mesh_basis()->index() );
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid dof type" );
                    }
                }
            }

            return aQ0 ;
        }


//------------------------------------------------------------------------------
    }
}