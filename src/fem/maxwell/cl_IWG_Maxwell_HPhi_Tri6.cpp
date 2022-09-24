//
// Created by christian on 12/22/21.
//

#include "cl_IWG_Maxwell_HPhi_Tri6.hpp"
#include "cl_FEM_MaxwellBoundaryConditionMagfield.hpp"
#include "fn_trans.hpp"
#include "fn_crossmat.hpp"
#include "fn_inv.hpp"
#include "fn_trans.hpp"
#include "cl_EF_LINE3.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_DofMgr_SolverData.hpp"
#include "fn_gesv.hpp"

#define HPHI_TRI6_LAMBDAN

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Maxwell_HPhi_Tri6::IWG_Maxwell_HPhi_Tri6() :
                IWG_Maxwell( ElementType::TRI6,
                             IwgType::MAXWELL_HPHI_TRI6,
                             IwgMode::Iterative,
                             SymmetryMode::Unsymmetric,
                             SideSetDofLinkMode::MasterAndSlave,
                             true )
        {
            // dof fields
            mFields.Superconductor = { "edge_h", "face_h" };
            mFields.Air = { "phi" };
            mFields.Ferro = { "az" };

            mFields.InterfaceScAir = { "lambda" };
            mFields.Cut = { "lambda" };

#ifdef HPHI_TRI6_LAMBDAN
            mFields.ThinShell = { "lambda_m0", "lambda_m1", "lambda_s0", "lambda_s1",  "lambda_tsn0", "lambda_tsn1" };
#else
            mFields.ThinShell = { "lambda_m0", "lambda_m1", "lambda_s0", "lambda_s1" };
#endif

            //mFields.ThinShell = { "lambda_m0", "lambda_m1", "lambda_s0", "lambda_s1", "lambda_n0", "lambda_n1" };
            //mFields.ThinShell = { "lambda_m0", "lambda_m1", "lambda_s0", "lambda_s1" };
            mFields.Farfield = { "lambda_n" };

            // non-dof fields
            mFields.MagneticFieldDensity     =  { "bx", "by", "bz" };
            mFields.CurrentDensity = {  "jz" };
            mFields.CurrentBC = { "lambda_I" };
            mFields.Ghost = { "elementEJ", "elementJ" };

            mEdgeDofMultiplicity = 2 ;
            mFaceDofMultiplicity = 2 ;
            mLambdaDofMultiplicity = 1 ;

            mNumberOfRhsDofsPerEdge = 2 ;
            mNumberOfRhsDofsPerFace = 2 ;

            // create the TS function
            mEdgeFunctionTS = new EF_LINE3() ;
            mHt.set_size( 3 );
            mC.set_size( 1, 6 );
            mH.set_size( 2 );
            mU.set_size( 2 );
            mE.set_size( 2 );

            mHn.set_size( 16 );
            mW.set_size( 16 );

            mAHn.set_size( 2, 2 );
            mBHn.set_size( 2 );

            mArho.set_size( 9, 9 );
            mBrho.set_size( 9 );
            mCrho.set_size( 9 );
            mPivot.set_size( 9 );
            mL.set_size( 2, 6 );
        }

//------------------------------------------------------------------------------

        IWG_Maxwell_HPhi_Tri6::~IWG_Maxwell_HPhi_Tri6()
        {

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tri6::link_jacobian_function( Group * aGroup )
        {
            // link edge function with group
            mEdgeFunction->link( aGroup );
            mEdgeFunctionTS->link( aGroup );

            switch( aGroup->domain_type() )
            {
                case( DomainType::SuperConductor ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_superconductor ;
                    break ;
                }
                case( DomainType::FerroMagnetic ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_ferro ;
                    break ;
                }
                case( DomainType::Air ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_air ;
                    break ;
                }
                case( DomainType::InterfaceScAir ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_scair ;
                    break ;
                }
                case( DomainType::InterfaceScFm ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_scfm ;
                    break ;
                }
                case( DomainType::InterfaceFmAir ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_fmair ;
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
                                            &IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_wave;
                                    break;
                                }
                                case ( MagfieldBcType::Farfied ) :
                                {
                                    mFunJacobian =
                                            &IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_farfield ;
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
                            & IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_cut ;
                    break ;
                }
                case( DomainType::ThinShell ) :
                {
                    mLayerData.set_size( mFields.Ghost.size() + 1 , mGroup->number_of_thin_shell_layers() );

                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_thinshell ;
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
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_scair(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // link edge function with element
            const IntegrationData * tNodeFunction = this->interface_data( aElement );

            // crossed expressions
            Matrix< real > & tB    = mGroup->work_B() ;

            // get weight functions
            const Vector< real > & tW = mGroup->integration_weights() ;

            Vector< real > & tNxE  = mGroup->work_tau() ;
            Vector< real > & tNxB  = mGroup->work_sigma() ;


            // now we can build the stiffness matrix
            Matrix< real > & tK = aJacobian ;

            // reset K-matrix
            tK.fill( 0.0 );

            // scaling parameter
            real tScale = constant::mu0 * this->timestep() ;

            if( aElement->master()->element()->is_curved() )
            {
                // reset components for parallel coupling
                tNxE.fill( 0.0 );
                tNxB.fill( 0.0 );

                // help parameter
                real tOmega ;

                // loop over all integration points
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // Jacobian for air element
                    mGroup->work_J() = tNodeFunction->dNdXi( k ) * mGroup->work_Xs() ;

                    // gradient operator for air element
                    tB = inv( mGroup->work_J() ) * tNodeFunction->dNdXi( k ) ;

                    // node function at integration point for air element
                    const Vector< real > & tN = tNodeFunction->phi( k );

                    // edge function at integration point
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // normal at integration point ( also computes mGroup->work_det_J() )
                    const Vector< real > & tn = this->normal_curved_2d( aElement, k );

                    // tOmega = tW( k ) * mGroup->work_det_J() ;
                    tOmega = tW( k ) * mGroup->work_det_J() ;

                    // add components to for parallel coupling
                    crossmat( tn, tE, tOmega, tNxE );
                    crossmat( tn, tB, tOmega, tNxB );

                    // add components for perpendicular coupling
                    for( uint j=0; j<8; ++j )
                    {
                        for( uint i=0; i<6; ++i )
                        {
                            tK( i+8, j ) -= tOmega * tN( i ) *
                                    ( tn( 0 ) * tE( 0, j )
                                    + tn( 1 ) * tE( 1, j ) );
                        }
                    }
                }
            } // end element with curved edges
            else
            {
                // help matrices for integration
                Matrix< real > & tIntNxi  = mGroup->work_Sigma() ;
                Matrix< real > & tIntE    = mGroup->work_Tau() ;
                tIntNxi.fill( 0.0 );
                tIntE.fill( 0.0 );

                // get normal vector ( is constant for linear element )
                const Vector< real > & tn = this->normal_straight_2d( aElement );

                // add length of edge to scaling parameter
                tScale *= mGroup->work_det_J() ;

                // integrate components for N_xi and E
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // edge function at integration point
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // node function at integration point
                    const Vector< real > & tN = tNodeFunction->phi( k );

                    // perpendicular coupling
                    for( uint j=0; j<8; ++j )
                    {
                        for( uint i=0; i<6; ++i )
                        {
                            tK( i+8, j ) -= tW( k ) * tN( i ) *
                                    ( tn( 0 ) * tE( 0, j )
                                    + tn( 1 ) * tE( 1, j ) );
                        }
                    }

                    tIntE   += tW( k ) * tE ;
                    tIntNxi += tW( k ) * tNodeFunction->dNdXi( k );
                }

                // compute gradient operator ( assuming constant geometry Jacobian  )
                tB = inv( tNodeFunction->dNdXi( 0 ) * mGroup->work_Xs() ) * tIntNxi ;

                // compute the cross products
                crossmat( tn, tB, tNxB );
                crossmat( tn, tIntE, tNxE );

            } // end element with straight edges

            // scale parallel coupling
            tNxB *= tScale ;
            tNxE *= tScale ;

            // scale perpendicular coupling

            for( uint j=0; j<8; ++j )
            {
                for( uint i=0; i<6; ++i )
                {
                    tK( i+8, j ) *= tScale ;
                }
            }

            // assemble parallel components
            tK(  0, 14 ) = tNxE( 0 );
            tK(  1, 14 ) = tNxE( 1 );
            tK(  2, 14 ) = tNxE( 2 );
            tK(  3, 14 ) = tNxE( 3 );
            tK(  4, 14 ) = tNxE( 4 );
            tK(  5, 14 ) = tNxE( 5 );
            tK(  6, 14 ) = tNxE( 6 );
            tK(  7, 14 ) = tNxE( 7 );
            tK(  8, 14 ) = tNxB( 0 );
            tK(  9, 14 ) = tNxB( 1 );
            tK( 10, 14 ) = tNxB( 2 );
            tK( 11, 14 ) = tNxB( 3 );
            tK( 12, 14 ) = tNxB( 4 );
            tK( 13, 14 ) = tNxB( 5 );
            tK.set_row( 14, tK.col( 14 ) );

            // no RHS, since M=0 and backward implicit
            aRHS.fill( 0.0 );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_fmair(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            const IntegrationData * tMaster =  mGroup->master_integration(
                    aElement->facet()->master_index() ) ;

            const IntegrationData * tSlave =  mGroup->slave_integration(
                    aElement->facet()->slave_index() ) ;

            // get node coords, needed to compute geometry Jacobian
            aElement->slave()->get_node_coors( mGroup->work_Xs() );
            // the B-Matrix
            Matrix< real > & tK = mGroup->work_K() ;
            Matrix< real > & tM = aJacobian ;
            //Matrix< real > & tJ = mGroup->work_J() ;
            Matrix< real > & tB = mGroup->work_B() ;
            Vector< real > & tnxB = mGroup->work_sigma() ;
            const Vector< real > & tW = mGroup->integration_weights() ;

            tK.fill( 0.0 );

            if( aElement->element()->is_curved() )
            {
                // grab node coords for normal vector
                aElement->get_node_coors( mGroup->work_X() );

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    const Vector< real > & tn = this->normal_curved_2d( aElement, k );

                    real tOmega = -tW( k ) * mGroup->work_det_J() * constant::mu0 ;
                    const Vector< real > & tN = tMaster->phi( k );
                    tB = inv( tSlave->dNdXi( k ) * mGroup->work_Xs() ) * tSlave->dNdXi( k );
                    crossmat( tn, tB, tnxB );

                    for( uint j=0; j<6; ++j )
                    {
                        for( uint i=0; i<6; ++i )
                        {
                            tK( i, j+6 ) += tOmega * tN( i ) * tnxB( j );
                        }
                    }
                }
            }
            else
            {
                // get the normal
                const Vector< real > & tn = this->normal_straight_2d( aElement );

                // compute the Jacobian
                Matrix< real > & tInvJ = mGroup->work_invJ();
                tInvJ = inv( tSlave->dNdXi( 0 ) * mGroup->work_Xs() );

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    const Vector< real > & tN = tMaster->phi( k );
                    tB = tInvJ * tSlave->dNdXi( k );
                    crossmat( tn, tB, tnxB );
                    real tOmega = -tW( k ) * mGroup->work_det_J() * constant::mu0 ;

                    for( uint j=0; j<6; ++j )
                    {
                        for( uint i=0; i<6; ++i )
                        {
                            tK( i, j+6 ) += tOmega * tN( i ) * tnxB( j );
                        }
                    }
                }
            }
            tM = - trans( tK );
            aRHS = tM * this->collect_q0_aphi_2d( aElement );
            aJacobian += tK * this->timestep() ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_wave(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {

            // reset the matrix
            aJacobian.fill( 0.0 );

            if ( aElement->element()->is_curved())
            {
                // integration data
                const IntegrationData * tMaster = mGroup->master_integration(
                        aElement->facet()->master_index());


                // grab node coords for normal vector
                aElement->get_node_coors( mGroup->work_X() );

                // integration weights
                const Vector< real > & tW = mGroup->integration_weights();

                for ( uint k = 0; k < mNumberOfIntegrationPoints; ++k )
                {
                    // vector for mass matrix
                    const Matrix< real > & tN = tMaster->N( k );

                    // we don't need the normal. But calling this function
                    // writes the integration element into mGroup->work_det_J()
                    this->normal_curved_2d( aElement, k );

                    aJacobian += ( tW( k ) * mGroup->work_det_J()) * trans( tN ) * tN;
                }
            }
            else
            {
                // get element length
                real tA = aElement->element()->node( 1 )->x() - aElement->element()->node( 0 )->x();
                real tB = aElement->element()->node( 1 )->y() - aElement->element()->node( 0 )->y();

                // element length
                real tL = std::sqrt( tA * tA + tB * tB );

                // helpers
                tA = tL / 15.;

                switch ( aElement->facet()->master_index() )
                {
                    case ( 0 ):
                    {
                        aJacobian( 0, 0 ) = tA + tA + tA + tA;
                        aJacobian( 1, 0 ) = -tA;
                        aJacobian( 3, 0 ) = tA + tA;

                        aJacobian( 0, 1 ) = -tA;
                        aJacobian( 1, 1 ) = tA + tA + tA + tA;
                        aJacobian( 3, 1 ) = tA + tA;

                        aJacobian( 0, 3 ) = tA + tA;
                        aJacobian( 1, 3 ) = tA + tA;
                        aJacobian( 3, 3 ) = 16. * tA;

                        break;

                    }
                    case ( 1 ) :
                    {
                        aJacobian( 1, 1 ) = tA + tA + tA + tA;
                        aJacobian( 2, 1 ) = -tA;
                        aJacobian( 4, 1 ) = tA + tA;

                        aJacobian( 1, 2 ) = -tA;
                        aJacobian( 2, 2 ) = tA + tA + tA + tA;
                        aJacobian( 4, 2 ) = tA + tA;

                        aJacobian( 1, 4 ) = tA + tA;
                        aJacobian( 2, 4 ) = tA + tA;
                        aJacobian( 4, 4 ) = 16. * tA;
                        break;

                    }
                    case ( 2 ) :
                    {
                        aJacobian( 0, 0 ) = tA + tA + tA + tA;
                        aJacobian( 2, 0 ) = -tA;
                        aJacobian( 5, 0 ) = tA + tA;

                        aJacobian( 0, 2 ) = -tA;
                        aJacobian( 2, 2 ) = tA + tA + tA + tA;
                        aJacobian( 5, 2 ) = tA + tA;

                        aJacobian( 0, 5 ) = tA + tA;
                        aJacobian( 2, 5 ) = tA + tA;
                        aJacobian( 5, 5 ) = 16. * tA;

                        break;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid side " );
                    }
                }
            }

            // apply penalty
            aJacobian *= mGroup->boundary_condition()->penalty() ;

            // get the field data ( actually negative h )
            real tHx = mGroup->boundary_condition()->data( 0 );
            real tHy = mGroup->boundary_condition()->data( 1 );


            Matrix< real > & tX = mGroup->node_coords();
            aElement->master()->get_node_coors( tX );

            aRHS = aJacobian * ( tHx * tX.col( 0 )
                               + tHy * tX.col( 1 ) );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_farfield(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {

            // the right hand side of the matrix
            aJacobian.fill( 0.0 );
            aRHS.fill( 0.0 );

            // get the node coordinates
            Matrix< real > & tXm = mGroup->work_Xm() ;
            this->collect_node_coords( aElement, mGroup->work_X() );
            this->collect_node_coords( aElement->master(), tXm );

            // get integration data from master
            const IntegrationData * tMaster =
                    mGroup->master_integration( aElement->facet()->master_index() );

            // get the integration weights
            const Vector< real > tW = tMaster->weights() ;

            // the gradient operator matrix
            Matrix< real > & tB = mGroup->work_B() ;

            // Edge Function
            mEdgeFunctionTS->link( aElement );

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // compute the gradient operator
                tB.matrix_data() = inv( tMaster->dNdXi( k ) * tXm ) * tMaster->dNdXi( k ) ;

                // compute the normal
                const Vector< real > & tn = this->normal_curved_2d( aElement, k );

                // loop over all nodes

                real tScale = tW( k ) * mGroup->work_det_J() * constant::mu0 ;

                for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                {
                    real tVal = tScale * (
                              tn( 0 ) * tB( 0, i )
                            + tn( 1 ) * tB( 1, i ) ) ;
                    aJacobian( mNumberOfNodesPerElement, i ) += tVal ;
                    aJacobian( i, mNumberOfNodesPerElement ) += tVal ;
                }
            }
        }


//------------------------------------------------------------------------------

#ifdef BELFEM_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

        void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_thinshell(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            //std::cout << "====== " << aElement->facet()->master_index() << " " << aElement->facet()->slave_index() << std::endl ;
            //aElement->element()->print() ;
            //aElement->master()->element()->print() ;
            //aElement->slave()->element()->print() ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // main containers
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // mass matrix
            Matrix< real > & tM = aJacobian;

            // stiffness matrix
            Matrix< real > & tK = mGroup->work_K();

            tM.fill( 0.0 );
            tK.fill( 0.0 );
            aRHS.fill( 0.0 );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Element Coordinates
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // element coordinates
            Matrix< real > & tX = mGroup->work_X();
            this->collect_node_coords( aElement, tX );

            // coordinates for master
            Matrix< real > & tXm = mGroup->work_Xm();
            this->collect_node_coords( aElement->master(), tXm );


            // coordinates for slave
            Matrix< real > & tXs = mGroup->work_Xs();
            this->collect_node_coords( aElement->slave(), tXs );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // functions and integration data
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // Edge Function
            mEdgeFunctionTS->link( aElement );

            // integration weights along edge
            const Vector< real > & tW = mGroup->integration_weights();

            // get integration data from master
            const IntegrationData * tIntMaster =
                    mGroup->master_integration( aElement->facet()->master_index() );

            // get integration data from slave
            const IntegrationData * tIntSlave =
                    mGroup->slave_integration( aElement->facet()->slave_index() );


            // field on master side
            Vector< real > & tPhiM = mGroup->work_phi();
            this->collect_node_data( aElement->master(), "phi", tPhiM );

            // todo: delete me
            Vector< real > tPhiS( 6 );
            this->collect_node_data( aElement->slave(), "phi", tPhiS );

            // data for one layer ( in-plane magnetic field )
            Vector< real > & tHt = mGroup->work_sigma();

            // container for gradient operator, master side
            Matrix< real > & tBm = mGroup->work_B();

            // container for gradient operator, slave side
            Matrix< real > & tBs = mGroup->work_D();

            // expression cross( n, B )
            Vector< real > tnxB = mGroup->work_psi();

            // mass matrix for individual layer
            Matrix< real > & tMlayer = mGroup->work_M();

            // stiffness matrix for individual layer
            Matrix< real > & tKlayer = mGroup->work_L();

            //Matrix< real > & tJ = mGroup->work_J();
            //Matrix< real > & tG = mGroup->work_G();
            //Matrix< real > & tH = mGroup->work_H();

            // reset the thin shell data
            mLayerData.fill( 0.0 );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // some sanity checks
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            real tSign = aElement->element()->physical_tag() == 1 ? 1.0 : -1.0 ;

            // make sure that sign is correct
            BELFEM_ASSERT( aElement->element()->physical_tag() == 1 ? 1.0 : -1.0
                        == ( aElement->edge_direction( 0 ) ? 1.0 : -1.0 ) ,
                           "Invalid Edge Direction for Element %lu",
                           ( long unsigned int ) aElement->id() );

            // make sure that integration weights make sense
            BELFEM_ASSERT( tIntMaster->weights().length() == mNumberOfIntegrationPoints,
                           "number of integraiton points on master does not match" );
            BELFEM_ASSERT( tIntSlave->weights().length() == mNumberOfIntegrationPoints,
                           "number of integraiton points on slave does not match" );

            BELFEM_ASSERT( mGroup->thinshell_integration()->weights().length() == mNumberOfIntegrationPoints,
                           "number of integraiton points on thin shell does not match" );

            if( aElement->id() == 17 )
            {
                std::cout << aElement->facet()->master_index() << " " << aElement->facet()->slave_index() << std::endl ;
                tXm.print("Xm");
                tXs.print("Xs");
                tPhiM.print("phiM");
                tPhiS.print("phiS");

            }
            // indices
            uint p ;
            uint q ;

            real tValue ;

            // element lengfth
            real tLength = 0 ;
            //std::cout << std::endl << "------- " << std::endl << "el " << aElement->id() << " " << tSign << " | " << aElement->master()->id() << " " << aElement->slave()->id() << " |  " << aElement->element()->node( 0 )->id() << " " << aElement->element()->node( 1 )->id() << std::endl ;

            // loop over all integration points
            for ( uint k = 0; k < mNumberOfIntegrationPoints; ++k )
            {
                // compute element normal, also writes integration increment
                // into mGroup->work_det_J()
                const Vector< real > & tn = this->normal_curved_2d( aElement, k );
                //const Vector< real > & tn = this->normal_straight_2d( aElement );

                // integration increment for this point
                real tWDetJ = tW( k ) * mGroup->work_det_J() ;
                mW( k ) = tWDetJ ;
                tLength += tWDetJ ;

                // evaluate points for edge function
                mE( 0 ) = mEdgeFunctionTS->E( k )( 0, 0 );
                mE( 1 ) = mEdgeFunctionTS->E( k )( 0, 1 );

                // gradient operator master
                tBm = inv( tIntMaster->dNdXi( k ) * tXm ) * tIntMaster->dNdXi( k );

                // gradient operator slave
                tBs = inv( tIntSlave->dNdXi( k ) * tXs ) * tIntSlave->dNdXi( k );

                // normal component of H ( needed for material properties)
                mHn( k ) = -dot( tn.vector_data(), tBm.matrix_data() * tPhiM.vector_data() );

                // Matrix< real > txm(  tIntMaster->N( k ).matrix_data() * tXm.matrix_data() ) ;
                // Matrix< real > txs(  tIntSlave->N( k ).matrix_data() * tXs.matrix_data() ) ;
                // txm.print("xm");
                // txs.print("xs");

                // std::cout << "n " << k << " " << mHn( k ) << " " <<  -dot( tn.vector_data(), tBs.matrix_data() * tPhiS.vector_data() ) << std::endl ;

                /*if( aElement->id() == 17 && k == 2 )
                {
                    tn.print("n");
                    Vector< real > bm( tBm.matrix_data() * tPhiM.vector_data() );
                    Vector< real > bs( tBs.matrix_data() * tPhiS.vector_data() ) ;
                    bm.print("bm");
                    bs.print("bs");

                }*/
                tWDetJ /= mDeltaTime ;

                crossmat( tn , tBm, tnxB );
                //real tH1 = dot( tnxB, tPhiM ) ;

                p = mNumberOfNodesPerElement+mNumberOfNodesPerElement;
                q = p + mNumberOfThinShellLayers * 4 + 2 ;

                for( uint j=0; j<mEdgeDofMultiplicity; ++j )
                {
                    for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                    {
                        tValue = tWDetJ * mE( j ) * tnxB( i );
                        tK( i, q+j ) += tValue ;
                        tK( q+j, i ) += tValue ;
                    }
                    for( uint i=0; i<mEdgeDofMultiplicity; ++i )
                    {
                        tValue = tWDetJ * mE( i ) * mE( j );
                        tK( p+i, q+j ) += tValue ;
                        tK( q+j, p+i ) += tValue ;
                    }
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // lambda bc for slave
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                p = q - 2 ;
                q += 2 ;
                crossmat( tn , tBs, tnxB );
                //real tH2 = dot( tnxB, tPhiS ) ;

                // hack x-coordinate for output
                //real tXi =  mGroup->integration_points()( 0, k ) ;

                //real tx = 0.5 * ( aElement->element()->node( 0 )->x() * ( 1 - tXi)
                 //       + aElement->element()->node( 1 )->x() * ( 1 + tXi) );

                //std::cout << aElement->id() << " j: " << tx << " " << tH1 << " " << tH2 << " " << ( tH2-tH1 ) / mGroup->thin_shell_thickness( 0 ) << std::endl ;

                for( uint j=0; j<mEdgeDofMultiplicity; ++j )
                {
                    for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                    {
                        tValue = tWDetJ * mE( j ) * tnxB( i );
                        tK( mNumberOfNodesPerElement+i, q+j ) += tValue ;
                        tK( q+j, mNumberOfNodesPerElement+i ) += tValue ;
                    }
                    for( uint i=0; i<mEdgeDofMultiplicity; ++i )
                    {
                        tValue = tWDetJ * mE( i ) * mE( j );
                        tK( p+i, q+j ) += tValue ;
                        tK( q+j, p+i ) += tValue ;
                    }
                }

#ifdef HPHI_TRI6_LAMBDAN
                q = mNumberOfDofsPerElement - 2 ;

                for( uint j=0; j<mEdgeDofMultiplicity; ++j )
                {
                    for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                    {
                        tValue = tWDetJ * mE( j ) * ( tn( 0 ) * tBm( 0, i ) + tn( 1 ) * tBm( 1, i ) );
                        tK( i, q+j ) += tValue ;
                        tK( q+j, i ) += tValue ;
                    }
                    for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                    {
                        tValue = tWDetJ * mE( j ) * ( tn( 0 ) * tBs( 0, i ) + tn( 1 ) * tBs( 1, i ) );
                        tK( mNumberOfNodesPerElement+i, q+j ) -= tValue ;
                        tK( q+j, mNumberOfNodesPerElement+i ) -= tValue ;
                    }
                }
#endif
            }

            /*mAHn.fill( 0.0 );
            mBHn.fill( 0.0 );
            // derivative of Hn to xi
            for ( uint k = 0; k < mNumberOfIntegrationPoints; ++k )
            {
                real tXi = mGroup->integration_points()( 0, k ) ;
                mAHn( 0, 0 ) += tW( k ) * tXi * tXi ;
                mAHn( 1, 0 ) += tW( k ) * tXi ;
                mAHn( 0, 1 ) += tW( k ) * tXi ;
                mAHn( 1, 1 ) += tW( k ) ;

                mBHn( 0 ) += tW( k ) * mHn( k ) * tXi ;
                mBHn( 1 ) += tW( k ) * tXi ;
            }
            gesv( mAHn, mBHn, mPivot );
            real tdHndX = -mBHn( 0 ) / tLength * 2.0 ; */

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Mass and Stiffness contributions
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            for( uint l=0; l<mNumberOfThinShellLayers; ++l)
            {
                mArho.fill( 0.0 );
                mBrho.fill( 0.0 );

                this->collect_edge_data_from_layer( aElement, "edge_h", l, tHt );

                // initial offset for node coordinate
                p = 2 * mNumberOfNodesPerElement ;

                for ( uint k = 0; k < mNumberOfIntegrationPoints; ++k )
                {
                    // evaluate points for edge function
                    mE( 0 ) = mEdgeFunctionTS->E( k )( 0, 0 );
                    mE( 1 ) = mEdgeFunctionTS->E( k )( 0, 1 );

                    this->compute_layer_mass( l, mE( 0 ), mE( 1 ), tMlayer );

                    this->compute_layer_stiffness( l, k, mE( 0 ), mE( 1 ), tHt, mHn( k ), 0.0, tLength, tKlayer );

                    for ( uint j = 0; j < 6; ++j )
                    {
                        for ( uint i = 0; i < 6; ++i )
                        {
                            tValue = mW( k ) * tMlayer( i, j );
                            tM( p + i, p + j ) += tValue;
                            tValue = mW( k ) * tKlayer( i, j );
                            tK( p + i, p + j ) += tValue;
                        }
                    }
                }

                // compute coeffs for rho
                gesv( mArho, mBrho, mPivot );

                // critical current
                const Material * tMaterial = mGroup->thin_shell_material( l );
                real tRhoCrit = tMaterial->type() == MaterialType::Maxwell ?
                        tMaterial->rho_el_crit( 0.0, 0.0, 0.0 ) : 0.0 ;


                this->compute_layer_stabilizer( tHt, tLength,
                                                mGroup->thin_shell_thickness( l ), tRhoCrit, tMlayer );

                for ( uint j = 0; j < 6; ++j )
                {
                    for ( uint i = 0; i < 6; ++i )
                    {
                        tM( p + i, p + j ) += tMlayer( i, j );
                    }
                }


                // jump to next layer
                p += 4 ;

            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Jacobian and RHS
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // compute the right hand side
            // const Vector< real > & tQ0 = this->collect_q0_thinshell( aElement ) ;

            aRHS = tM *  this->collect_q0_thinshell( aElement ) ;

            // finalize the Jacobian
            aJacobian += mDeltaTime * tK ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Layer Data
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // write layer data into mesh
            mLayerData *= 0.25 ;
            for( uint l=0; l<mNumberOfThinShellLayers; ++l )
            {
                // compute the field index
                index_t tIndex = mGhostElementMap(
                        aElement->id() * mGroup->number_of_thin_shell_layers() + l )->index();
                real tJz = mLayerData(0, l );

                // real tRho = mGroup->thin_shell_material( l )->rho_el( tJz );

                mMesh->field_data( "elementJ" )( tIndex ) =  tJz ;
                //mMesh->field_data( "elementEJ" )( tIndex ) = tRho * tJz * tJz ;
                mMesh->field_data( "elementEJ" )( tIndex ) = mLayerData( 1, l );

                // std::cout << "check " << aElement->id() << " " << tJz * tLength * mGroup->thin_shell_thickness( l ) << " " << mLayerData( 2, l ) << std::endl ;
            }

        }
#ifdef BELFEM_GCC
#pragma GCC diagnostic pop
#endif

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tri6::compute_layer_mass(
                const uint aLayer,
                const real aE0,
                const real aE1,
                Matrix< real > & aM )
        {
            mWork[ 0 ] =  constant::mu0 * mGroup->thin_shell_thickness( aLayer ) / 30.0 ;
            mWork[ 1 ] = aE0 * aE0 * mWork[ 0 ] ;
            mWork[ 2 ] = aE0 * aE1 * mWork[ 0 ] ;
            mWork[ 3 ] = aE1 * aE1 * mWork[ 0 ] ;
            mWork[ 4 ] = mWork[ 1 ] + mWork[ 1 ] ; // 2* e0 * e0 * t/30
            mWork[ 5 ] = mWork[ 2 ] + mWork[ 2 ] ; // 2* e0 * e1 * t/30
            mWork[ 6 ] = mWork[ 3 ] + mWork[ 3 ] ; // 2* e1 * e1 * t/30
            mWork[ 7 ] = mWork[ 4 ] + mWork[ 4 ] ; // 4* e0 * e0 * t/30
            mWork[ 8 ] = mWork[ 5 ] + mWork[ 5 ] ; // 4* e0 * e1 * t/30
            mWork[ 9 ] = mWork[ 6 ] + mWork[ 6 ] ; // 4* e1 * e1 * t/30

            aM( 0, 0 ) =  mWork[ 7 ];
            aM( 1, 0 ) =  mWork[ 8 ];
            aM( 2, 0 ) =  mWork[ 4 ];
            aM( 3, 0 ) =  mWork[ 5 ];
            aM( 4, 0 ) = -mWork[ 1 ];
            aM( 5, 0 ) = -mWork[ 2 ];
            aM( 0, 1 ) =  mWork[ 8 ];
            aM( 1, 1 ) =  mWork[ 9 ];
            aM( 2, 1 ) =  mWork[ 5 ];
            aM( 3, 1 ) =  mWork[ 6 ];
            aM( 4, 1 ) = -mWork[ 2 ];
            aM( 5, 1 ) = -mWork[ 3 ];
            aM( 0, 2 ) =  mWork[ 4 ];
            aM( 1, 2 ) =  mWork[ 5 ];
            aM( 2, 2 ) =  mWork[ 1 ]*16.;
            aM( 3, 2 ) =  mWork[ 2 ]*16.;
            aM( 4, 2 ) =  mWork[ 4 ];
            aM( 5, 2 ) =  mWork[ 5 ];
            aM( 0, 3 ) =  mWork[ 5 ];
            aM( 1, 3 ) =  mWork[ 6 ];
            aM( 2, 3 ) =  mWork[ 2 ]*16.;
            aM( 3, 3 ) =  mWork[ 3 ]*16.;
            aM( 4, 3 ) =  mWork[ 5 ];
            aM( 5, 3 ) =  mWork[ 6 ];
            aM( 0, 4 ) = -mWork[ 1 ];
            aM( 1, 4 ) = -mWork[ 2 ];
            aM( 2, 4 ) =  mWork[ 4 ];
            aM( 3, 4 ) =  mWork[ 5 ];
            aM( 4, 4 ) =  mWork[ 7 ];
            aM( 5, 4 ) =  mWork[ 8 ];
            aM( 0, 5 ) = -mWork[ 2 ];
            aM( 1, 5 ) = -mWork[ 3 ];
            aM( 2, 5 ) =  mWork[ 5 ];
            aM( 3, 5 ) =  mWork[ 6 ];
            aM( 4, 5 ) =  mWork[ 8 ];
            aM( 5, 5 ) =  mWork[ 9 ];
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tri6::compute_layer_stiffness(
                const uint aLayer,
                const uint aIntPoint,
                const real aE0,
                const real aE1,
                const Vector< real > & aHt,
                const real aHn,
                const real adHndX,
                const real  aXLength,
                Matrix< real > & aK )
        {
            // reset matrix
            aK.fill( 0.0 );

            // grab material
            const Material * tMaterial = mGroup->thin_shell_material( aLayer );

            // integration data
            const IntegrationData * tInteg = mGroup->thinshell_integration();

            const Vector< real > & tW = tInteg->weights() ;

            // get material
            // temperature, todo: make not constant
            real tT = 4.0 ;

            // grab layer thickness
            const real tThickness = mGroup->thin_shell_thickness( aLayer );
            real tDetJ = tThickness * 0.5 ;

            // element-wise curernt
            real tJz_el = 0 ;

            // value for ac losses (element-wise)
            real tEJ_El = 0.0 ;

            const Matrix< real > & tXi = tInteg->points() ;

            // note:: usually, line elements are numbered like this
            //        0 --- 2 --- 1
            //        but here, we use the scheme
            //        0 --- 1 --- 2
            //        that's why we need to swap the indices in mHt and tPhiXi!

            // edge-wise components of H
            mHt( 0 ) = aE0 * aHt( 0 ) + aE1 * aHt( 1 ) ; // field on master side
            mHt( 2 ) = aE0 * aHt( 2 ) + aE1 * aHt( 3 ) ; // middle field
            mHt( 1 ) = aE0 * aHt( 4 ) + aE1 * aHt( 5 ) ; // field on slave side

            // integrate over thickness
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {

                // derivative of shape function along thickness
                const Vector< real > & tPhiXi = tInteg->dphidxi( k ) ;

                // curl operator
                mC( 0, 0 ) = aE0 * tPhiXi( 0 );
                mC( 0, 1 ) = aE1 * tPhiXi( 0 );
                mC( 0, 2 ) = aE0 * tPhiXi( 2 );
                mC( 0, 3 ) = aE1 * tPhiXi( 2 );
                mC( 0, 4 ) = aE0 * tPhiXi( 1 );
                mC( 0, 5 ) = aE1 * tPhiXi( 1 );
                mC /= tDetJ ;


                // mHt.print("Ht");

                // current : tHt derived along thickness
                // real tJz = mC * aHt ;
                //std::cout << "x y " << adHndx << " " <<  dot( tPhiXi, mHt ) / tDetJ << std::endl ;

                //real tJz = -adHndx + dot( tPhiXi, mHt ) / tDetJ ;
                real tJz = dot( tPhiXi, mHt ) / tDetJ - adHndX ;

                // tangential component of h
                real tHt =   dot( tInteg->phi( k ), mHt );

                // std::cout << "   " << k << " " << tHt << " " << tJz << std::endl ;

                // resistivity
                real tRho = tMaterial->rho_el( tJz, tT, constant::mu0 * std::sqrt( tHt*tHt + aHn*aHn )  );
                tRho = tRho < BELFEM_EPSILON ? BELFEM_EPSILON : tRho ;

                real tX = 0.5 * ( tXi( 0, aIntPoint ) + 1.0 ) * aXLength ;
                real tY = 0.5 * ( tXi( 0, k )  + 1.0 ) * tThickness ;

                const Vector< real > & tQ = this->eval_quad9( tX, tY );
                for( uint j=0; j<9; ++j )
                {
                    for ( uint i = 0; i < 9; ++i )
                    {
                        mArho( i, j ) += tQ( i ) * tQ( j );
                    }
                    mBrho( j ) += std::log( tRho ) * tQ( j );
                }

                // contribution to stiffness matrix
                aK += tW( k ) * trans( mC ) * tRho * mC * tDetJ ;

                tJz_el += tW( k ) * tJz ;
                tEJ_El += tW( k ) *  tRho * tJz * tJz ;
            }

            mLayerData( 0, aLayer ) += tW( aIntPoint ) * tJz_el ;
            mLayerData( 1, aLayer ) += tW( aIntPoint ) * tEJ_El ;

            mLayerData( 2, aLayer ) += tW( aIntPoint ) *
                    (    aE0 * ( aHt( 0 ) - aHt( 4 ) )
                            + aE1 * ( aHt( 1 ) - aHt( 5 ) ) )  * aXLength * 2 ;


        }

//------------------------------------------------------------------------------

#ifdef BELFEM_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif


        void
        IWG_Maxwell_HPhi_Tri6::compute_layer_stabilizer(
                const Vector< real > & aHt,
                const real aXLength,
                const real aYLength,
                const real aRhoCrit,
                Matrix< real > & aG )
        {
            /*real tC = 0 ;
            for( uint k=0; k<6; ++k )
            {
                tC += aL( k, k )*aL( k, k ) ;
            }
            tC = std::sqrt( tC ); */

            // integration data
            const IntegrationData * tInteg = mGroup->thinshell_integration();

            const Vector< real > & tW = tInteg->weights() ;

            const Matrix< real > & tXi = tInteg->points() ;

            aG.fill( 0.0 );

            real tC = ( aXLength * aYLength ) ;


            for( uint l=0; l<mNumberOfIntegrationPoints; ++l )
            {
                real xi = tXi( 0, l );
                real tX = 0.5 * ( xi + 1.0 ) * aXLength ;


                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {

                    real eta = tXi( 0, k ) ;
                    real tY = 0.5 * ( eta + 1.0 ) * aYLength ;

                    // derivatives of rho
                    real tRho = std::exp( dot( this->eval_quad9( tX, tY ), mBrho ) );
                    real tdRhodX = tRho * dot( this->eval_quad9_dx( tX, tY ), mBrho );
                    real tdRhodY = tRho * dot( this->eval_quad9_dy( tX, tY ), mBrho );

                    real tVal = 2.0 / aYLength ;

                    mCrho( 0 ) =  ( 0.25*(1.0 - 2.0*eta)*(3.0*xi - 1.0) ) * tVal ;
                    mCrho( 1 ) =  ( 0.25*(2.0*eta - 1.0)*(3.0*xi + 1.0) ) * tVal ;
                    mCrho( 2 ) =  ( eta*(3.0*xi - 1.0)  ) * tVal ;
                    mCrho( 3 ) =  ( -eta*(3.0*xi + 1.0)  ) * tVal ;
                    mCrho( 4 ) =  ( -0.25*(2.0*eta + 1.0)*(3.0*xi - 1.0) ) * tVal ;
                    mCrho( 5 ) =  ( 0.25*(2.0*eta + 1.0)*(3.0*xi + 1.0)  ) * tVal ;

                    mL( 0, 0 ) = tdRhodY * mCrho( 0 );
                    mL( 0, 1 ) = tdRhodY * mCrho( 1 );
                    mL( 0, 2 ) = tdRhodY * mCrho( 2 );
                    mL( 0, 3 ) = tdRhodY * mCrho( 3 );
                    mL( 0, 4 ) = tdRhodY * mCrho( 4 );
                    mL( 0, 5 ) = tdRhodY * mCrho( 5 );

                    mL( 1, 0 ) = -tdRhodX * mCrho( 0 );
                    mL( 1, 1 ) = -tdRhodX * mCrho( 1 );
                    mL( 1, 2 ) = -tdRhodX * mCrho( 2 );
                    mL( 1, 3 ) = -tdRhodX * mCrho( 3 );
                    mL( 1, 4 ) = -tdRhodX * mCrho( 4 );
                    mL( 1, 5 ) = -tdRhodX * mCrho( 5 );


                    tVal = 4.0 / ( aYLength * aYLength ) * tRho ;

                    mL( 0, 0 ) += ( 0.5 - 1.5*xi )  * tVal ;
                    mL( 0, 1 ) += ( 1.5*xi + 0.5 )  * tVal ;
                    mL( 0, 2 ) += ( 3.0*xi - 1.0 )  * tVal ;
                    mL( 0, 3 ) += ( -3.0*xi - 1.0 ) * tVal ;
                    mL( 0, 4 ) += ( 0.5 - 1.5*xi )  * tVal ;
                    mL( 0, 5 ) += ( 1.5*xi + 0.5 )  * tVal ;

                    tVal = -4.0 / ( aXLength * aYLength ) * tRho ;

                    mL( 1, 0 ) += ( 0.75 - 1.5*eta )  * tVal ;
                    mL( 1, 1 ) += ( 1.5*eta - 0.75 )  * tVal ;
                    mL( 1, 2 ) += ( 3.0*eta )  * tVal ;
                    mL( 1, 3 ) += ( -3.0*eta )  * tVal ;
                    mL( 1, 4 ) += ( -1.5*eta - 0.75 ) * tVal ;
                    mL( 1, 5 ) += ( 1.5*eta + 0.75 )  * tVal ;

                        // std::cout << "delta " << tDelta << std::endl ;
                    real tDelta = aRhoCrit == 0 ? tC : tC * std::pow( tRho / aRhoCrit, 2 );
                    aG += tW( k ) * tDelta * trans( mL ) * mL ;

                    /*mL( 0, 0 ) = 0.5 * ( 1. - 3. * xi );
                    mL( 0, 1 ) = 0.5 * ( 1. + 3. * xi );
                    mL( 1, 0 ) = 0.5 * eta * ( eta - 1.0 );
                    mL( 1, 1 ) = 1.0 - eta * eta ;
                    mL( 1, 2 ) = 0.5 * eta * ( eta +  1.0 );

                    mC( 0, 0 ) = mL( 1, 0 ) * mL( 0, 0 ) ;
                    mC( 0, 1 ) = mL( 1, 0 ) * mL( 0, 1 ) ;
                    mC( 0, 2 ) = mL( 1, 1 ) * mL( 0, 0 ) ;
                    mC( 0, 3 ) = mL( 1, 1 ) * mL( 0, 1 ) ;
                    mC( 0, 4 ) = mL( 1, 2 ) * mL( 0, 0 ) ;
                    mC( 0, 5 ) = mL( 1, 2 ) * mL( 0, 1 ) ; */
                }
            }



            /*real tD = 0 ;
            for( uint k=0; k<6; ++k )
            {
                tD += aL( k, k )*tDelta*aL( k, k ) ;
            }
            tD = std::sqrt( tD ); */

        }

#ifdef BELFEM_GCC
#pragma GCC diagnostic pop
#endif

//------------------------------------------------------------------------------

        const Vector< real > &
        IWG_Maxwell_HPhi_Tri6::collect_q0_thinshell( Element * aElement )
        {
            // grab the output vector
            Vector< real > & aQ0 = mGroup->work_nedelec() ;

            // aQ0.fill( 0.0 );

            BELFEM_ASSERT( mGroup->domain_type() == DomainType::ThinShell,
                           "function IWG_Maxwell_HPhi_Tri6::collect_q0_thinshell can only be applied to a thin shell" );

            // grab field data from mesh
            const Vector< real > & tPhi  = mMesh->field_data( "phi0" );

            uint tCount = 0 ;

            // get the node dofs
            for( uint k=0; k<mNumberOfNodesPerElement ; ++k )
            {
                aQ0( tCount++ ) = tPhi( aElement->master()->element()->node( k )->index() );
            }
            for( uint k=0; k<mNumberOfNodesPerElement ; ++k )
            {
                aQ0( tCount++ ) = tPhi( aElement->slave()->element()->node( k )->index() );
            }


            // get the field
            Vector< real > & tData = mMesh->field_data( "edge_h0" );

            uint n = 2*mNumberOfThinShellLayers+1 ;

            // get the edge dofs
            for( uint l=0; l<n ; ++l )
            {

                // get the edge
                mesh::Edge * tEdge = mMesh->ghost_facet( aElement->id(), l )->edge( 0 );

                // get the data
                if( aElement->edge_direction( 0 ) )
                {
                    for ( uint k = 0; k < mEdgeDofMultiplicity; ++k )
                    {
                        aQ0( tCount++ ) = tData( mEdgeDofMultiplicity * tEdge->index() + k );
                    }
                }
                else
                {
                    for ( int k = mEdgeDofMultiplicity-1; k >= 0; k-- )
                    {
                        aQ0( tCount++ ) = tData( mEdgeDofMultiplicity * tEdge->index() + k );
                    }
                }
            }

            //todo: need to introduce edge-multiplicity here
            for( const string & tLabel : mFields.ThinShellLast )
            {
                // get lambda field
                const Vector< real > & tL = mMesh->field_data( tLabel );
                aQ0( tCount ) = tL( aElement->dof( tCount )->dof_index_on_field() );
                ++tCount ;
            }
            return aQ0 ;
        }

//------------------------------------------------------------------------------

        real
        IWG_Maxwell_HPhi_Tri6::compute_element_current( Element * aElement )
        {
            // link edge funcition with element
            mEdgeFunction->link( aElement );

            Matrix< real > & tX = mGroup->work_X() ;
            aElement->get_node_coors( tX );

            const Vector< real > & tW = mGroup->integration_weights() ;

            real aI = 0.0 ;
            real tV = 0.0 ;

            // grab edge data
            this->collect_nedelec_data( aElement,
                                        NedelecField::H,
                                        mGroup->work_nedelec() );

            Vector< real > tJz( 1 );

            // loop over all integration points
            for ( uint k = 0; k < mNumberOfIntegrationPoints; ++k )
            {
                // get operator for curl function
                const Matrix< real > & tC = mEdgeFunction->C( k );

                // compute current
                tJz = tC * mGroup->work_nedelec() ;

                // add current to integral
                aI += tW( k ) * tJz( 0 ) * mEdgeFunction->abs_det_J();

                // compute volume
                tV += tW( k ) *  mEdgeFunction->abs_det_J();
            }

            mMesh->field_data( "elementJ")( aElement->element()->index() ) = aI / tV ;

            return aI ;
        }

//------------------------------------------------------------------------------
    }
}