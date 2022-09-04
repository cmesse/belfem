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
            mFields.ThinShell = { "lambda_m0", "lambda_m1", "lambda_s0", "lambda_s1" };
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
                    mLayerData.set_size( mFields.Ghost.size(), mGroup->number_of_thin_shell_layers() );

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
                    mGroup->work_J() = tNodeFunction->dNdXi( k ) * mGroup->node_coords() ;

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
                tB = inv( tNodeFunction->dNdXi( 0 ) * mGroup->node_coords() ) * tIntNxi ;

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
            aElement->slave()->get_node_coors( mGroup->node_coords() );
            // the B-Matrix
            Matrix< real > & tK = mGroup->work_K() ;
            Matrix< real > & tM = aJacobian ;
            Matrix< real > & tJ = mGroup->work_J() ;
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

                    tJ = tSlave->dNdXi( k ) * mGroup->node_coords() ;
                    tB = inv( tJ ) * tSlave->dNdXi( k );
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
                tJ = tSlave->dNdXi( 0 ) * mGroup->node_coords() ;
                Matrix< real > & tInvJ = mGroup->work_invJ();
                tInvJ = inv( tJ );

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
                aElement->get_node_coors( mGroup->work_X());

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
            Matrix< real > & tX = mGroup->work_X() ;
            this->collect_node_coords( aElement->master(), tX );

            // get integration data from master
            const IntegrationData * tMaster =
                    mGroup->master_integration( aElement->facet()->master_index() );

            // get the integration weights
            const Vector< real > tW = tMaster->weights() ;

            // the gradient operator matrix
            Matrix< real > & tB = mGroup->work_B() ;

            if( aElement->master()->element()->is_curved() )
            {
                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // compute the gradient operator
                    tB.matrix_data() = inv( tMaster->dNdXi( k ) * tX ) * inv( tMaster->dNdXi( k ) * tX ) ;

                    // compute the normal
                    const Vector< real > & tn = this->normal_curved_2d( aElement, k );

                    real tScale = tW( k ) * mGroup->work_det_J() * constant::mu0 ;

                    // loop over all nodes
                    for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                    {
                        aJacobian( 0, i ) += tScale * (
                                  tn( 0 ) * tB( 0, i )
                                + tn( 1 ) * tB( 1, i ) ) ;
                    }
                }

                // make matrix symmetric
                for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                {
                    aJacobian( i, 0 ) =  aJacobian( 0, i ) ;
                }

            }
            else
            {
                // grab the inverse jacobian
                Matrix< real > & tInvJ = mGroup->work_J() ;
                tInvJ =  inv( tMaster->dNdXi( 0 ) * tX ) ;

                // compute the normal
                const Vector< real > & tn = this->normal_straight_2d( aElement );

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // compute the gradient operator
                    tB = tInvJ * tMaster->dNdXi( k );


                    // loop over all nodes
                    for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                    {
                        aJacobian( 0, i ) += tW( k ) * (
                                  tn( 0 ) * tB( 0, i )
                                + tn( 1 ) * tB( 1, i ) ) ;
                    }
                }

                // make matrix symmetric
                for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                {
                    aJacobian( i, 0 ) =  aJacobian( 0, i ) ;
                }

                // scale matrix ( det J has been computed by normal )
                aJacobian *= mGroup->work_det_J() * constant::mu0 ;
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
            mEdgeFunctionTS->link( aElement );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // containers and symbols
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // mass matrix
            Matrix< real > & tM = aJacobian;

            // stiffness matrix
            Matrix< real > & tK = mGroup->work_K();

            // the sign of the edge, in this case, we use the physical tag
            const real tSign = aElement->element()->physical_tag() == 1 ? 1.0 : -1.0;

            BELFEM_ASSERT( tSign == ( aElement->edge_direction( 0 ) ? 1.0 : -1.0 ),
                           "Invalid Edge Direction for Element %lu",
                           ( long unsigned int ) aElement->id());

            // container for geometry jacobian for master and slave
            Matrix< real > & tJ = mGroup->work_J();

            // container for gradient operator, master side
            Matrix< real > & tBm = mGroup->work_B();

            // container for gradient operator, slave side
            Matrix< real > & tBs = mGroup->work_D();

            // coordinates for master
            Matrix< real > & tX = mGroup->work_X();

            // coordinates for master
            Matrix< real > & tXm = mGroup->work_Xi();

            // coordinates for slave
            Matrix< real > & tXs = mGroup->work_Eta();

            // expression cross( n, B )
            Vector< real > tnxB = mGroup->work_psi();


            // mass matrix for individual layer
            Matrix< real > & tMlayer = mGroup->work_M();

            // stiffness matrix for individual layer
            Matrix< real > & tKlayer = mGroup->work_L();

            // additional indices
            uint o;
            uint p;
            uint q;
            uint r;
            uint s;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // reset matrices
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            tM.fill( 0.0 );
            tK.fill( 0.0 );
            aRHS.fill( 0.0 );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // geometry considerations and integration data
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // collect coordinates form element
            this->collect_node_coords( aElement, tX );

            // collect coordinates from master
            this->collect_node_coords( aElement->master(), tXm );

            // collect coordinates from slave
            this->collect_node_coords( aElement->slave(), tXs );

            // get integration data from master
            const IntegrationData * tMaster =
                    mGroup->master_integration( aElement->facet()->master_index());

            // get integration data from slave
            const IntegrationData * tSlave =
                    mGroup->slave_integration( aElement->facet()->slave_index());

            // integration weights along edge
            const Vector< real > & tW = mGroup->integration_weights();

            BELFEM_ASSERT( tMaster->weights().length() == mNumberOfIntegrationPoints,
                           "number of integraiton points on master does not match" );

            BELFEM_ASSERT( tSlave->weights().length() == mNumberOfIntegrationPoints,
                           "number of integraiton points on slave does not match" );

            // data for magnetic field (H-data)
            Vector< real > & tHt = mGroup->work_sigma();
            Vector< real > & tPhi = mGroup->work_phi();
            this->collect_node_data( aElement->master(), "phi", tPhi );

            // reset the thin shell data
            mLayerData.fill( 0.0 );

            //if ( aElement->element()->is_curved() )

            // needed for normal computation
            aElement->get_node_coors( mGroup->work_X());

            real tdM;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // lambda-condition master
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // get phi field
            this->collect_node_data( aElement->master(), "phi", tPhi );

            const IntegrationData * tIntMaster = mGroup->master_integration( aElement->facet()->master_index());
            const IntegrationData * tIntSlave = mGroup->slave_integration( aElement->facet()->slave_index() );

            r = 0 ;
            s = mNumberOfIntegrationPoints - 1 ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // lambda conditions for first dof
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // air side
            p = 0 ;
            q = 2 * mNumberOfNodesPerElement + ( 2 * mNumberOfThinShellLayers + 1 ) * mEdgeDofMultiplicity ;

            // first point  master
            tBm = inv( tIntMaster->dNdXi( r ) * tXm ) * tIntMaster->dNdXi( r );
            crossmat( this->normal_curved_2d( aElement, r ) , tBm, tnxB );

            for( uint i=0; i<mNumberOfNodesPerElement; ++i )
            {
                tM( p, q ) += tnxB( i ) ;
                tM( q, p ) += tnxB( i ) ;
                ++p ;
            }

            // edge side
            o = 2 * mNumberOfNodesPerElement ;
            tM( o, q ) = 1.0 ;
            tM( q, o ) = 1.0 ;

            // first point slave
            p = mNumberOfNodesPerElement ;
            q += 2 ;

            tBs = inv( tIntSlave->dNdXi( r ) * tXs ) * tIntSlave->dNdXi( r );
            crossmat( mNormal2D , tBs, tnxB );

            for( uint i=0; i<mNumberOfNodesPerElement; ++i )
            {
                tM( p, q ) += tnxB( i ) ;
                tM( q, p ) += tnxB( i ) ;
                ++p ;
            }

            o = 2 * mNumberOfNodesPerElement + ( 2 * mNumberOfThinShellLayers ) * mEdgeDofMultiplicity ;
            tM( o, q ) =  1.0 ;
            tM( q, o ) =  1.0 ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // lambda conditions for second dof, air side
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // second point master
            p = 0 ;
            q = 2 * mNumberOfNodesPerElement + ( 2 * mNumberOfThinShellLayers + 1 ) * mEdgeDofMultiplicity + 1 ;

            tBm = inv( tIntMaster->dNdXi( s ) * tXm ) * tIntMaster->dNdXi( s );
            crossmat( this->normal_curved_2d( aElement, s ) , tBm, tnxB );
            for( uint i=0; i<mNumberOfNodesPerElement; ++i )
            {
                tM( p, q ) += tnxB( i ) ;
                tM( q, p ) += tnxB( i ) ;
                ++p ;
            }

            // edge side
            o = 2 * mNumberOfNodesPerElement + 1 ;
            tM( o, q ) =  mGroup->work_det_J() ;
            tM( q, o ) =  mGroup->work_det_J() ;

            // second point slave
            p = mNumberOfNodesPerElement ;
            q += 2 ;
            tBs = inv( tIntSlave->dNdXi( s ) * tXs ) * tIntSlave->dNdXi( s );
            crossmat( mNormal2D , tBs, tnxB );

            for( uint i=0; i<mNumberOfNodesPerElement; ++i )
            {
                tM( p, q ) += tnxB( i ) ;
                tM( q, p ) += tnxB( i ) ;
                ++p ;
            }

            // edge side
            o =  2 * mNumberOfNodesPerElement + ( 2 * mNumberOfThinShellLayers ) * mEdgeDofMultiplicity + 1 ;
            tM( o, q ) = 1.0 ;
            tM( q, o ) = 1.0 ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // contribution for mass and stiffness matrix
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            for ( uint k = 0; k < mNumberOfIntegrationPoints; ++k )
            {
                // compute length increment
                const Vector< real > & tn = this->normal_curved_2d( aElement, k );

                // integration weight
                real tScale = tW( k ) * mGroup->work_det_J() ;

                mE( 0 ) = mEdgeFunctionTS->E( k )( 0, 0 );
                mE( 1 ) = mEdgeFunctionTS->E( k )( 0, 1 );

                // normal component of H
                mH.vector_data() = inv( tIntMaster->dNdXi( k ) * tXm ) * tIntMaster->dNdXi( k ) * tPhi ;
                real tHn = -dot( tn, mH );

                // DELETE ME
                if( aElement->id() == 75  )
                {
                    this->collect_node_data( aElement->master(), "phi", tPhi );
                    tBm = inv( tIntMaster->dNdXi( k ) * tXm ) * tIntMaster->dNdXi( k );

                    crossmat( this->normal_curved_2d( aElement, k ) , tBm, tnxB );

                    this->collect_edge_data_from_layer( aElement, "edge_h", 0, tHt );

                    //real tA0 = -dot( tnxB, tPhi );

                    /*std::cout << "--- " << k << " " << aElement->facet()->master_index() << " " << aElement->facet()->slave_index() << std::endl ;
                    tIntMaster->N( k ).print("N1");
                    tXm.print("X1");
                    tIntSlave->N( k ).print("N2");
                    tXs.print("X2");
                    mEdgeFunctionTS->N( k ).print("N");
                    tX.print("X");*/


                    //tIntMaster->dNdXi( k ).print("N_xi");
                    //tXm.print("X");

                     //mNormal2D.print("n");
                    //tBm.print("B");
                    //tPhi.print("Phi");


                    this->collect_node_data( aElement->slave(), "phi", tPhi );
                    tBm = inv( tIntSlave->dNdXi( k ) * tXm ) * tIntSlave->dNdXi( k );
                    crossmat( this->normal_curved_2d( aElement, k ) , tBm, tnxB );
                    this->collect_edge_data_from_layer( aElement, "edge_h", 0, tHt );
                    real tA1 = -dot( tnxB, tPhi );

                    real tA2 = mE( 0 ) * tHt( 0 ) + mE( 1 ) * tHt( 1 );

                    std::cout << " check " << k << " " << tA1 << " " << tA2 << " " << tA1/tA2 << std::endl;
                }

                p = 2 * mNumberOfNodesPerElement ;

                for( uint l=0; l<mNumberOfThinShellLayers; ++l )
                {
                    this->collect_edge_data_from_layer( aElement, "edge_h", l, tHt );

                    this->compute_layer_mass( l, mE( 0 ), mE( 1 ),  tMlayer );

                    this->compute_layer_stiffness( l, k, mE( 0 ), mE( 1 ),  tHt, tHn, tKlayer );

                    // add submatrix to macro element
                    for( uint j=0; j<6; ++j )
                    {
                        for( uint i=0; i<6; ++i )
                        {
                            tM( p+i, p+j ) += tScale * tMlayer( i, j );

                            tK( p+i, p+j)  += tScale * tKlayer( i, j );
                        }
                    }

                    // increment offset for the next layer
                    p+= 2 ;
                }

            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Jacobian and RHS
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // compute the right hand side
            aRHS = tM *  this->collect_q0_thinshell( aElement ) ;


            // hack
            q = 2 * mNumberOfNodesPerElement + ( 2 * mNumberOfThinShellLayers + 1 ) * mEdgeDofMultiplicity ;
            for( uint k=q; k<mNumberOfDofsPerElement; ++k )
            {
                aRHS( k ) = 0.0 ;
            }

            // finalize the Jacobian
            aJacobian += mDeltaTime * mGroup->work_K() ;


            if( aElement->id() == 75 )
            {
                this->print_dofs( aElement );
                const Vector< real > & tLHS =  this->collect_q0_thinshell( aElement );


                Vector< real > tRHS( aJacobian * tLHS );
                tLHS.print("LHS");
                aRHS.print("RHS0");
                tRHS.print("RHS1");
            }

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

                mMesh->field_data( "elementJ" )( tIndex ) =  mLayerData(0, l );
                //mMesh->field_data( "elementEJ" )( tIndex ) = mLayerData( 1, l );
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

            // integrate over thickness
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // note:: usually, line elements are numbered like this
                //        0 --- 2 --- 1
                //        but here, we use the scheme
                //        0 --- 1 --- 2
                //        that's why we need to swap the indices in mHt and tPhiXi!

                // edge-wise components of H
                mHt( 0 ) = aE0 * aHt( 0 ) + aE1 * aHt( 1 ) ;
                mHt( 2 ) = aE0 * aHt( 2 ) + aE1 * aHt( 3 ) ;
                mHt( 1 ) = aE0 * aHt( 4 ) + aE1 * aHt( 5 ) ;

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

                // current : tHt derived along thickness
                // real tJz = mC * aHt ;
                real tJz = dot( tPhiXi, mHt ) / tDetJ ;

                // tangential component of h
                real tHt =   dot( tInteg->phi( k ), mHt );

                // resistivity
                real tRho = tMaterial->rho_el( tJz, tT, constant::mu0 * std::sqrt( tHt*tHt + aHn*aHn ) );

                // contribution to stiffness matrix
                aK += tW( k ) * trans( mC ) * tRho * mC * tDetJ ;

                tJz_el += tW( k ) * tJz ;
                tEJ_El += tW( k ) *  tRho * tJz * tJz ;
            }

            mLayerData( 0, aLayer ) += tW( aIntPoint ) * tJz_el ;
            mLayerData( 1, aLayer ) += tW( aIntPoint ) * tEJ_El ;
        }

//------------------------------------------------------------------------------

        const Vector< real > &
        IWG_Maxwell_HPhi_Tri6::collect_q0_thinshell( Element * aElement )
        {
            // grab the output vector
            Vector< real > & aQ0 = mGroup->work_nedelec() ;

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

            // get the edge dofs
            for( uint l=0; l<mNumberOfThinShellLayers ; ++l )
            {
                tCount = this->collect_edge_data_from_layer( aElement, "edge_h0", l, aQ0, tCount );
            }

            // get facet index
            index_t tIndex = aElement->element()->edge( 0 )->index() ;

            //todo: need to introduce edge-multiplicity here
            for( const string & tLabel : mFields.ThinShellLast )
            {
                // get lambda field
                const Vector< real > & tL = mMesh->field_data( tLabel );
                aQ0( tCount++ ) = tL( tIndex );
                if( aElement->id() == 75 )
                {
                    std::cout << "#last " << tLabel << " " << tL( tIndex ) << std::endl ;
                }
            }

            return aQ0 ;
        }

//------------------------------------------------------------------------------
    }
}