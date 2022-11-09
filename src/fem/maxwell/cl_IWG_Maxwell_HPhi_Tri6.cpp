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
#include "fn_eigen.hpp"
#include "fn_max.hpp"
#include "fn_det.hpp"
#include "cl_IF_InterpolationFunction.hpp"
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


            mFields.InterfaceScAir = { "lambda_t0", "lambda_t1" };
            //mFields.InterfaceScFm = { "lambda_n0", "lambda_n1"  };

#ifdef BELFEM_FERROAIR_ENRICHED
            mFields.InterfaceFmAir = { "lambda_t0", "lambda_t1", "lambda_n0", "lambda_n1" };
            //mFields.InterfaceFmAir = { "lambda_t0", "lambda_t1" };

#endif

            mFields.Cut = { "lambda" };

            mFields.SymmetryAir = { "lambda_n0", "lambda_n1" };
            mFields.AntiSymmetryAir = { "lambda_n0", "lambda_n1" };

            mFields.InterfaceFmFm = { "_lambda" };

#ifdef BELFEM_FERRO_HPHIA
            mFields.Ferro = { "az" };
#ifdef BELFEM_FERRO_LINEAR
            mFields.SymmetryFerro = { "lambda" };
            mFields.AntiSymmetryFerro = { "lambda" };

#else
            mFields.SymmetryFerro = { "lambda_n0", "lambda_n1" };
            mFields.AntiSymmetryFerro = { "lambda_n0", "lambda_n1" };
#endif
#else
            mFields.Ferro = { "phi" };
            mFields.SymmetryFerro = { "lambda_n0", "lambda_n1" };
            mFields.AntiSymmetryFerro = { "lambda_n0", "lambda_n1" };
#endif

#ifdef HPHI_TRI6_LAMBDAN
            mFields.ThinShell = { "lambda_m0", "lambda_m1", "lambda_s0", "lambda_s1",  "lambda_t0", "lambda_t1" };
#else
            mFields.ThinShell = { "lambda_m0", "lambda_m1", "lambda_s0", "lambda_s1" };
#endif

            // non-dof fields
            mFields.MagneticFieldDensity     =  { "bx", "by", "bz" };
            mFields.CurrentDensity = {  "jz" };
            mFields.CurrentBC = { "lambda_I" };
            mFields.Ghost = { "elementEJ", "elementJ", "elementJJc", "elementRho" };

            mEdgeDofMultiplicity = 2 ;
            mFaceDofMultiplicity = 2 ;
            mLambdaDofMultiplicity = 1 ;

            mNumberOfRhsDofsPerEdge = 2 ;
            mNumberOfRhsDofsPerFace = 2 ;

            // create the TS function todo: tell if we want to use bernstein
            mEdgeFunctionTS = new EF_LINE3() ;
            mHt.set_size( 3 );
            mC.set_size( 1, 6 );
            mdCdx.set_size( 6 );
            mdCdy.set_size(  6 );

            mH.set_size( 2 );
            mU.set_size( 2 );
            mE.set_size( 2 );

            mHn.set_size( 16 );
            mW.set_size( 16 );
            mRho.set_size( 16, 16 );
            mJz.set_size( 16, 16 );
            mJc.set_size( 16, 16 );
            mdRhodx.set_size( 16, 16 );
            mdRhody.set_size( 16, 16 );
            mAHn.set_size( 2, 2 );
            mBHn.set_size( 2 );
            mPivot.set_size( 2 );


            // for stabilzation
            mL.set_size( 2, 6 );

            mXi.set_size( 5 );
            mEta.set_size( 5 );
            mPhi.set_size( 3 );
            mSigma.set_size( 5 );
            mdEdxi.set_size( 2 );
            mdEdy.set_size( 6 );
            md2Edy2.set_size( 6 );
            md2Edxdy.set_size( 6 );

            mN.set_size( 3 );
            mdNdeta.set_size( 3 );
            md2Ndeta2.set_size( 3 );
            mPsi.set_size( 6 );

            mX.set_size( 5 );
            mY.set_size( 5 );
            mT.set_size( 5, 9, 0.0 );

            mK.set_size( 5, 5 );
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
                case( DomainType::Conductor ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_superconductor ;
                    break ;
                }
                case( DomainType::Ferro ) :
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
                case( DomainType::InterfaceFmFm ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_fmfm ;
                    break ;
                }
                case( DomainType::InterfaceFmAir ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_fmair ;
                    break ;
                }
                case( DomainType::SymmetryAir ) :
                {
                    mFunJacobian =
                            &IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_symmetry_air ;
                    break;
                }
                case( DomainType::AntiSymmetryAir ) :
                {
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_antisymmetry_air ;
                    break ;
                }
                case( DomainType::SymmetryFerro ) :
                {
#ifdef BELFEM_FERRO_HPHIA
                    mFunJacobian =
                            &IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_symmetry_ferro ;
#else
                    mFunJacobian =
                            &IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_symmetry_air ;
#endif
                    break;
                }
                case( DomainType::AntiSymmetryFerro ) :
                {
#ifdef BELFEM_FERRO_HPHIA
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_antisymmetry_ferro ;
#else
                    mFunJacobian =
                            & IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_antisymmetry_air ;
#endif
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
                                    mGroup->boundary_condition())->subtype() )
                            {
                                case ( MagfieldBcType::Wave ) :
                                {
                                    mFunJacobian =
                                            &IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_wave;
                                    break;
                                }
                                case ( MagfieldBcType::Symmetry ) :
                                {
                                    mFunJacobian =
                                            &IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_symmetry_air ;
                                    break;
                                }
                                case ( MagfieldBcType::AntiSymmetry ) :
                                {
                                    mFunJacobian =
                                            &IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_antisymmetry_air ;
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
            // make sure that the facet is correctly oriented
            BELFEM_ASSERT( static_cast< DomainType >( aElement->master()->element()->physical_tag() ) == DomainType::Conductor,
                           "Master element of facet %lu must be of DomainType::Conductor", ( long unsigned int ) aElement->id() );

            BELFEM_ASSERT( static_cast< DomainType >( aElement->slave()->element()->physical_tag() ) == DomainType::Air,
                           "Master element of facet %lu must be of DomainType::Air", ( long unsigned int ) aElement->id() );

            // link edge function with element
            const IntegrationData * tNodeFunction = this->interface_data( aElement );

            // gradient operator
            Matrix< real > & tB    = mGroup->work_B() ;

            // get weight functions
            const Vector< real > & tW = mGroup->integration_weights() ;

            Vector< real > & tnxE  = mGroup->work_tau() ;
            Vector< real > & tnxB  = mGroup->work_sigma() ;

            // now we can build the stiffness matrix
            Matrix< real > & tK = aJacobian ;

            // reset K-matrix
            tK.fill( 0.0 );

            // scaling parameter
            real tScale = constant::mu0 * this->timestep() ;

            // Edge Function for interface
            mEdgeFunctionTS->link( aElement );;

            // help parameter
            real tWeight ;
            real tValue ;

            const uint p = 14 ;
            const uint q = 15 ;
            uint c ;

            real tPsi0 ;
            real tPsi1 ;

            const real & tDetJ = mGroup->work_det_J() ;

            if( aElement->master()->element()->is_curved() )
            {
                // needed for normal, but is also done by this->interface_data( aElement );
                // this->collect_node_coords( aElement->master(), mGroup->work_Xm() );
                // this->collect_node_coords( aElement->slave(), mGroup->work_Xs() );

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // gradient operator for air element
                    tB = inv( tNodeFunction->dNdXi( k ) * mGroup->work_Xs() )
                            * tNodeFunction->dNdXi( k ) ;

                    // node function at integration point for air element
                    const Vector< real > & tN = tNodeFunction->phi( k );

                    // edge function at integration point
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    // normal at integration point ( also computes mGroup->work_det_J() )
                    const Vector< real > & tn = this->normal_curved_2d( aElement, k );

                    crossmat( tn, tE, tnxE );
                    crossmat( tn, tB, tnxB );

                    // lambda multiplicators
                    tPsi0 = mEdgeFunctionTS->E( k )( 0, 0 );
                    tPsi1 = mEdgeFunctionTS->E( k )( 0, 1 );

                    // scaling value
                    tWeight = tW( k ) * tDetJ * tScale ;

                    // add components for parallel coupling
                    c = 0 ;

                    // edge dofs
                    for( uint e=0; e<8; ++e )
                    {
                        tValue = tWeight * tnxE( e );
                        tK( p, c ) += tPsi0 * tValue ;
                        tK( q, c ) += tPsi1 * tValue ;
                        tK( c, p ) += tPsi0 * tValue ;
                        tK( c, q ) += tPsi1 * tValue ;
                        ++c ;
                    }

                    // node dofs
                    for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                    {
                        tValue = tWeight * tnxB( i );
                        tK( p, c ) += tPsi0 * tValue ;
                        tK( q, c ) += tPsi1 * tValue ;
                        tK( c, p ) += tPsi0 * tValue ;
                        tK( c, q ) += tPsi1 * tValue ;
                        ++c ;
                    }

                    // add components for perpendicular coupling
                    for( uint j=0; j<8; ++j )
                    {
                        for( uint i=0; i<6; ++i )
                        {
                            tK( i+8, j ) -= tWeight * tN( i ) *
                                            ( tn( 0 ) * tE( 0, j )
                                              + tn( 1 ) * tE( 1, j ) );
                        }
                    }

                }
            }
            else
            {
                Matrix< real > tInvJ = mGroup->work_J() ;
                tInvJ = inv( tNodeFunction->dNdXi( 0 ) * mGroup->work_Xs() ) ;

                // normal at integration point ( also computes mGroup->work_det_J() )
                const Vector< real > & tn = this->normal_straight_2d( aElement );

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {
                    // gradient operator for air element
                    tB = tInvJ * tNodeFunction->dNdXi( k ) ;

                    // node function at integration point for air element
                    const Vector< real > & tN = tNodeFunction->phi( k );

                    // edge function at integration point
                    const Matrix< real > & tE = mEdgeFunction->E( k );

                    crossmat( tn, tE, tnxE );
                    crossmat( tn, tB, tnxB );

                    // lambda multiplicators
                    tPsi0 = mEdgeFunctionTS->E( k )( 0, 0 );
                    tPsi1 = mEdgeFunctionTS->E( k )( 0, 1 );

                    // scaling value
                    tWeight = tW( k ) * tDetJ * tScale ;

                    // add components for parallel coupling
                    c = 0 ;

                    // edge dofs
                    for( uint e=0; e<8; ++e )
                    {
                        tValue = tWeight * tnxE( e );
                        tK( p, c ) += tPsi0 * tValue ;
                        tK( q, c ) += tPsi1 * tValue ;
                        tK( c, p ) += tPsi0 * tValue ;
                        tK( c, q ) += tPsi1 * tValue ;
                        ++c ;
                    }

                    // node dofs
                    for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                    {
                        tValue = tWeight * tnxB( i );
                        tK( p, c ) += tPsi0 * tValue ;
                        tK( q, c ) += tPsi1 * tValue ;
                        tK( c, p ) += tPsi0 * tValue ;
                        tK( c, q ) += tPsi1 * tValue ;
                        ++c ;
                    }

                    // add components for perpendicular coupling
                    for( uint j=0; j<8; ++j )
                    {
                        for( uint i=0; i<6; ++i )
                        {
                            tK( i+8, j ) -= tWeight * tN( i ) *
                                            ( tn( 0 ) * tE( 0, j )
                                              + tn( 1 ) * tE( 1, j ) );
                        }
                    }

                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_fmair(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // make sure that the facet is correctly oriented
            BELFEM_ASSERT( static_cast< DomainType >( aElement->slave()->element()->physical_tag() ) == DomainType::Ferro,
                           "Slave element of facet %lu must be of DomainType::Ferro", ( long unsigned int ) aElement->id() );

            BELFEM_ASSERT( static_cast< DomainType >( aElement->master()->element()->physical_tag() ) == DomainType::Air,
                           "Master element of facet %lu must be of DomainType::Air", ( long unsigned int ) aElement->id() );


            const IntegrationData * tMaster =  mGroup->master_integration(
                    aElement->facet()->master_index() ) ;

            const IntegrationData * tSlave =  mGroup->slave_integration(
                    aElement->facet()->slave_index() ) ;


            // get node coords, needed to compute geometry Jacobian

            const Matrix< real > & tXi = mGroup->integration_points() ;

            // the B-Matrix
            Matrix< real > & tK = mGroup->work_K() ;
            Matrix< real > & tM = aJacobian ;

            const Vector< real > & tW = mGroup->integration_weights() ;

            tK.fill( 0.0 );
            tM.fill( 0.0 );


            // collect node coords, needed below and to compo
            aElement->master()->get_node_coors( mGroup->work_Xm() );
            aElement->slave()->get_node_coors( mGroup->work_Xs() );

            uint tNumNodesPerAir   = aElement->master()->element()->number_of_nodes() ;
            uint tNumNodesPerFerro = aElement->slave()->element()->number_of_nodes() ;


            Matrix< real > & tBm = mGroup->work_B() ;
            Matrix< real > & tBs = mGroup->work_C() ;




            uint p =  aElement->master()->number_of_dofs() + aElement->slave()->number_of_dofs() ;
            uint q = p + 1 ;
            uint r = q + 1 ;
            uint s = r + 1 ;

            // normal of air side
            const Vector< real > & tn = this->normal_straight_2d( aElement );

            Vector< real > tAz( mNumberOfNodesPerSlave );
            this->collect_node_data( aElement->slave(), "az", tAz );

            Vector< real > & tPhi = mPsi ;
            this->collect_node_data( aElement->master(), "phi", tPhi );

            real tValue ;

            const Material * tFerro = aElement->slave()->material() ;

            for ( uint k = 0; k < mNumberOfIntegrationPoints; ++k )
            {

                real tOmega = tW( k ) * mGroup->work_det_J();

                tBm = inv( tMaster->dNdXi( k ) * mGroup->work_Xm()) * tMaster->dNdXi( k );


#ifdef BELFEM_FERROAIR_ENRICHED
                const Vector< real > & tNm = tMaster->phi( k );
                const Vector< real > & tNs = tSlave->phi( k );

                tBs = inv( tSlave->dNdXi( k ) * mGroup->work_Xs()) * tSlave->dNdXi( k );




                real bsx = 0 ;
                real bsy = 0 ;

                for ( uint i = 0; i < tNumNodesPerFerro ; ++i )
                {
                    bsx += tBs( 1, i ) * tAz( i );
                    bsy -= tBs( 0, i ) * tAz( i );
                }

                // magnetic resistivity
                real tNu = aElement->slave()->material()->nu_s( std::sqrt( bsx*bsx + bsy*bsy) ) ;

                // begin debug printf
                /*if( aElement->id() == 330 ||  aElement->id() == 200 )
                {

                    real hmx = 0 ;
                    real hmy = 0 ;
                    real bmx = 0 ;
                    real bmy = 0 ;

                    for ( uint j = 0; j < tNumNodesPerAir ; ++j )
                    {
                        hmx -= tBm( 0, j ) * tPhi( j );
                        hmy -= tBm( 1, j ) * tPhi( j );
                    }
                    bmx = constant::mu0 * hmx ;
                    bmy = constant::mu0 * hmy ;

                    real hsx = tNu * bsx ;
                    real hsy = tNu * bsy ;

                    std::cout << aElement->id()  << " " << k << " bx: " << bsx/ bmx << " hy : " << hsy / hmy << " | " << constant::nu0 / tNu << " | " << mGroup->work_det_J() << std::endl;
                }*/


                // lambda multiplicators
                real tPsi0 = mEdgeFunctionTS->E( k )( 0, 0 ) ;
                real tPsi1 = mEdgeFunctionTS->E( k )( 0, 1 ) ;

                // in-plane condition: h != h
                for ( uint j = 0; j < tNumNodesPerAir ; ++j )
                {
                    tValue = tOmega * ( tn( 1 ) * tBm( 0, j ) - tn( 0 ) * tBm( 1, j ) );
                    tK( p, j ) += tPsi0 * tValue ;
                    tK( j, p ) += tPsi0 * tValue ;
                    tK( q, j ) += tPsi1 * tValue ;
                    tK( j, q ) += tPsi1 * tValue ;
                }
                for( uint i = 0; i<tNumNodesPerFerro ; ++i )
                {
                    tValue = tNu * tOmega * ( tn( 0 ) * tBs( 0, i ) + tn( 1 ) * tBs( 1, i ) );
                    tK( p, i +tNumNodesPerAir  ) += tPsi0 * tValue ;
                    tK( i+tNumNodesPerAir, p ) += tPsi0 * tValue ;
                    tK( q, i +tNumNodesPerAir  ) += tPsi1 * tValue ;
                    tK( i+tNumNodesPerAir, q ) += tPsi1 * tValue ;
                }

                // normal condition : b != b
                for ( uint j = 0; j < tNumNodesPerAir ; ++j )
                {
                    tValue = tOmega * ( tn( 0 ) * tBm( 0, j ) + tn( 1 ) * tBm( 1, j ) );
                    tK( r, j ) += tPsi0 * tValue ;
                    tK( j, r ) += tPsi0 * tValue ;
                    tK( s, j ) += tPsi1 * tValue ;
                    tK( j, s ) += tPsi1 * tValue ;
                }

                for( uint i = 0; i<tNumNodesPerFerro ; ++i )
                {
                    tValue = constant::nu0 * tOmega * ( tn( 0 ) * tBs( 1, i ) - tn( 1 ) * tBs( 0, i ) );
                    tK( r, i +tNumNodesPerAir  ) += tPsi0 * tValue ;
                    tK( i+tNumNodesPerAir, r ) += tPsi0 * tValue ;
                    tK( s, i +tNumNodesPerAir  ) += tPsi1 * tValue ;
                    tK( i+tNumNodesPerAir, s ) += tPsi1 * tValue ;
                }

                // flux condition
                /*for ( uint i = 0; i < tNumNodesPerAir ; ++i)
                {
                    for( uint j = 0; j<tNumNodesPerFerro ; ++j )
                    {
                        tValue = tOmega *  tNm( i ) * ( tn( 0 ) * tBs( 1, j ) - tn( 1 ) * tBs( 0, j ) );
                        tK( i, j+tNumNodesPerAir ) += tValue ;
                    }
                }*/

                // contributions to stiffness
                for ( uint j = 0; j < tNumNodesPerAir ; ++j )
                {
                    for ( uint i = 0; i < tNumNodesPerFerro; ++i )
                    {
                        real tValue =  tOmega * tNs( i ) *
                                      ( tn( 1 ) * tBm( 0, j ) - tn( 0 ) * tBm( 1, j ));

                        tK( i + tNumNodesPerAir, j ) -= tValue ;
                        tM( j, i + tNumNodesPerAir ) += tValue ;
                    }
                }


#else
                const Vector< real > & tNs = tSlave->phi( k );

                for ( uint j = 0; j < tNumNodesPerAir ; ++j )
                {
                    for ( uint i = 0; i < tNumNodesPerFerro; ++i )
                    {
                        real tValue =  tOmega * tNs( i ) *
                                      ( tn( 1 ) * tBm( 0, j ) - tn( 0 ) * tBm( 1, j ));

                        tK( i + tNumNodesPerAir, j ) -= tValue ;
                        tM( j, i + tNumNodesPerAir ) += tValue ;
                    }
                }
#endif
            }

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
                this->collect_node_coords( aElement->master(), mGroup->work_Xm() );

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
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_symmetry_air(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // that this works is a coincidence that we shamelessly exploit
            this->compute_jacobian_and_rhs_antisymmmetry_a_tri6(
                    aElement, aJacobian, aRHS );
            aJacobian *= constant::mu0 ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_antisymmetry_air(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // that this works is a coincidence that we shamelessly exploit
            this->compute_jacobian_and_rhs_symmmetry_a_tri6(
                    aElement, aJacobian, aRHS );
            aJacobian *= constant::mu0 ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_symmetry_ferro(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
#ifdef BELFEM_FERRO_LINEAR
            this->compute_jacobian_and_rhs_symmmetry_a_tri3(
                    aElement, aJacobian, aRHS );
#else
            this->compute_jacobian_and_rhs_symmmetry_a_tri6(
                    aElement, aJacobian, aRHS );
#endif
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi_Tri6::compute_jacobian_and_rhs_antisymmetry_ferro(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
#ifdef BELFEM_FERRO_LINEAR
            this->compute_jacobian_and_rhs_antisymmmetry_a_tri3(
                    aElement, aJacobian, aRHS );
#else
            this->compute_jacobian_and_rhs_antisymmmetry_a_tri6(
                    aElement, aJacobian, aRHS );
#endif
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

            // stiffness matrix
            Matrix< real > & tA = tK;


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

            // real tSign = aElement->element()->physical_tag() == 1 ? 1.0 : -1.0 ;

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

            // indices
            uint p ;
            uint q ;

            real tValue ;

            // element length
            real tLength = 0 ;

            //real tIz = 0.0 ;

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
                //tBm = inv( mGroup->work_J() ) * tIntMaster->dNdXi( k );
                tBm = inv( tIntMaster->dNdXi( k ) * tXm ) * tIntMaster->dNdXi( k );

                // gradient operator slave
                tBs = inv( tIntSlave->dNdXi( k ) * tXs ) * tIntSlave->dNdXi( k );

                // normal component of H ( needed for material properties)
                mHn( k ) = -dot( tn.vector_data(), tBm.matrix_data() * tPhiM.vector_data() );


                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // lambda bc for master
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                {
                    tnxB( i ) = tn( 0 ) * tBm( 1, i )
                                      -tn( 1 ) * tBm( 0, i );
                }

                real tHtm = -dot( tnxB, tPhiM );

                p = mNumberOfNodesPerElement+mNumberOfNodesPerElement;
                q = p + mNumberOfThinShellLayers * 4 + 2 ;

                for( uint j=0; j<mEdgeDofMultiplicity; ++j )
                {
                    for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                    {
                        tValue = tWDetJ * mE( j ) * tnxB( i ) ;
                        tA( i, q+j ) += tValue ;
                        tA( q+j, i ) += tValue ;
                    }
                    for( uint i=0; i<mEdgeDofMultiplicity; ++i )
                    {
                        tValue = tWDetJ * mE( i ) * mE( j )  ;
                        tA( p+i, q+j ) += tValue ;
                        tA( q+j, p+i ) += tValue ;
                    }
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // lambda bc for slave
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                p = q - 2 ;
                q += 2 ;
                for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                {
                    tnxB( i ) = tn( 0 ) * tBs( 1, i )
                                     - tn( 1 ) * tBs( 0, i );
                }
                real tHts = -dot( tnxB, tPhiS );

                for( uint j=0; j<mEdgeDofMultiplicity; ++j )
                {
                    for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                    {
                        tValue = tWDetJ * mE( j ) * tnxB( i )  ;
                        tA( mNumberOfNodesPerElement+i, q+j ) += tValue ;
                        tA( q+j, mNumberOfNodesPerElement+i ) += tValue ;
                    }
                    for( uint i=0; i<mEdgeDofMultiplicity; ++i )
                    {
                        tValue = tWDetJ * mE( i ) * mE( j ) ;
                        tA( p+i, q+j ) += tValue ;
                        tA( q+j, p+i ) += tValue ;
                    }
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // DEBUG printout
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                /*Matrix< real > & tX = mGroup->work_X() ;
                this->collect_node_coords( aElement, tX );

                this->collect_edge_data_from_layer( aElement, "edge_h", 0, tHt );


                std::cout << "el " << aElement->id() << " k " << k
                    << " x= " <<mE( 0 ) * tX( 0, 0 ) + mE( 1 ) * tX( 1, 0  )
                    << " y= " <<mE( 0 ) * tX( 0,  1) + mE( 1 ) * tX( 1, 1  )
                    << " m " << tHtm << " " <<  mE( 0 ) * tHt( 0 ) + mE( 1 ) * tHt( 1 )
                    << " s " << tHts << " " <<  mE( 0 ) * tHt( 4 ) + mE( 1 ) * tHt( 5 )
                    << std::endl ; */

#ifdef HPHI_TRI6_LAMBDAN
                q = mNumberOfDofsPerElement - 2 ;

                for( uint j=0; j<mEdgeDofMultiplicity; ++j )
                {
                    for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                    {
                        tValue = tWDetJ * mE( j ) * ( tn( 0 ) * tBm( 0, i ) + tn( 1 ) * tBm( 1, i ) ) ;
                        tA( i, q+j ) += tValue ;
                        tA( q+j, i ) += tValue ;
                    }
                    for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                    {
                        tValue = tWDetJ * mE( j ) * ( tn( 0 ) * tBs( 0, i ) + tn( 1 ) * tBs( 1, i ) ) ;
                        tA( mNumberOfNodesPerElement+i, q+j ) -= tValue ;
                        tA( q+j, mNumberOfNodesPerElement+i ) -= tValue ;
                    }
                }
#else
                /*q = mNumberOfNodesPerElement ;
                const Vector< real > & tNm = tIntMaster->phi( k );
                const Vector< real > & tNs = tIntSlave->phi( k );

                // flux from slave to master
                for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                {
                    for( uint j=0; j<mNumberOfNodesPerElement; ++j )
                    {
                        tK( i, q + j ) -= tNm( i ) *
                                ( tn( 0 ) * tBs( 0, j )
                                + tn( 1 ) * tBs( 1, j ) ) * tWDetJ ;
                    }
                }

                // flux from master to slave
                for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                {
                    for( uint j=0; j<mNumberOfNodesPerElement; ++j )
                    {
                        tK( i + q, j ) += tNs( i ) *
                                            ( tn( 0 ) * tBm( 0, j )
                                            + tn( 1 ) * tBm( 1, j ) ) * tWDetJ ;
                    }
                }*/
#endif
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Derivative of dHn/dx, careful with the sign!
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            mAHn.fill( 0.0 );
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
            real tdHndX = - mBHn( 0 ) / tLength * 2.0  ; // negative, because of element orientation

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Mass and Stiffness contributions
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            for( uint l=0; l<mNumberOfThinShellLayers; ++l)
            {
                //mArho.fill( 0.0 );
                //mBrho.fill( 0.0 );

                this->collect_edge_data_from_layer( aElement, "edge_h", l, tHt );

                // initial offset for node coo
                p = 2 * mNumberOfNodesPerElement ;

                for ( uint k = 0; k < mNumberOfIntegrationPoints; ++k )
                {

                    // evaluate points for edge function
                    mE( 0 ) = mEdgeFunctionTS->E( k )( 0, 0 );
                    mE( 1 ) = mEdgeFunctionTS->E( k )( 0, 1 );

                    // gradient operator master
                    tBm = inv( tIntMaster->dNdXi( k ) * tXm ) * tIntMaster->dNdXi( k );

                    // gradient operator slave
                    tBs = inv( tIntSlave->dNdXi( k ) * tXs ) * tIntSlave->dNdXi( k );

                    Vector< real > tHm( tBm * tPhiM );
                    Vector< real > tHs( tBs * tPhiS );

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


                // critical current
                const Material * tMaterial = mGroup->thin_shell_material( l );

                this->compute_layer_stabilizer( aElement, tHt, tLength,
                                                mGroup->thin_shell_thickness( l ), tMaterial->creep_expinent_minus_1(), tMlayer, tKlayer );

                for ( uint j = 0; j < 6; ++j )
                {
                    for ( uint i = 0; i < 6; ++i )
                    {
                        tM( p + i, p + j ) += tMlayer( i, j );
                        tK( p + i, p + j ) += tKlayer( i, j );
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
                mMesh->field_data( "elementJ" )( tIndex )   = mLayerData( 0, l );
                mMesh->field_data( "elementEJ" )( tIndex )  = mLayerData( 1, l );
                mMesh->field_data( "elementJJc" )( tIndex ) = mLayerData( 2, l );
                mMesh->field_data( "elementRho" )( tIndex ) = mLayerData( 3, l );

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
#ifdef BELFEM_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

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
            real tT = 12.0 ;

            // grab layer thickness
            const real tThickness = mGroup->thin_shell_thickness( aLayer );
            real tDetJ = tThickness * 0.5 ;

            // element-wise curernt
            real tJz_el = 0 ;
            real tJJc_el = 0 ;

            // value for ac losses (element-wise)
            real tEJ_el = 0.0 ;
            real tRho_el = 0.0 ;

            const Matrix< real > & tXi = tInteg->points() ;

            // note:: usually, line elements are numbered like this
            //        0 --- 2 --- 1
            //        but here, we use the scheme
            //        0 --- 1 --- 2
            //        that's why we need to swap the indices in mHt and tPhiXi!

            real xi = tXi( 0, aIntPoint );

            real tdxi  = 0.001 ;
            real tdeta = 0.001 ;

            mXi( 0 ) = xi - tdxi ;
            mXi( 1 ) = xi + tdxi ;
            mXi( 1 ) = xi ;
            mXi( 2 ) = xi ;
            mXi( 3 ) = xi ;
            mXi( 4 ) = xi ;


            // integrate over thickness
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // derivative of shape function along thickness
                // const Vector< real > & tPhiXi  = tInteg->dphidxi( k ) ;

                real eta = tXi( 0, k );

                mEta( 0 ) = eta ;
                mEta( 1 ) = eta ;
                mEta( 2 ) = eta-tdeta  ;
                mEta( 3 ) = eta+tdeta  ;
                mEta( 4 ) = eta ;

                real tJz = 0 ;
                real tRho = 0 ;
                real tB = 0 ;
                real tAlpha = 0 ;

                real tJc= 0 ;

                // comopute rho and gradients
                for( uint l=0; l<5; ++l )
                {
                    xi = mXi( l );
                    eta = mEta( l );

                    real tE0 = 0.5 * ( 1 - 3 * xi );
                    real tE1 = 0.5 * ( 1 + 3 * xi );
                    mPhi( 0 ) = eta - 0.5;
                    mPhi( 1 ) = -2 * eta;
                    mPhi( 2 ) = eta + 0.5;

                    // curl operator
                    mC( 0, 0 ) = tE0 * mPhi( 0 );
                    mC( 0, 1 ) = tE1 * mPhi( 0 );
                    mC( 0, 2 ) = tE0 * mPhi( 1 );
                    mC( 0, 3 ) = tE1 * mPhi( 1 );
                    mC( 0, 4 ) = tE0 * mPhi( 2 );
                    mC( 0, 5 ) = tE1 * mPhi( 2 );
                    mC /= tDetJ;

                    // edge-wise components of H
                    mHt( 0 ) = tE0 * aHt( 0 ) + tE1 * aHt( 1 ) ; // field on master side
                    mHt( 1 ) = tE0 * aHt( 4 ) + tE1 * aHt( 5 ) ; // field on slave side
                    mHt( 2 ) = tE0 * aHt( 2 ) + tE1 * aHt( 3 ) ; // middle field

                    // mHt.print("Ht");
                    // current : tHt derived along thickness
                    real tHt =   dot( tInteg->phi( k ), mHt );
                    tJz = dot( mC.row( 0 ), aHt );
                    tB = std::sqrt( tHt*tHt + aHn*aHn ) ;
                    tAlpha = tB < BELFEM_EPSILON ? 0 : std::abs( std::asin( tHt / tB ) );
                    tB *= constant::mu0 ;

                    tRho = tMaterial->rho_el( tJz, tT, tB, tAlpha );
                    mSigma( l ) = tRho ;
                }

                // remember values
                mRho( k, aIntPoint ) = tRho ;
                mdRhodx( k, aIntPoint ) = ( mSigma( 0 )-mSigma( 1 ) ) / ( tdxi * aXLength )  ;
                mdRhody( k, aIntPoint ) = ( mSigma( 3)-mSigma( 2 ) ) / ( tdeta * tThickness ) ;
                mJz( k, aIntPoint ) = tJz ;

                // std::cout << "    " << k << " j " << tJz << std::endl ;

                tJc = tMaterial->j_crit( tT, tB , tAlpha );
                mJc( k, aIntPoint ) = tJc ;

                //real tdJzdx = dot( mdCdx, aHt );
                //real tdJzdy = dot( mdCdy , aHt ); // <-- is  this wrong?


                //real tJz = -adHndx + dot( tPhiXi, mHt ) / tDetJ ;
               // real tJz = dot( tPhiXi, mHt ) / tDetJ - adHndX ;

                // tangential component of h
                // real tHt =   dot( tInteg->phi( k ), mHt );

                // std::cout << "   " << k << " " << tHt << " " << tJz << std::endl ;

                // resistivity
                //real tRho = std::max( tMaterial->rho_el( tJz, tT, constant::mu0 * std::sqrt( tHt*tHt + aHn*aHn )  ), BELFEM_EPS );

                // remember values
                mRho( k, aIntPoint ) = tRho ;

                // todo: #hack
                /*if( std::abs( tJz ) > BELFEM_EPSILON )
                {
                    mdRhodx( k, aIntPoint ) = tRho * tdJzdx / tJz * 34;
                    mdRhody( k, aIntPoint ) = tRho * tdJzdy / tJz * 34;
                }
                else
                {
                    mdRhodx( k, aIntPoint ) = 0.0 ;
                    mdRhody( k, aIntPoint ) = 0.0 ;
                }*/

                /*real tX = 0.5 * ( tXi( 0, aIntPoint ) + 1.0 ) * aXLength ;
                real tY = 0.5 * ( tXi( 0, k )  + 1.0 ) * tThickness ;

                const Vector< real > & tQ = this->eval_quad9( tX, tY );
                for( uint j=0; j<9; ++j )
                {
                    for ( uint i = 0; i < 9; ++i )
                    {
                        mArho( i, j ) += tW( aIntPoint ) * tW( k ) * tQ( i ) * tQ( j );
                    }
                    mBrho( j ) += tW( aIntPoint ) * tW( k ) * std::log( tRho ) * tQ( j );
                }*/

                // contribution to stiffness matrix
                aK += tW( k ) * trans( mC ) * tRho * mC * tDetJ ;

                tJz_el += tW( k ) * tJz ;
                tEJ_el += tW( k ) *  tRho * tJz * tJz ;
                tJJc_el += tW( k ) * tJz / tJc ;
                tRho_el += tW( k ) * tRho ;
            }

            mLayerData( 0, aLayer ) += tW( aIntPoint ) * tJz_el ;
            mLayerData( 1, aLayer ) += tW( aIntPoint ) * tEJ_el ;
            mLayerData( 2, aLayer ) += tW( aIntPoint ) * tJJc_el ;
            mLayerData( 3, aLayer ) += tW( aIntPoint ) * tRho_el ;
        }

//------------------------------------------------------------------------------



        void
        IWG_Maxwell_HPhi_Tri6::compute_layer_stabilizer(
                Element * aElement,
                const Vector< real > & aHt,
                const real aXLength,
                const real aYLength,
                const real aNminus1,
                Matrix< real > & aM,
                Matrix< real > & aK )
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

            aM.fill( 0.0 );
            aK.fill( 0.0 );

            real tDetJ = 0.25 * ( aXLength * aYLength ) ;
            real tR = 0.0 ;

            // real tC = 25000 ; // for single tape
            real tC = this->tau() ;

            for( uint l=0; l<mNumberOfIntegrationPoints; ++l )
            {
                real xi = tXi( 0, l );
                real tX = 0.5 * ( xi + 1.0 ) * aXLength ;


                // evaluate the functions
                mE( 0 ) = 0.5 - 1.5 * xi ;
                mE( 1 ) = 0.5 + 1.5 * xi ;

                mdEdxi( 0 ) = - 1.5 ;
                mdEdxi( 1 ) =   1.5 ;

                for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
                {

                    const Vector< real > & tPhi  = tInteg->phi( k ) ;

                    real eta = tXi( 0, k ) ;
                    //real tY = 0.5 * ( eta + 1.0 ) * aYLength ;

                    // derivatives of rho
                    real tRho = mRho( k, l ) ;
                    real tdRhodX = mdRhodx( k, l ) ;
                    real tdRhodY = mdRhody( k, l ) ;

                    // function in thickness direction
                    mN( 0 ) = 0.5*eta*( eta - 1.0 ) ;
                    mN( 1 ) = 1.0 - eta*eta ;
                    mN( 2 ) = 0.5*eta*( eta + 1.0 ) ;

                    mdNdeta( 0 ) = eta - 0.5 ;
                    mdNdeta( 1 ) = - 2. * eta ;
                    mdNdeta( 2 ) = eta + 0.5 ;

                    md2Ndeta2( 0 ) = 1.0 ;
                    md2Ndeta2( 1 ) = -2.0 ;
                    md2Ndeta2( 2 ) = 1.0 ;

                    mdEdy( 0 ) = mdNdeta( 0 ) * mE( 0 );
                    mdEdy( 1 ) = mdNdeta( 0 ) * mE( 1 );
                    mdEdy( 2 ) = mdNdeta( 1 ) * mE( 0 );
                    mdEdy( 3 ) = mdNdeta( 1 ) * mE( 1 );
                    mdEdy( 4 ) = mdNdeta( 2 ) * mE( 0 );
                    mdEdy( 5 ) = mdNdeta( 2 ) * mE( 1 );
                    mdEdy *= 2.0 / aYLength ;

                    md2Edy2( 0 ) = md2Ndeta2( 0 ) * mE( 0 );
                    md2Edy2( 1 ) = md2Ndeta2( 0 ) * mE( 1 );
                    md2Edy2( 2 ) = md2Ndeta2( 1 ) * mE( 0 );
                    md2Edy2( 3 ) = md2Ndeta2( 1 ) * mE( 1 );
                    md2Edy2( 4 ) = md2Ndeta2( 2 ) * mE( 0 );
                    md2Edy2( 5 ) = md2Ndeta2( 2 ) * mE( 1 );
                    md2Edy2 *= 4.0 / ( aYLength * aYLength ) ;

                    md2Edxdy( 0 ) = mdNdeta( 0 ) * mdEdxi( 0 );
                    md2Edxdy( 1 ) = mdNdeta( 0 ) * mdEdxi( 1 );
                    md2Edxdy( 2 ) = mdNdeta( 1 ) * mdEdxi( 0 );
                    md2Edxdy( 3 ) = mdNdeta( 1 ) * mdEdxi( 1 );
                    md2Edxdy( 4 ) = mdNdeta( 2 ) * mdEdxi( 0 );
                    md2Edxdy( 5 ) = mdNdeta( 2 ) * mdEdxi( 1 );

                    md2Edxdy *= - 4.0 / ( aXLength * aYLength );


                    for( uint i=0; i<6; ++i )
                    {
                        mL( 0, i ) = -tdRhodY * mdEdy( i ) - tRho * md2Edy2( i );
                        mL( 1, i ) = -tdRhodX * mdEdy( i ) + tRho * md2Edxdy( i );
                    }



                    // std::cout << "delta " << tDelta << std::endl ;
                    //real tTau = aRhoCrit == 0 ? tC : tC * std::pow( tRho / aRhoCrit, 2 ) ;
                    real tTau =  tC ; // aNminus1 == 0 ? tC : tC * std::pow( mJz( k, l ) / mJc( k, k ), aNminus1 ) ;

                    aM += tW( l ) * tW( k ) * trans( mL ) * mL * tDetJ * tTau ;


                    //aM += tW( l ) * tW( k ) * trans( mL ) * mN * tDetJ * tTau ;

                    //aK += tW( l ) * tW( k ) * trans( mL ) * mL * tDetJ * tTau ;


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


            //aM *= aXLength * aYLength * constant::mu0 / max( mCrho )  * 1e10 ;


            aK.fill( 0.0 );

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
