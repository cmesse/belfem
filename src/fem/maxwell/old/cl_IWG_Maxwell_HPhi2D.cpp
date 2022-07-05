//
// Created by christian on 10/18/21.
//

#include "commtools.hpp"
#include "constants.hpp"
#include "cl_IWG.hpp"
#include "cl_IWG_Maxwell_HPhi2D.hpp"
#include "fn_det.hpp"
#include "fn_trans.hpp"
#include "fn_dot.hpp"
#include "fn_inv.hpp"

#include "cl_FEM_DofManagerBase.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_MaxwellMaterial.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Maxwell_HPhi2D::IWG_Maxwell_HPhi2D() :
                IWG_Maxwell_Old( 2,
                             IwgType::MAXWELL_HPHI2D,
                             IwgMode::Iterative,
                             SymmetryMode::Unsymmetric,
                             SideSetDofLinkMode::MasterAndSlave,
                             true )
        {
            mScDofs = { "edge_h" };
            mAirDofs = { "phi" };
            mCoilDofs = { "phi" };
            mFerroDofs = { "az" };
            mFluxFields = { "jz" };
            mInterfaceDofs = { "lambda_t" };
            mFarfieldDofs = {"lambda"};
            mCutDofs = {"lambda"};
            mBfields = { "bx", "by" };
            mBoundaryDofs = {"lambda"};

            mOtherFields = { "a0z", "phi0", "edge_h0", "bx", "by", "lambda_I" };

            if( mTheta != 1.0 )
            {
                mOtherFields.push("lambda0");
                mOtherFields.push("lambda_I0");
                mOtherFields.push("jz0");
            }

            mHiddenFields = { "a0z", "phi0", "edge_h0", "lambda_I" };

            mNumberOfRhsDofsPerEdge = 1;
        }

//------------------------------------------------------------------------------

        IWG_Maxwell_HPhi2D::~IWG_Maxwell_HPhi2D()
        {
            //this->delete_boundary_conditions();
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi2D::shift_fields()
        {
            mMesh->field_data("a0z") = mMesh->field_data("az");
            mMesh->field_data("phi0") = mMesh->field_data("phi");
            mMesh->field_data("edge_h0")  = mMesh->field_data("edge_h");

            if( mTheta != 1.0 )
            {
                mMesh->field_data("jz0") = mMesh->field_data("jz");
                mMesh->field_data("lambda0") = mMesh->field_data("lambda");
                mMesh->field_data("lambda_I0") = mMesh->field_data("lambda_I");
            }

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi2D::allocate_work_matrices( Group * aGroup )
        {
            IWG_Maxwell_Old::allocate_work_matrices( aGroup );

            // get type of group
            DomainType tType = aGroup->type() == GroupType::BLOCK ?
                               mBlockTypes( aGroup->id() ) : mSideSetTypes( aGroup->id() );


            switch( tType )
            {
                case( DomainType::SuperConductor ) :
                {
                    // edge matrix
                    aGroup->work_E().set_size( 2,
                                               mNumberOfDofsPerElement, 0.0 );

                    // curl h
                    aGroup->work_C().set_size( 1,
                                               mNumberOfDofsPerElement, 0.0 );

                    // edge data form current timestep
                    aGroup->work_phi().set_size( mNumberOfDofsPerElement, 0.0 );

                    // edge data form last timestep
                    aGroup->work_psi().set_size( mNumberOfDofsPerElement, 0.0 );

                    // stiffness matrix
                    aGroup->work_K().set_size( mNumberOfDofsPerElement, mNumberOfDofsPerElement, 0.0 );

                    break ;

                }
                case( DomainType::Air ) :
                case( DomainType::Coil ) :
                case( DomainType::FerroMagnetic ) :
                {
                    // node data form last timestep
                    aGroup->work_phi().set_size( mNumberOfDofsPerElement, 0.0 );

                    // current data
                    aGroup->work_chi().set_size( mNumberOfNodesPerElement, 0.0 );

                    aGroup->work_B().set_size( 2, mNumberOfNodesPerElement );

                    // curl a
                    //aGroup->work_D().set_size( 2,
                    //                           mNumberOfNodesPerElement, 0.0 );
                    break ;
                }
                case( DomainType::InterfaceScAir ) :
                {
                    aGroup->work_B().set_size( 2, mNumberOfDofsPerElement, 0.0 );
                    aGroup->work_E().set_size( 2, mNumberOfDofsPerElement, 0.0 );

                    aGroup->work_phi().set_size( mNumberOfDofsPerElement, 0.0 );
                    aGroup->work_chi().set_size( mNumberOfEdgesPerElement, 0.0 );
                    aGroup->work_sigma().set_size( mNumberOfNodesPerElement );
                    aGroup->work_C().set_size( 1, mNumberOfNodesPerElement );

                    break ;
                }
                case( DomainType::InterfaceFmAir ) :
                {

                    aGroup->work_N().set_size( 1, mNumberOfDofsPerElement, 0.0 );
                    aGroup->work_B().set_size( 2, mNumberOfDofsPerElement, 0.0 );
                    aGroup->work_Psi().set_size(2, mNumberOfDofsPerElement, 0.0 );

                    aGroup->work_sigma().set_size( mNumberOfNodesPerElement );
                    aGroup->work_tau().set_size( mNumberOfNodesPerElement );
                    aGroup->work_psi().set_size(  mNumberOfDofsPerElement );

                    break ;
                }
                case( DomainType::Boundary ) :
                {
                    aGroup->work_phi().set_size( mNumberOfNodesPerElement );

                    aGroup->work_sigma().set_size( 1, mNumberOfNodesPerElement );

                    // for grad
                    aGroup->work_B().set_size( 2, mNumberOfDofsPerElement, 0.0 );

                    aGroup->work_chi().set_size( mNumberOfNodesPerElement, 0.0 );

                    aGroup->work_Chi().set_size( mNumberOfNodesPerElement, 2, 0.0 );

                    break ;
                }
                case( DomainType::Cut ) :
                {
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Invalid Group type");
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_superconductor(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // data from last timestep
            Vector< real > & tH0 = mGroup->work_psi() ;
            this->collect_edge_data( aElement, "edge_h0", tH0 );

            // data from current timestep
            Vector< real > & tH = mGroup->work_phi() ;
            this->collect_edge_data( aElement, "edge_h", tH );

            // mass matrix
            Matrix< real > & tM = aJacobian ;

            // stiffnes matrix
            Matrix< real > & tK = mGroup->work_K();

            // compute new balance if this is not a backwards Euler
            if( this->theta() != 1.0 )
            {
                tH *= this->theta();
                tH += ( 1.0 - this->theta()) * tH0;

                // call compute function for edge elements
                this->compute_sc_matrices( aElement, tM, tK );

                // compute residual vector
                aRHS = ( tM + mDeltaTime * ( mTheta-1.0 ) * tK ) * tH0;

                //  finalize Jacobian ( M + delta_t * theta * K )
                aJacobian += mDeltaTime * mTheta * tK;
            }
            else
            {
                // call compute function for edge elements
                this->compute_sc_matrices( aElement, tM, tK );

                // compute residual vector
                aRHS = tM * tH0;

                //  finalize Jacobian ( M + delta_t * theta * K )
                aJacobian += mDeltaTime * tK;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_air(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {

            // data from last timestep
            Vector< real > & tPhi0 = mGroup->work_psi() ;
            this->collect_node_data( aElement, "phi0", tPhi0 );

            // stiffmess matrix
            Matrix< real > & tM = aJacobian ;

            tM.fill( 0.0 );

            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights() ;

            // derivative operator
            Matrix< real > & tB = mGroup->work_B();

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // update nabla functions
                this->update_nabla( aElement, k );

                // compute derivative operator
                tB = mGroup->work_invJ() * mGroup->dNdXi( k );

                // compute mass matrix
                tM += tW( k ) * trans( tB ) * tB * mDomainIncrement ;
            }

            tM *= constant::mu0 ;
            aRHS = tM * tPhi0;

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_ferro(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            this->ferro_stiffness_2d( aElement, aJacobian );

            if( this->theta() != 1.0 )
            {
                // finish Jacobian step
                aRHS = ( mDeltaTime * ( mTheta -1.0 ) * aJacobian ) * mGroup->work_psi() ;

                //  finalize Jacobian
                aJacobian *= mDeltaTime * mTheta;
            }
            else
            {
                aRHS.fill( 0.0 );
                aJacobian *= mDeltaTime ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_interface_hphi(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            this->interface_hphi_2d( aElement, aJacobian, aRHS );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_interface_aphi(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Step 1: compute Jacobian etc
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            this->update_nabla( aElement->master(), 0 );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Step 2: collect dof data from last timestep
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // node dofs
            Vector< real > & tAz0   = mGroup->work_tau() ;
            Vector< real > & tPhi0  = mGroup->work_sigma() ;

            // dof vector
            Vector< real > & tQ0    = mGroup->work_psi() ;

            this->collect_node_data( aElement->master(), "az", tAz0 );
            this->collect_node_data( aElement->slave(), "phi0", tPhi0 );

            Matrix< real > & tM = aJacobian ;
            Matrix< real > & tK = mGroup->work_K() ;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Step 3: assemble dof vector
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            uint tCount = 0 ;
            for( uint i=0; i<mNumberOfNodesPerElement; ++i )
            {
                tQ0( tCount++ ) = tAz0( i );
            }
            for( uint i=0; i<mNumberOfNodesPerElement; ++i )
            {
                tQ0( tCount++ ) = tPhi0( i );
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Step 4: compute matrices
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            this->interface_phia_2d( aElement, tM, tK );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Step 4: compute timesetp
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            if( this->theta() != 1.0 )
            {
                // compute residual vector
                aRHS = ( tM + mDeltaTime * ( mTheta-1.0 ) * tK ) * tQ0;

                //  finalize Jacobian ( M + delta_t * theta * K )
                aJacobian += mDeltaTime * mTheta * tK;
            }
            else
            {
                // compute residual vector
                aRHS = tM * tQ0;

                //  finalize Jacobian ( M + delta_t * theta * K )
                aJacobian += mDeltaTime * tK;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_interface_ha(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Step 1: compute Jacobian etc
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            this->update_nabla( aElement->master(), 0 );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Step 2: collect dof data from last timestep
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // edge dofs
            Vector< real > & tH  = mGroup->work_chi() ;
            this->collect_edge_data( aElement->master(), "edge_h0", tH );

            // node dofs
            Vector< real > & tAz = mGroup->work_sigma() ;
            this->collect_node_data( aElement->slave(), "a0z", tAz );

            // assemble vector
            Vector< real > & tQ0 = mGroup->work_phi();
            uint tCount = 0 ;
            for( real h : tH )
            {
                tQ0( tCount++ ) = h ;
            }
            for( real az : tAz )
            {
                tQ0( tCount++ ) = az ;
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Step 3: get work matrices
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // mass matrix
            Matrix< real > & tM = aJacobian ;

            // stiffness matrix
            Matrix< real > & tK = mGroup->work_K();

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Step 4: compute matrices
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            this->interface_ha_2d( aElement, tM, tK );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Step 5: perform timestepping matrices
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            if( mTheta == 1.0 )
            {
                // compute right hand side vector
                aRHS = tM * tQ0 ;

                //  finalize Jacobian ( M + delta_t * theta * K )
                aJacobian += mDeltaTime * tK ;
            }
            else
            {
                // compute right hand side vector
                aRHS = ( tM + mDeltaTime * ( mTheta -1.0 ) * tK ) * tQ0 ;

                //  finalize Jacobian ( M + delta_t * theta * K )
                aJacobian += mDeltaTime  * mTheta * tK;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_cut(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {


            aJacobian( 0, 0 ) =  0.0 ;
            aJacobian( 1, 0 ) =  0.0 ;
            aJacobian( 2, 0 ) =  1.0 ;
            aJacobian( 0, 1 ) =  0.0 ;
            aJacobian( 1, 1 ) =  0.0 ;
            aJacobian( 2, 1 ) = -1.0 ;
            aJacobian( 0, 2 ) =  1.0 ;
            aJacobian( 1, 2 ) = -1.0 ;
            aJacobian( 2, 2 ) =  0.0 ;


            aRHS( 0 ) = 0.0 ;
            aRHS( 1 ) = 0.0 ;
            this->collect_lambda_data( aElement, "lambda_I", aRHS( 2 ) );

            if( mTheta != 1.0 )
            {
                BELFEM_ERROR( false, "fix this");

                real tI0 ;
                this->collect_lambda_data( aElement, "lambda_I0", tI0 );

                aRHS( 2 ) *= mTheta ;
                aRHS( 2 ) += ( mTheta - 1.0 ) * tI0  ;

                // compute residual vector

                // aRHS +=  ( mTheta-1.0 ) * aJacobian ) * a;

                aRHS*= mDeltaTime ;

                //  finalize Jacobian ( M + delta_t * theta * K )
                aJacobian *= mDeltaTime * mTheta ;
            }
            else
            {
                // compute residual vector
                aRHS *= mDeltaTime ;

                //  finalize Jacobian ( M + delta_t * theta * K )
                aJacobian *= mDeltaTime ;
            }

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi2D::link_jacobian_function( Group * aGroup )
        {
            if( aGroup->type() == GroupType::BLOCK )
            {
                switch( mBlockTypes( aGroup->id() ) )
                {
                    case( DomainType::SuperConductor ) :
                    {
                        mFunJacobian = & IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_superconductor ;
                        break ;
                    }
                    case( DomainType::FerroMagnetic ) :
                    {
                        mFunJacobian = & IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_ferro ;
                        break ;
                    }
                    case( DomainType::Air ) :
                    case( DomainType::Coil ) :
                    {
                        mFunJacobian = & IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_air ;
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "unsupported block type for block %lu: %s",
                                     ( long unsigned int ) aGroup->id(),
                                     to_string( mBlockTypes( aGroup->id() ) ).c_str() );
                    }
                }
            }
            else
            {
                switch( mSideSetTypes( aGroup->id() ) )
                {
                    case( DomainType::InterfaceScAir ) :
                    {
                        mFunJacobian = & IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_interface_hphi ;
                        break ;
                    }
                    case( DomainType::InterfaceFmAir ) :
                    {
                        mFunJacobian = & IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_interface_aphi ;
                        break ;
                    }
                    case( DomainType::InterfaceScFm ) :
                    {
                        mFunJacobian = & IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_interface_ha ;
                        break ;
                    }
                    case( DomainType::Cut ) :
                    {
                        mFunJacobian =  & IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_cut ;
                        break ;
                    }
                    case( DomainType::Boundary ) :
                    {
                        switch( mMagfieldTypeMap( aGroup->id() ) )
                        {
                            case( MagfieldBcType::Farfied ) :
                            {
                                mFunJacobian = & IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_farfield ;
                                break ;
                            }
                            case( MagfieldBcType::Symmetry ) :
                            {
                                mFunJacobian = & IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_symmetry ;
                                break ;
                            }
                            case( MagfieldBcType::Wave ) :
                            {
                                mFunJacobian = & IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_boundary ;
                                break ;
                            }
                            default :
                            {
                                BELFEM_ERROR( false, "Invalid Boundary Type for sideset %lu",
                                             ( long unsigned int ) aGroup->id() );
                            }
                        }

                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "unsupported sideset type for sideset %lu: %s",
                                     ( long unsigned int ) aGroup->id(),
                                     to_string( mSideSetTypes( aGroup->id() ) ).c_str());
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_interface(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            BELFEM_ERROR( false, "This function is not implemented");
        }

//-----------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_farfield(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            BELFEM_ASSERT(
                    aElement->master()->element()->type() == ElementType::TRI3,
                    "Element must be a linear triangle" );

            // pass
            aRHS.fill( 0.0 );
            aJacobian.fill( 0.0 );

            // integration points on master element ( air )
            /*Matrix< real > & tXi = mGroup->work_Tau() ;

            this->update_nabla( aElement->master(), 0.0 );

            Matrix< real > & tB = mGroup->work_B() ;
            this->compute_B_2d( 0, tB );

            this->project_intpoints( aElement, tXi, tXi );

            const Vector< real > & tn = this->normal( aElement, 0 );
            const real nx = tn( 0 );
            const real ny = tn( 1 );

            const Vector< real > & tW = mGroup->integration_weights() ;

            for( uint k=0; k<tW.length(); ++k )
            {
                for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                {
                    aJacobian( i, mNumberOfNodesPerElement ) += tW( k ) * ( nx*tB( 1, i ) -ny*tB( 0, i ) ) * this->domain_increment() ;
                }
            }
            aJacobian.set_row( mNumberOfNodesPerElement, aJacobian.col( mNumberOfNodesPerElement )) ; */
        }


//-----------------------------------------------------------------------------


        void
        IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_boundary(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            BELFEM_ERROR( false, "this sideset is blacklisted");

            BELFEM_ASSERT(
                    aElement->master()->element()->type() == ElementType::TRI3,
                    "Element must be a linear triangle" );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HPhi2D::compute_jacobian_and_rhs_symmetry(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            BELFEM_ERROR( false, "This function is not implemented");
        }

//------------------------------------------------------------------------------
    } /* end namespace fem */
}  /* end namespace belfem */