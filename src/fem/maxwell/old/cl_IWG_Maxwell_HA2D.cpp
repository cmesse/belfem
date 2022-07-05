//
// Created by christian on 10/4/21.
//

#include "constants.hpp"
#include "cl_IWG.hpp"
#include "cl_IWG_Maxwell_HA2D.hpp"
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

        IWG_Maxwell_HA2D::IWG_Maxwell_HA2D() :
                IWG_Maxwell_Old( 2,
                             IwgType::MAXWELL_HA2D,
                             IwgMode::Iterative,
                             SymmetryMode::Unsymmetric,
                             SideSetDofLinkMode::MasterAndSlave,
                             true )
        {
            mScDofs = { "edge_h" };
            mAirDofs = { "az" };
            mCoilDofs = { "az" };
            mFerroDofs = { "az" };
            mFluxFields = { "jz" };
            mFarfieldDofs = {"lambda_n"};
            mSymmetryDofs = {"lambda_p"};
            // mBoundaryDofs = { "az" };

            mOtherFields = { "a0z", "edge_h0", "bx", "by" };

            mHiddenFields = { "a0z", "edge_h0" };

            mNumberOfRhsDofsPerEdge = 1;

            mBfields = {"bx", "by"};
        }

//------------------------------------------------------------------------------

        IWG_Maxwell_HA2D::~IWG_Maxwell_HA2D()
        {
            //this->delete_boundary_conditions();
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HA2D::shift_fields()
        {
            mMesh->field_data("a0z") = mMesh->field_data("az");
            mMesh->field_data("edge_h0")   = mMesh->field_data("edge_h");
        }

            //------------------------------------------------------------------------------

        void
        IWG_Maxwell_HA2D::assign_dofs_per_sideset( const Vector< id_t > & aSideSetIDs )
        {
            // assume that all dofs sit on all sidesets
            mDofsPerSideSet.set_size( aSideSetIDs.length(), mDefaultDofTypes );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HA2D::allocate_work_matrices( Group * aGroup )
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

                    // curl a
                    aGroup->work_D().set_size( 2,
                                               mNumberOfNodesPerElement, 0.0 );
                    break ;
                }
                case( DomainType::InterfaceScAir ) :
                {
                    aGroup->work_N().set_size( 1,  mNumberOfDofsPerElement, 0.0 );
                    aGroup->work_E().set_size( 2, mNumberOfDofsPerElement, 0.0 );
                    aGroup->work_Phi().set_size( 1, mNumberOfDofsPerElement, 0.0 );
                    aGroup->work_Psi().set_size( 2, mNumberOfDofsPerElement, 0.0 );
                    break ;
                }
                case( DomainType::Boundary ) :
                {
                    // aGroup->work_N().set_size( 1, mNumberOfNodesPerElement, 0.0 );

                    aGroup->work_sigma().set_size( 1, mNumberOfNodesPerElement );

                    // for curl
                    aGroup->work_D().set_size( 2, mNumberOfDofsPerElement, 0.0 );

                    aGroup->work_chi().set_size( mNumberOfNodesPerElement, 0.0 );

                    aGroup->work_Chi().set_size( mNumberOfNodesPerElement, 2, 0.0 );
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Invalid Group type for group %lu", (long unsigned int ) aGroup->id() );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HA2D::compute_jacobian_and_rhs_superconductor(
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
        IWG_Maxwell_HA2D::compute_jacobian_and_rhs_air(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            Matrix< real > & tK      = aJacobian ;

            // reset result vectors
            tK.fill( 0.0 );

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // compute the new jacobian etc
                this->update_nabla( aElement, k );

                // compute curl operator
                const Matrix< real > & tC = this->curl_a( aElement, k ) ;

                tK += tW( k ) * trans( tC ) * tC * this->domain_increment() ;
            }

            tK *= constant::nu0 ;

            if( this->theta() != 1.0 )
            {
                // old potential
                Vector< real > & tQ0     = mGroup->work_psi() ;

                this->collect_node_data( aElement, "a0z", tQ0 );

                // finish Jacobian step
                aRHS = ( mDeltaTime * ( mTheta -1.0 ) * tK ) * tQ0;

                //  finalize Jacobian
                aJacobian *= mDeltaTime * mTheta;
            }
            else
            {
                // finish Jacobian step
                aRHS.fill( 0.0 );

                //  finalize Jacobian
                aJacobian *= mDeltaTime ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HA2D::compute_jacobian_and_rhs_coil(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            Matrix< real > & tK      = aJacobian ;

            // currents
            Vector< real > & tJz     = mGroup->work_chi() ;

            // reset result vectors
            tK.fill( 0.0 );

            // initialize right hand side
            aRHS.fill( 0.0 );


            this->collect_node_data( aElement, "jz", tJz );

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // compute the new jacobian etc
                this->update_nabla( aElement, k );

                // compute curl operator
                const Matrix< real > & tC = this->curl_a( aElement, k ) ;

                tK += tW( k ) * trans( tC ) * tC * this->domain_increment() ;

                aRHS += tW( k ) * mGroup->n( k ) * dot( mGroup->n( k ), tJz ) * this->domain_increment() ;
            }

            tK *= constant::nu0 ;

            aRHS *= mDeltaTime ;


            if( this->theta() != 1.0 )
            {
                // old potential
                Vector< real > & tQ0     = mGroup->work_psi() ;

                this->collect_node_data( aElement, "a0z", tQ0 );

                // finish Jacobian step
                aRHS += ( mDeltaTime * ( mTheta -1.0 ) * tK ) * tQ0;

                //  finalize Jacobian
                aJacobian *= mDeltaTime * mTheta;
            }
            else
            {

                //  finalize Jacobian
                aJacobian *= mDeltaTime ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HA2D::compute_jacobian_and_rhs_ferro(
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
        IWG_Maxwell_HA2D::compute_jacobian_and_rhs_interface(
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
        IWG_Maxwell_HA2D::compute_jacobian_and_rhs_boundary(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // we can do this  only once since the element is a linear triangle
            BELFEM_ASSERT(
                    aElement->master()->element()->type() == ElementType::TRI3,
                    "Element must be a linear triangle" );

            this->update_nabla( aElement->master(), 0 );

            // integration points on master element ( air )
            Matrix< real > & tXi = mGroup->work_Tau() ;

            // project the integration points on the edge
            this->project_intpoints( aElement, tXi, tXi );

            // const Matrix< real > & tC = this->curl_a( aElement->master(), 0 );
            const Vector< real > & tW = mGroup->integration_weights() ;
            Matrix< real > & tB = mGroup->work_Chi();
            Vector< real > & tN = mGroup->work_chi();

            // dot product
            aJacobian.fill( 0.0 );
            aRHS.fill( 0.0 );

            real tBx ;
            real tBy ;

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                tN( 0 ) = tXi( 0, k );
                tN( 1 ) = tXi( 1, k );
                tN( 2 ) = 1-  tXi( 0, k ) -  tXi( 1, k );

                const Vector< real > & tn = this->normal( aElement, k );

                // interpolate b-field
                tBx = dot( tN, tB.col( 0 ) );
                tBy = dot( tN, tB.col( 1 ) );

                aRHS += tW( k ) *
                        tN * ( tn( 1 ) * tBx - tn( 0 ) * tBy )
                        * this->domain_increment() ;
            }

            aRHS*= mDeltaTime * constant::nu0 ;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HA2D::compute_jacobian_and_rhs_farfield(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // dot product
            aJacobian.fill( 0.0 );
            aRHS.fill( 0.0 );

            // we can do this  only once since the element is a linear triangle
            BELFEM_ASSERT(
                    aElement->master()->element()->type() == ElementType::TRI3,
                    "Element must be a linear triangle" );


            this->update_nabla( aElement->master(), 0 );

            // integration points on master element ( air )
            Matrix< real > & tXi = mGroup->work_Tau() ;

            // project the integration points on the edge
            this->project_intpoints( aElement, tXi, tXi );

            const Matrix< real > & tC = this->curl_a( aElement->master(), 0 );
            const Vector< real > & tn = this->normal( aElement->master(), 0.0 );

            aJacobian.fill( 0.0 );
            aRHS.fill( 0.0 );

            for( uint i=0; i<mNumberOfNodesPerElement; ++i )
            {
                aJacobian( mNumberOfNodesPerElement, i ) = tn( 0 ) * tC( 0, i ) + tn( 1 ) * tC( 1, i ) ;
                aJacobian( i, mNumberOfNodesPerElement ) = aJacobian( mNumberOfNodesPerElement, i );
            }


            aJacobian *= 2.0 * this->domain_increment() * mDeltaTime * constant::nu0 ;

        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HA2D::compute_jacobian_and_rhs_symmetry(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {

            // dot product
            aJacobian.fill( 0.0 );
            aRHS.fill( 0.0 );

            /*
            // we can do this  only once since the element is a linear triangle
            BELFEM_ASSERT(
                    aElement->master()->element()->type() == ElementType::TRI3,
                    "Element must be a linear triangle" );

            this->update_nabla( aElement->master(), 0 );

            // integration points on master element ( air )
            Matrix< real > & tXi = mGroup->work_Tau() ;

            // project the integration points on the edge
            this->project_intpoints( aElement, tXi, tXi );

            const Matrix< real > & tC = this->curl_a( aElement->master(), 0 );
            const Vector< real > & tn = this->normal( aElement->master(), 0.0 );



            for( uint i=0; i<mNumberOfNodesPerElement; ++i )
            {
                aJacobian( mNumberOfNodesPerElement, i ) =
                        tn( 0 ) * tC( 1, i ) - tn( 1 ) * tC( 0, i ) ;
                aJacobian( i, mNumberOfNodesPerElement ) = aJacobian( mNumberOfNodesPerElement, i );
            }


            aJacobian *= 2.0 * this->domain_increment() * mDeltaTime * constant::nu0 ; */
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HA2D::link_jacobian_function( Group * aGroup )
        {
            if( aGroup->type() == GroupType::BLOCK )
            {
                switch( mBlockTypes( aGroup->id() ) )
                {
                    case( DomainType::SuperConductor ) :
                    {
                        mFunJacobian = & IWG_Maxwell_HA2D::compute_jacobian_and_rhs_superconductor ;
                        break ;
                    }
                    case( DomainType::Air ) :
                    {
                        mFunJacobian = & IWG_Maxwell_HA2D::compute_jacobian_and_rhs_air ;
                        break ;
                    }
                    case( DomainType::Coil ) :
                    {
                        mFunJacobian = & IWG_Maxwell_HA2D::compute_jacobian_and_rhs_coil ;
                        break ;
                    }
                    case( DomainType::FerroMagnetic ) :
                    {
                        mFunJacobian = & IWG_Maxwell_HA2D::compute_jacobian_and_rhs_ferro ;
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
                        mFunJacobian = & IWG_Maxwell_HA2D::compute_jacobian_and_rhs_interface ;
                        break ;
                    }
                    case( DomainType::Boundary ) :
                    {
                        switch( mMagfieldTypeMap( aGroup->id() ) )
                        {
                            case( MagfieldBcType::Farfied ) :
                            {
                                mFunJacobian = & IWG_Maxwell_HA2D::compute_jacobian_and_rhs_farfield ;
                                break ;
                            }
                            case( MagfieldBcType::Symmetry ) :
                            {
                                mFunJacobian = & IWG_Maxwell_HA2D::compute_jacobian_and_rhs_symmetry ;
                                break ;
                            }
                            case( MagfieldBcType::Wave ) :
                            {
                                mFunJacobian = & IWG_Maxwell_HA2D::compute_jacobian_and_rhs_boundary ;
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
                        BELFEM_ERROR( false, "unsupported sideset type for group %lu",
                            ( long unsigned int )  aGroup->id() );
                    }
                }
            }
        }

//------------------------------------------------------------------------------
    } /* end namespace fem */
}  /* end namespace belfem */