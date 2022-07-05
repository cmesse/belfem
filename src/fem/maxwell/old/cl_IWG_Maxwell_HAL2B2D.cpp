//
// Created by christian on 10/14/21.
//

#include "cl_IWG_Maxwell_HAL2B2D.hpp"
#include "fn_trans.hpp"

namespace belfem
{

    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Maxwell_HAL2B2D::IWG_Maxwell_HAL2B2D() :
                IWG_Maxwell_Old( 2,
                             IwgType::MAXWELL_HAL2B2D,
                             IwgMode::Direct,
                             SymmetryMode::PositiveDefiniteSymmetric,
                             SideSetDofLinkMode::FacetOnly,
                             true )
        {
            mScDofs =  { "bx", "by" };
            mAirDofs =  { "bx", "by" };
            mCoilDofs =  { "bx", "by" };
            mFerroDofs =  { "bx", "by" };

            mOtherFields = { "az", "edge_h" };
            mNumberOfRhsDofsPerEdge = 1;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HAL2B2D::allocate_work_matrices( Group * aGroup )
        {
            IWG_Maxwell_Old::allocate_work_matrices( aGroup );

            // get type of group
            DomainType tType = aGroup->type() == GroupType::BLOCK ?
                               mBlockTypes( aGroup->id() ) : mSideSetTypes( aGroup->id() );

            aGroup->work_N().set_size( 2,  mNumberOfNodesPerElement, 0.0 );

            switch( tType )
            {
                case( DomainType::SuperConductor ) :
                {
                    // node shape function
                    aGroup->work_N().set_size( 2, mNumberOfDofsPerElement, 0.0 );

                    // edge function
                    aGroup->work_E().set_size( 2, mNumberOfEdgesPerElement, 0.0 );

                    // curl h
                    aGroup->work_C().set_size( 1,
                                               mNumberOfEdgesPerElement, 0.0 );

                    // edge data form current timestep
                    // change this for #facedof
                    aGroup->work_phi().set_size( mNumberOfEdgesPerElement, 0.0 );

                    break ;
                }
                case( DomainType::Air ) :
                case( DomainType::Coil ) :
                case( DomainType::FerroMagnetic ) :
                {
                    // node shape function
                    aGroup->work_N().set_size( 2, mNumberOfDofsPerElement, 0.0 );

                    // node data form last timestep
                    aGroup->work_phi().set_size( mNumberOfNodesPerElement, 0.0 );

                    // curl a
                    aGroup->work_D().set_size( 2,
                                               mNumberOfNodesPerElement, 0.0 );
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
        IWG_Maxwell_HAL2B2D::compute_jacobian(
                Element        * aElement,
                Matrix< real > & aJacobian )
        {
            this->projection_jacobian_b2d( aElement, aJacobian );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HAL2B2D::compute_rhs_superconductor(
                Element        * aElement,
                Vector< real > & aRHS )
        {
            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            Vector< real > & tH = mGroup->work_phi();
            this->collect_edge_data( aElement, "edge_h", tH );

            aRHS.fill( 0.0 );

            uint tCount = 0 ;

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // compute the new jacobian etc
                this->update_nabla( aElement, k );

                // compute geometry jacobian
                // edge function
                const Matrix< real > & tE = this->E( aElement, k );

                // compute shape
                const Vector< real > & tn = mGroup->n( k );
                Matrix< real > & tN = mGroup->work_N() ;

                tCount = 0 ;
                for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                {
                    tN( 0, tCount++ ) = tn( i );
                    tN( 1, tCount++ ) = tn( i );
                }

                aRHS += tW( k ) * trans( tN ) * ( tE * tH )
                        * this->domain_increment() ;
            }

            aRHS *= mGroup->material()->mu_r() * constant::mu0 ;
        }

//------------------------------------------------------------------------------


        void
        IWG_Maxwell_HAL2B2D::compute_rhs_air(
                Element        * aElement,
                Vector< real > & aRHS )
        {
            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            Vector< real > & tAz = mGroup->work_phi();
            this->collect_node_data( aElement, "az", tAz );

            aRHS.fill( 0.0 );

            uint tCount ;

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // compute the new jacobian etc
                this->update_nabla( aElement, k );

                // compute geometry jacobian
                const Matrix< real > & tC = this->curl_a( aElement, k );

                // compute shape
                const Vector< real > & tn = mGroup->n( k );
                Matrix< real > & tN = mGroup->work_N() ;

                tCount = 0 ;
                for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                {
                    tN( 0, tCount++ ) = tn( i );
                    tN( 1, tCount++ ) = tn( i );
                }

                aRHS += tW( k ) * trans( tN ) * ( tC * tAz ) * this->domain_increment() ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HAL2B2D::link_jacobian_function( Group * aGroup )
        {
            if( aGroup->type() == GroupType::BLOCK )
            {
                switch( mBlockTypes( aGroup->id() ) )
                {
                    case( DomainType::SuperConductor ) :
                    {
                        mFunRHS = & IWG_Maxwell_HAL2B2D::compute_rhs_superconductor ;
                        break ;
                    }
                    case( DomainType::Air ) :
                    case( DomainType::Coil ) :
                    case( DomainType::FerroMagnetic ) :
                    {
                        mFunRHS = & IWG_Maxwell_HAL2B2D::compute_rhs_air ;
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "unsupported block type");
                    }
                }
            }
            else
            {
                BELFEM_ERROR( false, "IWG_Maxwell_HAL2B2D can only be linked to blocks");
            }
        }

//------------------------------------------------------------------------------
    } /* end namespace fem */
}  /* end namespace belfem */