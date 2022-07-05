//
// Created by christian on 10/14/21.
//

#include "cl_IWG_Maxwell_HAL2J2D.hpp"
#include "fn_trans.hpp"

namespace belfem
{

    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Maxwell_HAL2J2D::IWG_Maxwell_HAL2J2D() :
                IWG_Maxwell_Old( 2,
                             IwgType::MAXWELL_HAL2J2D,
                             IwgMode::Direct,
                             SymmetryMode::PositiveDefiniteSymmetric,
                             SideSetDofLinkMode::FacetOnly,
                             true )
        {
            mScDofs = { "jz" };
            mOtherFields = { "edge_h" };
            mNumberOfRhsDofsPerEdge = 1;
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HAL2J2D::allocate_work_matrices( Group * aGroup )
        {
            IWG_Maxwell_Old::allocate_work_matrices( aGroup );

            // get type of group
            DomainType tType = aGroup->type() == GroupType::BLOCK ?
                               mBlockTypes( aGroup->id() ) : mSideSetTypes( aGroup->id() );

            aGroup->work_N().set_size( 1,  mNumberOfNodesPerElement, 0.0 );

            switch( tType )
            {
                case( DomainType::SuperConductor ) :
                {
                    // node shape function
                    aGroup->work_N().set_size( 1, mNumberOfDofsPerElement, 0.0 );

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
                default:
                {
                    BELFEM_ERROR( false, "Invalid Group type");
                }
            }
        }
//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HAL2J2D::link_jacobian_function( Group * aGroup )
        {
            BELFEM_ASSERT( aGroup->type() == GroupType::BLOCK &&
                mBlockTypes( aGroup->id() ) == DomainType::SuperConductor,
                "L2 projection for j works only on superconducting parts" );
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HAL2J2D::compute_jacobian(
                Element        * aElement,
                Matrix< real > & aJacobian )
        {
            // reset result vectors
            aJacobian.fill( 0.0 );

            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // loop over all integration points
            for ( uint k = 0; k < mNumberOfIntegrationPoints; ++k )
            {
                // compute the new jacobian etc
                this->update_nabla( aElement, k );

                // interpolation function
                const Matrix< real > & tN = mGroup->N( k );

                aJacobian += tW( k ) * trans( tN  ) * tN * this->domain_increment();
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_Maxwell_HAL2J2D::compute_rhs(
                Element        * aElement,
                Vector< real > & aRHS )
        {
            // collect edge values for H
            Vector< real > & tH = mGroup->work_phi() ;
            this->collect_edge_data( aElement, "edge_h", tH );

            aRHS.fill( 0.0 );

            const Vector< real > & tW = mGroup->integration_weights();

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // compute the new jacobian etc
                this->update_nabla( aElement, k );

                // compute curl
                const Matrix< real > & tC = this->curl_h( aElement, k );

                aRHS += tW( k ) * mGroup->n( k ) * dot( tC.row( 0 ), tH.vector_data() )
                        * this->domain_increment();
            }
        }

//------------------------------------------------------------------------------
    } /* end namespace fem */
}  /* end namespace belfem */