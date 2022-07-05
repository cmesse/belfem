//
// Created by Christian Messe on 03.08.20.
//

#include "cl_IWG_TransientHeatConduction.hpp"

#include "commtools.hpp"
#include "cl_FEM_Field.hpp"
#include "cl_FEM_Kernel.hpp"

#include "fn_det.hpp"
#include "fn_inv.hpp"
#include "fn_trans.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_TransientHeatConduction::IWG_TransientHeatConduction(
                const uint aNumberOfDimensions,
                const IwgType aType ) :
                IWG_StationaryHeatConduction( aNumberOfDimensions, aType )
        {
            // T0 is the temperature from the last timestep
            this->add_fields( { "T0" } );

            // set default value foe euler method
            this->set_euler_method( EulerMethod::BackwardImplicit );
        }


//------------------------------------------------------------------------------
// Functions called by Field during assembly
//------------------------------------------------------------------------------

        void
        IWG_TransientHeatConduction::compute_jacobian_and_rhs(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
#ifdef BELFEM_DEBUG
            // get the number of nodes
            const uint tNumberOfNodes = mGroup->number_of_nodes_per_Element();

            // make sure that Jacobian has the right dimensions
            BELFEM_ASSERT( aJacobian.n_cols() == tNumberOfNodes &&
                          aJacobian.n_rows() == tNumberOfNodes,
                          "Jacobian has wrong dimensions ( is [ %u x %u ] but should be [ %u x %u ] )",
                          ( unsigned int ) aJacobian.n_cols(),
                          ( unsigned int ) aJacobian.n_rows(),
                          ( unsigned int ) tNumberOfNodes,
                          ( unsigned int ) tNumberOfNodes );


            // make sure that the residual has the right dimensions
            BELFEM_ASSERT( aRHS.length() == tNumberOfNodes,
                          "Resudual has wrong length ( is %u x %u )",
                          ( unsigned int ) aRHS.length(),
                          ( unsigned int ) tNumberOfNodes );
#endif
            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // nodal temperatures
            Vector< real > & tThat = mGroup->work_phi();

            // temperatures from last timestep
            Vector< real > tT0hat = mGroup->work_psi();

            // heat capacity matrix
            Matrix< real > & tC = aJacobian;

            // conductivity matrix
            Matrix< real > & tK = mGroup->work_H();

            // thermal conductivity matrix
            Matrix< real > tLambda = mGroup->work_C();

            // B-Matrix
            Matrix< real > & tB = mGroup->work_dNdX();

            // reset matrices
            tC.fill( 0.0 );
            tK.fill( 0.0 );

            this->collect_node_data( aElement, "T", tThat );
            this->collect_node_data( aElement, "T0", tT0hat );

            // get the density
            const real tRho = mMaterial->rho();

            for ( uint k = 0; k < mNumberOfIntegrationPoints; ++k )
            {
                // compute geometry Jacobian
                Matrix< real > & tJ = aElement->J( k );

                // determinant
                real tDetJ = std::abs( det( tJ ) );

                // get shape function
                const Matrix< real > & tN = mGroup->N( k );

                // derivative matrix
                tB = inv( tJ ) * mGroup->dNdXi( k );

                // interpolate temperature for this point
                real tT = this->compute_T( k );

                // contribution to heat capacity matrix
                tC += tW( k ) * trans( tN ) * tRho * mMaterial->c( tT ) * tN * tDetJ;

                // compute thermal conductivity
                mGroup->thermal_conductivity( tLambda, tT ) ;

                // contribution to conductivity matrix
                tK += tW( k ) * trans( tB ) * tLambda * tB * tDetJ;
            }

            aRHS = ( tC + mDeltaTime * ( mTheta - 1.0 ) * tK ) * tT0hat ;

            //  finalize Jacobian
            aJacobian += mDeltaTime  * mTheta * tK;

        }


//------------------------------------------------------------------------------

        void
        IWG_TransientHeatConduction::shift_fields()
        {
            if( mGroup->parent()->is_master() )
            {
                // get mesh
                Mesh * tMesh = mGroup->parent()->mesh() ;

                // shift T field
                tMesh->field_data("T0")
                        = tMesh->field_data("T");
            }

            // wait
            comm_barrier() ;

            // synchronize data over all procs
            mGroup->parent()->distribute_fields( {"T0"} );
        }

//------------------------------------------------------------------------------
// Private
//------------------------------------------------------------------------------

        void
        IWG_TransientHeatConduction::allocate_work_matrices( Group    * aGroup )
        {
            // call function from parent
            IWG_StationaryHeatConduction::allocate_work_matrices( aGroup );

            // change size of Psi Matrix
            aGroup->work_Psi().set_size( aGroup->number_of_nodes_per_element(),
                    mMesh->number_of_dimensions() );


            // Help Matrix containing K
            aGroup->work_H().set_size( mNumberOfDofsPerNode * mNumberOfNodesPerElement,
                                       mNumberOfDofsPerNode * mNumberOfNodesPerElement,
                                       0.0 );

            // needed for dotq at T0
            aGroup->work_psi().set_size( mNumberOfNodesPerElement );
        }


    } /* end namespace fem */
}  /* end namespace belfem */