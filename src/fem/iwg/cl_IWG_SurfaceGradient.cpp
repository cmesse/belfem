//
// Created by Christian Messe on 11.11.19.
//

#include "stringtools.hpp"

#include "cl_IWG_SurfaceGradient.hpp"
#include "cl_FEM_Group.hpp"
#include "cl_FEM_Field.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Dof.hpp"

#include "fn_det.hpp"
#include "fn_inv.hpp"
#include "fn_trans.hpp"
#include "fn_rotation_matrix.hpp"
#include "fn_norm.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_SurfaceGradient::IWG_SurfaceGradient() :
            IWG( IwgType::SurfaceGradient, IwgMode::Direct, SymmetryMode::PositiveDefiniteSymmetric  )
        {
            mNumberOfDofsPerNode = 1;
            mNumberOfDerivativeDimensions = 3;
            mNumberOfRhsCols = 3;

            // set the names for the fields
            mDofFields = { "grad_Tx", "grad_Ty", "grad_Tz" };

            mOtherFields = { "T" };

            this->initialize() ;
        }

//------------------------------------------------------------------------------

        void
        IWG_SurfaceGradient::compute_jacobian(
                Element        * aElement,
                Matrix <real> & aJacobian )
        {

            // reset result vectors
            aJacobian.fill( 0.0 );

            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // compute geometry Jacobian
                Matrix< real > & tJ = aElement->J( k );

                // determinant
                real tDetJ = std::abs( det( tJ ) );

                // get computen N
                Matrix< real > const & tN = mGroup->N( k );

                // add to integration
                aJacobian += tW( k ) * trans( tN ) * tN * tDetJ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_SurfaceGradient::compute_rhs(
                Element        * aElement,
                Matrix< real > & aRHS )
        {
            // get field from mesh
            // container for DOF data
            Matrix< real > & tPhiHat = mGroup->work_Phi();

            // collect element data
            this->collect_node_data( aElement, mDofFields, tPhiHat );

            // get rotation data
            Matrix< real > & tAxes = mGroup->work_Chi();
            Vector< real > & tAngles = mGroup->work_chi() ;
            aElement->get_rotation_data( tAxes, tAngles );

            aRHS.set_size( mNumberOfNodesPerElement, 3, 0.0 );

            // get the B-Matrix
            Matrix< real > & tB   = mGroup->work_dNdX();

            Matrix< real > tGrad( 3, 1 );

            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // rotation matrix
            Matrix< real > tT( 2, 3 );

            // rotation axis
            Vector< real > tR( 3 );

            // rotation angle
            real tAlpha;

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // compute geometry Jacobian
                Matrix< real > & tJ = aElement->J( k );

                // determinant
                real tDetJ = std::abs( det( tJ ) );

                // compute B
                tB = inv( tJ ) * mGroup->dNdXi( k );

                // get shape function
                const Matrix< real > & tN = mGroup->N( k );

                // compute rotation axis and angle
                tR.fill( 0.0 );
                tAlpha = 0.0;

                // interpolate rotation axis and angle
                for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                {
                    tR( 0 ) += tN( 0, i ) * tAxes( 0, i );
                    tR( 1 ) += tN( 0, i ) * tAxes( 1, i );
                    tR( 2 ) += tN( 0, i ) * tAxes( 2, i );
                    tAlpha  -= tN( 0, i ) * tAngles( i );
                }
                tR /= norm( tR );

                // compute rotation matrix
                rotation_matrix_strip( tR, tAlpha, tT );

                // compute gradient in global coordinate system
                tGrad        = trans( tT ) * tB * tPhiHat;

                // ( 3 x 1 ) =      ( 3 x 2 )      * ( 2 * n ) * ( n )

                // add to integration
                aRHS += tW( k ) * trans( tN ) * trans( tGrad ) * tDetJ;
            }

        }

//------------------------------------------------------------------------------

        void
        IWG_SurfaceGradient::allocate_work_matrices( Group * aGroup )
        {
            IWG::allocate_work_matrices( aGroup );

            uint tNumNodes = mesh::number_of_nodes( aGroup->element_type() );

            // work array for rotation axes
            aGroup->work_Chi().set_size( 3, tNumNodes );

            // work array for rotation angles
            aGroup->work_chi().set_size( tNumNodes );

        }

//------------------------------------------------------------------------------
    } /* end namespace fem */
} /* end namespace belfem */


