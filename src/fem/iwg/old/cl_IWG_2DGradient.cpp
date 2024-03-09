//
// Created by Christian Messe on 16.11.19.
//

#include "cl_IWG_2DGradient.hpp"
#include "cl_FEM_Group.hpp"
#include "cl_FEM_Field.hpp"
#include "cl_FEM_Element.hpp"

#include "fn_det.hpp"
#include "fn_inv.hpp"
#include "fn_trans.hpp"

#include "cl_Cell.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_2DGradient::IWG_2DGradient() :
            IWG( IwgType::Gradient2D, IwgMode::Direct, SymmetryMode::PositiveDefiniteSymmetric )
        {
            mNumberOfDofsPerNode = 2;
            mNumberOfDerivativeDimensions = 2;
            mNumberOfRhsCols = 2;

            mDofFields = { "grad_uxx", "grad_uxy" };
            mOtherFields = { "ux", "uy" };
            mTensorFields = { "grad_uxx", "grad_uxy", "grad_uyx", "grad_uyy" };
            this->initialize() ;
        }

//------------------------------------------------------------------------------

        void
        IWG_2DGradient::compute_jacobian(
                Element        * aElement,
                Matrix <real> & aJacobian )
        {

            // reset result vectors
            aJacobian.fill( 0.0 );

            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // get the N-Matrix
            Matrix< real > & tN      = mGroup->work_N();
            Matrix< real > & tB      = mGroup->work_dNdX();

            // get number of nodes from this element
            uint tNumNodes = mGroup->number_of_nodes_per_element() ;

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // compute geometry Jacobian
                Matrix< real > & tJ = aElement->J( k );

                const Matrix< real > & tNk = mGroup->N( k );

                // determinant
                real tDetJ = std::abs( det( tJ ) );

                // compute B
                tB = inv( tJ ) * mGroup->dNdXi( k );

                // populate N
                uint tCount = 0;
                for( uint i=0; i<tNumNodes; ++i )
                {
                    tN( 0, tCount ) = tNk( 0, i );
                    ++tCount;
                    tN( 1, tCount ) = tNk( 0, i );
                    ++tCount;
                }

                // add to integration
                aJacobian += tW( k ) * trans( tN ) * tN * tDetJ;
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_2DGradient::compute_rhs(
                Element        * aElement,
                Matrix< real > & aRHS )
        {
            aRHS.fill( 0.0 );

            // get number of nodes from this element
            uint tNumNodes = mGroup->number_of_nodes_per_element() ;

            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // get the N-Matrix
            Matrix< real > & tN   = mGroup->work_N();

            // loop over all integration points
            tN.fill( 0.0 );

            // get the B-Matrix
            Matrix< real > & tB = mGroup->work_dNdX();


            // container for DOF data
            Matrix< real > & tPhiHat = mGroup->work_Phi();

            // collect element data
            this->collect_node_data( aElement, mOtherFields, tPhiHat );

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // compute geometry Jacobian
                Matrix< real > & tJ = aElement->J( k );

                const Matrix< real > & tNk = mGroup->N( k );

                // determinant
                real tDetJ = std::abs( det( tJ ) );

                tB = inv( tJ ) * mGroup->dNdXi( k );

                // populate N
                uint tCount = 0;
                for( uint i=0; i<tNumNodes; ++i )
                {
                    tN( 0, tCount ) = tNk( 0, i );
                    ++tCount;
                    tN( 1, tCount ) = tNk( 0, i );
                    ++tCount;
                }

                // add to integration
                aRHS += tW( k ) * trans( tN ) * tB * tDetJ * tPhiHat ;
            }
        }

//------------------------------------------------------------------------------
    }
}
