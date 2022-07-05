//
// Created by Christian Messe on 08.11.19.
//

#include "cl_IWG_SimpleDiffusion.hpp"
#include "cl_FEM_Group.hpp"
#include "cl_FEM_Field.hpp"
#include "cl_FEM_Element.hpp"

#include "fn_det.hpp"
#include "fn_inv.hpp"
#include "fn_trans.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_SimpleDiffusion::IWG_SimpleDiffusion( const uint aNumberOfDimensions ) :
            IWG( IwgType::SimpleDiffusion, IwgMode::Direct, SymmetryMode::PositiveDefiniteSymmetric  )
        {
            mNumberOfDofsPerNode = 1;
            mNumberOfDerivativeDimensions = aNumberOfDimensions ;
            mNumberOfRhsCols = 1;

            // set the names for the fields
            mDofFields = { "T" };
            mFluxFields = { "dotQ" };

            this->initialize() ;
        }

//------------------------------------------------------------------------------

        void
        IWG_SimpleDiffusion::compute_jacobian(
                Element        * aElement,
                Matrix <real> & aJacobian )
        {

            // reset result vectors
            aJacobian.fill( 0.0 );

            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // get the B-Matrix
            Matrix< real > & tB   = mGroup->work_dNdX();

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // compute geometry Jacobian
                Matrix< real > & tJ = aElement->J( k );

                // determinant
                real tDetJ = std::abs( det( tJ ) );

                // compute B
                tB = inv( tJ ) * mGroup->dNdXi( k );

                // add to integration
                aJacobian += tW( k ) * trans( tB ) * tB * tDetJ;
            }
        }

//------------------------------------------------------------------------------
    }
}


