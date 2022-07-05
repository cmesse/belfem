//
// Created by Christian Messe on 15.11.19.
//

#include "cl_IWG_PlaneStress.hpp"
#include "cl_FEM_Group.hpp"
#include "cl_FEM_Field.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_Material.hpp"

#include "fn_det.hpp"
#include "fn_inv.hpp"
#include "fn_trans.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_PlaneStress::IWG_PlaneStress() :
            IWG( IwgType::PlaneStress, IwgMode::Direct, SymmetryMode::PositiveDefiniteSymmetric  )
        {
            mNumberOfDofsPerNode = 2;
            mNumberOfDerivativeDimensions = 3;
            mNumberOfRhsCols = 1;

            // set the names for the fields
            mDofFields   = { "ux", "uy" };
            mFluxFields  = { "Fx", "Fy" };

            this->initialize() ;
        }

//------------------------------------------------------------------------------

        void
        IWG_PlaneStress::compute_jacobian(
                Element        * aElement,
                Matrix <real> & aJacobian )
        {

            // reset result vectors
            aJacobian.fill( 0.0 );

            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // get the B-Matrix
            Matrix< real > & tdNdX   = mGroup->work_dNdX();
            Matrix< real > & tB      = mGroup->work_B();

            // elasticity matrix
            Matrix< real > & tC      = mGroup->work_C();

            // loop over all integration points
            tB.fill( 0.0 );

            // compute c
            mMaterial->C_ps( tC );

            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // compute geometry Jacobian
                Matrix< real > & tJ = aElement->J( k );

                // determinant
                real tDetJ = std::abs( det( tJ ) );

                // compute B
                tdNdX = inv( tJ ) * mGroup->dNdXi( k );

                // populate B
                uint tCount = 0;
                for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                {
                    tB( 0, tCount ) = tdNdX( 0, i );
                    tB( 2, tCount ) = tdNdX( 1, i );
                    ++tCount;
                    tB( 1, tCount ) = tdNdX( 1, i );
                    tB( 2, tCount ) = tdNdX( 0, i );
                    ++tCount;
                }

                // add to integration
                aJacobian += tW( k ) * trans( tB ) * tC * tB * tDetJ;
            }
        }

//------------------------------------------------------------------------------
    }
}
