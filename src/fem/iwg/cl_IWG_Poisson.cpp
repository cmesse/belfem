/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include "cl_IWG_Poisson.hpp"
#include "cl_FEM_Calculator.hpp"
#include "fn_trans.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Poisson::IWG_Poisson (
                const ModelDimensionality  aModelDimensionality,
                const IwgType aType,
                const IwgMode aMode ) :
                IWG_Timestep( aType, aModelDimensionality, aMode )
        {
            mNumberOfDofsPerNode = 1;

            switch ( aModelDimensionality )
            {
                case( ModelDimensionality::TwoD ) :
                case( ModelDimensionality::AxSymmX ) :
                case( ModelDimensionality::AxSymmY ) :
                {
                    mNumberOfSpatialDimensions = 2 ;
                    mNumberOfDerivativeDimensions = 2 ;
                    break ;
                }
                case( ModelDimensionality::ThreeD ) :
                {
                    mNumberOfSpatialDimensions = 3 ;
                    mNumberOfDerivativeDimensions = 3 ;
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Invalid Model Dimensionality");
                }
            }

            // timestep must be 1.0 for the stationary class
            // do not touch, only children should overwrite this
            mDeltaTime = 1.0 ;

            // set the names for all the fields
            mDofFields = { "phi" };

            mHasConvection = false ;

            // conditional initialization, skipped by children
            this->initialize( IwgType::Poisson ) ;

            mTimeStepMatrices->set_flag( MatrixFlag::K ) ;
            mTimeStepMatrices->set_flag( MatrixFlag::F ) ;

        }

//------------------------------------------------------------------------------

        void
        IWG_Poisson::compute_jacobian(
                Element        * aElement,
                Matrix< real > & aJacobian )
        {
            Calculator * tCalc = this->calc();

            aJacobian.fill( 0.0 );

            // link calculator
            tCalc->link( aElement );

            // get integration weights
            const Vector< real > & tW =  tCalc->integration()->weights();

            // loop over all integration points
            for( uint k=0; k<tCalc->num_intpoints(); ++k )
            {
                // get the gradient operator
                const Matrix< real > & tB = tCalc->B( k );

                // add to integration
                aJacobian += tW( k ) * trans( tB ) *
                      tB * tCalc->dV( k );
            }
        }

//------------------------------------------------------------------------------
    }
}