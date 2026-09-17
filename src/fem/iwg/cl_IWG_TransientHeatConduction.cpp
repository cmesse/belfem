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

#include "cl_IWG_TransientHeatConduction.hpp"
#include "cl_FEM_Dof.hpp"
#include "cl_FEM_DofManagerBase.hpp"
#include "cl_FEM_Kernel.hpp"
#include "fn_trans.hpp"

namespace belfem
{
    namespace fem
    {
        IWG_TransientHeatConduction::IWG_TransientHeatConduction(
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

            // set the names for all the fields
            mDofFields = { "T" };

            mHasConvection = false ;

            // conditional initialization, skipped by children
            this->initialize( IwgType::TransientHeatConduction ) ;
        }

//------------------------------------------------------------------------------

        void
        IWG_TransientHeatConduction::compute_mkf(
            Element * aElement )
        {
            // link calculator with element
            mGroup->calculator()->link( aElement );

            mTimeStepMatrices->reset() ;

            Matrix< real > & M = mTimeStepMatrices->M();
            Matrix< real > & K = mTimeStepMatrices->K();
            Vector< real > & f = mTimeStepMatrices->f();

            // reset matrices
            M.fill( 0.0 );
            K.fill( 0.0 );
            f.fill( 0.0 );

            // get integration weights
            Calculator * aCalc = mGroup->calculator();
            const Vector< real > & w =  aCalc->integration()->weights();
            //const Vector< real > tTnodes = aCalc->node_data( "T") ;

            // loop over all integration points
            for( uint k=0; k<aCalc->num_intpoints(); ++k )
            {

                //real tT = aCalc->node_interp( k, tTnodes );

                // get the gradient operator
                const Matrix< real > & B = aCalc->B( k );

                // get the gradient operator
                const Matrix< real > & N = aCalc->N(k) ;

                // mass matrix
                // TODO(cm): the mass matrix uses a unit coefficient; use the
                // material's cp once this IWG reads thermal properties
                M += w( k ) * trans( N ) * 1.0 *  N * aCalc->dV( k );

                //stiffness matrix
                //aK += w( k ) * trans( B ) * aCalc->material()->lambda( tT ) * B * aCalc->dV( k );
                K += w( k ) * trans( B ) * 111.0 * B * aCalc->dV( k );
            }
        }

        //------------------------------------------------------------------------------

    }

}