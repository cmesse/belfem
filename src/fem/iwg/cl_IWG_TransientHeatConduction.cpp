//
// Created by gregorygiard on 10/23/25.
//

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
                // todo:: replace by cp and lambda when we have thermal properties of material
                //aM += w( k ) * trans( N ) * aCalc->material()->cp( tT ) *  N * aCalc->dV( k );
                M += w( k ) * trans( N ) * 1.0 *  N * aCalc->dV( k );

                //stiffness matrix
                //aK += w( k ) * trans( B ) * aCalc->material()->lambda( tT ) * B * aCalc->dV( k );
                K += w( k ) * trans( B ) * 111.0 * B * aCalc->dV( k );
            }
        }

        //------------------------------------------------------------------------------

    }

}