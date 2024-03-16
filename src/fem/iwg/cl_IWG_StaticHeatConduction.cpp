   //
// Created by Christian Messe on 24.07.20.
//

#include "cl_IWG_StaticHeatConduction.hpp"
#include "cl_FEM_Group.hpp"
#include "cl_FEM_Element.hpp"
#include "fn_det.hpp"
#include "fn_trans.hpp"
#include "fn_inv.hpp"
#include "fn_dot.hpp"
#include "geometrytools.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_StaticHeatConduction::IWG_StaticHeatConduction(
                const ModelDimensionality  aModelDimensionality,
                const IwgType aType ) :
                IWG_Timestep( aType, aModelDimensionality )
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
            mDofFields = { "T" };

            mFluxFields = { "dotQ" };

            mOtherFields = { "alpha", "Tinf" } ;

            // heat conduction can have a convection term
            mHasConvection = true ;

            this->initialize() ;

            // allocate lambda matrix
            mLambda.set_size( mNumberOfSpatialDimensions, mNumberOfSpatialDimensions );
        }
//------------------------------------------------------------------------------
// Functions called by Field during assembly
//------------------------------------------------------------------------------

        void
        IWG_StaticHeatConduction::compute_conduction(
                Element        * aElement,
                Matrix< real > & aK,
                Vector< real > & aQ )
        {
            // reset result vectors
            aK.fill( 0.0 );
            aQ.fill( 0.0 );

            // link calculator
            mCalc->link( aElement );

            // get integration weights
            const Vector< real > & tW =  mCalc->integration()->weights();

            // nodal temperatures
            const Vector< real > & tTnodes   = mCalc->node_data("T");

            // loop over all integration points
            for( uint k=0; k<mCalc->num_intpoints(); ++k )
            {
                // get the B-Operator
                const Matrix< real > & tB = mCalc->B( k );

                real tT = mCalc->node_interp( k, tTnodes );

                // add to integration
                aK += tW( k ) * trans( tB ) *
                        mMaterial->lambda( tT ) *
                        tB * mCalc->dV( k );
            }
        }

//------------------------------------------------------------------------------

        void
        IWG_StaticHeatConduction::compute_alpha_boundary_condition(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {

            // reset result vectors
            aJacobian.fill( 0.0 );
            aRHS.fill( 0.0 );

            // link calculator
            mCalc->link( aElement );

            // get integration weights
            const Vector< real > & tW = mCalc->integration()->weights();

            // collect temperatures from last iteration
            const Vector< real > & tAlphaNodes = mCalc->node_data("alpha");
            const Vector< real > & tTinfNodes  = mCalc->node_data("Tinf");

            // loop over all integration points
            for( uint k=0; k<mCalc->num_intpoints(); ++k )
            {
                // interpolate alpha for this point
                real tAlpha = mCalc->node_interp( k , tAlphaNodes );

                // interpolate T inf for this point
                real tTinf = mCalc->node_interp( k, tTinfNodes );


                // add to integration
                aJacobian += tW( k ) * trans( mCalc->N( k ) ) * tAlpha * mCalc->N( k ) * mCalc->dS( k );


                aRHS += tW( k ) * trans( mCalc->N( k ).row( 0 ) ) * tAlpha
                               * tTinf * mCalc->dS( k );

            }
        }

//------------------------------------------------------------------------------

    void
    IWG_StaticHeatConduction::compute_convection(
            Element * aElement, Vector< real > & aConvection )
    {
        // Check input
        BELFEM_ASSERT( aConvection.length() ==
                       mesh::number_of_nodes( aElement->element()->type() ),
                      "load vector has wrong length" );

        BELFEM_ASSERT( mMesh->number_of_dimensions() == 2, "the convection routine has been written for the 2D case" );

        // reset the load vector
        aConvection.fill( 0.0 );

        // link element
        mCalc->link( aElement );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // get data containers
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // integration weights
        const Vector< real > & tW = mCalc->integration()->weights();
        const Vector < real > & tdotq = mCalc->node_data("dotQ");

        real tDotQ ;

        // loop over all integration points
        for( uint k=0; k<mCalc->num_intpoints(); ++k )
        {
            // compute heatflux
            tDotQ = mCalc->node_interp( k, tdotq );

            aConvection += tW( k ) * trans( mCalc->N( k ).row( 0 ) ) * tDotQ
                           * mCalc->dS( k );
        }
    }

//------------------------------------------------------------------------------

        void
        IWG_StaticHeatConduction::link_to_group( Group * aGroup )
        {
            IWG::link_to_group( aGroup );

            switch ( aGroup->domain_type() )
            {

                case( DomainType::Default ) :
                {
                    mFunComputeJacobian = & IWG_StaticHeatConduction::compute_conduction ;
                    break ;
                }
                case( DomainType::Neumann ) :
                {
                   //  Neumann blocks don't have a Jacobian
                   mFunComputeJacobian = nullptr ;
                   break ;
                }
                case( DomainType::ThermalAlpha ) :
                {
                    mFunComputeJacobian = & IWG_StaticHeatConduction::compute_alpha_boundary_condition ;

                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Illegal Domain type." );
                    mFunComputeJacobian = nullptr ;
                }
            }
        }


//------------------------------------------------------------------------------
    } /* end namespace fem */
} /* end namespace belfem */