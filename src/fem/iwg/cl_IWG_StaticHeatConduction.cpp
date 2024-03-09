   //
// Created by Christian Messe on 24.07.20.
//

#include "cl_IWG_StaticHeatConduction.hpp"
#include "cl_FEM_Group.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Field.hpp"
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

            mNumberOfRhsCols = 1;

            // timestep must be 1.0 for the stationary class
            // do not touch, only children should overwrite this
            mDeltaTime = 1.0 ;

            // set the names for all the fields
            mDofFields = { "T" };

            mFluxFields = { "dotQ" };

            mBcFields = { "alpha", "Tinf" } ;

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

            mCalc->link( aElement );


            // get integration weights
            const Vector< real > & tW =  mCalc->weights();

            // nodal temperatures
            Vector< real > & tTnodes   = mCalc->vector("T");

            // collect temperatures from last iteration
            this->collect_node_data( aElement, "T", tTnodes );

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // get the B-Operator
                const Matrix< real > & tB = mCalc->B( k );

                real tT = mCalc->node_interp( k, tTnodes );

                // add to integration
                aK += tW( k ) * trans( tB ) *
                        mMaterial->lambda( tT ) *
                        tB * det( mCalc->J( k ) );
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

            // get integration weights
            //const Vector< real > & tW = mGroup->integration_weights();


            // node coordinates
           // Matrix< real > & tX = mCalc->matrix("X");
            //this->collect_node_coords( aElement, tX );

            // collect temperatures from last iteration
            //this->collect_node_data( aElement, mBcFields );


            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // interpolate alpha for this point
                /*real tAlpha = dot( mGroup->n( k ).vector_data(), mCalc->vector( "alpha").vector_data() );

                // interpolate T inf for this point
                real tTinf = dot( mGroup->n( k ).vector_data(), mCalc->vector( "T").vector_data() );

                // determinant along surface
                real tDetJ =  mesh::compute_surface_increment( mGroup->dNdXi( k ), tX );

                // add to integration
                aJacobian += tW( k ) * trans( mGroup->N( k ).row( 0 ) ) * tAlpha
                             * mGroup->N( k ).row( 0 ) * tDetJ ;


                aRHS += tW( k ) * trans( mGroup->N( k ).row( 0 ) ) * tAlpha
                               * tTinf * tDetJ ; */

            }




            //aJacobian *= mDeltaTime ;
            //aRHS *= mDeltaTime ;
        }

//------------------------------------------------------------------------------

    void
    IWG_StaticHeatConduction::compute_convection(
            Element * aElement, Vector< real > & aConvection )
    {
        // Check input
        BELFEM_ASSERT( aConvection.length() ==
                      mGroup->parent()->enforce_linear_interpolation() ?
                      mesh::number_of_corner_nodes( aElement->element()->type() )
                      : mesh::number_of_nodes( aElement->element()->type() ),
                      "load vector has wrong length" );

        BELFEM_ASSERT( mMesh->number_of_dimensions() == 2, "the convection routine has been written for the 2D case" );

        BELFEM_ERROR( false, "this must be redone");

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // get data containers
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // integration weights
        const Vector< real > & tW = mGroup->integration()->weights();

        // nodal information for source fields
        Vector < real > & tdotq = mCalc->vector("dotq");


        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // create shortcuts for better readability
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // reset the load vector
        aConvection.fill( 0.0 );

        this->collect_node_data( aElement, "dotQ", tdotq );

        // collect node coords from master
        uint tNumDim = mField->mesh()->number_of_dimensions() ;
        Matrix< real > & tX = mGroup->calculator()->matrix("X");

        mesh::collect_node_coords( aElement->element(),
                                   tX,
                                   tNumDim,
                                   mGroup->parent()->enforce_linear_interpolation() );

        real tDotQ ;

        // loop over all integration points
        for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
        {
            // compute heatflux
            tDotQ = dot( mGroup->integration()->phi( k ), tdotq );

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