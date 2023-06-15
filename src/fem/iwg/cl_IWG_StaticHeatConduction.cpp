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
                const uint aNumberOfDimensions,
                const IwgType aType ) :
                IWG_Timestep( aType )
        {
            mNumberOfDofsPerNode = 1;
            mNumberOfDerivativeDimensions = aNumberOfDimensions ;
            mNumberOfRhsCols = 1;

            // timestep must be 1.0 for the stationary class
            // do not touch, only children should overwrite this
            mDeltaTime = 1.0 ;

            // set the names for all the fields
            mDofFields = { "T" };

            mFluxFields = { "dotQ" };

            this->initialize() ;
        }

//------------------------------------------------------------------------------
// Functions called by Field during assembly
//------------------------------------------------------------------------------

        void
        IWG_StaticHeatConduction::compute_jacobian_and_rhs(
                Element        * aElement,
                Matrix< real > & aJacobian,
                Vector< real > & aRHS )
        {
            // reset result vectors
            aJacobian.fill( 0.0 );
            aRHS.fill( 0.0 );

            // get integration weights
            const Vector< real > & tW = mGroup->integration_weights();

            // get the B-Matrix
            Matrix< real > & tB      = mGroup->work_dNdX();

            // nodal temperatures
            Vector< real > & tThat   = mGroup->work_phi();

            // matrix for thermal conduction
            Matrix< real > & tLambda = mGroup->work_C() ;

            // collect temperatures from last iteration
            this->collect_node_data( aElement, "T", tThat );

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {
                // compute geometry Jacobian
                Matrix< real > & tJ = aElement->J( k );

                // interpolate temperature for this point
                real tT = dot( mGroup->n( k ), tThat );

                // determinant
                real tDetJ = std::abs( det( tJ ) );

                // compute B
                tB = inv( tJ ) * mGroup->dNdXi( k );

                // compute thermal conductivity
                mGroup->thermal_conductivity( tLambda, tT );

                // add to integration
                aJacobian += tW( k ) * trans( tB ) * tLambda * tB * tDetJ;

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
            const Vector< real > & tW = mGroup->integration_weights();

            // matrix for boundary conditions
            Matrix< real > & tPsi = mGroup->work_Psi() ;

            // node coordinates
            Matrix< real > & tX = mGroup->node_coords() ;

            Cell< string > tFieldLabels = { "alpha", "Tinf" } ;

            // collect temperatures from last iteration
            this->collect_node_data( aElement, tFieldLabels, tPsi );

            // loop over all integration points
            for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
            {

                // interpolate alpha for this point
                real tAlpha = dot( mGroup->n( k ).vector_data(), tPsi.col( 0 ) );

                // interpolate T inf for this point
                real tTinf = dot( mGroup->n( k ).vector_data(), tPsi.col( 1 ) );

                // determinant along surface
                real tDetJ =  mesh::compute_surface_increment( mGroup->dNdXi( k ), tX );

                // add to integration
                aJacobian += tW( k ) * trans( mGroup->N( k ).row( 0 ) ) * tAlpha
                             * mGroup->N( k ).row( 0 ) * tDetJ ;


                aRHS += tW( k ) * trans( mGroup->N( k ).row( 0 ) ) * tAlpha
                               * tTinf * tDetJ ;

            }

            aJacobian *= mDeltaTime ;
            aRHS *= mDeltaTime ;
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

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // get data containers
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // integration weights
        const Vector< real > & tW = mGroup->integration_weights();

        // nodal information for source fields
        Vector < real > & tpsi = mGroup->work_psi();



        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // create shortcuts for better readability
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // reset the load vector
        aConvection.fill( 0.0 );

        this->collect_node_data( aElement, "dotQ", tpsi );

        // collect node coords from master
        uint tNumDim = mField->mesh()->number_of_dimensions() ;
        Matrix< real > & tX = mGroup->node_coords() ;

        mesh::collect_node_coords( aElement->element(),
                                   tX,
                                   tNumDim,
                                   mGroup->parent()->enforce_linear_interpolation() );

        real tDotQ ;

        // loop over all integration points
        for( uint k=0; k<mNumberOfIntegrationPoints; ++k )
        {
            // compute heatflux
            tDotQ = dot( mGroup->N( k ).row( 0 ), tpsi.vector_data() );

            aConvection += tW( k ) * trans( mGroup->N( k ).row( 0 ) ) * tDotQ
                           *  mesh::compute_surface_increment( mGroup->dNdXi( k ), tX );
        }

    }

//------------------------------------------------------------------------------

        void
        IWG_StaticHeatConduction::allocate_work_matrices( Group    * aGroup )
        {
            // call function from parent
            IWG::allocate_work_matrices( aGroup );

            // change size of Psi Matrix, needed for alpha and T_inf
            // and alpha and T_inf boundary condition
            aGroup->work_Psi().set_size(
                    mNumberOfNodesPerElement,
                    2 );

            // needed for dotq
            aGroup->work_psi().set_size( mNumberOfNodesPerElement );

            aGroup->work_C().set_size( mNumberOfSpatialDimensions,
                                       mNumberOfSpatialDimensions );

            /*
            // Help Matrix
            aGroup->work_H().set_size( mNumberOfDofsPerNode * mNumberOfNodesPerElement,
                                       mNumberOfDofsPerNode * mNumberOfNodesPerElement,
                                       0.0 );

            // change size of psi vector
             ) ; */
        }

//------------------------------------------------------------------------------
    } /* end namespace fem */
} /* end namespace belfem */