//
// Created by Christian Messe on 16.12.20.
//


#include "fn_Mesh_integrate_scalar_over_sidesets.hpp"
#include "cl_Cell.hpp"
#include "cl_Element.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"
#include "fn_intpoints.hpp"
#include "meshtools.hpp"
#include "fn_intpoints_auto_integration_order.hpp"
#include "fn_intpoints.hpp"
#include "commtools.hpp"
#include "fn_det.hpp"
#include "fn_sum.hpp"
#include "fn_inv.hpp"
#include "fn_dot.hpp"
#include "geometrytools.hpp"
#include "fn_trans.hpp"
#include "assert.hpp"
#include "fn_max.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        real
        integrate_scalar_over_sidesets( Mesh * aMesh,
                                       const string       & aFieldLabel,
                                       const Vector<id_t> & aSidesetIDs ,
                                       const proc_t         aMasterRank  )
        {
            // reset the value
            real aValue = 0.0;

            // this routine only works on the master
            if ( comm_rank() == aMasterRank )
            {

                Vector< real > & tValues = aMesh->field_data( aFieldLabel );

                for ( id_t tSideset : aSidesetIDs )
                {
                    if ( aMesh->sideset_exists( tSideset ) )
                    {
                        // number of dimensions on this mesh
                        uint tNumDim = aMesh->number_of_dimensions();

                        // get pointer to facets
                        Cell< mesh::Facet * > & tFacets = aMesh->sideset( tSideset )->facets();

                        // skip this sideset if it is empty
                        if ( tFacets.size() > 0 )
                        {
                            // create the factory
                            fem::InterpolationFunctionFactory tFactory;

                            // get the type of the surface element
                            ElementType tElemType = tFacets( 0 )->element()->type();

                            // create the shape function
                            fem::InterpolationFunction * tShape = tFactory.create_lagrange_function( tElemType );

                            // number of nodes for the surface element
                            uint tNumNodes = number_of_nodes( tElemType );

                            // integration points for facet
                            Matrix< real > tXi;

                            // integreation weights for facet
                            Vector< real > tW;

                            uint tOrder = auto_integration_order( tElemType );

                            // populate data
                            intpoints(
                                    IntegrationScheme::GAUSS,
                                    geometry_type( tElemType ),
                                    tOrder,
                                    tW,
                                    tXi );

                            // number of integration points
                            uint tNumPoints = tW.length();

                            // number of surface dimensions
                            uint tNumSurfaceDim = dimension( tElemType );

                            // container with shape function
                            Cell< Vector< real > > tN( tNumPoints, Vector< real >( tNumNodes ));


                            // container with derivatives
                            Cell< Matrix< real > > tdNdXi( tNumPoints, Matrix< real >( tNumSurfaceDim, tNumNodes ));

                            // help matrix
                            Matrix< real > tNasMatrix( 1, tNumNodes );

                            // populate data
                            for ( uint k = 0; k < tNumPoints; ++k )
                            {
                                Vector< real > & tNasVector = tN( k );

                                // per definition, the shape funciton is a matrix
                                tShape->N( tXi.col( k ), tNasMatrix );

                                // convert matrix to vector
                                for ( uint i = 0; i < tNumNodes; ++i )
                                {
                                    tNasVector( i ) = tNasMatrix( 0, i );
                                }

                                // compute derivative
                                tShape->dNdXi( tXi.col( k ), tdNdXi( k ));
                            }

                            // delete the shape function
                            delete tShape;

                            // container for node coordinates
                            Matrix< real > tNodeCoords( tNumNodes, tNumDim );

                            // container for nodal values of field
                            Vector< real > tNodeValues( tNumNodes );

                            // loop over all elements
                            for ( mesh::Facet * tFacet : tFacets )
                            {
                                // check of this element is owned
                                //if ( tFacet->master()->owner() == tMyRank )
                                    real tS = 0.0;

                                    // collect node coordinates
                                    collect_node_coords( tFacet->element(), tNodeCoords, tNumDim );

                                    // collect the node values on this element
                                    collect_node_data( tFacet->element(), tValues, tNodeValues );

                                    // loop over all integration points
                                    for ( uint k = 0; k < tNumPoints; ++k )
                                    {
                                        // compute surface increment
                                        real tDetJ = mesh::compute_surface_increment(
                                                tdNdXi( k ), tNodeCoords );

                                        // add contribution to surface
                                        tS += tW( k ) * dot( tN( k ), tNodeValues ) * tDetJ;
                                    }

                                    aValue += tS;
                                //}// end if own element
                            } // end element loop
                        } // end number of facets > 0
                    } // end sideset exists
                } // end sideset loop
            } // end master rank

            comm_barrier() ;
            broadcast( aMasterRank, aValue );

            return aValue;
        }
//------------------------------------------------------------------------------
    }
}