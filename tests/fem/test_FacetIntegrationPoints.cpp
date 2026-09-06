/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include <gtest/gtest.h>
#include <algorithm>

#include "typedefs.hpp"
#include "cl_Mesh.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"
#include "fn_intpoints.hpp"
#include "en_IntegrationScheme.hpp"
#include "fn_IF_initialize_integration_points_on_facet.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"

using namespace belfem;

namespace
{
    const real tEpsilon = 1.0e-12 ;

    /**
     * Test that integration points mapped from volume parametric space
     * to a facet produce the same physical coordinates as direct surface
     * integration, for a given element type.
     *
     * @param aVolType         the volume element type
     * @param aTestOrientations  if true, also test slave (orientation) mappings
     */
    void
    test_facet_integration_points(
        const ElementType aVolType,
        const bool        aTestOrientations = true,
        const uint        aMaxFacets = BELFEM_UINT_MAX )
    {
        mesh::ElementFactory              tElFactory ;
        fem::InterpolationFunctionFactory tIfFactory ;

        auto tRef = tElFactory.create_reference_element( aVolType );

        fem::InterpolationFunction * tVolFun =
            tIfFactory.create_lagrange_function( aVolType );

        uint tOrder = 2 * mesh::interpolation_order_numeric( aVolType ) + 1 ;

        // spatial dimension (2 or 3)
        uint tSpaceDim = 3 ;

        // surface element dimension = volume dimension - 1
        uint tSurfDim = mesh::dimension( aVolType ) - 1 ;

        // volume node coordinates (always stored as 3D)
        uint nnv = tRef->element()->number_of_nodes();
        Matrix< real > Xv( tSpaceDim, nnv );

        for ( uint k = 0; k < nnv; ++k )
        {
            Xv( 0, k ) = tRef->element()->node( k )->x();
            Xv( 1, k ) = tRef->element()->node( k )->y();
            Xv( 2, k ) = tRef->element()->node( k )->z();
        }

        Vector< real > w_v ;
        Vector< real > w_s ;
        Matrix< real > xi_v ;
        Matrix< real > xi_s ;
        Cell< mesh::Node * > tNodesS ;

        uint nf = std::min( tRef->element()->number_of_facets(),
                            aMaxFacets );

        Matrix< real > Ns ;
        Matrix< real > Nv ;
        Matrix< real > Xs ;
        Vector< real > p( tSpaceDim ) ;
        Vector< real > q( tSpaceDim ) ;

        // orientation table (only when testing orientations)
        Matrix< uint > Table ;
        uint tOffset = 0 ;
        if ( aTestOrientations )
        {
            tElFactory.create_orientation_table( tRef->element()->type(), Table );
        }

        for ( uint f = 0; f < nf; ++f )
        {
            tRef->element()->get_nodes_of_facet( f, tNodesS );
            uint nns = tNodesS.size();

            ElementType tSurfType =
                mesh::element_type_from_numnodes( tSurfDim, nns );

            // surface node coordinates (canonical ordering)
            Xs.set_size( tSpaceDim, nns );
            for ( uint k = 0; k < nns; ++k )
            {
                Xs( 0, k ) = tNodesS( k )->x();
                Xs( 1, k ) = tNodesS( k )->y();
                Xs( 2, k ) = tNodesS( k )->z();
            }

            fem::InterpolationFunction * tSurfFun =
                tIfFactory.create_lagrange_function( tSurfType );

            intpoints( IntegrationScheme::GAUSS,
                       mesh::geometry_type( tSurfType ),
                       tOrder, w_s, xi_s );

            fem::initialize_integration_points_on_facet(
                aVolType, f, w_v, xi_v, tOrder );

            uint ng = w_v.length();

            // --- master (canonical) side ---

            EXPECT_EQ( w_v.length(), w_s.length() )
                << "facet " << f << " weight count mismatch" ;

            for ( uint k = 0; k < ng; ++k )
            {
                EXPECT_NEAR( w_v( k ), w_s( k ), tEpsilon )
                    << "facet " << f << " weight " << k ;
            }

            for ( uint k = 0; k < ng; ++k )
            {
                tVolFun->N( xi_v.col( k ), Nv );
                tSurfFun->N( xi_s.col( k ), Ns );

                for ( uint d = 0; d < tSpaceDim; ++d )
                {
                    p( d ) = dot( Nv.row( 0 ), Xv.row( d ) );
                    q( d ) = dot( Ns.row( 0 ), Xs.row( d ) );
                }

                EXPECT_NEAR( norm( p - q ), 0.0, tEpsilon )
                    << "facet " << f << " master point " << k ;
            }

            // --- slave (orientation) side ---

            if ( aTestOrientations )
            {
                uint tNumOrientations =
                    mesh::number_of_corner_nodes( tSurfType );

                for ( uint o = 0; o < tNumOrientations; ++o )
                {
                    // permute surface node coordinates via orientation table
                    {
                        uint tCount = 0 ;
                        for ( uint r = 0; r < Table.n_rows(); ++r )
                        {
                            uint i = Table( r, tOffset );
                            if ( i != BELFEM_UINT_MAX )
                            {
                                Xs( 0, tCount ) = tRef->node( i )->x();
                                Xs( 1, tCount ) = tRef->node( i )->y();
                                Xs( 2, tCount ) = tRef->node( i )->z();
                                ++tCount ;
                            }
                        }
                    }
                    ++tOffset ;

                    fem::initialize_integration_points_on_facet(
                        aVolType, f, o, w_v, xi_v, tOrder );

                    for ( uint k = 0; k < ng; ++k )
                    {
                        EXPECT_NEAR( w_v( k ), w_s( k ), tEpsilon )
                            << "facet " << f << " orient " << o << " weight " << k ;
                    }

                    for ( uint k = 0; k < ng; ++k )
                    {
                        tVolFun->N( xi_v.col( k ), Nv );
                        tSurfFun->N( xi_s.col( k ), Ns );

                        for ( uint d = 0; d < tSpaceDim; ++d )
                        {
                            p( d ) = dot( Nv.row( 0 ), Xv.row( d ) );
                            q( d ) = dot( Ns.row( 0 ), Xs.row( d ) );
                        }

                        EXPECT_NEAR( norm( p - q ), 0.0, tEpsilon )
                            << "facet " << f << " orient " << o << " point " << k ;
                    }
                }
            }

            delete tSurfFun ;
        }

        delete tVolFun ;
        delete tRef ;
    }
}

// ---------------------------------------------------------------------------
//  Tetrahedra
// ---------------------------------------------------------------------------

TEST( FacetIntegrationPoints, TET4  ) { test_facet_integration_points( ElementType::TET4  ); }
TEST( FacetIntegrationPoints, TET10 ) { test_facet_integration_points( ElementType::TET10 ); }
TEST( FacetIntegrationPoints, TET20 ) { test_facet_integration_points( ElementType::TET20 ); }
TEST( FacetIntegrationPoints, TET35 ) { test_facet_integration_points( ElementType::TET35 ); }

// ---------------------------------------------------------------------------
//  Hexahedra
// ---------------------------------------------------------------------------

TEST( FacetIntegrationPoints, HEX8  ) { test_facet_integration_points( ElementType::HEX8  ); }
TEST( FacetIntegrationPoints, HEX20 ) { test_facet_integration_points( ElementType::HEX20 ); }
TEST( FacetIntegrationPoints, HEX27 ) { test_facet_integration_points( ElementType::HEX27 ); }
TEST( FacetIntegrationPoints, HEX64 ) { test_facet_integration_points( ElementType::HEX64 ); }

// ---------------------------------------------------------------------------
//  Pentahedra (wedges)
// ---------------------------------------------------------------------------

TEST( FacetIntegrationPoints, PENTA6  ) { test_facet_integration_points( ElementType::PENTA6  ); }
TEST( FacetIntegrationPoints, PENTA15 ) { test_facet_integration_points( ElementType::PENTA15 ); }
TEST( FacetIntegrationPoints, PENTA18 ) { test_facet_integration_points( ElementType::PENTA18 ); }

// ---------------------------------------------------------------------------
//  Pyramids
// ---------------------------------------------------------------------------

TEST( FacetIntegrationPoints, PYRA5  ) { test_facet_integration_points( ElementType::PYRA5  ); }
TEST( FacetIntegrationPoints, PYRA13 ) { test_facet_integration_points( ElementType::PYRA13 ); }
TEST( FacetIntegrationPoints, PYRA14 ) { test_facet_integration_points( ElementType::PYRA14 ); }

// ---------------------------------------------------------------------------
//  Thin-shell QUAD (4 facets, 1D LINE facets with 2 orientations each)
// ---------------------------------------------------------------------------

TEST( FacetIntegrationPoints, QUAD4TS )
{
    test_facet_integration_points( ElementType::QUAD4TS, false );
}

TEST( FacetIntegrationPoints, QUAD9TS )
{
    test_facet_integration_points( ElementType::QUAD9TS, false );
}

// ---------------------------------------------------------------------------
//  Thin-shell PENTA (5 facets, matches volume PENTA topology).
//  Since Option-B canonicalization, PENTA*TS routes through the volume
//  PENTA integration-point / orientation dispatch.
// ---------------------------------------------------------------------------

TEST( FacetIntegrationPoints, PENTA6TS )
{
    test_facet_integration_points( ElementType::PENTA6TS );
}

TEST( FacetIntegrationPoints, PENTA18TS )
{
    test_facet_integration_points( ElementType::PENTA18TS );
}
