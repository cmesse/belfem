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
#include <cmath>
#include <memory>
#include "cl_ReferenceElement.hpp"
#include "cl_IF_IntegrationData.hpp"
#include "meshtools.hpp"

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

// ---------------------------------------------------------------------------
//  2-D volumes, master side: every edge facet against direct LINE integration
// ---------------------------------------------------------------------------

TEST( FacetIntegrationPoints, TRI3  ) { test_facet_integration_points( ElementType::TRI3,  false ); }
TEST( FacetIntegrationPoints, TRI6  ) { test_facet_integration_points( ElementType::TRI6,  false ); }
TEST( FacetIntegrationPoints, QUAD4 ) { test_facet_integration_points( ElementType::QUAD4, false ); }
TEST( FacetIntegrationPoints, QUAD9 ) { test_facet_integration_points( ElementType::QUAD9, false ); }

// ---------------------------------------------------------------------------
//  2-D slave side, through IntegrationData as the calculators use it.
//
//  The neighbor B is the point reflection of the reference element A through
//  the midpoint of edge f: B is a valid, counterclockwise element whose edge
//  f covers the same segment as A's, traversed the other way. The master
//  points of A on f and the slave points of B on f must then coincide one by
//  one, with equal weights. Volume TRI and QUAD reverse the columns ( one
//  implicit orientation ); QUAD*TS carries two explicit orientations, and
//  exactly one of them must coincide while the other runs backwards.
// ---------------------------------------------------------------------------

namespace
{
    //! physical coordinates of the points of aData on an element with nodes aX
    void
    map_points( fem::InterpolationFunction * aFunction,
                const Matrix< real > & aX,
                const fem::IntegrationData & aData,
                Matrix< real > & aPhysical )
    {
        const Matrix< real > & tXi = aData.points() ;
        const uint tNumPoints = tXi.n_cols() ;
        aPhysical.set_size( 2, tNumPoints, 0.0 );

        Matrix< real > tN ;
        for ( uint k = 0; k < tNumPoints; ++k )
        {
            aFunction->N( tXi.col( k ), tN );
            for ( uint d = 0; d < 2; ++d )
            {
                aPhysical( d, k ) = dot( tN.row( 0 ), aX.row( d ) );
            }
        }
    }

    //! largest distance between column k of aA and column k ( or n-1-k if
    //! aReversed ) of aB; NaN if any distance is NaN
    real
    point_distance( const Matrix< real > & aA, const Matrix< real > & aB, const bool aReversed )
    {
        const uint n = aA.n_cols() ;
        real tMax = 0.0 ;
        for ( uint k = 0; k < n; ++k )
        {
            const uint j = aReversed ? n - 1 - k : k ;
            const real tDistance = std::sqrt(
                    ( aA( 0, k ) - aB( 0, j ) ) * ( aA( 0, k ) - aB( 0, j ) )
                  + ( aA( 1, k ) - aB( 1, j ) ) * ( aA( 1, k ) - aB( 1, j ) ) );

            // written so that a NaN distance makes the result NaN, which
            // no tolerance accepts; std::max would drop it
            if ( ! ( tDistance <= tMax ) ) tMax = tDistance ;
        }
        return tMax ;
    }

    void
    test_slave_points_2d( const ElementType aType )
    {
        mesh::ElementFactory              tElFactory ;
        fem::InterpolationFunctionFactory tIfFactory ;

        std::unique_ptr< mesh::ReferenceElement > tRef(
                tElFactory.create_reference_element( aType ) );
        std::unique_ptr< fem::InterpolationFunction > tFun(
                tIfFactory.create_lagrange_function( aType ) );

        const uint tOrder = 2 * mesh::interpolation_order_numeric( aType ) + 1 ;
        const bool tIsThinShell = aType == ElementType::QUAD4TS || aType == ElementType::QUAD9TS ;

        const uint n = tRef->element()->number_of_nodes() ;
        Matrix< real > Xa( 2, n );
        for ( uint k = 0; k < n; ++k )
        {
            Xa( 0, k ) = tRef->element()->node( k )->x() ;
            Xa( 1, k ) = tRef->element()->node( k )->y() ;
        }

        Cell< mesh::Node * > tFacetNodes ;
        Matrix< real > tMasterX ;
        Matrix< real > tSlaveX ;

        for ( uint f = 0; f < tRef->element()->number_of_facets(); ++f )
        {
            // B = point reflection of A through the midpoint of edge f
            tRef->element()->get_nodes_of_facet( f, tFacetNodes );
            const real tMx = 0.5 * ( tFacetNodes( 0 )->x() + tFacetNodes( 1 )->x() );
            const real tMy = 0.5 * ( tFacetNodes( 0 )->y() + tFacetNodes( 1 )->y() );

            Matrix< real > Xb( 2, n );
            for ( uint k = 0; k < n; ++k )
            {
                Xb( 0, k ) = 2.0 * tMx - Xa( 0, k );
                Xb( 1, k ) = 2.0 * tMy - Xa( 1, k );
            }

            fem::IntegrationData tMaster( aType );
            tMaster.populate_for_master( f, tOrder );
            map_points( tFun.get(), Xa, tMaster, tMasterX );

            ASSERT_GT( tMaster.weights().length(), 1u ) << "facet " << f ;

            const uint tNumOrientations = tIsThinShell ? 2 : 1 ;
            uint tNumMatches = 0 ;

            for ( uint o = 0; o < tNumOrientations; ++o )
            {
                fem::IntegrationData tSlave( aType );
                tSlave.populate_for_slave( f, o, tOrder );
                map_points( tFun.get(), Xb, tSlave, tSlaveX );

                ASSERT_EQ( tSlave.weights().length(), tMaster.weights().length() )
                    << "facet " << f << " orientation " << o ;

                const real tForward  = point_distance( tMasterX, tSlaveX, false );
                const real tBackward = point_distance( tMasterX, tSlaveX, true );

                if ( tForward < tEpsilon )
                {
                    ++tNumMatches ;
                    for ( uint k = 0; k < tMaster.weights().length(); ++k )
                    {
                        EXPECT_NEAR( tSlave.weights()( k ), tMaster.weights()( k ), tEpsilon )
                            << "facet " << f << " orientation " << o << " weight " << k ;
                    }
                }
                else
                {
                    // the other orientation covers the same segment backwards
                    EXPECT_TRUE( tIsThinShell )
                        << "facet " << f << ": slave points do not coincide, distance " << tForward ;
                    EXPECT_LT( tBackward, tEpsilon )
                        << "facet " << f << " orientation " << o << " is not the reversal" ;
                }
            }

            EXPECT_EQ( tNumMatches, 1u ) << "facet " << f ;
        }
    }
}

TEST( FacetSlavePoints, TRI3  ) { test_slave_points_2d( ElementType::TRI3  ); }
TEST( FacetSlavePoints, TRI6  ) { test_slave_points_2d( ElementType::TRI6  ); }
TEST( FacetSlavePoints, QUAD4 ) { test_slave_points_2d( ElementType::QUAD4 ); }
TEST( FacetSlavePoints, QUAD9 ) { test_slave_points_2d( ElementType::QUAD9 ); }
TEST( FacetSlavePoints, QUAD4TS_Orientations ) { test_slave_points_2d( ElementType::QUAD4TS ); }
TEST( FacetSlavePoints, QUAD9TS_Orientations ) { test_slave_points_2d( ElementType::QUAD9TS ); }

// the pairing the thin-shell factory writes for a ghost facet: the lower
// layer's top facet is the master, the upper layer's bottom facet the slave
// with orientation 1. The layers are stacked, not reflected
TEST( FacetSlavePoints, QUAD4TS_GhostPairing )
{
    const ElementType tType = ElementType::QUAD4TS ;

    mesh::ElementFactory              tElFactory ;
    fem::InterpolationFunctionFactory tIfFactory ;

    std::unique_ptr< mesh::ReferenceElement > tRef(
            tElFactory.create_reference_element( tType ) );
    std::unique_ptr< fem::InterpolationFunction > tFun(
            tIfFactory.create_lagrange_function( tType ) );

    real tYMin = BELFEM_REAL_MAX ;
    real tYMax = -BELFEM_REAL_MAX ;
    Matrix< real > Xa( 2, 4 );
    for ( uint k = 0; k < 4; ++k )
    {
        Xa( 0, k ) = tRef->element()->node( k )->x() ;
        Xa( 1, k ) = tRef->element()->node( k )->y() ;
        tYMin = std::min( tYMin, Xa( 1, k ) );
        tYMax = std::max( tYMax, Xa( 1, k ) );
    }

    // the upper layer sits directly on top of the lower one
    Matrix< real > Xb = Xa ;
    for ( uint k = 0; k < 4; ++k )
    {
        Xb( 1, k ) += tYMax - tYMin ;
    }

    fem::IntegrationData tMaster( tType );
    tMaster.populate_for_master( mesh::top_facet_index( tType ), 3 );

    fem::IntegrationData tSlave( tType );
    tSlave.populate_for_slave( mesh::bottom_facet_index( tType ), 1, 3 );

    Matrix< real > tMasterX ;
    Matrix< real > tSlaveX ;
    map_points( tFun.get(), Xa, tMaster, tMasterX );
    map_points( tFun.get(), Xb, tSlave,  tSlaveX );

    ASSERT_GT( tMasterX.n_cols(), 1u );
    EXPECT_LT( point_distance( tMasterX, tSlaveX, false ), tEpsilon );
}
