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

#include "typedefs.hpp"
#include "cl_Mesh.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"
#include "fn_intpoints.hpp"
#include "en_IntegrationScheme.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"

using namespace belfem;

namespace
{
    const real tEpsilon   = 1.0e-12 ;  // exact comparisons
    const real tEpsilonFD = 1.0e-4 ;   // finite difference tolerance
    const real tDxi       = 1.0e-5 ;   // finite difference step

    /**
     * Test Lagrange interpolation for a given element type:
     *   1. Kronecker delta: N_i(xi_j) = delta_ij
     *   2. Coordinate interpolation at nodes
     *   3. Partition of unity at integration points
     *   4. First derivative via central difference
     *   5. Second derivative via central difference
     */
    void
    test_lagrange_interpolation( const ElementType aType )
    {
        mesh::ElementFactory              tElFactory ;
        fem::InterpolationFunctionFactory tIfFactory ;

        auto tRef = tElFactory.create_reference_element( aType );
        fem::InterpolationFunction * tFun = tIfFactory.create_lagrange_function( aType );

        uint tNumNodes = tRef->element()->number_of_nodes();
        uint tNumDim   = mesh::dimension( aType );

        // parametric node coordinates
        Matrix< real > Xi ;
        tFun->param_coords( Xi );

        // physical node coordinates (dim x nodes)
        Matrix< real > X( tNumDim, tNumNodes );
        for ( uint k = 0; k < tNumNodes; ++k )
        {
            for ( uint d = 0; d < tNumDim; ++d )
            {
                X( d, k ) = tRef->element()->node( k )->x( d );
            }
        }

        // --- Test 1: Kronecker delta ---

        Matrix< real > N ;
        for ( uint j = 0; j < tNumNodes; ++j )
        {
            tFun->N( Xi.col( j ), N );
            for ( uint i = 0; i < tNumNodes; ++i )
            {
                EXPECT_NEAR( N( 0, i ), static_cast< real >( i == j ), tEpsilon )
                    << "Kronecker delta N(" << i << ", xi_" << j << ")" ;
            }
        }

        // --- Test 2: coordinate interpolation ---

        for ( uint j = 0; j < tNumNodes; ++j )
        {
            tFun->N( Xi.col( j ), N );
            for ( uint d = 0; d < tNumDim; ++d )
            {
                real x = dot( N.row( 0 ), X.row( d ) );
                EXPECT_NEAR( x, tRef->element()->node( j )->x( d ), tEpsilon )
                    << "coordinate interpolation node " << j << " dim " << d ;
            }
        }

        // --- Test 3: partition of unity at integration points ---

        uint tOrder = 2 * mesh::interpolation_order_numeric( tRef->element()->type() ) + 1 ;

        Vector< real > w ;
        Matrix< real > XiInt ;
        intpoints( IntegrationScheme::GAUSS,
                   mesh::geometry_type( tRef->element()->type() ),
                   tOrder, w, XiInt );

        uint tNumIntPoints = w.length();

        for ( uint k = 0; k < tNumIntPoints; ++k )
        {
            tFun->N( XiInt.col( k ), N );
            real tSum = 0.0 ;
            for ( uint i = 0; i < tNumNodes; ++i )
            {
                tSum += N( 0, i );
            }
            EXPECT_NEAR( tSum, 1.0, tEpsilon )
                << "partition of unity at integration point " << k ;
        }

        // --- Test 4: first derivative via central difference ---

        Matrix< real > dNdXi ;

        for ( uint k = 0; k < tNumIntPoints; ++k )
        {
            tFun->dNdXi( XiInt.col( k ), dNdXi );

            for ( uint d = 0; d < tNumDim; ++d )
            {
                Vector< real > xi0( XiInt.col( k ) );
                Vector< real > xi1( XiInt.col( k ) );
                xi0( d ) -= tDxi ;
                xi1( d ) += tDxi ;

                Matrix< real > N0 ;
                Matrix< real > N1 ;
                tFun->N( xi0, N0 );
                tFun->N( xi1, N1 );

                for ( uint i = 0; i < tNumNodes; ++i )
                {
                    real tFD  = ( N1( 0, i ) - N0( 0, i ) ) / ( 2.0 * tDxi );
                    real tAna = dNdXi( d, i );
                    EXPECT_NEAR( tAna, tFD, tEpsilonFD )
                        << "dNdXi(" << d << "," << i << ") at point " << k ;
                }
            }
        }

        // --- Test 5: second derivative via central difference on dNdXi ---

        Matrix< real > d2NdXi2 ;

        for ( uint k = 0; k < tNumIntPoints; ++k )
        {
            tFun->d2NdXi2( XiInt.col( k ), d2NdXi2 );

            // evaluate dNdXi at perturbed positions
            Cell< Matrix< real > > dN_minus ;
            Cell< Matrix< real > > dN_plus ;
            dN_minus.set_size( tNumDim );
            dN_plus.set_size( tNumDim );

            for ( uint d = 0; d < tNumDim; ++d )
            {
                Vector< real > xi_m( XiInt.col( k ) );
                Vector< real > xi_p( XiInt.col( k ) );
                xi_m( d ) -= tDxi ;
                xi_p( d ) += tDxi ;
                tFun->dNdXi( xi_m, dN_minus( d ) );
                tFun->dNdXi( xi_p, dN_plus( d ) );
            }

            // diagonal: d2N/dxi_d^2  (Voigt rows 0..tNumDim-1)
            for ( uint d = 0; d < tNumDim; ++d )
            {
                for ( uint i = 0; i < tNumNodes; ++i )
                {
                    real tFD = ( dN_plus( d )( d, i ) - dN_minus( d )( d, i ) )
                               / ( 2.0 * tDxi );
                    EXPECT_NEAR( d2NdXi2( d, i ), tFD, tEpsilonFD )
                        << "d2NdXi2(" << d << "," << i << ") at point " << k ;
                }
            }

            // cross: d2N/(dxi*deta)  (Voigt 2 in 2D, 5 in 3D)
            if ( tNumDim >= 2 )
            {
                uint tVoigt = ( tNumDim == 2 ) ? 2 : 5 ;
                for ( uint i = 0; i < tNumNodes; ++i )
                {
                    real tFD = 0.5 * (
                        ( dN_plus( 0 )( 1, i ) - dN_minus( 0 )( 1, i ) )
                      + ( dN_plus( 1 )( 0, i ) - dN_minus( 1 )( 0, i ) )
                    ) / ( 2.0 * tDxi );
                    EXPECT_NEAR( d2NdXi2( tVoigt, i ), tFD, tEpsilonFD )
                        << "d2N/dxideta(" << i << ") at point " << k ;
                }
            }

            if ( tNumDim == 3 )
            {
                // d2N/(deta*dzeta)  (Voigt 3)
                for ( uint i = 0; i < tNumNodes; ++i )
                {
                    real tFD = 0.5 * (
                        ( dN_plus( 1 )( 2, i ) - dN_minus( 1 )( 2, i ) )
                      + ( dN_plus( 2 )( 1, i ) - dN_minus( 2 )( 1, i ) )
                    ) / ( 2.0 * tDxi );
                    EXPECT_NEAR( d2NdXi2( 3, i ), tFD, tEpsilonFD )
                        << "d2N/detadzeta(" << i << ") at point " << k ;
                }

                // d2N/(dxi*dzeta)  (Voigt 4)
                for ( uint i = 0; i < tNumNodes; ++i )
                {
                    real tFD = 0.5 * (
                        ( dN_plus( 0 )( 2, i ) - dN_minus( 0 )( 2, i ) )
                      + ( dN_plus( 2 )( 0, i ) - dN_minus( 2 )( 0, i ) )
                    ) / ( 2.0 * tDxi );
                    EXPECT_NEAR( d2NdXi2( 4, i ), tFD, tEpsilonFD )
                        << "d2N/dxidzeta(" << i << ") at point " << k ;
                }
            }
        }

        delete tFun ;
        delete tRef ;
    }
}

// ---------------------------------------------------------------------------
//  1D elements
// ---------------------------------------------------------------------------

TEST( LagrangeInterpolation, LINE2 )  { test_lagrange_interpolation( ElementType::LINE2 ); }
TEST( LagrangeInterpolation, LINE3 )  { test_lagrange_interpolation( ElementType::LINE3 ); }
TEST( LagrangeInterpolation, LINE4 )  { test_lagrange_interpolation( ElementType::LINE4 ); }

// ---------------------------------------------------------------------------
//  2D triangles
// ---------------------------------------------------------------------------

TEST( LagrangeInterpolation, TRI3  )  { test_lagrange_interpolation( ElementType::TRI3  ); }
TEST( LagrangeInterpolation, TRI6  )  { test_lagrange_interpolation( ElementType::TRI6  ); }
TEST( LagrangeInterpolation, TRI10 )  { test_lagrange_interpolation( ElementType::TRI10 ); }

// ---------------------------------------------------------------------------
//  2D quadrilaterals
// ---------------------------------------------------------------------------

TEST( LagrangeInterpolation, QUAD4  ) { test_lagrange_interpolation( ElementType::QUAD4  ); }
TEST( LagrangeInterpolation, QUAD8  ) { test_lagrange_interpolation( ElementType::QUAD8  ); }
TEST( LagrangeInterpolation, QUAD9  ) { test_lagrange_interpolation( ElementType::QUAD9  ); }
TEST( LagrangeInterpolation, QUAD16 ) { test_lagrange_interpolation( ElementType::QUAD16 ); }

// ---------------------------------------------------------------------------
//  Tetrahedra
// ---------------------------------------------------------------------------

TEST( LagrangeInterpolation, TET4  )  { test_lagrange_interpolation( ElementType::TET4  ); }
TEST( LagrangeInterpolation, TET10 )  { test_lagrange_interpolation( ElementType::TET10 ); }
TEST( LagrangeInterpolation, TET20 )  { test_lagrange_interpolation( ElementType::TET20 ); }
TEST( LagrangeInterpolation, TET35 )  { test_lagrange_interpolation( ElementType::TET35 ); }

// ---------------------------------------------------------------------------
//  Hexahedra
// ---------------------------------------------------------------------------

TEST( LagrangeInterpolation, HEX8  )  { test_lagrange_interpolation( ElementType::HEX8  ); }
TEST( LagrangeInterpolation, HEX20 )  { test_lagrange_interpolation( ElementType::HEX20 ); }
TEST( LagrangeInterpolation, HEX27 )  { test_lagrange_interpolation( ElementType::HEX27 ); }
TEST( LagrangeInterpolation, HEX64 )  { test_lagrange_interpolation( ElementType::HEX64 ); }

// ---------------------------------------------------------------------------
//  Pentahedra (wedges)
// ---------------------------------------------------------------------------

TEST( LagrangeInterpolation, PENTA6  ) { test_lagrange_interpolation( ElementType::PENTA6  ); }
TEST( LagrangeInterpolation, PENTA15 ) { test_lagrange_interpolation( ElementType::PENTA15 ); }
TEST( LagrangeInterpolation, PENTA18 ) { test_lagrange_interpolation( ElementType::PENTA18 ); }

// ---------------------------------------------------------------------------
//  Pyramids
// ---------------------------------------------------------------------------

TEST( LagrangeInterpolation, PYRA5  )  { test_lagrange_interpolation( ElementType::PYRA5  ); }
TEST( LagrangeInterpolation, PYRA13 )  { test_lagrange_interpolation( ElementType::PYRA13 ); }
TEST( LagrangeInterpolation, PYRA14 )  { test_lagrange_interpolation( ElementType::PYRA14 ); }
