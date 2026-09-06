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
#include "fn_det.hpp"
#include "fn_trans.hpp"

using namespace belfem;

namespace
{
    const real tEpsilon = 1.0e-12 ;

    /**
     * Verify that the reference element has unit volume and that its
     * centroid is at the origin by numerical integration using Gauss
     * quadrature with the Lagrange shape functions.
     */
    void
    test_volume_and_centroid( const ElementType aType )
    {
        mesh::ElementFactory              tElFactory ;
        fem::InterpolationFunctionFactory tIfFactory ;

        auto tRef = tElFactory.create_reference_element( aType );
        fem::InterpolationFunction * tFun = tIfFactory.create_lagrange_function( aType );

        uint tNumNodes = tRef->element()->number_of_nodes();
        uint tNumDim   = mesh::dimension( aType );

        // physical node coordinates (nodes x dim)
        Matrix< real > X( tNumNodes, tNumDim );
        for ( uint i = 0; i < tNumNodes; ++i )
        {
            for ( uint j = 0; j < tNumDim; ++j )
            {
                X( i, j ) = tRef->element()->node( i )->x( j );
            }
        }

        // integration points
        Vector< real > w ;
        Matrix< real > Xi ;
        intpoints( IntegrationScheme::GAUSS,
                   mesh::geometry_type( tRef->element()->type() ),
                   0, w, Xi );

        uint tNumPoints = w.length();

        real V = 0.0 ;
        Vector< real > Xs( tNumDim, 0.0 );

        Matrix< real > J ;
        Matrix< real > N ;
        Matrix< real > dNdxi ;
        Matrix< real > Xk ;

        for ( uint k = 0; k < tNumPoints; ++k )
        {
            tFun->N( Xi.col( k ), N );
            tFun->dNdXi( Xi.col( k ), dNdxi );

            J = dNdxi * X ;
            real tDetJ = det( J );

            V += w( k ) * tDetJ ;

            Xk = trans( N * X );
            Xs += w( k ) * Xk.col( 0 ) * tDetJ ;
        }

        Xs /= V ;

        EXPECT_NEAR( V, 1.0, tEpsilon ) << "volume" ;

        for ( uint j = 0; j < tNumDim; ++j )
        {
            EXPECT_NEAR( Xs( j ), 0.0, tEpsilon )
                << "centroid component " << j ;
        }

        delete tFun ;
        delete tRef ;
    }
}

// ---------------------------------------------------------------------------
//  2D triangles
// ---------------------------------------------------------------------------

TEST( Integration, TRI3  )  { test_volume_and_centroid( ElementType::TRI3  ); }
TEST( Integration, TRI6  )  { test_volume_and_centroid( ElementType::TRI6  ); }
TEST( Integration, TRI10 )  { test_volume_and_centroid( ElementType::TRI10 ); }
TEST( Integration, TRI15 )  { test_volume_and_centroid( ElementType::TRI15 ); }

// ---------------------------------------------------------------------------
//  2D quadrilaterals
// ---------------------------------------------------------------------------

TEST( Integration, QUAD4  ) { test_volume_and_centroid( ElementType::QUAD4  ); }
TEST( Integration, QUAD8  ) { test_volume_and_centroid( ElementType::QUAD8  ); }
TEST( Integration, QUAD9  ) { test_volume_and_centroid( ElementType::QUAD9  ); }
TEST( Integration, QUAD16 ) { test_volume_and_centroid( ElementType::QUAD16 ); }

// ---------------------------------------------------------------------------
//  Tetrahedra
// ---------------------------------------------------------------------------

TEST( Integration, TET4  )  { test_volume_and_centroid( ElementType::TET4  ); }
TEST( Integration, TET10 )  { test_volume_and_centroid( ElementType::TET10 ); }
TEST( Integration, TET20 )  { test_volume_and_centroid( ElementType::TET20 ); }
TEST( Integration, TET35 )  { test_volume_and_centroid( ElementType::TET35 ); }

// ---------------------------------------------------------------------------
//  Hexahedra
// ---------------------------------------------------------------------------

TEST( Integration, HEX8  )  { test_volume_and_centroid( ElementType::HEX8  ); }
TEST( Integration, HEX20 )  { test_volume_and_centroid( ElementType::HEX20 ); }
TEST( Integration, HEX27 )  { test_volume_and_centroid( ElementType::HEX27 ); }
TEST( Integration, HEX64 )  { test_volume_and_centroid( ElementType::HEX64 ); }

// ---------------------------------------------------------------------------
//  Pentahedra (wedges)
// ---------------------------------------------------------------------------

TEST( Integration, PENTA6  ) { test_volume_and_centroid( ElementType::PENTA6  ); }
TEST( Integration, PENTA15 ) { test_volume_and_centroid( ElementType::PENTA15 ); }
TEST( Integration, PENTA18 ) { test_volume_and_centroid( ElementType::PENTA18 ); }

// ---------------------------------------------------------------------------
//  Pyramids
// ---------------------------------------------------------------------------

TEST( Integration, PYRA5  )  { test_volume_and_centroid( ElementType::PYRA5  ); }
TEST( Integration, PYRA13 )  { test_volume_and_centroid( ElementType::PYRA13 ); }
TEST( Integration, PYRA14 )  { test_volume_and_centroid( ElementType::PYRA14 ); }
