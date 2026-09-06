/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 *
 * Circulation and curl battery for the volume Nedelec elements
 * ( TRI3, TRI6, TET4, TET10 ). See todo/falsification_tooling.md ( D1 ).
 *
 * This battery pins the class of defect found 2026-08-14 ( EF_TET4 edge 2,
 * EF_TET10 table-wide eta/zeta exchange, todo/nedelec_edge_function_defects.md ):
 * a wrong gradient or a wrong scalar factor breaks the circulation identity
 * within machine precision, on the reference element and on a distorted one.
 *
 * Conventions pinned here ( independent of the code under test ):
 *  - the EXODUS TET node map carries the booby trap: node 1 holds zeta and
 *    node 2 holds eta; the naive triangle-extension map has detJ < 0
 *    ( see src/fem/interpolation/doc/nedelec_derivation.md, section 4 )
 *  - edge k runs between the EXODUS corner pairs
 *    TRI { 01, 12, 20 },  TET { 01, 12, 20, 03, 13, 23 }
 *  - expected circulations along the canonical ( unflipped ) edge direction:
 *    first order 1 per own edge, TRI6 pair ( 1/2, 1/2 ), TET10 pair ( 1, 1 ),
 *    zero on every foreign edge, zero for every face dof
 *    ( TRI6 value derived symbolically 2026-08-14; TET10 values verified
 *    against the corrected tables and the notes' polynomial set )
 *
 * The distorted-element coordinates are the exact vertices used by the
 * symbolic probes ( tmp/tet10/tet4_circulation_probe.py and friends ), so the
 * two verification chains meet on the same numbers.
 *
 * Not covered here: LINE3 ( needs the thin-shell facet chain, see
 * test_InterfaceOrientation.cpp ) and the TS/TB family ( ditto ).
 */

#include <gtest/gtest.h>
#include <cmath>
#include <iostream>

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "cl_Cell.hpp"
#include "support/cl_EF_TestVolume.hpp"
#include "cl_EdgeFunctionFactory.hpp"

using namespace belfem;

namespace
{
    const real tEpsilon   = 1.0e-11 ;  // exact identities
    const real tEpsilonFD = 1.0e-8 ;   // central differences ( exact for
                                       // polynomials of degree <= 2, so this
                                       // absorbs roundoff only )
    const real tDxi       = 1.0e-4 ;   // finite difference step

    /**
     * everything the battery needs to know about one element type
     */
    struct ElementSpec
    {
        ElementType    mType ;
        uint           mDim ;
        uint           mNumCorners ;
        uint           mNumEdges ;
        uint           mDofsPerEdge ;
        uint           mNumFaceDofs ;
        real           mOwnEdgeCirculation ;
    };

    ElementSpec
    spec( const ElementType aType )
    {
        switch ( aType )
        {
            case ElementType::TRI3  : return { aType, 2, 3, 3, 1, 0, 1.0 };
            case ElementType::TRI6  : return { aType, 2, 3, 3, 2, 2, 0.5 };
            case ElementType::TET4  : return { aType, 3, 4, 6, 1, 0, 1.0 };
            case ElementType::TET10 : return { aType, 3, 4, 6, 2, 8, 1.0 };
            default :
            {
                BELFEM_ERROR( false, "unsupported element type" );
                return { aType, 0, 0, 0, 0, 0, 0.0 };
            }
        }
    }

    /**
     * pinned parameter coordinates of the corner nodes ( dim x corners ).
     * TET: node 1 carries zeta, node 2 carries eta — THE trap.
     */
    Matrix< real >
    param_corners( const ElementSpec & aSpec )
    {
        Matrix< real > aXi( aSpec.mDim, aSpec.mNumCorners, 0.0 );
        if ( aSpec.mDim == 2 )
        {
            aXi( 0, 0 ) = 1.0 ;   // node 0 : xi
            aXi( 1, 1 ) = 1.0 ;   // node 1 : eta
        }
        else
        {
            aXi( 0, 0 ) = 1.0 ;   // node 0 : xi
            aXi( 2, 1 ) = 1.0 ;   // node 1 : ZETA
            aXi( 1, 2 ) = 1.0 ;   // node 2 : ETA
        }
        return aXi ;
    }

    const uint gEdgeCorners[ 6 ][ 2 ] =
            { { 0, 1 }, { 1, 2 }, { 2, 0 }, { 0, 3 }, { 1, 3 }, { 2, 3 } };

    /**
     * reference corners: physical == parameter coordinates, J = I
     */
    Matrix< real >
    reference_corners( const ElementSpec & aSpec )
    {
        return param_corners( aSpec );
    }

    /**
     * distorted corners, identical to the symbolic probes
     */
    Matrix< real >
    distorted_corners( const ElementSpec & aSpec )
    {
        Matrix< real > aX( aSpec.mDim, aSpec.mNumCorners );
        if ( aSpec.mDim == 2 )
        {
            aX( 0, 0 ) = 1.0 ; aX( 1, 0 ) = 1.0 ;
            aX( 0, 1 ) = 4.0 ; aX( 1, 1 ) = 2.0 ;
            aX( 0, 2 ) = 2.0 ; aX( 1, 2 ) = 4.0 ;
        }
        else
        {
            aX( 0, 0 ) = 1.0 ; aX( 1, 0 ) = 2.0 ; aX( 2, 0 ) = -1.0 ;
            aX( 0, 1 ) = 4.0 ; aX( 1, 1 ) = 1.0 ; aX( 2, 1 ) =  0.0 ;
            aX( 0, 2 ) = 2.0 ; aX( 1, 2 ) = 5.0 ; aX( 2, 2 ) =  1.0 ;
            aX( 0, 3 ) = 2.0 ; aX( 1, 3 ) = 2.0 ; aX( 2, 3 ) =  3.0 ;
        }
        return aX ;
    }

    /**
     * evaluation points. Column 0 is the INTERIOR base point on purpose:
     * link() evaluates its Jacobian at point 0, and a boundary point there
     * ( eta = 0 on the first edge ) masked the eta/zeta shape-derivative
     * defect D3 of 2026-08-14. Then the finite difference cluster
     * ( minus/plus per parameter direction ), then 3-point Gauss per edge.
     */
    Matrix< real >
    evaluation_points( const ElementSpec & aSpec )
    {
        const real tS = 0.5 * std::sqrt( 0.6 );
        const real tT[ 3 ] = { 0.5 - tS, 0.5, 0.5 + tS };

        Matrix< real > tXiCorner = param_corners( aSpec );

        const uint tNumCols = 1 + 2 * aSpec.mDim + 3 * aSpec.mNumEdges ;
        Matrix< real > aXi( aSpec.mDim, tNumCols );

        // interior base point, column 0
        const real tBase[ 3 ] = { 0.3, 0.25, 0.2 };
        uint tCount = 0 ;
        for ( uint d = 0; d < aSpec.mDim; ++d )
        {
            aXi( d, tCount ) = tBase[ d ];
        }
        ++tCount ;

        // finite difference cluster
        for ( uint d = 0; d < aSpec.mDim; ++d )
        {
            for ( uint p = 0; p < 2; ++p )
            {
                for ( uint i = 0; i < aSpec.mDim; ++i )
                {
                    aXi( i, tCount ) = tBase[ i ];
                }
                aXi( d, tCount ) += ( p == 0 ? -tDxi : tDxi );
                ++tCount ;
            }
        }

        // edge Gauss points
        for ( uint e = 0; e < aSpec.mNumEdges; ++e )
        {
            const uint a = gEdgeCorners[ e ][ 0 ];
            const uint b = gEdgeCorners[ e ][ 1 ];
            for ( uint g = 0; g < 3; ++g )
            {
                for ( uint d = 0; d < aSpec.mDim; ++d )
                {
                    aXi( d, tCount ) = ( 1.0 - tT[ g ] ) * tXiCorner( d, a )
                                     + tT[ g ] * tXiCorner( d, b );
                }
                ++tCount ;
            }
        }

        return aXi ;
    }

    /**
     * circulation of dof column aDof along edge aEdge, canonical direction,
     * 3-point Gauss with the constant physical tangent X_b - X_a
     */
    real
    edge_circulation(
            fem::EdgeFunction    * aEF,
            const ElementSpec    & aSpec,
            const Matrix< real > & aCorners,
            const uint             aEdge,
            const uint             aDof )
    {
        const real tW[ 3 ] = { 5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0 };

        const uint a = gEdgeCorners[ aEdge ][ 0 ];
        const uint b = gEdgeCorners[ aEdge ][ 1 ];

        real tCirc = 0.0 ;
        for ( uint g = 0; g < 3; ++g )
        {
            const Matrix< real > & tE =
                    aEF->E( 1 + 2 * aSpec.mDim + 3 * aEdge + g );
            real tDot = 0.0 ;
            for ( uint d = 0; d < aSpec.mDim; ++d )
            {
                tDot += tE( d, aDof )
                        * ( aCorners( d, b ) - aCorners( d, a ) );
            }
            tCirc += tW[ g ] * tDot ;
        }
        return tCirc ;
    }

    /**
     * test-side Jacobian J( m, d ) = dx_m / dxi_d from the pinned barycentric
     * map ( TRI: xi, eta on nodes 0, 1; TET: xi, ZETA, ETA on nodes 0, 1, 2;
     * the last node carries the dependent coordinate )
     */
    void
    test_jacobian(
            const ElementSpec    & aSpec,
            const Matrix< real > & aCorners,
            Matrix< real >       & aJ,
            Matrix< real >       & aInvJ,
            real                 & aDetJ )
    {
        const uint n = aSpec.mDim ;
        aJ.set_size( n, n );
        aInvJ.set_size( n, n );

        // column d of J: corner carrying parameter d minus the last corner
        // 2D: ( n0 - n2, n1 - n2 ) ; 3D: ( n0 - n3, n2 - n3, n1 - n3 )
        const uint tCarrier2D[ 2 ] = { 0, 1 };
        const uint tCarrier3D[ 3 ] = { 0, 2, 1 };
        const uint * tCarrier = ( n == 2 ) ? tCarrier2D : tCarrier3D ;
        const uint tLast = aSpec.mNumCorners - 1 ;

        for ( uint d = 0; d < n; ++d )
        {
            for ( uint m = 0; m < n; ++m )
            {
                aJ( m, d ) = aCorners( m, tCarrier[ d ] )
                           - aCorners( m, tLast );
            }
        }

        if ( n == 2 )
        {
            aDetJ = aJ( 0, 0 ) * aJ( 1, 1 ) - aJ( 0, 1 ) * aJ( 1, 0 );
            aInvJ( 0, 0 ) =  aJ( 1, 1 ) / aDetJ ;
            aInvJ( 0, 1 ) = -aJ( 0, 1 ) / aDetJ ;
            aInvJ( 1, 0 ) = -aJ( 1, 0 ) / aDetJ ;
            aInvJ( 1, 1 ) =  aJ( 0, 0 ) / aDetJ ;
        }
        else
        {
            aDetJ = aJ( 0, 0 ) * ( aJ( 1, 1 ) * aJ( 2, 2 ) - aJ( 1, 2 ) * aJ( 2, 1 ) )
                  - aJ( 0, 1 ) * ( aJ( 1, 0 ) * aJ( 2, 2 ) - aJ( 1, 2 ) * aJ( 2, 0 ) )
                  + aJ( 0, 2 ) * ( aJ( 1, 0 ) * aJ( 2, 1 ) - aJ( 1, 1 ) * aJ( 2, 0 ) );

            aInvJ( 0, 0 ) = ( aJ( 1, 1 ) * aJ( 2, 2 ) - aJ( 1, 2 ) * aJ( 2, 1 ) ) / aDetJ ;
            aInvJ( 0, 1 ) = ( aJ( 0, 2 ) * aJ( 2, 1 ) - aJ( 0, 1 ) * aJ( 2, 2 ) ) / aDetJ ;
            aInvJ( 0, 2 ) = ( aJ( 0, 1 ) * aJ( 1, 2 ) - aJ( 0, 2 ) * aJ( 1, 1 ) ) / aDetJ ;
            aInvJ( 1, 0 ) = ( aJ( 1, 2 ) * aJ( 2, 0 ) - aJ( 1, 0 ) * aJ( 2, 2 ) ) / aDetJ ;
            aInvJ( 1, 1 ) = ( aJ( 0, 0 ) * aJ( 2, 2 ) - aJ( 0, 2 ) * aJ( 2, 0 ) ) / aDetJ ;
            aInvJ( 1, 2 ) = ( aJ( 0, 2 ) * aJ( 1, 0 ) - aJ( 0, 0 ) * aJ( 1, 2 ) ) / aDetJ ;
            aInvJ( 2, 0 ) = ( aJ( 1, 0 ) * aJ( 2, 1 ) - aJ( 1, 1 ) * aJ( 2, 0 ) ) / aDetJ ;
            aInvJ( 2, 1 ) = ( aJ( 0, 1 ) * aJ( 2, 0 ) - aJ( 0, 0 ) * aJ( 2, 1 ) ) / aDetJ ;
            aInvJ( 2, 2 ) = ( aJ( 0, 0 ) * aJ( 1, 1 ) - aJ( 0, 1 ) * aJ( 1, 0 ) ) / aDetJ ;
        }
    }

    /**
     * link an edge function to a fresh fixture and return both;
     * caller owns the edge function
     */
    fem::EdgeFunction *
    make_linked_ef(
            const ElementSpec       & aSpec,
            fem::test::EF_TestVolume & aFixture )
    {
        fem::EdgeFunctionFactory tFactory ;
        fem::EdgeFunction * aEF = tFactory.create_edge_function( aSpec.mType );
        aEF->precompute( evaluation_points( aSpec ) );
        aEF->link( aFixture.element() );
        return aEF ;
    }

//------------------------------------------------------------------------------

    /**
     * Test 1: circulation identity. Every dof against every edge:
     * own-edge value for edge dofs, zero for everything else.
     * This is the test that catches both 2026-08-14 defects.
     */
    void
    test_circulation(
            const ElementSpec    & aSpec,
            const Matrix< real > & aCorners )
    {
        fem::test::EF_TestVolume tFixture( aSpec.mType, aCorners );
        fem::EdgeFunction * tEF = make_linked_ef( aSpec, tFixture );

        const uint tNumDofs = aSpec.mNumEdges * aSpec.mDofsPerEdge
                            + aSpec.mNumFaceDofs ;

        for ( uint e = 0; e < aSpec.mNumEdges; ++e )
        {
            for ( uint k = 0; k < tNumDofs; ++k )
            {
                const bool tOwn =
                        k < aSpec.mNumEdges * aSpec.mDofsPerEdge
                        && k / aSpec.mDofsPerEdge == e ;

                const real tExpect = tOwn ? aSpec.mOwnEdgeCirculation : 0.0 ;

                EXPECT_NEAR(
                        edge_circulation( tEF, aSpec, aCorners, e, k ),
                        tExpect, tEpsilon )
                        << "circulation of dof " << k << " on edge " << e ;
            }
        }

        delete tEF ;
    }

//------------------------------------------------------------------------------

    /**
     * Test 2: the C operator equals the curl of E, checked against central
     * differences of E in physical space ( exact for these polynomial
     * degrees ). Catches E-vs-C inconsistencies like the TET4 defect from
     * the other side.
     */
    void
    test_curl(
            const ElementSpec    & aSpec,
            const Matrix< real > & aCorners,
            const bool             aCurved = false )
    {
        fem::test::EF_TestVolume tFixture( aSpec.mType, aCorners, -1, aCurved );
        fem::EdgeFunction * tEF = make_linked_ef( aSpec, tFixture );

        Matrix< real > tJ ;
        Matrix< real > tInvJ ;
        real tDetJ ;
        test_jacobian( aSpec, aCorners, tJ, tInvJ, tDetJ );

        const uint tNumDofs = aSpec.mNumEdges * aSpec.mDofsPerEdge
                            + aSpec.mNumFaceDofs ;
        const uint tBase = 0 ;

        // dE_c / dx_m = sum_d dE_c / dxi_d * invJ( d, m )
        // dE[ d ] = ( E( plus_d ) - E( minus_d ) ) / ( 2 dxi ), evaluated
        // before use because E() returns a reference to a work matrix
        Cell< Matrix< real > > tDE ;
        tDE.set_size( aSpec.mDim, {} );
        for ( uint d = 0; d < aSpec.mDim; ++d )
        {
            Matrix< real > tMinus( tEF->E( 1 + 2 * d ) );
            const Matrix< real > & tPlus = tEF->E( 2 + 2 * d );
            tDE( d ).set_size( aSpec.mDim, tNumDofs );
            for ( uint c = 0; c < aSpec.mDim; ++c )
            {
                for ( uint k = 0; k < tNumDofs; ++k )
                {
                    tDE( d )( c, k ) =
                            ( tPlus( c, k ) - tMinus( c, k ) )
                            / ( 2.0 * tDxi );
                }
            }
        }

        const Matrix< real > & tC = tEF->C( tBase );

        for ( uint k = 0; k < tNumDofs; ++k )
        {
            // physical gradient g( c, m ) = dE_c / dx_m
            real g[ 3 ][ 3 ] = { { 0.0 } };
            for ( uint c = 0; c < aSpec.mDim; ++c )
            {
                for ( uint m = 0; m < aSpec.mDim; ++m )
                {
                    for ( uint d = 0; d < aSpec.mDim; ++d )
                    {
                        g[ c ][ m ] += tDE( d )( c, k ) * tInvJ( d, m );
                    }
                }
            }

            if ( aSpec.mDim == 2 )
            {
                // scalar curl: dEy/dx - dEx/dy
                EXPECT_NEAR( tC( 0, k ), g[ 1 ][ 0 ] - g[ 0 ][ 1 ], tEpsilonFD )
                        << "curl of dof " << k ;
            }
            else
            {
                EXPECT_NEAR( tC( 0, k ), g[ 2 ][ 1 ] - g[ 1 ][ 2 ], tEpsilonFD )
                        << "curl x of dof " << k ;
                EXPECT_NEAR( tC( 1, k ), g[ 0 ][ 2 ] - g[ 2 ][ 0 ], tEpsilonFD )
                        << "curl y of dof " << k ;
                EXPECT_NEAR( tC( 2, k ), g[ 1 ][ 0 ] - g[ 0 ][ 1 ], tEpsilonFD )
                        << "curl z of dof " << k ;
            }
        }

        delete tEF ;
    }

//------------------------------------------------------------------------------

    /**
     * Test 3: EXODUS orientation gives a positive Jacobian ( the node-map
     * booby trap turns this negative ), and the edge function agrees with
     * the test-side determinant.
     */
    void
    test_det_j(
            const ElementSpec    & aSpec,
            const Matrix< real > & aCorners )
    {
        fem::test::EF_TestVolume tFixture( aSpec.mType, aCorners );
        fem::EdgeFunction * tEF = make_linked_ef( aSpec, tFixture );

        Matrix< real > tJ ;
        Matrix< real > tInvJ ;
        real tDetJ ;
        test_jacobian( aSpec, aCorners, tJ, tInvJ, tDetJ );

        EXPECT_GT( tDetJ, 0.0 ) << "test-side Jacobian must be positive" ;
        EXPECT_GT( tEF->det_J(), 0.0 ) << "EXODUS orientation gives detJ > 0" ;
        EXPECT_NEAR( tEF->det_J(), tDetJ, tEpsilon * std::abs( tDetJ ) )
                << "det_J against the pinned node map" ;

        delete tEF ;
    }

//------------------------------------------------------------------------------

    /**
     * Test 4: flipping the stored corner order of one mesh edge negates
     * exactly that edge's dof column(s) and nothing else. ( In the global
     * edge direction the circulation is still positive; the sign lives in
     * the local view. The swap of the two SECOND-ORDER dofs of a negative
     * edge happens at dof-numbering level, not inside the edge function. )
     */
    void
    test_edge_flip(
            const ElementSpec    & aSpec,
            const Matrix< real > & aCorners )
    {
        const uint tFlip = 1 ;   // an edge that exists for TRI and TET

        fem::test::EF_TestVolume tFixture( aSpec.mType, aCorners, tFlip );
        fem::EdgeFunction * tEF = make_linked_ef( aSpec, tFixture );

        for ( uint e = 0; e < aSpec.mNumEdges; ++e )
        {
            for ( uint i = 0; i < aSpec.mDofsPerEdge; ++i )
            {
                const uint k = e * aSpec.mDofsPerEdge + i ;
                const real tExpect = ( e == tFlip )
                        ? -aSpec.mOwnEdgeCirculation
                        :  aSpec.mOwnEdgeCirculation ;

                EXPECT_NEAR(
                        edge_circulation( tEF, aSpec, aCorners, e, k ),
                        tExpect, tEpsilon )
                        << "own-edge circulation of dof " << k
                        << " with edge " << tFlip << " flipped" ;
            }
        }

        delete tEF ;
    }

//------------------------------------------------------------------------------

    /**
     * Test 5 ( quadratic elements ): the curved-path evaluation must
     * reproduce the straight-path E and C on straight geometry.
     */
    void
    test_curved_path(
            const ElementSpec    & aSpec,
            const Matrix< real > & aCorners )
    {
        fem::test::EF_TestVolume tStraight( aSpec.mType, aCorners, -1, false );
        fem::test::EF_TestVolume tCurved(   aSpec.mType, aCorners, -1, true );

        fem::EdgeFunction * tEFS = make_linked_ef( aSpec, tStraight );
        fem::EdgeFunction * tEFC = make_linked_ef( aSpec, tCurved );

        const uint tNumDofs = aSpec.mNumEdges * aSpec.mDofsPerEdge
                            + aSpec.mNumFaceDofs ;
        const uint tBase = 0 ;

        Matrix< real > tES( tEFS->E( tBase ) );
        const Matrix< real > & tEC = tEFC->E( tBase );
        for ( uint c = 0; c < aSpec.mDim; ++c )
        {
            for ( uint k = 0; k < tNumDofs; ++k )
            {
                EXPECT_NEAR( tES( c, k ), tEC( c, k ), tEpsilon )
                        << "E straight vs curved path, entry ("
                        << c << "," << k << ")" ;
            }
        }

        Matrix< real > tCS( tEFS->C( tBase ) );
        const Matrix< real > & tCC = tEFC->C( tBase );
        const uint tNumCurl = ( aSpec.mDim == 2 ) ? 1 : 3 ;
        for ( uint c = 0; c < tNumCurl; ++c )
        {
            for ( uint k = 0; k < tNumDofs; ++k )
            {
                EXPECT_NEAR( tCS( c, k ), tCC( c, k ), tEpsilon )
                        << "C straight vs curved path, entry ("
                        << c << "," << k << ")" ;
            }
        }

        delete tEFS ;
        delete tEFC ;
    }

//------------------------------------------------------------------------------

    /**
     * Test 6 ( G-operator, first-order elements ): the gradient operator
     * against the layout contract documented at EdgeFunction::mGrad.
     *   1. FD      : G equals the physical gradient of E, via the same
     *                central-difference stencil as test_curl ( exact for
     *                polynomials of degree <= 2, roundoff-limited here )
     *   2. curl tie: the antisymmetric rows of G reproduce C to round-off
     *   3. Frobenius: ||G_e||^2 = 0.5 * ||C_e||^2 per dof ( antisymmetry;
     *                catches a wrong scalar factor, transpose-invariant, so
     *                the curl tie stays the first gate )
     *   4. trace   : sum_i G( i + d*i, e ) = 0 ( Whitney divergence-free —
     *                weak falsifier, kept as the elementwise-div pin )
     *   5. gradient mode: for the edge dofs of the interpolant of
     *                grad( phi ), phi = x + 2y ( + 3z ), G * q is the zero
     *                vector — the curl-null-space blindness at class level.
     *                Signs come from element()->edge_directions.
     *   6. shape   : n_rows == d*d, n_cols == ndofs ( a return-mG typo
     *                returns the barycentric coefficient table instead )
     * aFlipEdge exercises the mS path: on the default fixtures every edge
     * direction is positive, so a missing mS[ e ] in the fill would pass
     * every check unflipped. On the flipped fixture the discriminators are
     * FD, curl tie, and gradient mode — Frobenius and trace are
     * sign-invariant and stay green either way.
     */
    void
    test_gradient(
            const ElementSpec    & aSpec,
            const Matrix< real > & aCorners,
            const int              aFlipEdge = -1 )
    {
        fem::test::EF_TestVolume tFixture( aSpec.mType, aCorners, aFlipEdge );
        fem::EdgeFunction * tEF = make_linked_ef( aSpec, tFixture );

        Matrix< real > tJ ;
        Matrix< real > tInvJ ;
        real tDetJ ;
        test_jacobian( aSpec, aCorners, tJ, tInvJ, tDetJ );

        const uint d = aSpec.mDim ;
        const uint tNumDofs = aSpec.mNumEdges * aSpec.mDofsPerEdge
                            + aSpec.mNumFaceDofs ;
        const uint tBase = 0 ;

        const Matrix< real > & tG = tEF->G( tBase );

        // 6. shape ( ASSERT: a return-mG typo would return the coefficient
        //    table, and the reads below would then be out of bounds )
        ASSERT_EQ( tG.n_rows(), d * d ) << "G row count" ;
        ASSERT_EQ( tG.n_cols(), tNumDofs ) << "G column count" ;

        // 1. FD of E, same stencil as test_curl; E() aliases mE, so the
        //    minus side must be copied before the plus side is evaluated
        Cell< Matrix< real > > tDE ;
        tDE.set_size( d, {} );
        for ( uint p = 0; p < d; ++p )
        {
            Matrix< real > tMinus( tEF->E( 1 + 2 * p ) );
            const Matrix< real > & tPlus = tEF->E( 2 + 2 * p );
            tDE( p ).set_size( d, tNumDofs );
            for ( uint c = 0; c < d; ++c )
            {
                for ( uint k = 0; k < tNumDofs; ++k )
                {
                    tDE( p )( c, k ) = ( tPlus( c, k ) - tMinus( c, k ) )
                                     / ( 2.0 * tDxi );
                }
            }
        }

        real tMaxFD = 0.0 ;
        for ( uint k = 0; k < tNumDofs; ++k )
        {
            for ( uint c = 0; c < d; ++c )        // field component
            {
                for ( uint m = 0; m < d; ++m )    // derivative direction
                {
                    real tFD = 0.0 ;
                    for ( uint p = 0; p < d; ++p )
                    {
                        tFD += tDE( p )( c, k ) * tInvJ( p, m );
                    }
                    EXPECT_NEAR( tG( m + d * c, k ), tFD, tEpsilonFD )
                            << "G vs FD, derivative " << m
                            << " of component " << c << ", dof " << k ;
                    const real tDev = std::abs( tG( m + d * c, k ) - tFD );
                    if ( tDev > tMaxFD ) { tMaxFD = tDev ; }
                }
            }
        }

        // 2. curl tie ( round-off epsilon, not the FD one )
        const Matrix< real > & tC = tEF->C( tBase );
        real tMaxTie = 0.0 ;
        for ( uint k = 0; k < tNumDofs; ++k )
        {
            real tDev ;
            if ( d == 2 )
            {
                tDev = std::abs( tC( 0, k ) - ( tG( 2, k ) - tG( 1, k ) ) );
                if ( tDev > tMaxTie ) { tMaxTie = tDev ; }
            }
            else
            {
                tDev = std::abs( tC( 0, k ) - ( tG( 7, k ) - tG( 5, k ) ) );
                if ( tDev > tMaxTie ) { tMaxTie = tDev ; }
                tDev = std::abs( tC( 1, k ) - ( tG( 2, k ) - tG( 6, k ) ) );
                if ( tDev > tMaxTie ) { tMaxTie = tDev ; }
                tDev = std::abs( tC( 2, k ) - ( tG( 3, k ) - tG( 1, k ) ) );
                if ( tDev > tMaxTie ) { tMaxTie = tDev ; }
            }
        }
        EXPECT_LT( tMaxTie, tEpsilon ) << "curl tie" ;

        // 3. Frobenius identity per dof
        const uint tNumCurl = ( d == 2 ) ? 1 : 3 ;
        for ( uint k = 0; k < tNumDofs; ++k )
        {
            real tG2 = 0.0 ;
            for ( uint r = 0; r < d * d; ++r )
            {
                tG2 += tG( r, k ) * tG( r, k );
            }
            real tC2 = 0.0 ;
            for ( uint c = 0; c < tNumCurl; ++c )
            {
                tC2 += tC( c, k ) * tC( c, k );
            }
            const real tScale = tC2 > 1.0 ? tC2 : 1.0 ;
            EXPECT_NEAR( tG2, 0.5 * tC2, tEpsilon * tScale )
                    << "Frobenius identity of dof " << k ;
        }

        // 4. trace ( absolute tolerance — comparing exact zeros )
        for ( uint k = 0; k < tNumDofs; ++k )
        {
            real tTr = 0.0 ;
            for ( uint i = 0; i < d; ++i )
            {
                tTr += tG( i + d * i, k );
            }
            EXPECT_NEAR( tTr, 0.0, tEpsilon ) << "trace of dof " << k ;
        }

        // 5. gradient mode ( first-order elements only: one dof per edge )
        real tS[ 12 ];
        tFixture.element()->edge_directions( tS );
        real tQ[ 12 ] = { 0.0 };
        for ( uint e = 0; e < aSpec.mNumEdges; ++e )
        {
            const uint a = gEdgeCorners[ e ][ 0 ];
            const uint b = gEdgeCorners[ e ][ 1 ];
            real tPhiA = 0.0 ;
            real tPhiB = 0.0 ;
            for ( uint c = 0; c < d; ++c )
            {
                tPhiA += ( c + 1.0 ) * aCorners( c, a );
                tPhiB += ( c + 1.0 ) * aCorners( c, b );
            }
            tQ[ e ] = tS[ e ] * ( tPhiB - tPhiA );
        }
        for ( uint r = 0; r < d * d; ++r )
        {
            real tGq = 0.0 ;
            for ( uint e = 0; e < aSpec.mNumEdges; ++e )
            {
                tGq += tG( r, e ) * tQ[ e ];
            }
            EXPECT_NEAR( tGq, 0.0, tEpsilon )
                    << "gradient-mode blindness, row " << r ;
        }

        std::cout << "    [ G battery ] max |G - FD(E)| = " << tMaxFD
                  << ", max curl-tie dev = " << tMaxTie << std::endl ;

        delete tEF ;
    }
}

//------------------------------------------------------------------------------
//  TRI3
//------------------------------------------------------------------------------

TEST( EdgeFunctions, Tri3CirculationReference )
{
    ElementSpec s = spec( ElementType::TRI3 );
    test_circulation( s, reference_corners( s ) );
}

TEST( EdgeFunctions, Tri3CirculationDistorted )
{
    ElementSpec s = spec( ElementType::TRI3 );
    test_circulation( s, distorted_corners( s ) );
}

TEST( EdgeFunctions, Tri3Curl )
{
    ElementSpec s = spec( ElementType::TRI3 );
    test_curl( s, distorted_corners( s ) );
}

TEST( EdgeFunctions, Tri3DetJ )
{
    ElementSpec s = spec( ElementType::TRI3 );
    test_det_j( s, distorted_corners( s ) );
}

TEST( EdgeFunctions, Tri3EdgeFlip )
{
    ElementSpec s = spec( ElementType::TRI3 );
    test_edge_flip( s, distorted_corners( s ) );
}

TEST( EdgeFunctions, Tri3GradientReference )
{
    ElementSpec s = spec( ElementType::TRI3 );
    test_gradient( s, reference_corners( s ) );
}

TEST( EdgeFunctions, Tri3GradientDistorted )
{
    ElementSpec s = spec( ElementType::TRI3 );
    test_gradient( s, distorted_corners( s ) );
}

TEST( EdgeFunctions, Tri3GradientDistortedFlipped )
{
    ElementSpec s = spec( ElementType::TRI3 );
    test_gradient( s, distorted_corners( s ), 1 );
}

//------------------------------------------------------------------------------
//  TRI6
//------------------------------------------------------------------------------

TEST( EdgeFunctions, Tri6CirculationReference )
{
    ElementSpec s = spec( ElementType::TRI6 );
    test_circulation( s, reference_corners( s ) );
}

TEST( EdgeFunctions, Tri6CirculationDistorted )
{
    ElementSpec s = spec( ElementType::TRI6 );
    test_circulation( s, distorted_corners( s ) );
}

TEST( EdgeFunctions, Tri6Curl )
{
    ElementSpec s = spec( ElementType::TRI6 );
    test_curl( s, distorted_corners( s ) );
}

TEST( EdgeFunctions, Tri6DetJ )
{
    ElementSpec s = spec( ElementType::TRI6 );
    test_det_j( s, distorted_corners( s ) );
}

TEST( EdgeFunctions, Tri6EdgeFlip )
{
    ElementSpec s = spec( ElementType::TRI6 );
    test_edge_flip( s, distorted_corners( s ) );
}

//------------------------------------------------------------------------------
//  TET4
//------------------------------------------------------------------------------

TEST( EdgeFunctions, Tet4CirculationReference )
{
    ElementSpec s = spec( ElementType::TET4 );
    test_circulation( s, reference_corners( s ) );
}

TEST( EdgeFunctions, Tet4CirculationDistorted )
{
    ElementSpec s = spec( ElementType::TET4 );
    test_circulation( s, distorted_corners( s ) );
}

TEST( EdgeFunctions, Tet4Curl )
{
    ElementSpec s = spec( ElementType::TET4 );
    test_curl( s, distorted_corners( s ) );
}

TEST( EdgeFunctions, Tet4DetJ )
{
    ElementSpec s = spec( ElementType::TET4 );
    test_det_j( s, distorted_corners( s ) );
}

TEST( EdgeFunctions, Tet4EdgeFlip )
{
    ElementSpec s = spec( ElementType::TET4 );
    test_edge_flip( s, distorted_corners( s ) );
}

TEST( EdgeFunctions, Tet4GradientReference )
{
    ElementSpec s = spec( ElementType::TET4 );
    test_gradient( s, reference_corners( s ) );
}

TEST( EdgeFunctions, Tet4GradientDistorted )
{
    ElementSpec s = spec( ElementType::TET4 );
    test_gradient( s, distorted_corners( s ) );
}

TEST( EdgeFunctions, Tet4GradientDistortedFlipped )
{
    ElementSpec s = spec( ElementType::TET4 );
    test_gradient( s, distorted_corners( s ), 1 );
}

//------------------------------------------------------------------------------
//  TET10
//------------------------------------------------------------------------------

TEST( EdgeFunctions, Tet10CirculationReference )
{
    ElementSpec s = spec( ElementType::TET10 );
    test_circulation( s, reference_corners( s ) );
}

TEST( EdgeFunctions, Tet10CirculationDistorted )
{
    ElementSpec s = spec( ElementType::TET10 );
    test_circulation( s, distorted_corners( s ) );
}

TEST( EdgeFunctions, Tet10Curl )
{
    ElementSpec s = spec( ElementType::TET10 );
    test_curl( s, distorted_corners( s ) );
}

TEST( EdgeFunctions, Tet10CurlCurvedPath )
{
    ElementSpec s = spec( ElementType::TET10 );
    test_curl( s, distorted_corners( s ), true );
}

TEST( EdgeFunctions, Tet10DetJ )
{
    ElementSpec s = spec( ElementType::TET10 );
    test_det_j( s, distorted_corners( s ) );
}

TEST( EdgeFunctions, Tet10EdgeFlip )
{
    ElementSpec s = spec( ElementType::TET10 );
    test_edge_flip( s, distorted_corners( s ) );
}

TEST( EdgeFunctions, Tet10CurvedPathEquivalence )
{
    ElementSpec s = spec( ElementType::TET10 );
    test_curved_path( s, distorted_corners( s ) );
}

//------------------------------------------------------------------------------
//  HEX8 ( first non-affine volume element; also the curl-sign gate )
//------------------------------------------------------------------------------

namespace
{
    // 12-edge corner table, read from cl_Element_HEX8.hpp get_nodes_of_edge:
    // bottom loop, top loop, verticals
    const uint gHexEdges[ 12 ][ 2 ] = {
            { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 },
            { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 },
            { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } };

    // reference corners of [-1,1]^3, EXODUS order
    const real gHexXi[ 8 ][ 3 ] = {
            { -1, -1, -1 }, { 1, -1, -1 }, { 1, 1, -1 }, { -1, 1, -1 },
            { -1, -1,  1 }, { 1, -1,  1 }, { 1, 1,  1 }, { -1, 1,  1 } };

    // pinned distortion offsets, O( 0.1 ), no RNG ( the truncation floor
    // and the H-controls need genuine trilinear distortion )
    const real gHexJitter[ 8 ][ 3 ] = {
            {  0.11, -0.07,  0.05 }, { -0.09,  0.12, -0.06 },
            {  0.08,  0.05,  0.13 }, { -0.12, -0.10, -0.04 },
            { -0.05,  0.09, -0.11 }, {  0.13, -0.08,  0.07 },
            { -0.06, -0.13,  0.09 }, {  0.10,  0.06, -0.12 } };

    Matrix< real >
    hex_corners( const uint aFlavor )   // 0 box, 1 sheared, 2 distorted
    {
        const real tScale[ 3 ] = { 1.2, 0.9, 1.1 };
        Matrix< real > aX( 3, 8 );
        for ( uint n = 0; n < 8; ++n )
        {
            real x = tScale[ 0 ] * gHexXi[ n ][ 0 ];
            real y = tScale[ 1 ] * gHexXi[ n ][ 1 ];
            real z = tScale[ 2 ] * gHexXi[ n ][ 2 ];
            if ( aFlavor == 1 )
            {
                // affine shear: nonzero divergence content, still H = 0
                aX( 0, n ) = x + 0.3 * y ;
                aX( 1, n ) = y + 0.2 * z ;
                aX( 2, n ) = z + 0.1 * x ;
            }
            else
            {
                aX( 0, n ) = x ;
                aX( 1, n ) = y ;
                aX( 2, n ) = z ;
                if ( aFlavor == 2 )
                {
                    aX( 0, n ) += gHexJitter[ n ][ 0 ];
                    aX( 1, n ) += gHexJitter[ n ][ 1 ];
                    aX( 2, n ) += gHexJitter[ n ][ 2 ];
                }
            }
        }
        return aX ;
    }

    // trilinear dN/dxi at a point, EXODUS order: tdN( d, n )
    void
    hex_dN( const real * aXi, real adN[ 3 ][ 8 ] )
    {
        for ( uint n = 0; n < 8; ++n )
        {
            const real a = gHexXi[ n ][ 0 ];
            const real b = gHexXi[ n ][ 1 ];
            const real c = gHexXi[ n ][ 2 ];
            adN[ 0 ][ n ] = a * ( 1. + b * aXi[ 1 ] ) * ( 1. + c * aXi[ 2 ] ) / 8. ;
            adN[ 1 ][ n ] = ( 1. + a * aXi[ 0 ] ) * b * ( 1. + c * aXi[ 2 ] ) / 8. ;
            adN[ 2 ][ n ] = ( 1. + a * aXi[ 0 ] ) * ( 1. + b * aXi[ 1 ] ) * c / 8. ;
        }
    }

    // test-side NON-transposed geometry Jacobian J( m, d ) = dx_m / dxi_d
    // at a reference point, from the TEST corners ( dim x nodes ), and its
    // cofactor inverse — never production mInvJ ( transposed convention )
    void
    hex_jacobian(
            const Matrix< real > & aX,
            const real           * aXi,
            Matrix< real >       & aInvJ )
    {
        real tdN[ 3 ][ 8 ];
        hex_dN( aXi, tdN );
        Matrix< real > tJ( 3, 3 );
        for ( uint m = 0; m < 3; ++m )
        {
            for ( uint d = 0; d < 3; ++d )
            {
                real tSum = 0.0 ;
                for ( uint n = 0; n < 8; ++n )
                {
                    tSum += tdN[ d ][ n ] * aX( m, n );
                }
                tJ( m, d ) = tSum ;
            }
        }
        const real tDet =
              tJ( 0, 0 ) * ( tJ( 1, 1 ) * tJ( 2, 2 ) - tJ( 1, 2 ) * tJ( 2, 1 ) )
            - tJ( 0, 1 ) * ( tJ( 1, 0 ) * tJ( 2, 2 ) - tJ( 1, 2 ) * tJ( 2, 0 ) )
            + tJ( 0, 2 ) * ( tJ( 1, 0 ) * tJ( 2, 1 ) - tJ( 1, 1 ) * tJ( 2, 0 ) );
        aInvJ.set_size( 3, 3 );
        aInvJ( 0, 0 ) = ( tJ( 1, 1 ) * tJ( 2, 2 ) - tJ( 1, 2 ) * tJ( 2, 1 ) ) / tDet ;
        aInvJ( 0, 1 ) = ( tJ( 0, 2 ) * tJ( 2, 1 ) - tJ( 0, 1 ) * tJ( 2, 2 ) ) / tDet ;
        aInvJ( 0, 2 ) = ( tJ( 0, 1 ) * tJ( 1, 2 ) - tJ( 0, 2 ) * tJ( 1, 1 ) ) / tDet ;
        aInvJ( 1, 0 ) = ( tJ( 1, 2 ) * tJ( 2, 0 ) - tJ( 1, 0 ) * tJ( 2, 2 ) ) / tDet ;
        aInvJ( 1, 1 ) = ( tJ( 0, 0 ) * tJ( 2, 2 ) - tJ( 0, 2 ) * tJ( 2, 0 ) ) / tDet ;
        aInvJ( 1, 2 ) = ( tJ( 0, 2 ) * tJ( 1, 0 ) - tJ( 0, 0 ) * tJ( 1, 2 ) ) / tDet ;
        aInvJ( 2, 0 ) = ( tJ( 1, 0 ) * tJ( 2, 1 ) - tJ( 1, 1 ) * tJ( 2, 0 ) ) / tDet ;
        aInvJ( 2, 1 ) = ( tJ( 0, 1 ) * tJ( 2, 0 ) - tJ( 0, 0 ) * tJ( 2, 1 ) ) / tDet ;
        aInvJ( 2, 2 ) = ( tJ( 0, 0 ) * tJ( 1, 1 ) - tJ( 0, 1 ) * tJ( 1, 0 ) ) / tDet ;
    }

    /**
     * circulation of every dof along every edge ( 2-pt Gauss, straight
     * edges ), expected s_k * delta_jk
     */
    void
    hex8_circulation_battery(
            const Matrix< real > & aX,
            const int              aFlipEdge = -1 )
    {
        fem::test::EF_TestVolume tFixture( ElementType::HEX8, aX, aFlipEdge );

        // 24-column pack: 2 Gauss points per edge, reference coordinates
        const real tG = 1.0 / std::sqrt( 3.0 );
        Matrix< real > tXi( 3, 24 );
        for ( uint e = 0; e < 12; ++e )
        {
            const uint a = gHexEdges[ e ][ 0 ];
            const uint b = gHexEdges[ e ][ 1 ];
            for ( uint g = 0; g < 2; ++g )
            {
                const real t = ( g == 0 ) ? -tG : tG ;
                for ( uint d = 0; d < 3; ++d )
                {
                    tXi( d, 2 * e + g ) =
                          0.5 * ( gHexXi[ a ][ d ] + gHexXi[ b ][ d ] )
                        + 0.5 * t * ( gHexXi[ b ][ d ] - gHexXi[ a ][ d ] );
                }
            }
        }

        fem::EdgeFunctionFactory tFactory ;
        fem::EdgeFunction * tEF =
                tFactory.create_edge_function( ElementType::HEX8 );
        tEF->precompute( tXi );
        tEF->link( tFixture.element() );

        real tS[ 12 ];
        tFixture.element()->edge_directions( tS );

        for ( uint e = 0; e < 12; ++e )
        {
            const uint a = gHexEdges[ e ][ 0 ];
            const uint b = gHexEdges[ e ][ 1 ];
            real tDx[ 3 ];
            for ( uint d = 0; d < 3; ++d )
            {
                tDx[ d ] = 0.5 * ( aX( d, b ) - aX( d, a ) );
            }
            for ( uint k = 0; k < 12; ++k )
            {
                real tCirc = 0.0 ;
                for ( uint g = 0; g < 2; ++g )
                {
                    const Matrix< real > & tE = tEF->E( 2 * e + g );
                    tCirc += tE( 0, k ) * tDx[ 0 ]
                           + tE( 1, k ) * tDx[ 1 ]
                           + tE( 2, k ) * tDx[ 2 ];
                }
                // dof circulation along the MESH edge direction: s_k on the
                // own edge, zero elsewhere
                const real tExpect = ( e == k ) ? tS[ k ] : 0.0 ;
                EXPECT_NEAR( tCirc, tExpect, tEpsilon )
                        << "edge " << e << " dof " << k ;
            }
        }
        delete tEF ;
    }

    /**
     * Stokes consistency on the box top face ( THE curl-sign gate ):
     * for every dof, the circulation around the top loop ( edges 4..7,
     * whose canonical directions already form the CCW-from-+z loop ) must
     * equal the 2x2-Gauss surface integral of C_z. BOTH sides integrated;
     * a curl tie alone passes when E and C are negated together.
     */
    void
    hex8_stokes_battery()
    {
        Matrix< real > tX = hex_corners( 0 );   // box, unflipped only
        fem::test::EF_TestVolume tFixture( ElementType::HEX8, tX );

        // pack: cols 0..7 = 2 Gauss points on each of edges 4..7,
        //       cols 8..11 = 2x2 face Gauss at zeta = +1
        const real tG = 1.0 / std::sqrt( 3.0 );
        Matrix< real > tXi( 3, 12 );
        for ( uint e = 0; e < 4; ++e )
        {
            const uint a = gHexEdges[ 4 + e ][ 0 ];
            const uint b = gHexEdges[ 4 + e ][ 1 ];
            for ( uint g = 0; g < 2; ++g )
            {
                const real t = ( g == 0 ) ? -tG : tG ;
                for ( uint d = 0; d < 3; ++d )
                {
                    tXi( d, 2 * e + g ) =
                          0.5 * ( gHexXi[ a ][ d ] + gHexXi[ b ][ d ] )
                        + 0.5 * t * ( gHexXi[ b ][ d ] - gHexXi[ a ][ d ] );
                }
            }
        }
        const real tFacePts[ 4 ][ 2 ] =
                { { -tG, -tG }, { tG, -tG }, { -tG, tG }, { tG, tG } };
        for ( uint g = 0; g < 4; ++g )
        {
            tXi( 0, 8 + g ) = tFacePts[ g ][ 0 ];
            tXi( 1, 8 + g ) = tFacePts[ g ][ 1 ];
            tXi( 2, 8 + g ) = 1.0 ;
        }

        fem::EdgeFunctionFactory tFactory ;
        fem::EdgeFunction * tEF =
                tFactory.create_edge_function( ElementType::HEX8 );
        tEF->precompute( tXi );
        tEF->link( tFixture.element() );

        // physical face area element of the box top face:
        // dS = ( Lx/2 )( Ly/2 ) dxi deta, from the same corner scales
        const real tJf = 1.2 * 0.9 ;

        for ( uint k = 0; k < 12; ++k )
        {
            // LHS: loop circulation, canonical directions = loop directions
            real tLoop = 0.0 ;
            for ( uint e = 0; e < 4; ++e )
            {
                const uint a = gHexEdges[ 4 + e ][ 0 ];
                const uint b = gHexEdges[ 4 + e ][ 1 ];
                real tDx[ 3 ];
                for ( uint d = 0; d < 3; ++d )
                {
                    tDx[ d ] = 0.5 * ( tX( d, b ) - tX( d, a ) );
                }
                for ( uint g = 0; g < 2; ++g )
                {
                    const Matrix< real > & tE = tEF->E( 2 * e + g );
                    tLoop += tE( 0, k ) * tDx[ 0 ]
                           + tE( 1, k ) * tDx[ 1 ]
                           + tE( 2, k ) * tDx[ 2 ];
                }
            }
            // RHS: surface integral of C_z over the top face
            real tFlux = 0.0 ;
            for ( uint g = 0; g < 4; ++g )
            {
                const Matrix< real > & tC = tEF->C( 8 + g );
                tFlux += tC( 2, k ) * tJf ;
            }
            EXPECT_NEAR( tLoop, tFlux, tEpsilon )
                    << "Stokes, dof " << k
                    << " ( loop " << tLoop << " flux " << tFlux << " )" ;
        }
        delete tEF ;
    }

    /**
     * gradient battery: shape, FD ( with the order window on the
     * distorted fixture ), curl tie at base AND non-base column,
     * G-must-move, box-only trace, and the three controls.
     * aFlavor: 0 box, 1 sheared, 2 distorted.
     */
    void
    hex8_gradient_battery(
            const uint aFlavor,
            const int  aFlipEdge = -1 )
    {
        Matrix< real > tX = hex_corners( aFlavor );
        fem::test::EF_TestVolume tFixture( ElementType::HEX8, tX, aFlipEdge );

        // 13-column pack: base + 6 at tHexDxi + 6 at tHexDxi/2.
        // The distorted fixture uses a COARSER step than the affine ones:
        // at tDxi = 1e-4 its truncation error ( ~5e-12, measured ) sits in
        // the roundoff regime and the order ratio is meaningless — the
        // audit's caveat, confirmed by the first run of this battery. At
        // 1e-2 truncation dominates and the O( h^2 ) window has teeth.
        // Affine fixtures are FD-exact at any step; they keep tDxi.
        const real tHexDxi = ( aFlavor == 2 ) ? 1.0e-2 : tDxi ;
        const real tBasePt[ 3 ] = { 0.3, -0.2, 0.45 };
        Matrix< real > tXi( 3, 13 );
        for ( uint d = 0; d < 3; ++d ) { tXi( d, 0 ) = tBasePt[ d ]; }
        uint tCount = 1 ;
        for ( uint h = 0; h < 2; ++h )
        {
            const real tH = ( h == 0 ) ? tHexDxi : 0.5 * tHexDxi ;
            for ( uint p = 0; p < 3; ++p )
            {
                for ( uint s = 0; s < 2; ++s )
                {
                    for ( uint d = 0; d < 3; ++d )
                    {
                        tXi( d, tCount ) = tBasePt[ d ];
                    }
                    tXi( p, tCount ) += ( s == 0 ? -tH : tH );
                    ++tCount ;
                }
            }
        }

        fem::EdgeFunctionFactory tFactory ;
        fem::EdgeFunction * tEF =
                tFactory.create_edge_function( ElementType::HEX8 );
        tEF->precompute( tXi );
        tEF->link( tFixture.element() );

        // shape first ( G() returns mGrad by reference — copy )
        Matrix< real > tGm( tEF->G( 0 ) );
        ASSERT_EQ( tGm.n_rows(), ( uint ) 9 ) << "G row count" ;
        ASSERT_EQ( tGm.n_cols(), ( uint ) 12 ) << "G column count" ;

        // curl tie at the base point ( fixed C, round-off )
        {
            Matrix< real > tC( tEF->C( 0 ) );
            for ( uint k = 0; k < 12; ++k )
            {
                EXPECT_NEAR( tC( 0, k ), tGm( 7, k ) - tGm( 5, k ), tEpsilon )
                        << "curl tie x, dof " << k ;
                EXPECT_NEAR( tC( 1, k ), tGm( 2, k ) - tGm( 6, k ), tEpsilon )
                        << "curl tie y, dof " << k ;
                EXPECT_NEAR( tC( 2, k ), tGm( 3, k ) - tGm( 1, k ), tEpsilon )
                        << "curl tie z, dof " << k ;
            }
        }

        // non-base column: curl tie there too, and G must move ( a G that
        // always reads column 0 passes every base-point check )
        {
            Matrix< real > tG1( tEF->G( 1 ) );
            const Matrix< real > & tC1 = tEF->C( 1 );
            real tMaxMove = 0.0 ;
            for ( uint k = 0; k < 12; ++k )
            {
                EXPECT_NEAR( tC1( 0, k ), tG1( 7, k ) - tG1( 5, k ), tEpsilon )
                        << "curl tie x at non-base, dof " << k ;
                EXPECT_NEAR( tC1( 1, k ), tG1( 2, k ) - tG1( 6, k ), tEpsilon )
                        << "curl tie y at non-base, dof " << k ;
                EXPECT_NEAR( tC1( 2, k ), tG1( 3, k ) - tG1( 1, k ), tEpsilon )
                        << "curl tie z at non-base, dof " << k ;
                for ( uint r = 0; r < 9; ++r )
                {
                    const real tV = std::abs( tG1( r, k ) - tGm( r, k ) );
                    if ( tV > tMaxMove ) { tMaxMove = tV ; }
                }
            }
            EXPECT_GT( tMaxMove, tEpsilon )
                    << "G must depend on the integration point" ;
        }

        // FD: test-side inverse Jacobian at the base point
        Matrix< real > tInvJ ;
        hex_jacobian( tX, tBasePt, tInvJ );

        real tMaxDevH  = 0.0 ;
        real tMaxDevH2 = 0.0 ;
        for ( uint h = 0; h < 2; ++h )
        {
            const real tH = ( h == 0 ) ? tHexDxi : 0.5 * tHexDxi ;
            const uint tOff = 1 + 6 * h ;
            Cell< Matrix< real > > tDE ;
            tDE.set_size( 3, {} );
            for ( uint p = 0; p < 3; ++p )
            {
                Matrix< real > tMinus( tEF->E( tOff + 2 * p ) );
                const Matrix< real > & tPlus = tEF->E( tOff + 2 * p + 1 );
                tDE( p ).set_size( 3, 12 );
                for ( uint c = 0; c < 3; ++c )
                {
                    for ( uint k = 0; k < 12; ++k )
                    {
                        tDE( p )( c, k ) =
                                ( tPlus( c, k ) - tMinus( c, k ) )
                                / ( 2.0 * tH );
                    }
                }
            }
            real tMaxDev = 0.0 ;
            for ( uint k = 0; k < 12; ++k )
            {
                for ( uint c = 0; c < 3; ++c )
                {
                    for ( uint m = 0; m < 3; ++m )
                    {
                        real tFD = 0.0 ;
                        for ( uint p = 0; p < 3; ++p )
                        {
                            tFD += tDE( p )( c, k ) * tInvJ( p, m );
                        }
                        const real tDev =
                                std::abs( tGm( m + 3 * c, k ) - tFD );
                        if ( tDev > tMaxDev ) { tMaxDev = tDev ; }
                        if ( aFlavor != 2 )
                        {
                            // affine fixtures: E is polynomial, FD exact
                            EXPECT_NEAR( tGm( m + 3 * c, k ), tFD,
                                    tEpsilonFD )
                                    << "G vs FD ( affine ), derivative "
                                    << m << " component " << c
                                    << " dof " << k ;
                        }
                    }
                }
            }
            if ( h == 0 ) { tMaxDevH = tMaxDev ; }
            else          { tMaxDevH2 = tMaxDev ; }
        }
        if ( aFlavor == 2 )
        {
            // distorted: no absolute tolerance — the expected-order window.
            // Central FD is O( h^2 ): halving h divides the deviation ~4
            EXPECT_GT( tMaxDevH, 1.0e-10 )
                    << "distorted FD must be truncation-dominated" ;
            const real tRatio = tMaxDevH / tMaxDevH2 ;
            EXPECT_GT( tRatio, 3.0 ) << "FD order window ( low )" ;
            EXPECT_LT( tRatio, 5.0 ) << "FD order window ( high )" ;
            std::cout << "    [ HEX8 G ] distorted FD dev( h ) = "
                      << tMaxDevH << ", dev( h/2 ) = " << tMaxDevH2
                      << ", ratio = " << tRatio << std::endl ;
        }

        // trace: asserted on the box ONLY ( sheared and distorted hexes
        // have genuinely nonzero basis divergence )
        if ( aFlavor == 0 )
        {
            for ( uint k = 0; k < 12; ++k )
            {
                EXPECT_NEAR(
                        tGm( 0, k ) + tGm( 4, k ) + tGm( 8, k ), 0.0,
                        tEpsilon ) << "box trace of dof " << k ;
            }
        }

        // controls. q_e = s_e * ( phi_b - phi_a ) on the local edge pairs
        real tS[ 12 ];
        tFixture.element()->edge_directions( tS );

        // negative: linear phi is interpolated exactly on EVERY trilinear
        // hex ( isoparametric identity ) -> G*q = 0
        {
            real tQ[ 12 ];
            for ( uint e = 0; e < 12; ++e )
            {
                const uint a = gHexEdges[ e ][ 0 ];
                const uint b = gHexEdges[ e ][ 1 ];
                real tPhiA = 0.0, tPhiB = 0.0 ;
                for ( uint d = 0; d < 3; ++d )
                {
                    tPhiA += ( d + 1.0 ) * tX( d, a );
                    tPhiB += ( d + 1.0 ) * tX( d, b );
                }
                tQ[ e ] = tS[ e ] * ( tPhiB - tPhiA );
            }
            for ( uint r = 0; r < 9; ++r )
            {
                real tGq = 0.0 ;
                for ( uint e = 0; e < 12; ++e )
                {
                    tGq += tGm( r, e ) * tQ[ e ];
                }
                EXPECT_NEAR( tGq, 0.0, tEpsilon )
                        << "linear-phi negative control, row " << r ;
            }
        }

        if ( aFlavor == 0 )
        {
            // positive-sym on the box: phi = x*y is in the interpolation
            // space exactly; G*q = vec( Hess ) with rows 1 and 3 equal 1
            real tQ[ 12 ];
            for ( uint e = 0; e < 12; ++e )
            {
                const uint a = gHexEdges[ e ][ 0 ];
                const uint b = gHexEdges[ e ][ 1 ];
                tQ[ e ] = tS[ e ] * ( tX( 0, b ) * tX( 1, b )
                                    - tX( 0, a ) * tX( 1, a ) );
            }
            for ( uint r = 0; r < 9; ++r )
            {
                real tGq = 0.0 ;
                for ( uint e = 0; e < 12; ++e )
                {
                    tGq += tGm( r, e ) * tQ[ e ];
                }
                const real tExpect = ( r == 1 || r == 3 ) ? 1.0 : 0.0 ;
                EXPECT_NEAR( tGq, tExpect, tEpsilon )
                        << "xy positive control, row " << r ;
            }
        }

        if ( aFlavor == 2 )
        {
            // positive-H on the distorted hex: nodal phi = xi ( parameter
            // values +-1 ) gives grad( phi_h ) = grad( xi ), so G*q is the
            // inverse-map Hessian H^(0): curl-free, symmetric, NONZERO —
            // the only control that sees the Hessian term without FD
            real tQ[ 12 ];
            for ( uint e = 0; e < 12; ++e )
            {
                const uint a = gHexEdges[ e ][ 0 ];
                const uint b = gHexEdges[ e ][ 1 ];
                tQ[ e ] = tS[ e ] *
                        ( gHexXi[ b ][ 0 ] - gHexXi[ a ][ 0 ] );
            }
            Matrix< real > tC( tEF->C( 0 ) );
            for ( uint r = 0; r < 3; ++r )
            {
                real tCq = 0.0 ;
                for ( uint e = 0; e < 12; ++e )
                {
                    tCq += tC( r, e ) * tQ[ e ];
                }
                EXPECT_NEAR( tCq, 0.0, tEpsilon )
                        << "xi control is curl-free, row " << r ;
            }
            real tMaxGq = 0.0 ;
            real tGq[ 9 ];
            for ( uint r = 0; r < 9; ++r )
            {
                tGq[ r ] = 0.0 ;
                for ( uint e = 0; e < 12; ++e )
                {
                    tGq[ r ] += tGm( r, e ) * tQ[ e ];
                }
                const real tV = std::abs( tGq[ r ] );
                if ( tV > tMaxGq ) { tMaxGq = tV ; }
            }
            EXPECT_GT( tMaxGq, tEpsilon )
                    << "xi control must see the Hessian term" ;
            // G*q is a Hessian: explicitly symmetric, so the pin does not
            // depend on the curl tie staying in this function ( audit
            // hardening, HEX8 C.3 round )
            EXPECT_NEAR( tGq[ 1 ], tGq[ 3 ], tEpsilon ) << "Hess sym xy" ;
            EXPECT_NEAR( tGq[ 2 ], tGq[ 6 ], tEpsilon ) << "Hess sym xz" ;
            EXPECT_NEAR( tGq[ 5 ], tGq[ 7 ], tEpsilon ) << "Hess sym yz" ;
        }

        delete tEF ;
    }
}

TEST( EdgeFunctions, Hex8CirculationBox )
{
    hex8_circulation_battery( hex_corners( 0 ) );
}

TEST( EdgeFunctions, Hex8CirculationDistorted )
{
    hex8_circulation_battery( hex_corners( 2 ) );
}

TEST( EdgeFunctions, Hex8StokesConsistency )
{
    hex8_stokes_battery();
}

TEST( EdgeFunctions, Hex8GradientBox )
{
    hex8_gradient_battery( 0 );
}

TEST( EdgeFunctions, Hex8GradientSheared )
{
    hex8_gradient_battery( 1 );
}

TEST( EdgeFunctions, Hex8GradientDistorted )
{
    hex8_gradient_battery( 2 );
}

TEST( EdgeFunctions, Hex8GradientDistortedFlipped )
{
    hex8_gradient_battery( 2, 1 );
}

//------------------------------------------------------------------------------
//  Higher-order gradient battery ( TRI6, TET10 — straight AND curved paths )
//------------------------------------------------------------------------------

#include "cl_IF_InterpolationFunctionFactory.hpp"

namespace
{
    // coarse FD step for genuinely curved quadratic geometry: at the affine
    // battery's tDxi = 1e-4 the truncation error is roundoff-dominated and
    // the order ratio is meaningless ( the HEX8 lesson, measured there )
    const real tDxiHo = 1.0e-2 ;

    /**
     * pinned midnode offsets ( dim x numEdges ), every axis nonzero
     * somewhere — a single-axis offset would leave Voigt rows untested
     */
    Matrix< real >
    ho_mid_offsets( const ElementSpec & aSpec )
    {
        Matrix< real > aOff( aSpec.mDim, aSpec.mNumEdges );
        if ( aSpec.mDim == 2 )
        {
            const real tOff[ 2 ][ 3 ] = {
                    {  0.09, -0.06,  0.07 },
                    { -0.07,  0.10,  0.05 } };
            for ( uint e = 0; e < 3; ++e )
            {
                aOff( 0, e ) = tOff[ 0 ][ e ];
                aOff( 1, e ) = tOff[ 1 ][ e ];
            }
        }
        else
        {
            const real tOff[ 3 ][ 6 ] = {
                    {  0.08, -0.06,  0.05,  0.07, -0.05,  0.04 },
                    { -0.05,  0.09, -0.07,  0.04,  0.06, -0.08 },
                    {  0.06,  0.04, -0.08, -0.06,  0.09,  0.05 } };
            for ( uint e = 0; e < 6; ++e )
            {
                aOff( 0, e ) = tOff[ 0 ][ e ];
                aOff( 1, e ) = tOff[ 1 ][ e ];
                aOff( 2, e ) = tOff[ 2 ][ e ];
            }
        }
        return aOff ;
    }

    /**
     * test-side Jacobian of the ( possibly curved ) quadratic map at a
     * reference point, from ALL nodes of the fixture element via the
     * LAGRANGE interpolation class ( independent of the Nedelec code under
     * test ): J( m, d ) = dx_m / dxi_d, inverted by cofactors.
     * Returns detJ so the caller can gate positivity.
     */
    real
    ho_jacobian(
            mesh::Element                * aElement,
            fem::InterpolationFunction   * aLagrange,
            const real                   * aXi,
            const uint                     aDim,
            Matrix< real >               & aInvJ )
    {
        const uint tNumNodes = aElement->number_of_nodes();

        Vector< real > tXi( aDim );
        for ( uint p = 0; p < aDim; ++p ) { tXi( p ) = aXi[ p ]; }

        Matrix< real > tdN ;
        aLagrange->dNdXi( tXi, tdN );          // dim x nodes

        Matrix< real > tJ( aDim, aDim, 0.0 );  // J( m, d ) = dx_m / dxi_d
        for ( uint n = 0; n < tNumNodes; ++n )
        {
            const mesh::Node * tNode = aElement->node( n );
            const real tXn[ 3 ] = { tNode->x(), tNode->y(), tNode->z() };
            for ( uint p = 0; p < aDim; ++p )
            {
                for ( uint m = 0; m < aDim; ++m )
                {
                    tJ( m, p ) += tdN( p, n ) * tXn[ m ];
                }
            }
        }

        aInvJ.set_size( aDim, aDim );
        real tDet ;
        if ( aDim == 2 )
        {
            tDet = tJ( 0, 0 ) * tJ( 1, 1 ) - tJ( 0, 1 ) * tJ( 1, 0 );
            aInvJ( 0, 0 ) =  tJ( 1, 1 ) / tDet ;
            aInvJ( 0, 1 ) = -tJ( 0, 1 ) / tDet ;
            aInvJ( 1, 0 ) = -tJ( 1, 0 ) / tDet ;
            aInvJ( 1, 1 ) =  tJ( 0, 0 ) / tDet ;
        }
        else
        {
            tDet = tJ( 0, 0 ) * ( tJ( 1, 1 ) * tJ( 2, 2 ) - tJ( 1, 2 ) * tJ( 2, 1 ) )
                 - tJ( 0, 1 ) * ( tJ( 1, 0 ) * tJ( 2, 2 ) - tJ( 1, 2 ) * tJ( 2, 0 ) )
                 + tJ( 0, 2 ) * ( tJ( 1, 0 ) * tJ( 2, 1 ) - tJ( 1, 1 ) * tJ( 2, 0 ) );
            aInvJ( 0, 0 ) = ( tJ( 1, 1 ) * tJ( 2, 2 ) - tJ( 1, 2 ) * tJ( 2, 1 ) ) / tDet ;
            aInvJ( 0, 1 ) = ( tJ( 0, 2 ) * tJ( 2, 1 ) - tJ( 0, 1 ) * tJ( 2, 2 ) ) / tDet ;
            aInvJ( 0, 2 ) = ( tJ( 0, 1 ) * tJ( 1, 2 ) - tJ( 0, 2 ) * tJ( 1, 1 ) ) / tDet ;
            aInvJ( 1, 0 ) = ( tJ( 1, 2 ) * tJ( 2, 0 ) - tJ( 1, 0 ) * tJ( 2, 2 ) ) / tDet ;
            aInvJ( 1, 1 ) = ( tJ( 0, 0 ) * tJ( 2, 2 ) - tJ( 0, 2 ) * tJ( 2, 0 ) ) / tDet ;
            aInvJ( 1, 2 ) = ( tJ( 0, 2 ) * tJ( 1, 0 ) - tJ( 0, 0 ) * tJ( 1, 2 ) ) / tDet ;
            aInvJ( 2, 0 ) = ( tJ( 1, 0 ) * tJ( 2, 1 ) - tJ( 1, 1 ) * tJ( 2, 0 ) ) / tDet ;
            aInvJ( 2, 1 ) = ( tJ( 0, 1 ) * tJ( 2, 0 ) - tJ( 0, 0 ) * tJ( 2, 1 ) ) / tDet ;
            aInvJ( 2, 2 ) = ( tJ( 0, 0 ) * tJ( 1, 1 ) - tJ( 0, 1 ) * tJ( 1, 0 ) ) / tDet ;
        }
        return tDet ;
    }

    /**
     * gradient battery for the higher-order pair. NOT test_gradient():
     * that helper is first-order only ( Frobenius, trace == 0 and the
     * gradient-mode control are all FALSE here ). Checks:
     *   shape ASSERT; curl tie ( base + non-base + G-must-move );
     *   FD vs E through the quadratic test-side Jacobian
     *     ( straight: exact at tEpsilonFD; curved: coarse step pair,
     *       truncation floor + O( h^2 ) ratio window, NO absolute tol );
     *   weak nonzero trace ( straight distorted only — the Q5 expectation );
     *   detJ > 0 at every cluster point.
     */
    void
    test_gradient_ho(
            const ElementSpec    & aSpec,
            const Matrix< real > & aCorners,
            const int              aFlipEdge,
            const bool             aCurvedGeometry,
            const bool             aCheckTrace = false )
    {
        const uint d = aSpec.mDim ;
        const uint tNumDofs = aSpec.mNumEdges * aSpec.mDofsPerEdge
                            + aSpec.mNumFaceDofs ;

        Matrix< real > tOffsets ;
        const Matrix< real > * tOffPtr = nullptr ;
        if ( aCurvedGeometry )
        {
            tOffsets = ho_mid_offsets( aSpec );
            tOffPtr  = & tOffsets ;
        }
        fem::test::EF_TestVolume tFixture( aSpec.mType, aCorners,
                aFlipEdge, aCurvedGeometry, tOffPtr );

        // pack: base + 2d at tH + 2d at tH/2
        const real tH0 = aCurvedGeometry ? tDxiHo : tDxi ;
        const real tBasePt[ 3 ] = { 0.3, 0.25, 0.2 };
        Matrix< real > tXi( d, 1 + 4 * d );
        for ( uint p = 0; p < d; ++p ) { tXi( p, 0 ) = tBasePt[ p ]; }
        uint tCount = 1 ;
        for ( uint h = 0; h < 2; ++h )
        {
            const real tH = ( h == 0 ) ? tH0 : 0.5 * tH0 ;
            for ( uint p = 0; p < d; ++p )
            {
                for ( uint s = 0; s < 2; ++s )
                {
                    for ( uint q = 0; q < d; ++q )
                    {
                        tXi( q, tCount ) = tBasePt[ q ];
                    }
                    tXi( p, tCount ) += ( s == 0 ? -tH : tH );
                    ++tCount ;
                }
            }
        }

        fem::EdgeFunctionFactory tFactory ;
        fem::EdgeFunction * tEF =
                tFactory.create_edge_function( aSpec.mType );
        tEF->precompute( tXi );
        tEF->link( tFixture.element() );

        // test-side Lagrange for the quadratic map
        fem::InterpolationFunctionFactory tIfFactory ;
        fem::InterpolationFunction * tLagrange =
                tIfFactory.create_lagrange_function( aSpec.mType );
        mesh::Element * tMeshElement = tFixture.element()->element();

        // detJ > 0 at every cluster point ( curved offsets must not fold )
        Matrix< real > tInvJ ;
        for ( uint c = 0; c < 1 + 4 * d; ++c )
        {
            real tPt[ 3 ] = { 0.0, 0.0, 0.0 };
            for ( uint p = 0; p < d; ++p ) { tPt[ p ] = tXi( p, c ); }
            const real tDet = ho_jacobian(
                    tMeshElement, tLagrange, tPt, d, tInvJ );
            ASSERT_GT( tDet, 0.0 ) << "detJ at cluster point " << c ;
        }

        // shape ( G returns mGrad by reference — copy )
        Matrix< real > tGm( tEF->G( 0 ) );
        ASSERT_EQ( tGm.n_rows(), d * d ) << "G row count" ;
        ASSERT_EQ( tGm.n_cols(), tNumDofs ) << "G column count" ;

        // curl tie at base + non-base + G-must-move
        {
            Matrix< real > tC( tEF->C( 0 ) );
            for ( uint k = 0; k < tNumDofs; ++k )
            {
                if ( d == 2 )
                {
                    EXPECT_NEAR( tC( 0, k ), tGm( 2, k ) - tGm( 1, k ),
                            tEpsilon ) << "curl tie, dof " << k ;
                }
                else
                {
                    EXPECT_NEAR( tC( 0, k ), tGm( 7, k ) - tGm( 5, k ),
                            tEpsilon ) << "curl tie x, dof " << k ;
                    EXPECT_NEAR( tC( 1, k ), tGm( 2, k ) - tGm( 6, k ),
                            tEpsilon ) << "curl tie y, dof " << k ;
                    EXPECT_NEAR( tC( 2, k ), tGm( 3, k ) - tGm( 1, k ),
                            tEpsilon ) << "curl tie z, dof " << k ;
                }
            }
            Matrix< real > tG1( tEF->G( 1 ) );
            const Matrix< real > & tC1 = tEF->C( 1 );
            real tMaxMove = 0.0 ;
            for ( uint k = 0; k < tNumDofs; ++k )
            {
                if ( d == 2 )
                {
                    EXPECT_NEAR( tC1( 0, k ), tG1( 2, k ) - tG1( 1, k ),
                            tEpsilon ) << "curl tie at non-base, dof " << k ;
                }
                else
                {
                    EXPECT_NEAR( tC1( 0, k ), tG1( 7, k ) - tG1( 5, k ),
                            tEpsilon ) << "curl tie x at non-base, dof " << k ;
                    EXPECT_NEAR( tC1( 1, k ), tG1( 2, k ) - tG1( 6, k ),
                            tEpsilon ) << "curl tie y at non-base, dof " << k ;
                    EXPECT_NEAR( tC1( 2, k ), tG1( 3, k ) - tG1( 1, k ),
                            tEpsilon ) << "curl tie z at non-base, dof " << k ;
                }
                for ( uint r = 0; r < d * d; ++r )
                {
                    const real tV = std::abs( tG1( r, k ) - tGm( r, k ) );
                    if ( tV > tMaxMove ) { tMaxMove = tV ; }
                }
            }
            EXPECT_GT( tMaxMove, tEpsilon )
                    << "G must depend on the integration point" ;
        }

        // FD vs E through the quadratic test-side Jacobian at the base
        real tBase3[ 3 ] = { 0.0, 0.0, 0.0 };
        for ( uint p = 0; p < d; ++p ) { tBase3[ p ] = tBasePt[ p ]; }
        ho_jacobian( tMeshElement, tLagrange, tBase3, d, tInvJ );

        real tMaxDevH  = 0.0 ;
        real tMaxDevH2 = 0.0 ;
        for ( uint h = 0; h < 2; ++h )
        {
            const real tH = ( h == 0 ) ? tH0 : 0.5 * tH0 ;
            const uint tOff = 1 + 2 * d * h ;
            Cell< Matrix< real > > tDE ;
            tDE.set_size( d, {} );
            for ( uint p = 0; p < d; ++p )
            {
                Matrix< real > tMinus( tEF->E( tOff + 2 * p ) );
                const Matrix< real > & tPlus = tEF->E( tOff + 2 * p + 1 );
                tDE( p ).set_size( d, tNumDofs );
                for ( uint c = 0; c < d; ++c )
                {
                    for ( uint k = 0; k < tNumDofs; ++k )
                    {
                        tDE( p )( c, k ) =
                                ( tPlus( c, k ) - tMinus( c, k ) )
                                / ( 2.0 * tH );
                    }
                }
            }
            real tMaxDev = 0.0 ;
            for ( uint k = 0; k < tNumDofs; ++k )
            {
                for ( uint c = 0; c < d; ++c )
                {
                    for ( uint m = 0; m < d; ++m )
                    {
                        real tFD = 0.0 ;
                        for ( uint p = 0; p < d; ++p )
                        {
                            tFD += tDE( p )( c, k ) * tInvJ( p, m );
                        }
                        const real tDev =
                                std::abs( tGm( m + d * c, k ) - tFD );
                        if ( tDev > tMaxDev ) { tMaxDev = tDev ; }
                        if ( ! aCurvedGeometry )
                        {
                            // straight: E is a quadratic polynomial in x,
                            // central FD is exact
                            EXPECT_NEAR( tGm( m + d * c, k ), tFD,
                                    tEpsilonFD )
                                    << "G vs FD ( straight ), derivative "
                                    << m << " component " << c
                                    << " dof " << k ;
                        }
                    }
                }
            }
            if ( h == 0 ) { tMaxDevH = tMaxDev ; }
            else          { tMaxDevH2 = tMaxDev ; }
        }
        if ( aCurvedGeometry )
        {
            EXPECT_GT( tMaxDevH, 1.0e-10 )
                    << "curved FD must be truncation-dominated" ;
            const real tRatio = tMaxDevH / tMaxDevH2 ;
            EXPECT_GT( tRatio, 3.0 ) << "FD order window ( low )" ;
            EXPECT_LT( tRatio, 5.0 ) << "FD order window ( high )" ;
            std::cout << "    [ HO G ] curved FD dev( h ) = " << tMaxDevH
                      << ", dev( h/2 ) = " << tMaxDevH2
                      << ", ratio = " << tRatio << std::endl ;
        }

        // weak nonzero trace ( the Q5 expectation ) — the CALLER opts in
        // on the straight distorted unflipped fixture only ( reference
        // symmetries can shrink traces; binding v2 delta 10 )
        if ( aCheckTrace )
        {
            real tMaxTr = 0.0 ;
            for ( uint k = 0; k < tNumDofs; ++k )
            {
                real tTr = 0.0 ;
                for ( uint i = 0; i < d; ++i )
                {
                    tTr += tGm( i + d * i, k );
                }
                const real tV = std::abs( tTr );
                if ( tV > tMaxTr ) { tMaxTr = tV ; }
            }
            EXPECT_GT( tMaxTr, tEpsilon )
                    << "higher-order basis divergence must be nonzero" ;
        }

        delete tLagrange ;
        delete tEF ;
    }

    /**
     * path equivalence on STRAIGHT geometry: the curved flag selects the
     * other code path but the map is affine, so G must agree at every
     * cluster point. Pins the fork and — with the NaN-filled mCurv of the
     * straight branch — a missed hess skip ( which would produce NaN ).
     * On straight geometry S = 0, so this test has NO power over the hess
     * algebra itself; the curved battery carries that.
     */
    void
    test_gradient_ho_path_equivalence(
            const ElementSpec    & aSpec,
            const Matrix< real > & aCorners )
    {
        const uint d = aSpec.mDim ;
        const uint tNumDofs = aSpec.mNumEdges * aSpec.mDofsPerEdge
                            + aSpec.mNumFaceDofs ;

        fem::test::EF_TestVolume tStraight( aSpec.mType, aCorners, -1, false );
        fem::test::EF_TestVolume tCurved(   aSpec.mType, aCorners, -1, true );

        Matrix< real > tXi = evaluation_points( aSpec );

        fem::EdgeFunctionFactory tFactory ;
        fem::EdgeFunction * tEFS =
                tFactory.create_edge_function( aSpec.mType );
        fem::EdgeFunction * tEFC =
                tFactory.create_edge_function( aSpec.mType );
        tEFS->precompute( tXi );
        tEFS->link( tStraight.element() );
        tEFC->precompute( tXi );
        tEFC->link( tCurved.element() );

        const uint tNumCols = tXi.n_cols();
        for ( uint c = 0; c < tNumCols; ++c )
        {
            Matrix< real > tGS( tEFS->G( c ) );
            const Matrix< real > & tGC = tEFC->G( c );
            for ( uint r = 0; r < d * d; ++r )
            {
                for ( uint k = 0; k < tNumDofs; ++k )
                {
                    EXPECT_NEAR( tGS( r, k ), tGC( r, k ), tEpsilon )
                            << "G straight vs curved path, point " << c
                            << " entry ( " << r << ", " << k << " )" ;
                }
            }
        }
        delete tEFS ;
        delete tEFC ;
    }
}

TEST( EdgeFunctions, Tri6GradientReference )
{
    ElementSpec s = spec( ElementType::TRI6 );
    test_gradient_ho( s, reference_corners( s ), -1, false );
}

TEST( EdgeFunctions, Tri6GradientDistorted )
{
    ElementSpec s = spec( ElementType::TRI6 );
    test_gradient_ho( s, distorted_corners( s ), -1, false, true );
}

TEST( EdgeFunctions, Tri6GradientDistortedFlipped )
{
    ElementSpec s = spec( ElementType::TRI6 );
    test_gradient_ho( s, distorted_corners( s ), 1, false );
}

TEST( EdgeFunctions, Tri6GradientPathEquivalence )
{
    ElementSpec s = spec( ElementType::TRI6 );
    test_gradient_ho_path_equivalence( s, distorted_corners( s ) );
}

TEST( EdgeFunctions, Tri6GradientCurved )
{
    ElementSpec s = spec( ElementType::TRI6 );
    test_gradient_ho( s, distorted_corners( s ), -1, true );
}

TEST( EdgeFunctions, Tet10GradientReference )
{
    ElementSpec s = spec( ElementType::TET10 );
    test_gradient_ho( s, reference_corners( s ), -1, false );
}

TEST( EdgeFunctions, Tet10GradientDistorted )
{
    ElementSpec s = spec( ElementType::TET10 );
    test_gradient_ho( s, distorted_corners( s ), -1, false, true );
}

TEST( EdgeFunctions, Tet10GradientDistortedFlipped )
{
    ElementSpec s = spec( ElementType::TET10 );
    test_gradient_ho( s, distorted_corners( s ), 1, false );
}

TEST( EdgeFunctions, Tet10GradientPathEquivalence )
{
    ElementSpec s = spec( ElementType::TET10 );
    test_gradient_ho_path_equivalence( s, distorted_corners( s ) );
}

TEST( EdgeFunctions, Tet10GradientCurved )
{
    ElementSpec s = spec( ElementType::TET10 );
    test_gradient_ho( s, distorted_corners( s ), -1, true );
}
