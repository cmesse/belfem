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
 * The inverse element map. The oracle of every round trip is the mapper's
 * own forward map: weights( xi ) gives the shape values, x = sum N_i X_i,
 * and evaluate( x ) must return xi. The elements are standalone mesh
 * elements over the reference element's nodes, moved by an affine or a
 * curved transform; no kernel is built.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Node.hpp"
#include "cl_Element.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_ReferenceElement.hpp"
#include "meshtools.hpp"
#include "cl_IF_ElementMapper.hpp"

using namespace belfem;

namespace
{
    using Transform = std::function< void( const real aRef[ 3 ], real aX[ 3 ] ) > ;

    /**
     * a mesh element whose nodes are the reference element's nodes moved by
     * aTransform. Owns the element and the nodes.
     */
    class MappedElement
    {
        mesh::Element        * mElement = nullptr ;
        Cell< mesh::Node * >   mNodes ;

    public:

        MappedElement( const ElementType aType, const Transform & aTransform )
        {
            mesh::ElementFactory tFactory ;
            std::unique_ptr< mesh::ReferenceElement > tRef(
                    tFactory.create_reference_element( aType ) );

            const uint n = tRef->element()->number_of_nodes() ;
            mElement = tFactory.create_element( aType, 1 );
            mNodes.set_size( n, nullptr );

            for ( uint k = 0; k < n; ++k )
            {
                const real tRefX[ 3 ] = {
                        tRef->element()->node( k )->x(),
                        tRef->element()->node( k )->y(),
                        tRef->element()->node( k )->z() };
                real tX[ 3 ] = { 0.0, 0.0, 0.0 };
                aTransform( tRefX, tX );

                mNodes( k ) = new mesh::Node( k + 1, tX[ 0 ], tX[ 1 ], tX[ 2 ] );
                mElement->insert_node( mNodes( k ), k );
            }
        }

        MappedElement( const MappedElement & ) = delete ;
        MappedElement & operator=( const MappedElement & ) = delete ;

        ~MappedElement()
        {
            delete mElement ;
            for ( mesh::Node * tNode : mNodes )
            {
                delete tNode ;
            }
        }

        mesh::Element *
        element() { return mElement ; }

        //! the forward map through the mapper's own shape functions
        void
        forward( fem::ElementMapper & aMapper, const Vector< real > & aXi, Vector< real > & aX )
        {
            const uint n = mElement->number_of_nodes() ;
            const uint d = aX.length() ;
            Vector< real > tN( n, 0.0 );
            aMapper.weights( aXi, tN );

            for ( uint i = 0; i < d; ++i )
            {
                real tSum = 0.0 ;
                for ( uint k = 0; k < n; ++k )
                {
                    tSum += tN( k ) * mNodes( k )->x( i ) ;
                }
                aX( i ) = tSum ;
            }
        }
    };

    void
    identity( const real aRef[ 3 ], real aX[ 3 ] )
    {
        aX[ 0 ] = aRef[ 0 ] ; aX[ 1 ] = aRef[ 1 ] ; aX[ 2 ] = aRef[ 2 ] ;
    }

    // rotation, anisotropic stretch and shift; no axis stays aligned
    void
    affine_2d( const real aRef[ 3 ], real aX[ 3 ] )
    {
        aX[ 0 ] =  1.7 * aRef[ 0 ] + 0.4 * aRef[ 1 ] + 3.0 ;
        aX[ 1 ] = -0.3 * aRef[ 0 ] + 1.2 * aRef[ 1 ] - 1.0 ;
        aX[ 2 ] = 0.0 ;
    }

    void
    affine_3d( const real aRef[ 3 ], real aX[ 3 ] )
    {
        aX[ 0 ] =  1.7 * aRef[ 0 ] + 0.4 * aRef[ 1 ] - 0.2 * aRef[ 2 ] + 3.0 ;
        aX[ 1 ] = -0.3 * aRef[ 0 ] + 1.2 * aRef[ 1 ] + 0.1 * aRef[ 2 ] - 1.0 ;
        aX[ 2 ] =  0.2 * aRef[ 0 ] - 0.1 * aRef[ 1 ] + 0.9 * aRef[ 2 ] + 0.5 ;
    }

    // a bilinear, non-affine quadrilateral: the right edge is tilted, so
    // the xi*eta coefficients do not vanish
    void
    trapezoid( const real aRef[ 3 ], real aX[ 3 ] )
    {
        aX[ 0 ] = aRef[ 0 ] + 0.3 * aRef[ 0 ] * aRef[ 1 ] + 0.1 * aRef[ 1 ] ;
        aX[ 1 ] = aRef[ 1 ] + 0.2 * aRef[ 0 ] * aRef[ 1 ] ;
        aX[ 2 ] = 0.0 ;
    }

    // a smooth curved distortion; every tested element stays valid
    void
    curved_2d( const real aRef[ 3 ], real aX[ 3 ] )
    {
        aX[ 0 ] = aRef[ 0 ] + 0.08 * aRef[ 1 ] * aRef[ 1 ] ;
        aX[ 1 ] = aRef[ 1 ] + 0.06 * aRef[ 0 ] * aRef[ 0 ] ;
        aX[ 2 ] = 0.0 ;
    }

    void
    curved_3d( const real aRef[ 3 ], real aX[ 3 ] )
    {
        aX[ 0 ] = aRef[ 0 ] + 0.05 * aRef[ 1 ] * aRef[ 2 ] ;
        aX[ 1 ] = aRef[ 1 ] + 0.04 * aRef[ 0 ] * aRef[ 0 ] ;
        aX[ 2 ] = aRef[ 2 ] + 0.03 * aRef[ 0 ] * aRef[ 1 ] ;
    }

    /**
     * forward map, inverse map, compare. The tolerance is absolute in the
     * natural coordinates, which are O(1)
     */
    void
    expect_round_trip(
            const ElementType aType,
            const Transform & aTransform,
            const Cell< Vector< real > > & aPoints,
            const real aTolerance )
    {
        MappedElement tElement( aType, aTransform );
        fem::ElementMapper tMapper ;
        tMapper.link( tElement.element() );

        const uint tNatDim = mesh::dimension( aType );
        Vector< real > tX( tNatDim, 0.0 );
        Vector< real > tXi( tNatDim, 0.0 );

        for ( const Vector< real > & tXiRef : aPoints )
        {
            tElement.forward( tMapper, tXiRef, tX );

            ASSERT_TRUE( tMapper.evaluate( tX, tXi ) )
                << "type " << ( int ) aType << ": interior point reported outside" ;

            for ( uint i = 0; i < tNatDim; ++i )
            {
                EXPECT_NEAR( tXi( i ), tXiRef( i ), aTolerance )
                    << "type " << ( int ) aType << " coordinate " << i ;
            }
        }
    }

    Cell< Vector< real > >
    tri_points()
    {
        return { { 0.2, 0.3 }, { 0.6, 0.1 }, { 0.05, 0.85 }, { 1.0 / 3.0, 1.0 / 3.0 } };
    }

    Cell< Vector< real > >
    quad_points()
    {
        return { { 0.0, 0.0 }, { 0.5, -0.4 }, { -0.9, 0.7 }, { 0.3, 0.95 } };
    }

    Cell< Vector< real > >
    tet_points()
    {
        return { { 0.1, 0.2, 0.3 }, { 0.25, 0.25, 0.25 }, { 0.6, 0.1, 0.05 } };
    }

    Cell< Vector< real > >
    hex_points()
    {
        return { { 0.0, 0.0, 0.0 }, { 0.5, -0.4, 0.2 }, { -0.9, 0.7, -0.6 } };
    }

    Cell< Vector< real > >
    pyra_points()
    {
        return { { 0.0, 0.0, 0.25 }, { 0.3, -0.2, 0.4 }, { -0.5, 0.4, 0.1 } };
    }
}

//------------------------------------------------------------------------------

TEST( ElementMapper, Tri3RoundTrip )
{
    expect_round_trip( ElementType::TRI3, identity,  tri_points(), 1e-12 );
    expect_round_trip( ElementType::TRI3, affine_2d, tri_points(), 1e-12 );
}

TEST( ElementMapper, Quad4RoundTripAffineAndTrapezoid )
{
    // the xi*eta == 0 branch, then the quadratic-in-eta branch
    expect_round_trip( ElementType::QUAD4, affine_2d, quad_points(), 1e-12 );
    expect_round_trip( ElementType::QUAD4, trapezoid, quad_points(), 1e-12 );
}

TEST( ElementMapper, Tet4RoundTrip )
{
    expect_round_trip( ElementType::TET4, identity,  tet_points(), 1e-12 );
    expect_round_trip( ElementType::TET4, affine_3d, tet_points(), 1e-12 );
}

TEST( ElementMapper, NewtonRoundTripOnCurvedElements )
{
    // the general path: seeded by the linear inverse for TRI and QUAD,
    // by the centroid for the solids
    expect_round_trip( ElementType::TRI6,  curved_2d, tri_points(),  1e-12 );
    expect_round_trip( ElementType::QUAD9, curved_2d, quad_points(), 1e-12 );
    expect_round_trip( ElementType::HEX8,  curved_3d, hex_points(),  1e-12 );
    expect_round_trip( ElementType::HEX27, curved_3d, hex_points(),  1e-12 );
    expect_round_trip( ElementType::PYRA5, affine_3d, pyra_points(), 1e-12 );
}

//------------------------------------------------------------------------------

TEST( ElementMapper, OutsidePointsAreRejected )
{
    Vector< real > tXi ;

    {
        // TRI3: below each leg, and beyond the hypotenuse xi + eta > 1
        MappedElement tElement( ElementType::TRI3, affine_2d );
        fem::ElementMapper tMapper ;
        tMapper.link( tElement.element() );
        tXi.set_size( 2 );

        Vector< real > tX( 2 );
        for ( const Vector< real > & tOut : Cell< Vector< real > >(
                { { 0.6, 0.6 }, { -0.1, 0.5 }, { 0.5, -0.1 } } ) )
        {
            tElement.forward( tMapper, tOut, tX );
            EXPECT_FALSE( tMapper.evaluate( tX, tXi ) )
                << "xi = " << tOut( 0 ) << ", eta = " << tOut( 1 ) ;
        }

        tElement.forward( tMapper, { 0.4, 0.4 }, tX );
        EXPECT_TRUE( tMapper.evaluate( tX, tXi ) );
    }
    {
        // QUAD4 on both branches: |xi| > 1 and |eta| > 1
        for ( const Transform & tTransform : { Transform( identity ), Transform( trapezoid ) } )
        {
            MappedElement tElement( ElementType::QUAD4, tTransform );
            fem::ElementMapper tMapper ;
            tMapper.link( tElement.element() );
            tXi.set_size( 2 );

            Vector< real > tX( 2 );
            for ( const Vector< real > & tOut : Cell< Vector< real > >(
                    { { 1.2, 0.0 }, { -1.2, 0.3 }, { 0.0, 1.2 }, { 0.2, -1.2 } } ) )
            {
                tElement.forward( tMapper, tOut, tX );
                EXPECT_FALSE( tMapper.evaluate( tX, tXi ) )
                    << "xi = " << tOut( 0 ) << ", eta = " << tOut( 1 ) ;
            }
        }
    }
    {
        // PYRA5 through the general path: inside_pyra must reject a point
        // below the base, one above the apex and one beside the base
        MappedElement tElement( ElementType::PYRA5, affine_3d );
        fem::ElementMapper tMapper ;
        tMapper.link( tElement.element() );
        tXi.set_size( 3 );

        Vector< real > tX( 3 );
        for ( const Vector< real > & tOut : Cell< Vector< real > >(
                { { 0.0, 0.0, -0.2 }, { 0.0, 0.0, 1.2 }, { 1.3, 0.0, 0.1 } } ) )
        {
            tElement.forward( tMapper, tOut, tX );
            EXPECT_FALSE( tMapper.evaluate( tX, tXi ) )
                << "xi = " << tOut( 0 ) << ", eta = " << tOut( 1 ) << ", zeta = " << tOut( 2 ) ;
        }

        tElement.forward( tMapper, { 0.1, -0.2, 0.3 }, tX );
        EXPECT_TRUE( tMapper.evaluate( tX, tXi ) );
    }
}

//------------------------------------------------------------------------------

TEST( ElementMapper, Tri3JacobianSignIsIndependentOfNodeOrder )
{
    // the same triangle and its mirror image: one counterclockwise, one
    // clockwise. Mirrored points must map to the same natural coordinates
    auto tMirror = []( const real aRef[ 3 ], real aX[ 3 ] )
    {
        affine_2d( aRef, aX );
        aX[ 0 ] = -aX[ 0 ] ;
    };

    MappedElement tA( ElementType::TRI3, affine_2d );
    MappedElement tB( ElementType::TRI3, tMirror );

    fem::ElementMapper tMapperA ;
    fem::ElementMapper tMapperB ;
    tMapperA.link( tA.element() );
    tMapperB.link( tB.element() );

    Vector< real > tX( 2 );
    Vector< real > tXiA( 2 );
    Vector< real > tXiB( 2 );

    for ( const Vector< real > & tXiRef : tri_points() )
    {
        tA.forward( tMapperA, tXiRef, tX );
        ASSERT_TRUE( tMapperA.evaluate( tX, tXiA ) );

        tX( 0 ) = -tX( 0 ) ;
        ASSERT_TRUE( tMapperB.evaluate( tX, tXiB ) );

        for ( uint i = 0; i < 2; ++i )
        {
            EXPECT_NEAR( tXiA( i ), tXiRef( i ), 1e-12 );
            EXPECT_NEAR( tXiB( i ), tXiRef( i ), 1e-12 );
        }
    }
}

//------------------------------------------------------------------------------

TEST( ElementMapperDebug, NonConvergentNewtonThrows )
{
    // a valid QUAD9 whose center node is lifted by c = 0.4, so that
    // y( xi, eta ) = eta + c ( 1 - xi^2 )( 1 - eta^2 ); dy/deta > 0 on the
    // element. On the line x = 0, y is bounded above by c + 1 / ( 4 c ),
    // so the point ( 0, 3 ) has no preimage anywhere and the Newton
    // iterate wanders without converging. This pins today's behavior, an
    // abort at the iteration cap, so that the day it becomes a status
    // return is a visible change; what it must never do is return a point
    auto tLifted = []( const real aRef[ 3 ], real aX[ 3 ] )
    {
        const real c = 0.4 ;
        aX[ 0 ] = aRef[ 0 ] ;
        aX[ 1 ] = aRef[ 1 ] + c * ( 1.0 - aRef[ 0 ] * aRef[ 0 ] ) * ( 1.0 - aRef[ 1 ] * aRef[ 1 ] ) ;
        aX[ 2 ] = 0.0 ;
    };

    MappedElement tElement( ElementType::QUAD9, tLifted );
    fem::ElementMapper tMapper ;
    tMapper.link( tElement.element() );

    // the refusal must be the iteration cap, not some other guard on the
    // way ( a singular Jacobian or a dimension check )
    Vector< real > tXi( 2 );
    bool tThrown = false ;
    try
    {
        tMapper.evaluate( { 0.0, 3.0 }, tXi );
    }
    catch( const std::runtime_error & aError )
    {
        tThrown = true ;
        EXPECT_NE( std::string( aError.what() ).find( "failed to converge" ), std::string::npos )
            << "message was: " << aError.what() ;
    }
    EXPECT_TRUE( tThrown ) << "returned xi = ( " << tXi( 0 ) << ", " << tXi( 1 ) << " )" ;
}
