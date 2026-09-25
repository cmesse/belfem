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
 * The facet and edge tables of every Lagrange element type, checked against
 * the geometry of the reference element rather than against a copy of the
 * table. The nodes are placed at their parameter coordinates, so the nodes
 * of facet f are exactly the element nodes inside the convex hull of the
 * facet's corner nodes, and the nodes of edge e exactly those on the
 * segment between its first two. A table that lists the right nodes in a
 * different order passes; a table that lists a wrong node, a node twice or
 * a null fails and names the type, the facet and both sets. Types the two
 * factories refuse are skipped and counted.
 */

#include <gtest/gtest.h>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "Mesh_Enums.hpp"
#include "cl_Node.hpp"
#include "cl_Element.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"
#include "cl_IF_InterpolationFunction.hpp"

using namespace belfem ;
using namespace belfem::mesh ;

namespace
{
    const real tHullTol = 1.0e-9 ;

    struct Point
    {
        real mX ;
        real mY ;
        real mZ ;
    };

    Point
    sub( const Point & aA, const Point & aB )
    {
        return { aA.mX - aB.mX, aA.mY - aB.mY, aA.mZ - aB.mZ };
    }

    real
    dot( const Point & aA, const Point & aB )
    {
        return aA.mX * aB.mX + aA.mY * aB.mY + aA.mZ * aB.mZ ;
    }

    Point
    cross( const Point & aA, const Point & aB )
    {
        return { aA.mY * aB.mZ - aA.mZ * aB.mY,
                 aA.mZ * aB.mX - aA.mX * aB.mZ,
                 aA.mX * aB.mY - aA.mY * aB.mX };
    }

    real
    norm( const Point & aA )
    {
        return std::sqrt( dot( aA, aA ) );
    }

    // aP on the closed segment aA-aB
    bool
    on_segment( const Point & aP, const Point & aA, const Point & aB )
    {
        const Point tD = sub( aB, aA );
        const Point tQ = sub( aP, aA );
        const real tT = dot( tQ, tD ) / dot( tD, tD );
        const Point tR = { tQ.mX - tT * tD.mX, tQ.mY - tT * tD.mY, tQ.mZ - tT * tD.mZ };
        return norm( tR ) < tHullTol && tT > -tHullTol && tT < 1.0 + tHullTol ;
    }

    // aP in the closed triangle aA-aB-aC ( coplanar and inside )
    bool
    in_triangle( const Point & aP, const Point & aA, const Point & aB, const Point & aC )
    {
        const Point tN = cross( sub( aB, aA ), sub( aC, aA ) );
        const real tNN = norm( tN );
        if ( tNN < tHullTol ) return false ;
        if ( std::abs( dot( sub( aP, aA ), tN ) ) / tNN > tHullTol ) return false ;
        const real tS0 = dot( cross( sub( aB, aA ), sub( aP, aA ) ), tN );
        const real tS1 = dot( cross( sub( aC, aB ), sub( aP, aB ) ), tN );
        const real tS2 = dot( cross( sub( aA, aC ), sub( aP, aC ) ), tN );
        const real tTol = -tHullTol * tNN * tNN ;
        return tS0 >= tTol && tS1 >= tTol && tS2 >= tTol ;
    }

    // the convex hull of two, three or four corner points ( a facet of a
    // reference element is a segment, a triangle or a planar quad )
    bool
    in_hull( const Point & aP, const std::vector< Point > & aCorners )
    {
        switch ( aCorners.size() )
        {
            case 2 : return on_segment( aP, aCorners[ 0 ], aCorners[ 1 ] );
            case 3 : return in_triangle( aP, aCorners[ 0 ], aCorners[ 1 ], aCorners[ 2 ] );
            case 4 : return in_triangle( aP, aCorners[ 0 ], aCorners[ 1 ], aCorners[ 2 ] )
                         || in_triangle( aP, aCorners[ 0 ], aCorners[ 2 ], aCorners[ 3 ] );
            default : return false ;
        }
    }

    std::string
    set_to_string( const std::set< uint > & aSet )
    {
        std::string tOut = "{" ;
        for ( uint k : aSet ) tOut += " " + std::to_string( k );
        return tOut + " }" ;
    }

    std::string
    type_name( const ElementType aType )
    {
        return to_string( aType ) + " ( " + std::to_string( static_cast< uint >( aType ) ) + " )" ;
    }

    // an element of aType with node k at its parameter coordinates and
    // id k + 1; the fixture owns nodes and element
    struct ReferenceElement
    {
        Element * mElement = nullptr ;
        Cell< Node * > mNodes ;
        std::vector< Point > mPoints ;

        ReferenceElement( const ElementType aType )
        {
            fem::InterpolationFunctionFactory tIfFactory ;
            fem::InterpolationFunction * tFun = tIfFactory.create_lagrange_function( aType );
            Matrix< real > tXi ;
            tFun->param_coords( tXi );
            delete tFun ;

            const uint tDim = tXi.n_rows();
            const uint tN   = tXi.n_cols();

            ElementFactory tFactory ;
            mElement = tFactory.create_element( aType, 1 );

            mPoints.resize( tN );
            mNodes.set_size( tN, nullptr );
            for ( uint k = 0; k < tN; ++k )
            {
                mPoints[ k ] = { tXi( 0, k ),
                                 tDim > 1 ? tXi( 1, k ) : 0.0,
                                 tDim > 2 ? tXi( 2, k ) : 0.0 };
                mNodes( k ) = new Node( k + 1, mPoints[ k ].mX, mPoints[ k ].mY, mPoints[ k ].mZ );
                mElement->insert_node( mNodes( k ), k );
            }
        }

        ~ReferenceElement()
        {
            delete mElement ;
            for ( Node * tNode : mNodes ) delete tNode ;
        }

        ReferenceElement( const ReferenceElement & ) = delete ;
        ReferenceElement & operator=( const ReferenceElement & ) = delete ;
        ReferenceElement( ReferenceElement && ) = delete ;
        ReferenceElement & operator=( ReferenceElement && ) = delete ;

        uint
        number_of_nodes() const
        {
            return mPoints.size();
        }

        // the indices of the nodes inside the hull of aCorners
        std::set< uint >
        expected( const Cell< Node * > & aCorners ) const
        {
            std::vector< Point > tCorners ;
            for ( Node * tCorner : aCorners )
            {
                tCorners.push_back( mPoints[ tCorner->id() - 1 ] );
            }
            std::set< uint > tOut ;
            for ( uint k = 0; k < mPoints.size(); ++k )
            {
                if ( in_hull( mPoints[ k ], tCorners ) ) tOut.insert( k );
            }
            return tOut ;
        }
    };

    // true when aTable holds exactly aExpected, each node once, no null;
    // aWhat names the table for the failure text
    void
    expect_table( const std::string & aWhat,
                  const Cell< Node * > & aTable,
                  const std::set< uint > & aExpected )
    {
        std::set< uint > tActual ;
        bool tDuplicate = false ;
        bool tNull = false ;
        for ( Node * tNode : aTable )
        {
            if ( tNode == nullptr )
            {
                tNull = true ;
                continue ;
            }
            if ( ! tActual.insert( tNode->id() - 1 ).second ) tDuplicate = true ;
        }
        EXPECT_FALSE( tNull ) << aWhat << ": null node in the table" ;
        EXPECT_FALSE( tDuplicate ) << aWhat << ": a node is listed twice" ;
        EXPECT_EQ( tActual, aExpected ) << aWhat << ": table " << set_to_string( tActual )
                                        << " but the hull holds " << set_to_string( aExpected );
    }

    // returns false when a factory refuses the type ( the type has no
    // Lagrange function or no element ), so a family test can count skips
    bool
    check_type( const ElementType aType )
    {
        try
        {
            ReferenceElement tRef( aType );
            Element * tElement = tRef.mElement ;
            const std::string tName = type_name( aType );

            EXPECT_EQ( tElement->number_of_nodes(), tRef.number_of_nodes() )
                << tName << ": element node count differs from the interpolation's" ;

            Cell< Node * > tCorners ;
            Cell< Node * > tTable ;

            for ( uint f = 0; f < tElement->number_of_facets(); ++f )
            {
                tElement->get_corner_nodes_of_facet( f, tCorners );
                tElement->get_nodes_of_facet( f, tTable );
                expect_table( tName + " facet " + std::to_string( f ),
                              tTable, tRef.expected( tCorners ) );
            }

            // an edge table lists its two end nodes first; a type without
            // edge tables throws on the first index and is done, a throw on
            // a later index is a broken table
            for ( uint e = 0; e < tElement->number_of_edges(); ++e )
            {
                try
                {
                    tElement->get_nodes_of_edge( e, tTable );
                }
                catch ( const std::runtime_error & aError )
                {
                    if ( e == 0 ) break ;
                    ADD_FAILURE() << tName << " edge " << e << " throws although edge 0 did not: " << aError.what();
                    continue ;
                }
                EXPECT_GE( tTable.size(), 2u ) << tName << " edge " << e ;
                if ( tTable.size() < 2 ) continue ;
                tCorners.set_size( 2, nullptr );
                tCorners( 0 ) = tTable( 0 );
                tCorners( 1 ) = tTable( 1 );
                expect_table( tName + " edge " + std::to_string( e ),
                              tTable, tRef.expected( tCorners ) );
            }
            return true ;
        }
        catch ( const std::runtime_error & )
        {
            return false ;
        }
    }

    // every type of a family must have both a Lagrange function and an
    // element; aSkipsAllowed types may be refused
    void
    check_family( const std::vector< ElementType > & aTypes, const uint aSkipsAllowed = 0 )
    {
        uint tSkipped = 0 ;
        for ( ElementType tType : aTypes )
        {
            if ( ! check_type( tType ) ) ++tSkipped ;
        }
        EXPECT_LE( tSkipped, aSkipsAllowed ) << tSkipped << " types of the family were refused by a factory" ;
    }
}

//------------------------------------------------------------------------------

TEST( FacetTables, Lines )
{
    check_family( { ElementType::LINE2, ElementType::LINE3, ElementType::LINE4, ElementType::LINE5 } );
}

TEST( FacetTables, Triangles )
{
    check_family( { ElementType::TRI3, ElementType::TRI6, ElementType::TRI10, ElementType::TRI15 } );
}

TEST( FacetTables, Quads )
{
    check_family( { ElementType::QUAD4, ElementType::QUAD8, ElementType::QUAD9, ElementType::QUAD16 } );
}

TEST( FacetTables, Tets )
{
    check_family( { ElementType::TET4, ElementType::TET10, ElementType::TET20, ElementType::TET35 } );
}

TEST( FacetTables, Prisms )
{
    check_family( { ElementType::PENTA6, ElementType::PENTA15, ElementType::PENTA18 } );
}

TEST( FacetTables, Pyramids )
{
    check_family( { ElementType::PYRA5, ElementType::PYRA13, ElementType::PYRA14 } );
}

TEST( FacetTables, Hexes )
{
    check_family( { ElementType::HEX8, ElementType::HEX20, ElementType::HEX27, ElementType::HEX64 } );
}

TEST( FacetTables, ThinShellTypes )
{
    check_family( { ElementType::QUAD4TS, ElementType::QUAD9TS,
                    ElementType::PENTA6TS, ElementType::PENTA18TS,
                    ElementType::HEX8TS, ElementType::HEX8TB } );
}
