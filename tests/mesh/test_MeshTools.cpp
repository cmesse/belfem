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
#include <cmath>
#include <stdexcept>
#include <string>
#include <memory>

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "cl_Node.hpp"
#include "cl_Element.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"
#include "meshtools.hpp"
#include "Mesh_Enums.hpp"

using namespace belfem;

namespace
{
    const real tEpsilon = 1.0e-12 ;

    //! an element whose nodes sit at the parametric node coordinates of its
    //! Lagrange interpolant, so the through-thickness coordinate of every
    //! node is known: eta in 2-D, zeta in 3-D, both in [ -1, 1 ]. Owns the
    //! element and the nodes; releases both on every exit, including a
    //! throw out of the test body
    class ParametricElement
    {
        public:

            mesh::Element        * mElement ;
            Cell< mesh::Node * >   mNodes ;
            uint                   mDimension ;

//------------------------------------------------------------------------------

            explicit ParametricElement( const ElementType aType )
            {
                mesh::ElementFactory              tElementFactory ;
                fem::InterpolationFunctionFactory tFunctionFactory ;

                mDimension = mesh::dimension( aType );

                // the coordinate table is column-per-node, rows xi/eta[/zeta];
                // the interpolant is read before anything is owned, so a
                // throwing factory leaks nothing
                Matrix< real > tXi ;
                {
                    std::unique_ptr< fem::InterpolationFunction > tFunction(
                            tFunctionFactory.create_lagrange_function( aType ) );
                    tFunction->param_coords( tXi );
                }

                mElement = tElementFactory.create_element( aType, 1 );

                const uint tNumNodes = mElement->number_of_nodes();
                EXPECT_EQ( tXi.n_cols(), tNumNodes );
                EXPECT_EQ( tXi.n_rows(), mDimension );

                mNodes.set_size( tNumNodes, nullptr );
                for ( uint k = 0; k < tNumNodes; ++k )
                {
                    mNodes( k ) = new mesh::Node(
                            k + 1,
                            tXi( 0, k ),
                            mDimension > 1 ? tXi( 1, k ) : 0.0,
                            mDimension > 2 ? tXi( 2, k ) : 0.0 );
                    mElement->insert_node( mNodes( k ), k );
                }
            }

//------------------------------------------------------------------------------

            ~ParametricElement()
            {
                delete mElement ;
                for ( uint k = 0; k < mNodes.size(); ++k )
                {
                    delete mNodes( k );
                }
            }

            // owning type: a copy would double-free
            ParametricElement( const ParametricElement & ) = delete ;
            ParametricElement & operator=( const ParametricElement & ) = delete ;

//------------------------------------------------------------------------------

            real
            thickness_coordinate( const mesh::Node * aNode ) const
            {
                return mDimension == 2 ? aNode->y() : aNode->z() ;
            }

//------------------------------------------------------------------------------

            bool
            owns( const mesh::Node * aNode ) const
            {
                for ( uint k = 0; k < mNodes.size(); ++k )
                {
                    if ( mNodes( k ) == aNode ) return true ;
                }
                return false ;
            }

//------------------------------------------------------------------------------

            uint
            count_at( const real aThicknessCoordinate ) const
            {
                uint tCount = 0 ;
                for ( uint k = 0; k < mNodes.size(); ++k )
                {
                    if ( std::abs( this->thickness_coordinate( mNodes( k ) )
                                   - aThicknessCoordinate ) < tEpsilon )
                    {
                        ++tCount ;
                    }
                }
                return tCount ;
            }
    };

//------------------------------------------------------------------------------

    //! every pointer in the list occurs once
    bool
    all_distinct( const Cell< mesh::Node * > & aNodes )
    {
        for ( uint i = 0; i < aNodes.size(); ++i )
        {
            for ( uint j = i + 1; j < aNodes.size(); ++j )
            {
                if ( aNodes( i ) == aNodes( j ) ) return false ;
            }
        }
        return true ;
    }

//------------------------------------------------------------------------------

    //! one test case per element type, so a throw out of one type cannot
    //! hide the verdict on the others
    class MeshToolsBottomTop : public ::testing::TestWithParam< ElementType >
    {
    };

    // to_string() folds several thin-shell variants onto their base type and
    // has no name for others, so the enum value is appended to keep the
    // parameter names unique
    std::string
    type_name( const ::testing::TestParamInfo< ElementType > & aInfo )
    {
        return to_string( aInfo.param ) + "_"
               + std::to_string( static_cast< unsigned >( aInfo.param ) );
    }
}

//------------------------------------------------------------------------------

// the two helpers return the bottom and top face of a stacked element as
// positionally paired lists: top( k ) sits straight above bottom( k ). The
// top list comes from the element's facet table, the bottom list is
// hard-coded per type, so every case label of get_bottom_nodes is checked
// against the geometry of a parametric element
TEST_P( MeshToolsBottomTop, PositionalPairing )
{
    ParametricElement tCase( GetParam() );

    Cell< mesh::Node * > tBottom ;
    Cell< mesh::Node * > tTop ;
    mesh::get_bottom_nodes( tCase.mElement, tBottom );
    mesh::get_top_nodes( tCase.mElement, tTop );

    // the bottom face is exactly the set of nodes at thickness coordinate -1
    const uint tExpected = tCase.count_at( -1.0 );
    ASSERT_GT( tExpected, 0u );
    ASSERT_EQ( tBottom.size(), tExpected );
    ASSERT_EQ( tTop.size(), tExpected );
    EXPECT_TRUE( all_distinct( tBottom ) );
    EXPECT_TRUE( all_distinct( tTop ) );

    for ( uint k = 0; k < tExpected; ++k )
    {
        ASSERT_NE( tBottom( k ), nullptr ) << "bottom " << k ;
        ASSERT_NE( tTop( k ), nullptr )    << "top " << k ;
        ASSERT_TRUE( tCase.owns( tBottom( k ) ) ) << "bottom " << k ;
        ASSERT_TRUE( tCase.owns( tTop( k ) ) )    << "top " << k ;

        EXPECT_NEAR( tCase.thickness_coordinate( tBottom( k ) ), -1.0, tEpsilon )
            << "bottom " << k ;
        EXPECT_NEAR( tCase.thickness_coordinate( tTop( k ) ), 1.0, tEpsilon )
            << "top " << k ;

        // positional pairing: same in-plane position, +2 through the thickness
        EXPECT_NEAR( tTop( k )->x() - tBottom( k )->x(), 0.0, tEpsilon )
            << "pair " << k ;
        if ( tCase.mDimension == 2 )
        {
            EXPECT_NEAR( tTop( k )->y() - tBottom( k )->y(), 2.0, tEpsilon )
                << "pair " << k ;
        }
        else
        {
            EXPECT_NEAR( tTop( k )->y() - tBottom( k )->y(), 0.0, tEpsilon )
                << "pair " << k ;
            EXPECT_NEAR( tTop( k )->z() - tBottom( k )->z(), 2.0, tEpsilon )
                << "pair " << k ;
        }
    }
}

// the complete case list of get_bottom_nodes
INSTANTIATE_TEST_SUITE_P(
        MeshTools,
        MeshToolsBottomTop,
        ::testing::Values(
                ElementType::QUAD4,
                ElementType::QUAD4TS,
                ElementType::QUAD9,
                ElementType::QUAD9TS,
                ElementType::PENTA6,
                ElementType::PENTA6TS,
                ElementType::PENTA15,
                ElementType::PENTA18,
                ElementType::PENTA18TS,
                ElementType::HEX8,
                ElementType::HEX8TS,
                ElementType::HEX8TB,
                ElementType::HEX20,
                ElementType::HEX27 ),
        type_name );

//------------------------------------------------------------------------------

// element types without a bottom/top face, and types the bottom helper does
// not carry, are refused by both helpers, not silently mapped
TEST( MeshTools, BottomTopNodesRejectUnsupported )
{
    for ( const ElementType tType : { ElementType::TRI3,
                                      ElementType::TET4,
                                      ElementType::QUAD8,
                                      ElementType::QUAD16,
                                      ElementType::HEX64 } )
    {
        SCOPED_TRACE( to_string( tType ) );

        ParametricElement tCase( tType );
        Cell< mesh::Node * > tNodes ;

        EXPECT_THROW( mesh::get_bottom_nodes( tCase.mElement, tNodes ), std::runtime_error );
        EXPECT_THROW( mesh::get_top_nodes( tCase.mElement, tNodes ), std::runtime_error );
    }
}
