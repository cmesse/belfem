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
 * Column bound of IWG::collect_node_coords: the output carries exactly the
 * problem dimension the IWG constructor set, never the element's own
 * dimension and never one column more. Free of Maxwell headers on purpose:
 * this translation unit is compiled unconditionally.
 */

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "commtools.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "cl_Node.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_IWG.hpp"
#include "cl_TimestepMatrices.hpp"
#include "cl_IWG_TransientHeatConduction.hpp"
#include "en_IWGs.hpp"

using namespace belfem;

namespace
{
    /**
     * collect_node_coords is protected. The transient heat conduction
     * constructor is the cheapest seam that sets mNumberOfSpatialDimensions:
     * it needs no kernel, no mesh and no field.
     */
    struct CoordProbe : public fem::IWG_TransientHeatConduction
    {
        explicit CoordProbe( const ModelDimensionality aDimensionality ) :
                fem::IWG_TransientHeatConduction( aDimensionality )
        {
        }

        using fem::IWG::collect_node_coords ;
    };

    /**
     * the element under test together with the nodes it owns. The fem
     * element is built through the test constructor, whose destructor frees
     * only its own buffers -- the mesh element and the nodes are this
     * fixture's to delete.
     */
    struct ProbeElement
    {
        mesh::Element        * mElement = nullptr ;
        fem::Element         * mFemElement = nullptr ;
        Cell< mesh::Node * >   mNodes ;

        ProbeElement(
                const ElementType aType,
                const uint        aNumberOfNodes,
                const real        aCoords[][ 3 ] )
        {
            mesh::ElementFactory tFactory ;

            mElement = tFactory.create_element( aType, 1 );
            mNodes.set_size( aNumberOfNodes, nullptr );

            for ( uint k = 0; k < aNumberOfNodes; ++k )
            {
                mesh::Node * tNode = new mesh::Node( k + 1,
                        aCoords[ k ][ 0 ], aCoords[ k ][ 1 ], aCoords[ k ][ 2 ] );

                mNodes( k ) = tNode ;
                mElement->insert_node( tNode, k );
            }

            mFemElement = new fem::Element( mElement );
        }

        ProbeElement( const ProbeElement & ) = delete ;
        ProbeElement & operator=( const ProbeElement & ) = delete ;

        ~ProbeElement()
        {
            delete mFemElement ;
            delete mElement ;
            for ( mesh::Node * tNode : mNodes )
            {
                delete tNode ;
            }
        }
    };

    // integer coordinates so the copy is bit exact, all distinct and
    // asymmetric so a transposed or shifted copy cannot match by accident
    const real tTet4Coords[ 4 ][ 3 ] = {
            {  1.0,  2.0,  3.0 },
            {  4.0,  5.0,  6.0 },
            {  7.0,  9.0, 11.0 },
            { 13.0, 17.0, 19.0 } };

    // the z values are nonzero on purpose: a third column written on a 2-D
    // problem would carry them
    const real tTri3Coords[ 3 ][ 3 ] = {
            {  1.0,  2.0,  7.0 },
            {  4.0,  5.0,  8.0 },
            {  9.0, 11.0, 13.0 } };
}

//------------------------------------------------------------------------------

TEST( IwgNodeCoords, ThreeDVolumeElementFillsThreeColumns )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    ProbeElement tProbe( ElementType::TET4, 4, tTet4Coords );
    CoordProbe tIwg( ModelDimensionality::ThreeD );

    Matrix< real > tX ;
    tIwg.collect_node_coords( tProbe.mFemElement, tX );

    EXPECT_EQ( tX.n_rows(), 4u );
    EXPECT_EQ( tX.n_cols(), 3u );

    for ( uint i = 0; i < 3; ++i )
    {
        for ( uint k = 0; k < 4; ++k )
        {
            EXPECT_EQ( tX( k, i ), tTet4Coords[ k ][ i ] );
        }
    }
}

//------------------------------------------------------------------------------

TEST( IwgNodeCoords, TwoDElementFillsTwoColumnsAndResizesAnOversizedOutput )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    ProbeElement tProbe( ElementType::TRI3, 3, tTri3Coords );
    CoordProbe tIwg( ModelDimensionality::TwoD );

    Matrix< real > tX ;
    tIwg.collect_node_coords( tProbe.mFemElement, tX );

    EXPECT_EQ( tX.n_rows(), 3u );
    EXPECT_EQ( tX.n_cols(), 2u );

    for ( uint i = 0; i < 2; ++i )
    {
        for ( uint k = 0; k < 3; ++k )
        {
            EXPECT_EQ( tX( k, i ), tTri3Coords[ k ][ i ] );
        }
    }

    // an output that is already larger is sized down, not filled in place.
    // set_size does not preserve, so only the written entries are asserted.
    Matrix< real > tBig( 5, 3, BELFEM_QUIET_NAN );
    tIwg.collect_node_coords( tProbe.mFemElement, tBig );

    EXPECT_EQ( tBig.n_rows(), 3u );
    EXPECT_EQ( tBig.n_cols(), 2u );

    for ( uint i = 0; i < 2; ++i )
    {
        for ( uint k = 0; k < 3; ++k )
        {
            EXPECT_EQ( tBig( k, i ), tTri3Coords[ k ][ i ] );
        }
    }
}

//------------------------------------------------------------------------------

TEST( IwgNodeCoords, ColumnCountIsTheProblemDimensionNotTheElementDimension )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    // a surface element in a three-dimensional problem: the column count
    // follows mNumberOfSpatialDimensions, not mesh::dimension( TRI3 )
    ProbeElement tProbe( ElementType::TRI3, 3, tTri3Coords );
    CoordProbe tIwg( ModelDimensionality::ThreeD );

    Matrix< real > tX ;
    tIwg.collect_node_coords( tProbe.mFemElement, tX );

    EXPECT_EQ( tX.n_rows(), 3u );
    EXPECT_EQ( tX.n_cols(), 3u );

    for ( uint i = 0; i < 3; ++i )
    {
        for ( uint k = 0; k < 3; ++k )
        {
            EXPECT_EQ( tX( k, i ), tTri3Coords[ k ][ i ] );
        }
    }
}
