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
 * Sign contract of mesh::Pipette, and the sign of the thermal volume weight
 * on a thin-shell layer element. Free of Maxwell headers on purpose: this
 * translation unit is compiled unconditionally, while the ghost-stack
 * batteries are not.
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <string>

#include "typedefs.hpp"
#include "commtools.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "fn_det.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_Block.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_ReferenceElement.hpp"
#include "cl_Pipette.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_KernelParameters.hpp"
#include "cl_FEM_DofManagerBase.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_IWG.hpp"
#include "cl_TimestepMatrices.hpp"
#include "en_IWGs.hpp"
#include "en_SolverEnums.hpp"

using namespace belfem;

namespace
{
//------------------------------------------------------------------------------
// C: mesh::Pipette is a SIGNED measure for every type except QUAD4TS
//------------------------------------------------------------------------------

    //! one row of the battery. mSwap holds the transpositions the MeshChecker
    //! applies to this type; its swap methods are private, so they are
    //! reproduced here and must match cl_MeshChecker.hpp. Neither PENTA6TS
    //! nor HEX8TS has a swap of its own there: the checker shares swap_penta6
    //! and swap_hex8 with PENTA6 and HEX8. The reference geometry of every
    //! row, the TS types included, comes from create_unity_nodes for the
    //! row's own type, with no remap
    struct PipetteCase
    {
        ElementType mType ;
        bool        mIs3D ;
        bool        mAbsolute ;
        uint        mNumSwaps ;
        uint        mSwap[ 2 ][ 2 ] ;
    };

    const PipetteCase tPipetteCases[ 8 ] = {
            { ElementType::TRI3,     false, false, 1, { { 1, 2 }, { 0, 0 } } },
            { ElementType::QUAD4,    false, false, 1, { { 1, 3 }, { 0, 0 } } },
            { ElementType::QUAD4TS,  false, true,  1, { { 1, 3 }, { 0, 0 } } },
            { ElementType::TET4,     true,  false, 1, { { 1, 2 }, { 0, 0 } } },
            { ElementType::PENTA6,   true,  false, 2, { { 1, 2 }, { 4, 5 } } },
            { ElementType::PENTA6TS, true,  false, 2, { { 1, 2 }, { 4, 5 } } },
            { ElementType::HEX8,     true,  false, 2, { { 1, 3 }, { 5, 7 } } },
            { ElementType::HEX8TS,   true,  false, 2, { { 1, 3 }, { 5, 7 } } } };

    //! the element under test together with the nodes it owns
    struct PipetteElement
    {
        mesh::Element        * mElement = nullptr ;
        Cell< mesh::Node * >   mNodes ;

        PipetteElement() = default ;

        PipetteElement( const PipetteElement & ) = delete ;
        PipetteElement & operator=( const PipetteElement & ) = delete ;

        ~PipetteElement()
        {
            delete mElement ;
            for ( mesh::Node * tNode : mNodes )
            {
                delete tNode ;
            }
        }
    };

    //! reference element of measure one from ElementFactory::create_unity_nodes
    void
    build_reference_element( const ElementType aType, PipetteElement & aOut )
    {
        mesh::ElementFactory tFactory ;

        mesh::ReferenceElement * tReference =
                tFactory.create_reference_element( aType );

        const uint tNumNodes = tReference->element()->number_of_nodes();

        aOut.mElement = tFactory.create_element( aType, 1 );
        aOut.mNodes.set_size( tNumNodes, nullptr );

        for ( uint k = 0; k < tNumNodes; ++k )
        {
            const mesh::Node * tSource = tReference->node( k );

            mesh::Node * tNode = new mesh::Node( k + 1,
                    tSource->x( 0 ), tSource->x( 1 ), tSource->x( 2 ) );

            aOut.mNodes( k ) = tNode ;
            aOut.mElement->insert_node( tNode, k );
        }

        delete tReference ;
    }

    //! affine map for the two-dimensional types, deliberately not axis
    //! aligned and not symmetric
    Matrix< real >
    pipette_affine_2d()
    {
        Matrix< real > A( 2, 2 );
        A( 0, 0 ) =  1.50 ; A( 0, 1 ) = 0.80 ;
        A( 1, 0 ) = -0.40 ; A( 1, 1 ) = 1.32 ;
        return A ;
    }

    //! the three-dimensional counterpart. The last row is scaled to put the
    //! determinant on 2.3 rather than tuning an entry by hand
    Matrix< real >
    pipette_affine_3d()
    {
        Matrix< real > A( 3, 3 );
        A( 0, 0 ) =  1.50 ; A( 0, 1 ) =  0.80 ; A( 0, 2 ) = 0.30 ;
        A( 1, 0 ) = -0.40 ; A( 1, 1 ) =  1.32 ; A( 1, 2 ) = 0.60 ;
        A( 2, 0 ) =  0.70 ; A( 2, 1 ) = -0.60 ; A( 2, 2 ) = 0.90 ;

        const real tScale = 2.3 / det( A );
        for ( uint j = 0; j < 3; ++j )
        {
            A( 2, j ) *= tScale ;
        }
        return A ;
    }

    real
    pipette_measure( const ElementType aType, mesh::Element * aElement )
    {
        mesh::Pipette tPipette ;
        tPipette.set_element_type( aType );
        return tPipette.measure( aElement );
    }

    void
    pipette_swap_nodes( PipetteElement & aElement, const PipetteCase & aCase )
    {
        for ( uint s = 0; s < aCase.mNumSwaps; ++s )
        {
            const uint i = aCase.mSwap[ s ][ 0 ];
            const uint j = aCase.mSwap[ s ][ 1 ];

            mesh::Node * tSwap = aElement.mElement->node( i );
            aElement.mElement->insert_node( aElement.mElement->node( j ), i );
            aElement.mElement->insert_node( tSwap, j );
        }
    }

//------------------------------------------------------------------------------
// D1: the thermal volume weight on a single thin-shell layer element
//------------------------------------------------------------------------------

    /**
     * one QUAD4TS layer wound the way ThinShellFactory winds them: the offset
     * is the tape direction rotated to ( dy, -dx ), and the node order is
     * bottom0, bottom1, top1, top0. That winding is clockwise, so the scalar
     * map's determinant is negative -- the premise the test establishes
     * before it reads a matrix.
     *
     * The block carries the layer thickness for documentary value only; the
     * thermal dV takes it from the element geometry, not from the block.
     */
    class ShellThermalProblem
    {
        public:

            Mesh                  * mMesh   = nullptr ;
            fem::KernelParameters * mParams = nullptr ;
            fem::Kernel           * mKernel = nullptr ;
            fem::DofManager       * mField  = nullptr ;
            fem::IWG              * mIWG    = nullptr ;   // owned by the kernel

//------------------------------------------------------------------------------

            ShellThermalProblem( const real aLength,
                                 const real aThickness,
                                 const real aAngle )
            {
                // third argument false: the MeshChecker in the Kernel ctor
                // refuses a mesh that already has edges, and connectivity
                // computation inside finalize() would create them
                mMesh = new Mesh( 2, 0, false );

                const real tDx = aLength * std::cos( aAngle );
                const real tDy = aLength * std::sin( aAngle );

                const real tNx =  tDy / aLength ;
                const real tNy = -tDx / aLength ;

                Cell< mesh::Node * > & tNodes = mMesh->nodes();
                tNodes.push( new mesh::Node( 1, 0.0, 0.0 ) );
                tNodes.push( new mesh::Node( 2, tDx, tDy ) );
                tNodes.push( new mesh::Node( 3, tDx + aThickness * tNx,
                                                tDy + aThickness * tNy ) );
                tNodes.push( new mesh::Node( 4, aThickness * tNx,
                                                aThickness * tNy ) );

                mesh::ElementFactory tFactory ;
                mesh::Element * tElement =
                        tFactory.create_element( ElementType::QUAD4TS, 1 );
                for ( uint k = 0; k < 4; ++k )
                {
                    tElement->insert_node( tNodes( k ), k );
                }
                tElement->set_block_id( 1 );
                mMesh->elements().push( tElement );

                // domain type stays Default: a ThinShell block would send
                // Kernel::compute_element_volumes down the facet-area path,
                // which needs a mesh::ThinShell container this fixture has
                // no use for
                mesh::Block * tBlock = new mesh::Block( 1, 1 );
                tBlock->set_thickness( aThickness );
                tBlock->insert_element( tElement );
                mMesh->blocks().push( tBlock );

                mMesh->finalize();

                mParams = new fem::KernelParameters( mMesh );
                mKernel = new fem::Kernel( mParams );

                Vector< id_t > tBlockIDs( 1, 1 );
                Vector< id_t > tSideSetIDs ;

                mIWG = mKernel->create_equation(
                        IwgType::TransientHeatConduction,
                        ModelDimensionality::TwoD,
                        tBlockIDs,
                        tSideSetIDs );

                mField = mKernel->create_field( mIWG );
                mField->set_solver( gDefaultSolver );
                mField->initialize();
            }

//------------------------------------------------------------------------------

            ShellThermalProblem( const ShellThermalProblem & ) = delete ;
            ShellThermalProblem &
            operator=( const ShellThermalProblem & ) = delete ;

            ~ShellThermalProblem()
            {
                delete mKernel ;
                delete mParams ;
                delete mMesh ;
            }
    };
}

//------------------------------------------------------------------------------

TEST( Pipette, MeasuresAreSignedExceptQuad4Ts )
{
    const Matrix< real > tA2 = pipette_affine_2d();
    const Matrix< real > tA3 = pipette_affine_3d();

    ASSERT_NEAR( det( tA2 ), 2.3, 1.0e-12 );
    ASSERT_NEAR( det( tA3 ), 2.3, 1.0e-12 );

    const real tB2[ 2 ] = { 0.37, -0.11 };
    const real tB3[ 3 ] = { 0.37, -0.11, 0.53 };

    for ( const PipetteCase & tCase : tPipetteCases )
    {
        SCOPED_TRACE( to_string( tCase.mType ) );

        PipetteElement tRef ;
        build_reference_element( tCase.mType, tRef );

        // the factory's unity nodes give every reference element measure one
        EXPECT_NEAR( pipette_measure( tCase.mType, tRef.mElement ),
                1.0, 1.0e-12 );

        // affine image: the measure scales with det A and with nothing else
        const real tDetA = tCase.mIs3D ? det( tA3 ) : det( tA2 );

        for ( mesh::Node * tNode : tRef.mNodes )
        {
            const real tX = tNode->x( 0 );
            const real tY = tNode->x( 1 );
            const real tZ = tNode->x( 2 );

            if ( tCase.mIs3D )
            {
                tNode->set_coords(
                        tA3( 0, 0 ) * tX + tA3( 0, 1 ) * tY
                            + tA3( 0, 2 ) * tZ + tB3[ 0 ],
                        tA3( 1, 0 ) * tX + tA3( 1, 1 ) * tY
                            + tA3( 1, 2 ) * tZ + tB3[ 1 ],
                        tA3( 2, 0 ) * tX + tA3( 2, 1 ) * tY
                            + tA3( 2, 2 ) * tZ + tB3[ 2 ] );
            }
            else
            {
                tNode->set_coords(
                        tA2( 0, 0 ) * tX + tA2( 0, 1 ) * tY + tB2[ 0 ],
                        tA2( 1, 0 ) * tX + tA2( 1, 1 ) * tY + tB2[ 1 ] );
            }
        }

        EXPECT_NEAR( pipette_measure( tCase.mType, tRef.mElement ),
                tDetA, 1.0e-12 * tDetA );

        // reversed winding
        pipette_swap_nodes( tRef, tCase );

        if ( tCase.mAbsolute )
        {
            EXPECT_NEAR( pipette_measure( tCase.mType, tRef.mElement ),
                    tDetA, 1.0e-12 * tDetA )
                << "QUAD4TS must report the magnitude on either winding" ;
        }
        else
        {
            EXPECT_NEAR( pipette_measure( tCase.mType, tRef.mElement ),
                    -tDetA, 1.0e-12 * tDetA )
                << "a reversed element must measure negative" ;
        }

        // the permutation is an involution: applied twice it restores the
        // node order, and with it the sign
        pipette_swap_nodes( tRef, tCase );

        for ( uint k = 0; k < tRef.mNodes.size(); ++k )
        {
            EXPECT_EQ( tRef.mElement->node( k )->id(), tRef.mNodes( k )->id() )
                << "node order not restored at slot " << k ;
        }

        EXPECT_NEAR( pipette_measure( tCase.mType, tRef.mElement ),
                tDetA, 1.0e-12 * tDetA );
    }
}

//------------------------------------------------------------------------------

// The oracle is the integrated volume weight read off the assembled mass
// matrix, not a temperature: with the unit coefficient and the partition of
// unity, sum_ij M_ij IS the measure. A uniformly negated system solves to
// the same temperatures but reports -0.2 here.
TEST( Pipette, AllShellThermalMassIsPositive )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    const real tLength    = 2.0 ;
    const real tThickness = 0.1 ;
    const real tAngle     = 0.7 ;
    const real tArea      = tLength * tThickness ;

    ShellThermalProblem tProblem( tLength, tThickness, tAngle );

    // premise: the production winding is clockwise. Signed area from the
    // four corners, the expression mesh::Pipette uses for a plain QUAD4
    mesh::Element * tQuad = tProblem.mMesh->elements()( 0 );
    real tX[ 4 ];
    real tY[ 4 ];
    for ( uint k = 0; k < 4; ++k )
    {
        tX[ k ] = tQuad->node( k )->x( 0 );
        tY[ k ] = tQuad->node( k )->x( 1 );
    }
    const real tSignedArea =
            0.5 * ( ( tY[ 1 ] - tY[ 3 ] ) * ( tX[ 0 ] - tX[ 2 ] )
                  + ( tY[ 0 ] - tY[ 2 ] ) * ( tX[ 3 ] - tX[ 1 ] ) );

    ASSERT_LT( tSignedArea, 0.0 )
        << "fixture no longer reproduces the production clockwise winding" ;
    ASSERT_NEAR( std::abs( tSignedArea ), tArea, 1.0e-12 * tArea );

    fem::Block * tBlock = tProblem.mField->block( 1 );
    ASSERT_TRUE( tBlock != nullptr );
    ASSERT_EQ( tBlock->number_of_elements(), 1u );

    // compute_mkf reads the calculator off the group, and the element
    // matrices are sized by the caller: neither the IWG constructor nor
    // link_to_group does it
    tProblem.mIWG->link_to_group( tBlock );

    fem::TimestepMatrices * tMatrices = tProblem.mIWG->matrices();
    tMatrices->set_flag( fem::MatrixFlag::M );
    tMatrices->set_flag( fem::MatrixFlag::K );
    tMatrices->set_flag( fem::MatrixFlag::F );
    tMatrices->initialize( mesh::number_of_nodes( ElementType::QUAD4TS ) );

    tProblem.mIWG->compute_mkf( tBlock->elements()( 0 ) );

    const Matrix< real > & M = tMatrices->M();
    const Matrix< real > & K = tMatrices->K();

    ASSERT_EQ( M.n_rows(), 4u );
    ASSERT_EQ( M.n_cols(), 4u );

    real tSumM = 0.0 ;
    for ( uint i = 0; i < 4; ++i )
    {
        EXPECT_GT( M( i, i ), 0.0 )
            << "mass diagonal " << i << " is not positive" ;

        for ( uint j = 0; j < 4; ++j )
        {
            tSumM += M( i, j );
        }
    }

    EXPECT_NEAR( tSumM, tArea, 1.0e-12 * tArea )
        << "integrated thermal volume weight " << tSumM
        << " vs thickness times length " << tArea ;

    real tMaxK = 0.0 ;
    for ( uint i = 0; i < 4; ++i )
    {
        for ( uint j = 0; j < 4; ++j )
        {
            tMaxK = std::max( tMaxK, std::abs( K( i, j ) ) );
        }
    }
    ASSERT_GT( tMaxK, 0.0 );

    // a constant field is in the Laplacian's kernel: the Cartesian gradient
    // operator's node contributions cancel row by row
    for ( uint i = 0; i < 4; ++i )
    {
        real tRowSum = 0.0 ;
        for ( uint j = 0; j < 4; ++j )
        {
            tRowSum += K( i, j );
        }
        EXPECT_LT( std::abs( tRowSum ), 1.0e-12 * tMaxK )
            << "stiffness row " << i << " does not annihilate a constant" ;
    }
}
