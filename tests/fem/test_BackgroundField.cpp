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
 * The background-field sideset kernel ( maxwell::background_phi ) must
 * weight every integration point with that point's own surface increment.
 * The fixture is one HEX8 air element whose bottom face is not a
 * parallelogram, so the increment of normal_hex varies across the face.
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>

#include "typedefs.hpp"
#include "commtools.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "fn_trans.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_Facet.hpp"
#include "cl_Block.hpp"
#include "cl_SideSet.hpp"
#include "cl_Element_Factory.hpp"
#include "meshtools.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_KernelParameters.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_SideSet.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Calculator.hpp"
#include "cl_FEM_PhysicalBoundaryCondition.hpp"
#include "cl_MeshChecker.hpp"
#include "cl_IWG_Maxwell.hpp"
#include "cl_TimestepMatrices.hpp"
#include "en_Maxwell_Formulations.hpp"
#include "en_SolverEnums.hpp"
#include "mt_maxwell_background.hpp"

using namespace belfem;

namespace belfem
{
    namespace fem
    {
        namespace
        {
            // bottom-face footprint, counterclockwise seen from +z; not a
            // parallelogram, since ( 0,0 ) + ( 0.9,1.1 ) != ( 1.3,0.1 ) + ( 0,0.8 )
            const real gFootX[ 4 ] = { 0.0, 1.3, 0.9, 0.0 };
            const real gFootY[ 4 ] = { 0.0, 0.1, 1.1, 0.8 };
            const real gHeight     = 0.5 ;

            const id_t gSideSetID  = 10 ;

//------------------------------------------------------------------------------

            /**
             * One HEX8 air element ( block 1 ) and one QUAD4 boundary facet on
             * its bottom face ( HEX8 facet 4, zeta = -1, no slave ) in a
             * BackgroundField sideset, linked through a real
             * Kernel / DofManager / IWG_Maxwell ( h-phi ) chain.
             *
             * Ownership: the kernel deletes the dof manager, the IWG and every
             * registered boundary condition; the fixture deletes the kernel,
             * the parameters and the mesh; the mesh owns nodes, blocks,
             * elements, sidesets and facets.
             *
             * No edge objects are created: an air block carries only nodal
             * phi dofs, so neither the dof count nor the facet-and-master link
             * reads a mesh edge, and the master edge function only reads the
             * ( default ) edge-direction bitset.
             */
            class BackgroundHexStack
            {
                Mesh                      * mMesh ;
                KernelParameters          * mParams ;
                Kernel                    * mKernel ;
                DofManager                * mField ;  // owned by the kernel
                IWG_Maxwell               * mIWG ;    // owned by the kernel
                PhysicalBoundaryCondition * mBC ;     // owned by the kernel
                SideSet                   * mGroup ;  // owned by the field

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                explicit
                BackgroundHexStack( const real aFieldValue )
                {
                    // no connectivities: the MeshChecker must run on a
                    // mesh without edges
                    mMesh = new Mesh( 3, 0, false );

                    mesh::ElementFactory tFactory ;

                    Cell< mesh::Node * > & tNodes = mMesh->nodes();
                    id_t tID = 1 ;
                    for ( uint r = 0; r < 2; ++r )
                    {
                        for ( uint k = 0; k < 4; ++k )
                        {
                            tNodes.push( new mesh::Node( tID++,
                                    gFootX[ k ], gFootY[ k ], gHeight * r ) );
                        }
                    }

                    mesh::Block * tBlock = new mesh::Block( 1, 1 );
                    tBlock->set_domain_type( DomainType::Air );
                    mesh::Element * tHex =
                            tFactory.create_element( ElementType::HEX8, 1 );
                    for ( uint k = 0; k < 8; ++k )
                    {
                        tHex->insert_node( tNodes( k ), k );
                    }
                    tHex->set_block_id( 1 );
                    mMesh->elements().push( tHex );
                    tBlock->insert_element( tHex );
                    mMesh->blocks().push( tBlock );

                    // facet nodes in HEX8 facet-4 order: 0, 3, 2, 1
                    const uint tFacetNodes[ 4 ] = { 0, 3, 2, 1 };
                    const uint tFacetIndex =
                            mesh::bottom_facet_index( ElementType::HEX8 );

                    mesh::Element * tSurf =
                            tFactory.create_element( ElementType::QUAD4, 100 );
                    for ( uint k = 0; k < 4; ++k )
                    {
                        tSurf->insert_node( tNodes( tFacetNodes[ k ] ), k );
                    }
                    mesh::Facet * tFacet = new mesh::Facet( tSurf );
                    tFacet->set_master( tHex, tFacetIndex );
                    tFacet->set_sideset_id( gSideSetID );

                    mesh::SideSet * tSideSet = new mesh::SideSet( gSideSetID, 1 );
                    tSideSet->set_domain_type( DomainType::BackgroundField );
                    tSideSet->insert_facet( tFacet );
                    mMesh->sidesets().push( tSideSet );

                    mMesh->finalize();

                    // sets the checker flag, so the Kernel ctor skips its run
                    {
                        MeshChecker tCheck( mMesh );
                    }

                    mParams = new KernelParameters( mMesh );
                    mKernel = new Kernel( mParams );

                    mIWG = new IWG_Maxwell(
                            maxwell::Formulation::HPhi,
                            ModelDimensionality::ThreeD,
                            false, false );

                    Vector< id_t > tBlockIDs( 1, 1 );
                    Cell< DomainType > tBlockTypes( 1, DomainType::Air );
                    mIWG->set_blocks( tBlockIDs, tBlockTypes );

                    Vector< id_t > tSideSetIDs( 1, gSideSetID );
                    Cell< DomainType > tSideSetTypes( 1, DomainType::BackgroundField );
                    mIWG->set_sidesets( tSideSetIDs, tSideSetTypes );

                    mKernel->add_equation( mIWG );

                    // H = value * direction, read by background_phi through
                    // Kernel::boundary_condition( sideset id )
                    mBC = new PhysicalBoundaryCondition( 0 );
                    mBC->set_type( BoundaryConditionType::Background );
                    mBC->set_direction( 0.0, 0.0, 1.0 );
                    mBC->set_domains( tSideSetIDs );
                    mBC->fix( aFieldValue );
                    mKernel->add_boundary_condition( mBC );

                    mField = mKernel->create_field( mIWG );
                    mField->create_fields( mField->iwg() );
                    mField->set_solver( gDefaultSolver );
                    mField->initialize();

                    mGroup = mField->sideset( gSideSetID );

                    // binds background_phi; the calculator's "h" vector was
                    // created with the sideset calculator
                    mIWG->link_to_group( mGroup );
                }

                BackgroundHexStack( const BackgroundHexStack & ) = delete ;
                BackgroundHexStack &
                operator=( const BackgroundHexStack & ) = delete ;
                BackgroundHexStack( BackgroundHexStack && ) = delete ;
                BackgroundHexStack &
                operator=( BackgroundHexStack && ) = delete ;

                ~BackgroundHexStack()
                {
                    delete mKernel ;
                    delete mParams ;
                    delete mMesh ;
                }

                Calculator *
                calculator() { return mGroup->calculator() ; }

                fem::Element *
                facet_element() { return mGroup->elements()( 0 ) ; }

                const Vector< real > &
                direction() { return mBC->direction() ; }

                real
                field_value() { return mBC->value() ; }
            };

//------------------------------------------------------------------------------

            real
            absolute_footprint_area()
            {
                real tTwiceArea = 0.0 ;
                for ( uint k = 0; k < 4; ++k )
                {
                    const uint l = ( k + 1 ) % 4 ;
                    tTwiceArea += gFootX[ k ] * gFootY[ l ] - gFootX[ l ] * gFootY[ k ];
                }
                return 0.5 * std::abs( tTwiceArea );
            }

//------------------------------------------------------------------------------
        }
    }
}

//------------------------------------------------------------------------------

// Part A: the fixture can see the defect. On a planar bilinear face the
// increment is linear in ( xi, eta ), so the quad rule ( weights summing to 4 )
// integrates it to the absolute face area; and the increments differ between
// integration points, so a kernel that reuses the point-0 increment is
// distinguishable from one that does not.
TEST( BackgroundField, SurfaceIncrementVariesOnHexFace )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    fem::BackgroundHexStack tStack( 1.0 );
    fem::Calculator * tCalc = tStack.calculator();
    tCalc->link( tStack.facet_element() );

    const Vector< real > & w = tCalc->integration()->weights();
    const uint tNumPoints = tCalc->num_intpoints();
    ASSERT_GT( tNumPoints, 1u );

    const real tdS0 = tCalc->dS( 0 );
    real tArea = 0.0 ;
    real tMaxDeviation = 0.0 ;
    for ( uint k = 0; k < tNumPoints; ++k )
    {
        const real tdS = tCalc->dS( k );
        ASSERT_TRUE( std::isfinite( tdS ) ) << "dS at point " << k ;
        tArea += w( k ) * tdS ;
        tMaxDeviation = std::max( tMaxDeviation, std::abs( tdS - tdS0 ) );
    }

    const real tExpected = fem::absolute_footprint_area();
    EXPECT_NEAR( tArea, tExpected, 1.0e-12 * tExpected )
        << "the weighted increments do not integrate to the face area" ;

    EXPECT_GT( tMaxDeviation, 1.0e-3 * tdS0 )
        << "the increment is constant over the face: the fixture is blind" ;
}

//------------------------------------------------------------------------------

// Part B: background_phi against a reference assembled from the same
// calculator with the per-point increment dS( k ).
TEST( BackgroundField, KernelUsesPerPointSurfaceIncrement )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    fem::BackgroundHexStack tStack( 1.0 );
    fem::Calculator * tCalc = tStack.calculator();

    const uint n = mesh::number_of_nodes( ElementType::HEX8 );
    ASSERT_EQ( tStack.facet_element()->number_of_local_dofs(), n );

    tCalc->link( tStack.facet_element() );

    fem::TimestepMatrices tMat ;
    tMat.set_flag( fem::MatrixFlag::K );
    tMat.set_flag( fem::MatrixFlag::F );
    tMat.initialize( n );
    fem::maxwell::background_phi( tCalc, & tMat );

    const Matrix< real > tK = tMat.K();
    const Vector< real > tF = tMat.f();

    // reference: relink so no cached normal from the kernel call survives
    tCalc->link( tStack.facet_element() );

    Vector< real > tH( 3 );
    for ( uint i = 0; i < 3; ++i )
    {
        tH( i ) = tStack.field_value() * tStack.direction()( i );
    }

    const Vector< real > & w = tCalc->integration()->weights();
    Matrix< real > tKref( n, n, 0.0 );
    Vector< real > tFref( n, 0.0 );
    for ( uint k = 0; k < tCalc->num_intpoints(); ++k )
    {
        const Matrix< real > & B = tCalc->Bm( k );
        const real tWdS = w( k ) * tCalc->dS( k );
        tKref += tWdS * trans( B ) * B ;
        tFref -= tWdS * trans( B ) * tH ;
    }

    ASSERT_EQ( tK.n_rows(), n );
    ASSERT_EQ( tK.n_cols(), n );
    ASSERT_EQ( tF.length(), n );

    real tMaxKref = 0.0 ;
    real tMaxDK   = 0.0 ;
    for ( uint j = 0; j < n; ++j )
    {
        for ( uint i = 0; i < n; ++i )
        {
            ASSERT_TRUE( std::isfinite( tK( i, j ) ) ) << "K( " << i << ", " << j << " )" ;
            ASSERT_TRUE( std::isfinite( tKref( i, j ) ) ) << "Kref( " << i << ", " << j << " )" ;
            tMaxKref = std::max( tMaxKref, std::abs( tKref( i, j ) ) );
            tMaxDK   = std::max( tMaxDK, std::abs( tK( i, j ) - tKref( i, j ) ) );
        }
    }

    real tMaxFref = 0.0 ;
    real tMaxDF   = 0.0 ;
    for ( uint i = 0; i < n; ++i )
    {
        ASSERT_TRUE( std::isfinite( tF( i ) ) ) << "f( " << i << " )" ;
        ASSERT_TRUE( std::isfinite( tFref( i ) ) ) << "fref( " << i << " )" ;
        tMaxFref = std::max( tMaxFref, std::abs( tFref( i ) ) );
        tMaxDF   = std::max( tMaxDF, std::abs( tF( i ) - tFref( i ) ) );
    }

    ASSERT_GT( tMaxKref, 0.0 );
    ASSERT_GT( tMaxFref, 0.0 );

    EXPECT_LT( tMaxDK, 1.0e-12 * tMaxKref )
        << "K: relative deviation " << tMaxDK / tMaxKref ;
    EXPECT_LT( tMaxDF, 1.0e-12 * tMaxFref )
        << "f: relative deviation " << tMaxDF / tMaxFref ;
}
