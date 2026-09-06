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
#include "commtools.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_Facet.hpp"
#include "cl_Block.hpp"
#include "cl_SideSet.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_KernelParameters.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_SideSet.hpp"
#include "cl_FEM_Dof.hpp"
#include "en_IWGs.hpp"
#include "en_SolverEnums.hpp"

using namespace belfem;

namespace
{
    // the seeding contract under test:
    // DofManager::initialize( true ) must seed the FREE dofs from the field
    // and must NOT write the fixed dofs' values back into the field — that
    // back-write is what clobbers restored fields at Dirichlet nodes after
    // a warm restart ( Controller::load_memdump ).

    constexpr real tDirichletValue = 300.0 ;

    inline real
    sentinel( const index_t aIndex )
    {
        return 100.0 + static_cast< real >( aIndex ) ;
    }

    //! two-square strip with a Dirichlet sideset on the left edge, mimicking
    //! the restore ordering: dofs fixed and field populated BEFORE the
    //! first initialize() of the manager. The block meshes as four TRI3
    //! or as two QUAD4 ( the latter gates the Calculator's plain-QUAD4
    //! dispatch cases )
    class SeedingProblem
    {
        public:

            Mesh                   * mMesh ;
            fem::KernelParameters  * mParams ;
            fem::Kernel            * mKernel ;
            fem::DofManager        * mField ;

//------------------------------------------------------------------------------

            SeedingProblem( const ElementType aElementType = ElementType::TRI3 )
            {
                // third argument must stay false: connectivity computation
                // would create edges in finalize() and the MeshChecker in
                // the Kernel ctor refuses pre-existing edges
                mMesh = new Mesh( 2, 0, false );

                mesh::ElementFactory tFactory ;

                // node grid
                //   4---5---6
                //   | 1 | 2 |
                //   1---2---3
                Cell< mesh::Node * > & tNodes = mMesh->nodes();
                tNodes.push( new mesh::Node( 1, 0.0, 0.0 ) );
                tNodes.push( new mesh::Node( 2, 1.0, 0.0 ) );
                tNodes.push( new mesh::Node( 3, 2.0, 0.0 ) );
                tNodes.push( new mesh::Node( 4, 0.0, 1.0 ) );
                tNodes.push( new mesh::Node( 5, 1.0, 1.0 ) );
                tNodes.push( new mesh::Node( 6, 2.0, 1.0 ) );

                Cell< mesh::Element * > & tElements = mMesh->elements();

                mesh::Block * tBlock = nullptr ;

                // the element carrying the left edge, and the local index
                // of that edge on it
                mesh::Element * tDirichletMaster = nullptr ;
                uint            tDirichletEdge   = 0 ;

                if ( aElementType == ElementType::QUAD4 )
                {
                    tBlock = new mesh::Block( 1, 2 );

                    // both quads counterclockwise
                    const uint tQuads[ 2 ][ 4 ] = {
                            { 0, 1, 4, 3 },     // ids 1, 2, 5, 4
                            { 1, 2, 5, 4 } };   // ids 2, 3, 6, 5

                    for ( uint e = 0; e < 2; ++e )
                    {
                        mesh::Element * tQuad =
                                tFactory.create_element( ElementType::QUAD4, e + 1 );
                        for ( uint k = 0; k < 4; ++k )
                        {
                            tQuad->insert_node( tNodes( tQuads[ e ][ k ] ), k );
                        }
                        tQuad->set_block_id( 1 );
                        tElements.push( tQuad );
                        tBlock->insert_element( tQuad );
                    }

                    // the left edge ( ids 4->1 ) is local edge 3
                    // ( nodes 3->0 ) of quad 1
                    tDirichletMaster = tElements( 0 );
                    tDirichletEdge   = 3 ;
                }
                else
                {
                    // the fixture knows exactly two topologies; anything
                    // else must not silently fall through to triangles
                    EXPECT_EQ( aElementType, ElementType::TRI3 );

                    tBlock = new mesh::Block( 1, 4 );

                    // each square split into two counterclockwise triangles;
                    // node triples are indices into tNodes
                    const uint tTriangles[ 4 ][ 3 ] = {
                            { 0, 1, 4 },     // ids 1, 2, 5
                            { 0, 4, 3 },     // ids 1, 5, 4 — carries the left edge
                            { 1, 2, 5 },     // ids 2, 3, 6
                            { 1, 5, 4 } };   // ids 2, 6, 5

                    for ( uint e = 0; e < 4; ++e )
                    {
                        mesh::Element * tTri =
                                tFactory.create_element( ElementType::TRI3, e + 1 );
                        for ( uint k = 0; k < 3; ++k )
                        {
                            tTri->insert_node( tNodes( tTriangles[ e ][ k ] ), k );
                        }
                        tTri->set_block_id( 1 );
                        tElements.push( tTri );
                        tBlock->insert_element( tTri );
                    }

                    // the left edge ( ids 4->1 ) is local edge 2
                    // ( nodes 2->0 ) of triangle 2
                    tDirichletMaster = tElements( 1 );
                    tDirichletEdge   = 2 ;
                }

                mMesh->blocks().push( tBlock );

                // set_master relinks the facet nodes to the canonical
                // face order
                mesh::Element * tLine =
                        tFactory.create_element( ElementType::LINE2, 101 );
                tLine->insert_node( tNodes( 3 ), 0 );
                tLine->insert_node( tNodes( 0 ), 1 );

                mesh::Facet * tFacet = new mesh::Facet( tLine );
                tFacet->set_master( tDirichletMaster, tDirichletEdge );

                mesh::SideSet * tSideSet = new mesh::SideSet( 10, 1 );
                tSideSet->insert_facet( tFacet );
                mMesh->sidesets().push( tSideSet );

                mMesh->finalize();

                mParams = new fem::KernelParameters( mMesh );
                mKernel = new fem::Kernel( mParams );

                fem::IWG * tIWG = mKernel->create_equation(
                        IwgType::Poisson,
                        ModelDimensionality::TwoD,
                        { 1 },
                        { 10 } );

                mField = mKernel->create_field( tIWG );
                mField->set_solver( gDefaultSolver );

                // fix the boundary dofs, as the factories do before a restore
                mField->sideset( 10 )->impose_dirichlet( tDirichletValue );

                // stand-in for load_fields: size the dof field, then plant
                // the "restored" state. create_fields() inside initialize()
                // resizes only on length mismatch, so the values survive
                mField->create_fields( mField->iwg() );

                Vector< real > & tPhi = mMesh->field_data( "phi" );
                for ( index_t i = 0; i < tPhi.length(); ++i )
                {
                    tPhi( i ) = sentinel( i );
                }
            }

//------------------------------------------------------------------------------

            ~SeedingProblem()
            {
                delete mKernel ;
                delete mParams ;
                delete mMesh ;
            }

//------------------------------------------------------------------------------

            fem::Dof *
            node_dof( mesh::Node * aNode )
            {
                return mField->dof( mField->calculate_dof_id( aNode, 0 ) );
            }
    };
}

//------------------------------------------------------------------------------

TEST( DofSeeding, SeedModePreservesRestoredFields )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    SeedingProblem tProblem ;

    tProblem.mField->initialize( true );

    Vector< real > & tPhi = tProblem.mMesh->field_data( "phi" );

    uint tNumFixed = 0 ;

    for ( mesh::Node * tNode : tProblem.mMesh->nodes() )
    {
        const index_t i = tNode->index() ;

        // the restored state must survive everywhere, Dirichlet nodes included
        EXPECT_NEAR( tPhi( i ), sentinel( i ), 1e-12 )
            << "field clobbered at node id " << tNode->id() ;

        fem::Dof * tDof = tProblem.node_dof( tNode );

        if ( tDof->is_fixed() )
        {
            // fixed dofs keep their imposed value: seeding is free-dofs-only
            EXPECT_NEAR( tDof->value(), tDirichletValue, 1e-12 );
            ++tNumFixed ;
        }
        else
        {
            // free dofs were seeded from the field
            EXPECT_NEAR( tDof->value(), sentinel( i ), 1e-12 );
        }
    }

    // the left edge holds exactly two Dirichlet nodes
    EXPECT_EQ( tNumFixed, 2u );
}

//------------------------------------------------------------------------------

// the same seeding contract on a plain-QUAD4 block: gates the Calculator's
// QUAD4 dispatch cases ( nedelec data + 2D dV ), which rejected quads under
// every physics before 2026-08-27
TEST( DofSeeding, SeedModeOnQuad4Block )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    SeedingProblem tProblem( ElementType::QUAD4 );

    tProblem.mField->initialize( true );

    Vector< real > & tPhi = tProblem.mMesh->field_data( "phi" );

    uint tNumFixed = 0 ;

    for ( mesh::Node * tNode : tProblem.mMesh->nodes() )
    {
        const index_t i = tNode->index() ;

        EXPECT_NEAR( tPhi( i ), sentinel( i ), 1e-12 )
            << "field clobbered at node id " << tNode->id() ;

        fem::Dof * tDof = tProblem.node_dof( tNode );

        if ( tDof->is_fixed() )
        {
            EXPECT_NEAR( tDof->value(), tDirichletValue, 1e-12 );
            ++tNumFixed ;
        }
        else
        {
            EXPECT_NEAR( tDof->value(), sentinel( i ), 1e-12 );
        }
    }

    EXPECT_EQ( tNumFixed, 2u );
}

//------------------------------------------------------------------------------

// negative control: full-mode initialize() must overwrite the Dirichlet
// nodes' field slots with the imposed dof value. If this stops firing, the
// suite has lost its ability to see the H-C defect and the seed-mode test
// above proves nothing
TEST( DofSeeding, FullModeWritesFixedDofsIntoField )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    SeedingProblem tProblem ;

    tProblem.mField->initialize();

    Vector< real > & tPhi = tProblem.mMesh->field_data( "phi" );

    uint tNumFixed = 0 ;

    for ( mesh::Node * tNode : tProblem.mMesh->nodes() )
    {
        const index_t i = tNode->index() ;

        if ( tProblem.node_dof( tNode )->is_fixed() )
        {
            // the back-write replaced the sentinel with the imposed value
            EXPECT_NEAR( tPhi( i ), tDirichletValue, 1e-12 );
            ++tNumFixed ;
        }
        else
        {
            EXPECT_NEAR( tPhi( i ), sentinel( i ), 1e-12 );
        }
    }

    EXPECT_EQ( tNumFixed, 2u );
}
