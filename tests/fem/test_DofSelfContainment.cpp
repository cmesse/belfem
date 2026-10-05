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
 * The flattening invariant of hanging dofs, as a property of the built dof
 * set: a hanging dof never lists itself among its sources, and no source
 * hangs in turn. Dof::set_sources only asserts the second half, and only in
 * an assertions-active build; at that moment the dof being set does not
 * hang yet, so a self-source would pass the assert. The walk below runs in
 * every build.
 *
 * DofData::reset deletes every dof exactly once, also in the states a refusal
 * inside the setup leaves: collect_hanging_dofs shares pointers between the
 * dof and hanging-dof containers until its end, and create_dofs leaves null
 * slots. It also releases the mesh dof containers that still hold its dofs,
 * the empty ones a sibling's connect left behind, and the counters of bases
 * that were never allocated. These tests build those states by hand; the
 * refusal path itself needs a hanging-edge fixture, because the mesh
 * flattens node source chains before the dof stage.
 */

#include <gtest/gtest.h>
#include <memory>

#include "typedefs.hpp"
#include "commtools.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_Block.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_KernelParameters.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_Dof.hpp"
#include "cl_FEM_DofMgr_DofData.hpp"
#include "en_IWGs.hpp"
#include "en_SolverEnums.hpp"

using namespace belfem;

namespace
{
    //! two unit squares that do NOT share their interface nodes; the right
    //! square's interface nodes hang on the left square's, one by a single
    //! unit source and one by two half-weight sources
    //!
    //!   4---3  6---8      node 5 ( 1, 0 ) hangs on node 2
    //!   | L |  | R |      node 6 ( 1, 1 ) hangs on nodes 3 and 2, 0.5 each
    //!   1---2  5---7      TRI3 {1,2,3} {1,3,4} {5,7,8} {5,8,6}, all CCW
    class HangingStrip
    {
        public:

            // owners in construction order: a throw in the constructor
            // destroys what was built, in reverse order
            std::unique_ptr< Mesh >                   mMesh ;
            std::unique_ptr< fem::KernelParameters >  mParams ;
            std::unique_ptr< fem::Kernel >            mKernel ;
            fem::DofManager        * mField ;

            HangingStrip()
            {
                // third argument must stay false: the MeshChecker in the
                // Kernel ctor refuses pre-existing edges
                mMesh.reset( new Mesh( 2, 0, false ) );

                mesh::ElementFactory tFactory ;

                Cell< mesh::Node * > & tNodes = mMesh->nodes();
                tNodes.push( new mesh::Node( 1, 0.0, 0.0 ) );
                tNodes.push( new mesh::Node( 2, 1.0, 0.0 ) );
                tNodes.push( new mesh::Node( 3, 1.0, 1.0 ) );
                tNodes.push( new mesh::Node( 4, 0.0, 1.0 ) );
                tNodes.push( new mesh::Node( 5, 1.0, 0.0 ) );
                tNodes.push( new mesh::Node( 6, 1.0, 1.0 ) );
                tNodes.push( new mesh::Node( 7, 2.0, 0.0 ) );
                tNodes.push( new mesh::Node( 8, 2.0, 1.0 ) );

                // the ties must exist before the Kernel is built: its ctor
                // collects the hanging bases from the mesh
                {
                    Cell< mesh::Node * > tSources = { tNodes( 1 ) };
                    Vector< real >       tWeights( 1, 1.0 );
                    tNodes( 4 )->set_sources( tSources, tWeights );
                }
                {
                    Cell< mesh::Node * > tSources = { tNodes( 2 ), tNodes( 1 ) };
                    Vector< real >       tWeights( 2, 0.5 );
                    tNodes( 5 )->set_sources( tSources, tWeights );
                }

                Cell< mesh::Element * > & tElements = mMesh->elements();
                mesh::Block * tBlock = new mesh::Block( 1, 4 );

                const uint tTriangles[ 4 ][ 3 ] = {
                        { 0, 1, 2 },
                        { 0, 2, 3 },
                        { 4, 6, 7 },
                        { 4, 7, 5 } };

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

                mMesh->blocks().push( tBlock );
                mMesh->finalize();

                mParams.reset( new fem::KernelParameters( mMesh.get() ) );
                mKernel.reset( new fem::Kernel( mParams.get() ) );

                fem::IWG * tIWG = mKernel->create_equation(
                        IwgType::Poisson,
                        ModelDimensionality::TwoD,
                        { 1 } );

                mField = mKernel->create_field( tIWG );
                mField->set_solver( gDefaultSolver );
            }

            HangingStrip( const HangingStrip & ) = delete ;
            HangingStrip( HangingStrip && ) = delete ;
            HangingStrip & operator=( const HangingStrip & ) = delete ;
            HangingStrip & operator=( HangingStrip && ) = delete ;

            //! the dof of every mesh node, hanging or not
            Cell< fem::Dof * >
            node_dofs()
            {
                Cell< fem::Dof * > tDofs ;
                for ( mesh::Node * tNode : mMesh->nodes() )
                {
                    tDofs.push( mField->dof( mField->calculate_dof_id( tNode, 0 ) ) );
                }
                return tDofs ;
            }
    };

    //! walks the dof set; aPhase names the moment for the failure message
    void
    expect_flat_and_self_free( HangingStrip & aStrip, const char * aPhase,
                               uint & aNumHanging )
    {
        aNumHanging = 0 ;

        for ( fem::Dof * tDof : aStrip.node_dofs() )
        {
            ASSERT_NE( tDof, nullptr ) << aPhase ;

            if ( ! tDof->is_hanging() ) continue ;
            ++aNumHanging ;

            for ( uint k = 0; k < tDof->number_of_sources(); ++k )
            {
                fem::Dof * tSource = tDof->source( k );

                ASSERT_NE( tSource, nullptr ) << aPhase << ": dof " << tDof->id() ;
                EXPECT_NE( tSource, tDof )
                    << aPhase << ": dof " << tDof->id() << " is its own source " << k ;
                EXPECT_FALSE( tSource->is_hanging() )
                    << aPhase << ": source " << tSource->id() << " of dof "
                    << tDof->id() << " hangs itself" ;
            }
        }
    }
}

//------------------------------------------------------------------------------

TEST( DofSelfContainment, NoHangingDofContainsItselfOrAHangingSource )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    HangingStrip tStrip ;

    // the property must hold where the sources are resolved ( inside
    // create_field ) and after the dof manager is frozen
    uint tNumHanging = 0 ;
    expect_flat_and_self_free( tStrip, "after create_field", tNumHanging );
    EXPECT_EQ( tNumHanging, 2u ) << "the fixture must produce both ties" ;

    tStrip.mField->initialize();

    expect_flat_and_self_free( tStrip, "after initialize", tNumHanging );
    EXPECT_EQ( tNumHanging, 2u );
}

TEST( DofSelfContainment, SolvedDofContainerHoldsNoHangingDof )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    HangingStrip tStrip ;
    tStrip.mField->initialize();

    // the container the solver numbers holds only independent dofs, and
    // every source of a hanging dof is one of them
    Cell< fem::Dof * > & tDofs = tStrip.mField->dofs();
    EXPECT_EQ( tDofs.size(), 6u );

    for ( fem::Dof * tDof : tDofs )
    {
        EXPECT_FALSE( tDof->is_hanging() ) << "dof " << tDof->id() ;
    }

    for ( fem::Dof * tDof : tStrip.node_dofs() )
    {
        for ( uint k = 0; k < tDof->number_of_sources(); ++k )
        {
            bool tFound = false ;
            for ( fem::Dof * tIndependent : tDofs )
            {
                tFound = tFound || tIndependent == tDof->source( k );
            }
            EXPECT_TRUE( tFound ) << "source " << k << " of dof " << tDof->id()
                                  << " is not an independent dof" ;
        }
    }
}

//------------------------------------------------------------------------------

namespace
{
    //! a dof that counts its destructions, so the oracle sees a double
    //! delete and a leak alike
    class CountingDof : public fem::Dof
    {
        uint & mDestroyed ;

    public:

        CountingDof( const id_t aID, mesh::Node * aNode, uint & aDestroyed ) :
            fem::Dof( aID, 0, aNode ),
            mDestroyed( aDestroyed )
        {
        }

        ~CountingDof() override
        {
            ++mDestroyed ;
        }
    };

    //! a dof data on the strip's field whose containers the test fills
    struct ResetCase
    {
        HangingStrip  mStrip ;
        uint          mDestroyed = 0 ;
        std::unique_ptr< fem::dofmgr::DofData > mData ;

        ResetCase() :
            mData( new fem::dofmgr::DofData( mStrip.mField, nullptr ) )
        {
        }

        fem::Dof *
        make( const id_t aID )
        {
            return new CountingDof( aID, mStrip.mMesh->node( 1 ), mDestroyed );
        }

        fem::Dof *
        make_on( const id_t aID, mesh::Node * aNode )
        {
            return new CountingDof( aID, aNode, mDestroyed );
        }
    };

    //! replaces the node's dof container by one that holds only aDof, the
    //! readable state a finished connect_dofs_to_mesh() leaves on a basis it
    //! fills; the strip's field has already filled every node
    void
    connect_alone( mesh::Node * aNode, fem::Dof * aDof )
    {
        aNode->increment_dof_counter();
        aNode->allocate_dof_container();
        aNode->insert_dof( aDof );
    }
}

TEST( DofDataReset, SharedDofIsDeletedOnce )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    ResetCase tCase ;
    fem::Dof * tShared = tCase.make( 1 );
    tCase.mData->dofs().push( tShared );
    tCase.mData->hanging_dofs().push( tShared );

    tCase.mData->reset();

    EXPECT_EQ( tCase.mDestroyed, 1u );
    EXPECT_EQ( tCase.mData->dofs().size(), 0u );
    EXPECT_EQ( tCase.mData->hanging_dofs().size(), 0u );
}

TEST( DofDataReset, DisjointDofsAreEachDeleted )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    ResetCase tCase ;
    tCase.mData->dofs().push( tCase.make( 1 ) );
    tCase.mData->dofs().push( tCase.make( 2 ) );
    tCase.mData->hanging_dofs().push( tCase.make( 3 ) );

    tCase.mData->reset();

    EXPECT_EQ( tCase.mDestroyed, 3u );
}

TEST( DofDataReset, NullSlotsAreSkipped )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    ResetCase tCase ;
    fem::Dof * tShared = tCase.make( 1 );

    tCase.mData->dofs().push( tShared );
    tCase.mData->dofs().push( nullptr );
    tCase.mData->hanging_dofs().push( nullptr );
    tCase.mData->hanging_dofs().push( tShared );
    tCase.mData->hanging_dofs().push( tCase.make( 2 ) );

    ASSERT_NO_THROW( tCase.mData->reset() );

    EXPECT_EQ( tCase.mDestroyed, 2u );
}

TEST( DofDataReset, StaleFlagsDoNotDecide )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    // the worker branch of collect_hanging_dofs flags the dofs it is
    // collecting, so a refusal there leaves flags behind
    ResetCase tCase ;
    fem::Dof * tShared  = tCase.make( 1 );
    fem::Dof * tHanging = tCase.make( 2 );
    tShared->flag();
    tHanging->flag();

    tCase.mData->dofs().push( tShared );
    tCase.mData->hanging_dofs().push( tShared );
    tCase.mData->hanging_dofs().push( tHanging );

    tCase.mData->reset();

    EXPECT_EQ( tCase.mDestroyed, 2u );
}

TEST( DofDataReset, MeshContainersAreReleased )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    // reset releases every container the setup filled, including a basis
    // reached only through the hanging list
    ResetCase tCase ;
    mesh::Node * tNode1 = tCase.mStrip.mMesh->node( 1 );
    mesh::Node * tNode2 = tCase.mStrip.mMesh->node( 2 );
    mesh::Node * tNode3 = tCase.mStrip.mMesh->node( 3 );

    fem::Dof * tShared  = tCase.make_on( 1, tNode1 );
    fem::Dof * tOwned   = tCase.make_on( 2, tNode2 );
    fem::Dof * tHanging = tCase.make_on( 3, tNode3 );

    tCase.mData->dofs().push( tShared );
    tCase.mData->dofs().push( tOwned );
    tCase.mData->hanging_dofs().push( tShared );
    tCase.mData->hanging_dofs().push( tHanging );

    connect_alone( tNode1, tShared );
    connect_alone( tNode2, tOwned );
    connect_alone( tNode3, tHanging );

    ASSERT_EQ( tNode1->number_of_dofs(), 1u );
    ASSERT_EQ( tNode2->number_of_dofs(), 1u );
    ASSERT_EQ( tNode3->number_of_dofs(), 1u );

    tCase.mData->reset();

    EXPECT_EQ( tCase.mDestroyed, 3u );
    EXPECT_EQ( tNode1->number_of_dofs(), 0u );
    EXPECT_EQ( tNode2->number_of_dofs(), 0u );
    EXPECT_EQ( tNode3->number_of_dofs(), 0u );
}

TEST( DofDataReset, AnotherManagersContainerIsKept )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    // a manager's teardown releases only the containers that still hold its
    // own dofs; a sibling that replaced them keeps its links
    ResetCase tA ;
    uint tDestroyedB = 0 ;
    std::unique_ptr< fem::dofmgr::DofData > tB(
            new fem::dofmgr::DofData( tA.mStrip.mField, nullptr ) );

    mesh::Node * tNode1 = tA.mStrip.mMesh->node( 1 );

    fem::Dof * tDofA = new CountingDof( 1, tNode1, tA.mDestroyed );
    tA.mData->dofs().push( tDofA );
    connect_alone( tNode1, tDofA );

    fem::Dof * tDofB = new CountingDof( 2, tNode1, tDestroyedB );
    tB->dofs().push( tDofB );
    connect_alone( tNode1, tDofB );

    tA.mData->reset();

    EXPECT_EQ( tA.mDestroyed, 1u );
    ASSERT_EQ( tNode1->number_of_dofs(), 1u );
    EXPECT_EQ( tNode1->dof( 0 ), static_cast< graph::Vertex * >( tDofB ) );

    tB->reset();

    EXPECT_EQ( tDestroyedB, 1u );
    EXPECT_EQ( tNode1->number_of_dofs(), 0u );
}

TEST( DofDataReset, ACountedBasisWithoutAnArrayIsCleared )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    // a basis that was counted but never allocated has nothing to free
    // and its counter is cleared
    ResetCase tCase ;
    mesh::Node * tNode = tCase.mStrip.mMesh->node( 2 );

    tNode->reset_dof_container();
    tNode->increment_dof_counter();

    tCase.mData->dofs().push( tCase.make_on( 1, tNode ) );

    ASSERT_NO_THROW( tCase.mData->reset() );

    EXPECT_EQ( tCase.mDestroyed, 1u );
    EXPECT_FALSE( tNode->has_dof_container() );
    EXPECT_EQ( tNode->number_of_dofs(), 0u );
}

TEST( DofDataReset, AnEmptyAllocationIsReleased )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    // an allocated container with no entries belongs to nobody and is
    // released
    ResetCase tCase ;
    mesh::Node * tNode = tCase.mStrip.mMesh->node( 3 );

    tNode->increment_dof_counter();
    tNode->allocate_dof_container();

    ASSERT_TRUE( tNode->has_dof_container() );
    ASSERT_EQ( tNode->number_of_dofs(), 0u );

    tCase.mData->dofs().push( tCase.make_on( 1, tNode ) );

    tCase.mData->reset();

    EXPECT_EQ( tCase.mDestroyed, 1u );
    EXPECT_FALSE( tNode->has_dof_container() );
    EXPECT_EQ( tNode->number_of_dofs(), 0u );
}

TEST( DofDataReset, AFullBasisIsReleasedBeforeItIsCountedAgain )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    // connect releases what a previous manager left on a basis before it
    // counts, so a basis at the counter's ceiling can be connected again
    ResetCase tFull ;
    mesh::Node * tNode1 = tFull.mStrip.mMesh->node( 1 );

    for ( uint k = 0; k < 255; ++k )
    {
        tFull.mData->dofs().push( tFull.make_on( k + 1, tNode1 ) );
    }
    tFull.mData->connect_dofs_to_mesh();

    ASSERT_EQ( tNode1->number_of_dofs(), 255u );

    uint tDestroyedOne = 0 ;
    std::unique_ptr< fem::dofmgr::DofData > tOne(
            new fem::dofmgr::DofData( tFull.mStrip.mField, nullptr ) );

    fem::Dof * tDof = new CountingDof( 1000, tNode1, tDestroyedOne );
    tOne->dofs().push( tDof );

    ASSERT_NO_THROW( tOne->connect_dofs_to_mesh() );

    EXPECT_EQ( tNode1->number_of_dofs(), 1u );
    EXPECT_TRUE( tNode1->holds_dof( tDof ) );

    tOne->reset();

    EXPECT_EQ( tDestroyedOne, 1u );
    EXPECT_EQ( tNode1->number_of_dofs(), 0u );

    tFull.mData->reset();

    EXPECT_EQ( tFull.mDestroyed, 255u );
}
