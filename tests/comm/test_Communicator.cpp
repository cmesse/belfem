/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Tier 1 tests for the Communicator class accessors and
 * CommunicationObject lifecycle. Runs without MPI.
 * See: tests_06_comm.md §2
 *
 * gComm is initialized in test_main.cpp via Communicator(argc, argv).
 */

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"

extern belfem::Communicator gComm;

// =============================================================================
// §2.1  Non-MPI Defaults  [semantic]
// =============================================================================

#ifndef BELFEM_MPI

TEST( Communicator, DefaultRankIsZero )
{
    EXPECT_EQ( gComm.rank(), 0 );
}

TEST( Communicator, DefaultSizeIsOne )
{
    EXPECT_EQ( gComm.size(), 1 );
}

TEST( Communicator, DefaultNodeSizeIsOne )
{
    EXPECT_EQ( gComm.node_size(), 1 );
}

TEST( Communicator, DefaultMaxTagIsZero )
{
    EXPECT_EQ( gComm.max_tag(), 0 );
}

#endif // BELFEM_MPI

// =============================================================================
// §2.2  Argument Handling  [semantic]
// =============================================================================

TEST( Communicator, ExecPathCaptured )
{
    // after init, exec_path() should be non-empty (argv[0] was captured)
    EXPECT_FALSE( gComm.exec_path().empty() );
}

TEST( Communicator, WorkdirNonEmpty )
{
    EXPECT_FALSE( gComm.workdir().empty() );
}

TEST( Communicator, SetArgumentsOverwrites )
{
    gComm.set_arguments( "custom_test_args" );
    EXPECT_EQ( gComm.argument_string(), "custom_test_args" );
}

// =============================================================================
// §2.3  Random Engine  [semantic]
// =============================================================================

TEST( Communicator, RandomEngineUsable )
{
    // random() returns a reference to std::mt19937; calling it should not throw
    auto tVal = gComm.random()();
    ( void ) tVal;   // just verify it doesn't crash
}

// =============================================================================
// §2.4  CommunicationObject Lifecycle  [semantic]
// =============================================================================

// Minimal concrete subclass for testing
namespace
{
    class TestCommObject : public belfem::CommunicationObject
    {
        int mFreeCount = 0;
    public:
        TestCommObject() : belfem::CommunicationObject() {}

        void free() override
        {
            ++mFreeCount;
            belfem::CommunicationObject::free();
        }

        int free_count() const { return mFreeCount; }
    };
}

TEST( CommunicationObject, ObjectRegistersOnConstruction )
{
    belfem::index_t tBefore = gComm.objects().size();
    TestCommObject tObj;
    EXPECT_EQ( gComm.objects().size(), tBefore + 1 );
    tObj.free();   // prevent dangling pointer in gComm.objects()
}

TEST( CommunicationObject, ObjectIndexMatchesRegistration )
{
    TestCommObject tObj;
    belfem::index_t tIdx = tObj.index();
    EXPECT_EQ( gComm.objects()( tIdx ), &tObj );
    tObj.free();
}

TEST( CommunicationObject, ObjectFreeNullsSlot )
{
    TestCommObject tObj;
    belfem::index_t tIdx = tObj.index();

    tObj.free();

    EXPECT_EQ( gComm.objects()( tIdx ), nullptr );
    EXPECT_EQ( tObj.free_count(), 1 );
}

TEST( CommunicationObject, MultipleObjectsRegistered )
{
    belfem::index_t tBefore = gComm.objects().size();

    TestCommObject tObj1;
    TestCommObject tObj2;
    TestCommObject tObj3;

    EXPECT_EQ( gComm.objects().size(), tBefore + 3 );

    // each has a unique index
    EXPECT_NE( tObj1.index(), tObj2.index() );
    EXPECT_NE( tObj2.index(), tObj3.index() );
    EXPECT_NE( tObj1.index(), tObj3.index() );

    tObj1.free();
    tObj2.free();
    tObj3.free();
}
