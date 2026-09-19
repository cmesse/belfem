/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * A block must free the last generation of its enrichment tables when it
 * dies. A parentless block enriches TRI and TET geometry unconditionally,
 * so it builds those tables without a kernel.
 */

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "commtools.hpp"
#include "cl_FEM_Block.hpp"

using namespace belfem;

// No leak assertion of its own: the test exists to be run under Valgrind,
// where the tables of the last generation must not appear in a loss record
TEST( BlockRebuild, EnrichmentTablesAreFreed )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only" ;
    }

    fem::Block * tBlock = new fem::Block( nullptr, ElementType::TRI3 );
    tBlock->initialize_lookup_tables( 4 );
    tBlock->initialize_lookup_tables( 5 );
    EXPECT_NE( tBlock->enrichment_data( 0 ), nullptr );
    delete tBlock ;
}
