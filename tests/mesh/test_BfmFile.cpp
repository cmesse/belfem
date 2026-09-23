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

/*
 * The provenance stamp of a .bfm file. The writer stores the git hash of the
 * writing build in meta/git; a file from another commit, or one without the
 * stamp, still loads. The loader warns about the foreign hash; reuse is
 * gated elsewhere, by the base-mesh checksum and the settings fingerprint.
 */

#include <cstdio>
#include <memory>
#include <string>
#include <unistd.h>

#include <gtest/gtest.h>

#include "typedefs.hpp"
#include "commtools.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "cl_Mesh_BfmFile.hpp"
#include "cl_Cell.hpp"
#include "cl_Node.hpp"
#include "cl_Block.hpp"
#include "cl_Facet.hpp"
#include "cl_SideSet.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_HDF5.hpp"
#include "hdf5_tools.hpp"
#include "banner.hpp"

#ifdef BELFEM_HDF5
#include <hdf5.h>

namespace
{
    //! two TRI3 on the unit square in one block, saved to a scratch .bfm
    //! and removed on scope exit
    class SavedMesh
    {
        const std::string mPath ;

    public:

        belfem::Mesh mMesh ;

        explicit SavedMesh( const std::string & aName ) :
            mPath( ::testing::TempDir() + "belfem_test_bfm_" + aName + "_"
                   + std::to_string( ( long ) getpid() ) + ".bfm" ),
            mMesh( 2, 0, false )
        {
            belfem::mesh::ElementFactory tFactory ;

            belfem::Cell< belfem::mesh::Node * > & tNodes = mMesh.nodes();
            tNodes.push( new belfem::mesh::Node( 1, 0.0, 0.0 ) );
            tNodes.push( new belfem::mesh::Node( 2, 1.0, 0.0 ) );
            tNodes.push( new belfem::mesh::Node( 3, 1.0, 1.0 ) );
            tNodes.push( new belfem::mesh::Node( 4, 0.0, 1.0 ) );

            belfem::mesh::Block * tBlock = new belfem::mesh::Block( 1, 2 );
            const belfem::uint tTriangles[ 2 ][ 3 ] = { { 0, 1, 2 }, { 0, 2, 3 } };
            for ( belfem::uint e = 0; e < 2; ++e )
            {
                belfem::mesh::Element * tTri =
                        tFactory.create_element( belfem::ElementType::TRI3, e + 1 );
                for ( belfem::uint k = 0; k < 3; ++k )
                {
                    tTri->insert_node( tNodes( tTriangles[ e ][ k ] ), k );
                }
                tTri->set_block_id( 1 );
                mMesh.elements().push( tTri );
                tBlock->insert_element( tTri );
            }
            mMesh.blocks().push( tBlock );

            // one boundary facet: the writer throws on a mesh without any
            belfem::mesh::Element * tLine =
                    tFactory.create_element( belfem::ElementType::LINE2, 100 );
            tLine->insert_node( tNodes( 0 ), 0 );
            tLine->insert_node( tNodes( 1 ), 1 );
            belfem::mesh::Facet * tFacet = new belfem::mesh::Facet( tLine );
            tFacet->set_master( mMesh.elements()( 0 ), 0 );
            tFacet->set_sideset_id( 10 );

            belfem::mesh::SideSet * tSideSet = new belfem::mesh::SideSet( 10, 1 );
            tSideSet->insert_facet( tFacet );
            mMesh.sidesets().push( tSideSet );

            mMesh.finalize();

            mMesh.save( mPath );
        }

        ~SavedMesh()
        {
            std::remove( mPath.c_str() );
        }

        const std::string &
        path() const { return mPath ; }
    };

    //! loads the file; the file object hands the mesh over, so the caller
    //! owns what comes back
    std::unique_ptr< belfem::Mesh >
    reload( const std::string & aPath )
    {
        belfem::mesh::BfmFile tFile( aPath );
        tFile.load();
        return std::unique_ptr< belfem::Mesh >( tFile.get() );
    }

    //! replaces meta/git by aStamp, or removes it when aStamp is empty. The
    //! handle is closed before the function returns, so a reload can open
    //! the file again
    void
    rewrite_git_stamp( const std::string & aPath, const std::string & aStamp )
    {
        belfem::HDF5 tFile( aPath, belfem::FileMode::OPEN_RDWR );
        tFile.select_group( "meta" );
        if ( belfem::hdf5::dataset_exists( tFile.active_group(), "git" ) )
        {
            H5Ldelete( tFile.active_group(), "git", H5P_DEFAULT );
        }
        if ( aStamp.size() > 0 )
        {
            tFile.save_data( "git", aStamp );
        }
        tFile.close_active_group();
        tFile.close();
    }
}

TEST( BfmFile, MeshRoundTrips )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    SavedMesh tSaved( "roundtrip" );
    std::unique_ptr< belfem::Mesh > tLoaded = reload( tSaved.path() );

    ASSERT_NE( tLoaded, nullptr );
    ASSERT_EQ( tLoaded->number_of_nodes(), tSaved.mMesh.number_of_nodes() );
    ASSERT_EQ( tLoaded->number_of_elements(), tSaved.mMesh.number_of_elements() );
    EXPECT_EQ( tLoaded->number_of_dimensions(), 2u );

    // the checksum is stored and read back, not recomputed, so the geometry
    // and the connectivity are compared directly
    for ( belfem::index_t k = 0; k < tSaved.mMesh.number_of_nodes(); ++k )
    {
        const belfem::mesh::Node * tA = tSaved.mMesh.nodes()( k );
        const belfem::mesh::Node * tB = tLoaded->node( tA->id() );
        ASSERT_NE( tB, nullptr ) << "node " << tA->id() ;
        EXPECT_EQ( tB->x(), tA->x() ) << "node " << tA->id() ;
        EXPECT_EQ( tB->y(), tA->y() ) << "node " << tA->id() ;
    }
    for ( belfem::index_t e = 0; e < tSaved.mMesh.number_of_elements(); ++e )
    {
        const belfem::mesh::Element * tA = tSaved.mMesh.elements()( e );
        const belfem::mesh::Element * tB = tLoaded->element( tA->id() );
        ASSERT_NE( tB, nullptr ) << "element " << tA->id() ;
        ASSERT_EQ( tB->number_of_nodes(), tA->number_of_nodes() );
        for ( belfem::uint k = 0; k < tA->number_of_nodes(); ++k )
        {
            EXPECT_EQ( tB->node( k )->id(), tA->node( k )->id() )
                << "element " << tA->id() << ", node " << k ;
        }
    }
    ASSERT_EQ( tLoaded->number_of_sidesets(), 1u );
    EXPECT_EQ( tLoaded->sideset( 10 )->number_of_facets(), 1u );
}

TEST( BfmFile, MissingGitStampLoads )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    // a file written before the stamp existed, or by a build not made from git
    SavedMesh tSaved( "nostamp" );
    rewrite_git_stamp( tSaved.path(), "" );

    std::unique_ptr< belfem::Mesh > tLoaded ;
    ::testing::internal::CaptureStdout();
    ASSERT_NO_THROW( tLoaded = reload( tSaved.path() ) );
    const std::string tLog = ::testing::internal::GetCapturedStdout();
    EXPECT_EQ( tLoaded->number_of_nodes(), 4u );

    // nothing to compare, nothing to flag
    EXPECT_EQ( tLog.find( "written by git" ), std::string::npos ) << tLog ;
}

TEST( BfmFile, ForeignGitHashStillLoads )
{
    if ( belfem::comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    // a file from another commit is flagged, never refused.
    // Mutation: turn the warning into a refusal, or drop it, and this goes red
    SavedMesh tSaved( "foreign" );
    rewrite_git_stamp( tSaved.path(), "0123456789abcdef0123456789abcdef01234567" );

    std::unique_ptr< belfem::Mesh > tLoaded ;
    ::testing::internal::CaptureStdout();
    ASSERT_NO_THROW( tLoaded = reload( tSaved.path() ) );
    const std::string tLog = ::testing::internal::GetCapturedStdout();
    EXPECT_EQ( tLoaded->number_of_nodes(), 4u );

    // flagged: a build made from git compares and names both stamps
    if ( belfem::is_built_from_git() )
    {
        EXPECT_NE( tLog.find( "written by git 0123456789abcdef0123456789abcdef01234567" ),
                   std::string::npos ) << tLog ;
        EXPECT_NE( tLog.find( "but this build is git " + belfem::git_commit_hash() ),
                   std::string::npos ) << tLog ;
    }
}

#endif // BELFEM_HDF5
