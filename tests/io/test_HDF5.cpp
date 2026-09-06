#include <gtest/gtest.h>

#include <cstdio>
#include <stdexcept>
#include <string>

#include "cl_HDF5.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"
#include "filetools.hpp"

//-----------------------------------------------------------------------------
// BELFEM I/O Module Tests -- HDF5 Focus
//
// This file tests the high-level HDF5 wrapper contract exposed by belfem::HDF5,
// not the raw HDF5 C API and not the internal hdf5:: helper templates.
//
// Scope covered here:
//   - file modes and file lifecycle
//   - group creation / selection / tree navigation
//   - round-trip save/load for supported public overloads
//   - representative error paths guarded by BELFEM_ERROR
//
// Deferred on purpose:
//   - parallel HDF5 / OPEN_RDONLY_PARALLEL / make_path_parallel semantics
//   - Ascii / CsvFile / InputFile / XML
//   - mesh-level HDF5 reader/writer integration
//
// NOTE:
// The uploaded implementation exposes raw-array load overloads for index_t*
// and real*. The draft therefore follows the implementation, even though an
// earlier plan variant still mentioned uint* for the integer raw-array path.
//-----------------------------------------------------------------------------

#ifdef BELFEM_HDF5

namespace
{
    const belfem::real tEps = 1.0e-12;

    void
    remove_file_if_exists( const std::string & aPath )
    {
        if ( belfem::file_exists( aPath ) )
        {
            std::remove( aPath.c_str() );
        }
    }

    template< typename T >
    void
    expect_vector_eq( const belfem::Vector< T > & aA,
                      const belfem::Vector< T > & aB )
    {
        ASSERT_EQ( aA.length(), aB.length() );
        for ( belfem::index_t k = 0; k < aA.length(); ++k )
        {
            EXPECT_EQ( aA( k ), aB( k ) );
        }
    }

    void
    expect_vector_near( const belfem::Vector< belfem::real > & aA,
                        const belfem::Vector< belfem::real > & aB,
                        const belfem::real aTol = tEps )
    {
        ASSERT_EQ( aA.length(), aB.length() );
        for ( belfem::index_t k = 0; k < aA.length(); ++k )
        {
            EXPECT_NEAR( aA( k ), aB( k ), aTol );
        }
    }

    template< typename T >
    void
    expect_matrix_eq( const belfem::Matrix< T > & aA,
                      const belfem::Matrix< T > & aB )
    {
        ASSERT_EQ( aA.n_rows(), aB.n_rows() );
        ASSERT_EQ( aA.n_cols(), aB.n_cols() );
        for ( belfem::index_t i = 0; i < aA.n_rows(); ++i )
        {
            for ( belfem::index_t j = 0; j < aA.n_cols(); ++j )
            {
                EXPECT_EQ( aA( i, j ), aB( i, j ) );
            }
        }
    }

    void
    expect_matrix_near( const belfem::Matrix< belfem::real > & aA,
                        const belfem::Matrix< belfem::real > & aB,
                        const belfem::real aTol = tEps )
    {
        ASSERT_EQ( aA.n_rows(), aB.n_rows() );
        ASSERT_EQ( aA.n_cols(), aB.n_cols() );
        for ( belfem::index_t i = 0; i < aA.n_rows(); ++i )
        {
            for ( belfem::index_t j = 0; j < aA.n_cols(); ++j )
            {
                EXPECT_NEAR( aA( i, j ), aB( i, j ), aTol );
            }
        }
    }

    void
    expect_cell_string_eq( const belfem::Cell< belfem::string > & aA,
                           const belfem::Cell< belfem::string > & aB )
    {
        ASSERT_EQ( aA.size(), aB.size() );
        for ( belfem::index_t k = 0; k < aA.size(); ++k )
        {
            EXPECT_EQ( aA( k ), aB( k ) );
        }
    }
}

class HDF5Test : public ::testing::Test
{
    protected:
        std::string mPath;

        std::string
        make_path( const std::string & aSuffix = "data" ) const
        {
            const ::testing::TestInfo * tInfo =
                    ::testing::UnitTest::GetInstance()->current_test_info();

            return std::string( "test_hdf5_" )
                 + tInfo->test_suite_name()
                 + "_"
                 + tInfo->name()
                 + "_"
                 + aSuffix
                 + ".hdf5";
        }

        void
        reset_path( const std::string & aSuffix = "data" )
        {
            mPath = this->make_path( aSuffix );
            remove_file_if_exists( mPath );
        }

        void
        TearDown() override
        {
            if ( ! mPath.empty() )
            {
                remove_file_if_exists( mPath );
            }
        }
};

//-----------------------------------------------------------------------------
// 1. File modes and file lifecycle
//-----------------------------------------------------------------------------

TEST_F( HDF5Test, FileModeNewCreatesFile )
{
    this->reset_path();

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
    }

    EXPECT_TRUE( belfem::file_exists( mPath ) );
}

TEST_F( HDF5Test, FileModeNewTruncatesExisting )
{
    this->reset_path();

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        const belfem::sint tValue = -42;
        tFile.save_data( "old_value", tValue );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        const belfem::uint tValue = 17;
        tFile.save_data( "new_value", tValue );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );

        belfem::uint tNewValue = 0;
        tFile.load_data( "new_value", tNewValue );
        EXPECT_EQ( tNewValue, 17u );

        belfem::sint tOldValue = 0;
        EXPECT_THROW( tFile.load_data( "old_value", tOldValue ), std::runtime_error );
    }
}

TEST_F( HDF5Test, FileModeOpenRdonly )
{
    this->reset_path();

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        const belfem::real tValue = 3.25;
        tFile.save_data( "x", tValue );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        belfem::real tValue = 0.0;
        tFile.load_data( "x", tValue );
        EXPECT_NEAR( tValue, 3.25, tEps );
    }
}

TEST_F( HDF5Test, FileModeOpenRdwr )
{
    this->reset_path();

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        const belfem::sint tValue = 7;
        tFile.save_data( "a", tValue );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDWR );
        const belfem::real tValue = 2.5;
        tFile.save_data( "b", tValue );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        belfem::sint tA = 0;
        belfem::real tB = 0.0;
        tFile.load_data( "a", tA );
        tFile.load_data( "b", tB );
        EXPECT_EQ( tA, 7 );
        EXPECT_NEAR( tB, 2.5, tEps );
    }
}

//-----------------------------------------------------------------------------
// 2. Scalar round-trips
//-----------------------------------------------------------------------------

TEST_F( HDF5Test, ScalarSintRoundTrip )
{
    this->reset_path();

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        const belfem::sint tExpected = -42;
        tFile.save_data( "value", tExpected );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        belfem::sint tActual = 0;
        tFile.load_data( "value", tActual );
        EXPECT_EQ( tActual, -42 );
    }
}

TEST_F( HDF5Test, ScalarUintRoundTrip )
{
    this->reset_path();

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        const belfem::uint tExpected = 12345u;
        tFile.save_data( "value", tExpected );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        belfem::uint tActual = 0;
        tFile.load_data( "value", tActual );
        EXPECT_EQ( tActual, 12345u );
    }
}

TEST_F( HDF5Test, ScalarLuintRoundTrip )
{
    this->reset_path();

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        const belfem::luint tExpected = 987654321ul;
        tFile.save_data( "value", tExpected );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        belfem::luint tActual = 0;
        tFile.load_data( "value", tActual );
        EXPECT_EQ( tActual, 987654321ul );
    }
}

TEST_F( HDF5Test, ScalarRealRoundTrip )
{
    this->reset_path();

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        const belfem::real tExpected = 3.14159265358979;
        tFile.save_data( "value", tExpected );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        belfem::real tActual = 0.0;
        tFile.load_data( "value", tActual );
        EXPECT_NEAR( tActual, 3.14159265358979, tEps );
    }
}

TEST_F( HDF5Test, ScalarBoolTrueRoundTrip )
{
    this->reset_path();

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.save_data( "value", true );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        bool tActual = false;
        tFile.load_data( "value", tActual );
        EXPECT_TRUE( tActual );
    }
}

TEST_F( HDF5Test, ScalarBoolFalseRoundTrip )
{
    this->reset_path();

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.save_data( "value", false );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        bool tActual = true;
        tFile.load_data( "value", tActual );
        EXPECT_FALSE( tActual );
    }
}

//-----------------------------------------------------------------------------
// 3. String round-trips
//-----------------------------------------------------------------------------

TEST_F( HDF5Test, StringRoundTrip )
{
    this->reset_path();

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        const belfem::string tExpected = "Hello BELFEM";
        tFile.save_data( "text", tExpected );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        belfem::string tActual;
        tFile.load_data( "text", tActual );
        EXPECT_EQ( tActual, "Hello BELFEM" );
    }
}

TEST_F( HDF5Test, EmptyStringThrows )
{
    this->reset_path();

    belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
    const belfem::string tEmpty = "";
    EXPECT_THROW( tFile.save_data( "text", tEmpty ), std::runtime_error );
}

TEST_F( HDF5Test, StringFromCharPointer )
{
    this->reset_path();

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.save_data( "text", "literal" );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        belfem::string tActual;
        tFile.load_data( "text", tActual );
        EXPECT_EQ( tActual, "literal" );
    }
}

//-----------------------------------------------------------------------------
// 4. Vector round-trips
//-----------------------------------------------------------------------------

TEST_F( HDF5Test, VectorSintRoundTrip )
{
    this->reset_path();

    const belfem::Vector< belfem::sint > tExpected = { -1, 0, 1, 42 };

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.save_data( "vec", tExpected );
    }

    {
        belfem::Vector< belfem::sint > tActual;
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        tFile.load_data( "vec", tActual );
        expect_vector_eq( tActual, tExpected );
    }
}

TEST_F( HDF5Test, VectorUintRoundTrip )
{
    this->reset_path();

    const belfem::Vector< belfem::uint > tExpected = { 1u, 2u, 3u, 4u, 5u };

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.save_data( "vec", tExpected );
    }

    {
        belfem::Vector< belfem::uint > tActual;
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        tFile.load_data( "vec", tActual );
        expect_vector_eq( tActual, tExpected );
    }
}

TEST_F( HDF5Test, VectorLuintRoundTrip )
{
    this->reset_path();

    const belfem::Vector< belfem::luint > tExpected = { 10ul, 100ul, 1000ul };

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.save_data( "vec", tExpected );
    }

    {
        belfem::Vector< belfem::luint > tActual;
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        tFile.load_data( "vec", tActual );
        expect_vector_eq( tActual, tExpected );
    }
}

TEST_F( HDF5Test, VectorRealRoundTrip )
{
    this->reset_path();

    const belfem::Vector< belfem::real > tExpected = { 0.5, -1.25, 3.75, 8.0 };

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.save_data( "vec", tExpected );
    }

    {
        belfem::Vector< belfem::real > tActual;
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        tFile.load_data( "vec", tActual );
        expect_vector_near( tActual, tExpected );
    }
}

TEST_F( HDF5Test, VectorEmptyRoundTrip )
{
    this->reset_path();

    const belfem::Vector< belfem::real > tExpected;

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.save_data( "vec", tExpected );
    }

    {
        belfem::Vector< belfem::real > tActual;
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        tFile.load_data( "vec", tActual );
        EXPECT_EQ( tActual.length(), 0u );
    }
}

//-----------------------------------------------------------------------------
// 5. Matrix round-trips
//-----------------------------------------------------------------------------

TEST_F( HDF5Test, MatrixSintRoundTrip )
{
    this->reset_path();

    belfem::Matrix< belfem::sint > tExpected( 3, 4 );
    belfem::sint tValue = -6;
    for ( belfem::index_t i = 0; i < tExpected.n_rows(); ++i )
    {
        for ( belfem::index_t j = 0; j < tExpected.n_cols(); ++j )
        {
            tExpected( i, j ) = tValue++;
        }
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.save_data( "mat", tExpected );
    }

    {
        belfem::Matrix< belfem::sint > tActual;
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        tFile.load_data( "mat", tActual );
        expect_matrix_eq( tActual, tExpected );
    }
}

TEST_F( HDF5Test, MatrixUintRoundTrip )
{
    this->reset_path();

    belfem::Matrix< belfem::uint > tExpected( 2, 3 );
    tExpected( 0, 0 ) = 1u;
    tExpected( 0, 1 ) = 2u;
    tExpected( 0, 2 ) = 3u;
    tExpected( 1, 0 ) = 4u;
    tExpected( 1, 1 ) = 5u;
    tExpected( 1, 2 ) = 6u;

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.save_data( "mat", tExpected );
    }

    {
        belfem::Matrix< belfem::uint > tActual;
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        tFile.load_data( "mat", tActual );
        expect_matrix_eq( tActual, tExpected );
    }
}

TEST_F( HDF5Test, MatrixRealRoundTrip )
{
    this->reset_path();

    belfem::Matrix< belfem::real > tExpected( 4, 5 );
    for ( belfem::index_t i = 0; i < tExpected.n_rows(); ++i )
    {
        for ( belfem::index_t j = 0; j < tExpected.n_cols(); ++j )
        {
            tExpected( i, j ) = 0.25 * static_cast< belfem::real >( 10 * i + j ) - 1.0;
        }
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.save_data( "mat", tExpected );
    }

    {
        belfem::Matrix< belfem::real > tActual;
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        tFile.load_data( "mat", tActual );
        expect_matrix_near( tActual, tExpected );
    }
}

TEST_F( HDF5Test, MatrixRealDimensionsPreserved )
{
    this->reset_path();

    belfem::Matrix< belfem::real > tExpected( 4, 5, 0.0 );
    tExpected( 2, 3 ) = 9.0;

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.save_data( "mat", tExpected );
    }

    {
        belfem::Matrix< belfem::real > tActual;
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        tFile.load_data( "mat", tActual );
        EXPECT_EQ( tActual.n_rows(), 4u );
        EXPECT_EQ( tActual.n_cols(), 5u );
    }
}

TEST_F( HDF5Test, MatrixSingleElement )
{
    this->reset_path();

    belfem::Matrix< belfem::real > tExpected( 1, 1 );
    tExpected( 0, 0 ) = 42.5;

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.save_data( "mat", tExpected );
    }

    {
        belfem::Matrix< belfem::real > tActual;
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        tFile.load_data( "mat", tActual );
        EXPECT_EQ( tActual.n_rows(), 1u );
        EXPECT_EQ( tActual.n_cols(), 1u );
        EXPECT_NEAR( tActual( 0, 0 ), 42.5, tEps );
    }
}

//-----------------------------------------------------------------------------
// 6. Cell<string> and raw-array round-trips
//-----------------------------------------------------------------------------

TEST_F( HDF5Test, CellStringRoundTrip )
{
    this->reset_path();

    const belfem::Cell< belfem::string > tExpected = { "alpha", "beta", "gamma" };

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.save_data( "cell", tExpected );
    }

    {
        belfem::Cell< belfem::string > tActual;
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        tFile.load_data( "cell", tActual );
        expect_cell_string_eq( tActual, tExpected );
    }
}

TEST_F( HDF5Test, CellStringEmpty )
{
    this->reset_path();

    const belfem::Cell< belfem::string > tExpected;

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.save_data( "cell", tExpected );
    }

    {
        belfem::Cell< belfem::string > tActual;
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        tFile.load_data( "cell", tActual );
        EXPECT_EQ( tActual.size(), 0u );
    }
}

TEST_F( HDF5Test, CellStringWithWhitespace )
{
    this->reset_path();

    const belfem::Cell< belfem::string > tExpected = {
        "hello world",
        "foo bar",
        "  padded text  "
    };

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.save_data( "cell", tExpected );
    }

    {
        belfem::Cell< belfem::string > tActual;
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        tFile.load_data( "cell", tActual );
        expect_cell_string_eq( tActual, tExpected );
    }
}

TEST_F( HDF5Test, RawArrayIndexRoundTrip )
{
    this->reset_path();

    const belfem::Vector< belfem::index_t > tExpected = { 1, 2, 3, 5, 8 };

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.save_data( "vec", tExpected );
    }

    {
        belfem::index_t tBuffer[ 5 ] = { 0, 0, 0, 0, 0 };
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        tFile.load_data( "vec", tBuffer, 5 );

        for ( belfem::index_t k = 0; k < 5; ++k )
        {
            EXPECT_EQ( tBuffer[ k ], tExpected( k ) );
        }
    }
}

TEST_F( HDF5Test, RawArrayRealRoundTrip )
{
    this->reset_path();

    const belfem::Vector< belfem::real > tExpected = { 0.25, 0.5, 1.0, 2.0 };

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.save_data( "vec", tExpected );
    }

    {
        belfem::real tBuffer[ 4 ] = { 0.0, 0.0, 0.0, 0.0 };
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        tFile.load_data( "vec", tBuffer, 4 );

        for ( belfem::index_t k = 0; k < 4; ++k )
        {
            EXPECT_NEAR( tBuffer[ k ], tExpected( k ), tEps );
        }
    }
}

//-----------------------------------------------------------------------------
// 7. Group management and tree navigation
//-----------------------------------------------------------------------------

TEST_F( HDF5Test, CreateAndSelectGroup )
{
    this->reset_path();

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.create_group( "grp" );
        tFile.save_data( "value", 17u );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        tFile.select_group( "grp" );
        belfem::uint tValue = 0;
        tFile.load_data( "value", tValue );
        EXPECT_EQ( tValue, 17u );
        EXPECT_EQ( tFile.tree(), "/grp" );
    }
}

TEST_F( HDF5Test, NestedGroups )
{
    this->reset_path();

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        const hid_t tOuter = tFile.create_group( "outer" );
        tFile.save_data( "outer_value", 11u );
        tFile.create_group( "inner", tOuter );
        tFile.save_data( "inner_value", 22u );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );

        tFile.select_group( "outer" );
        belfem::uint tOuterValue = 0;
        tFile.load_data( "outer_value", tOuterValue );
        EXPECT_EQ( tOuterValue, 11u );
        EXPECT_EQ( tFile.tree(), "/outer" );

        tFile.select_group( "inner" );
        belfem::uint tInnerValue = 0;
        tFile.load_data( "inner_value", tInnerValue );
        EXPECT_EQ( tInnerValue, 22u );
        EXPECT_EQ( tFile.tree(), "/outer/inner" );
    }
}

TEST_F( HDF5Test, MultipleGroupsSameLevel )
{
    this->reset_path();

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.create_group( "grp1" );
        tFile.save_data( "value", 1u );
        tFile.close_active_group();

        tFile.create_group( "grp2" );
        tFile.save_data( "value", 2u );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );

        tFile.select_group( "grp1" );
        belfem::uint tA = 0;
        tFile.load_data( "value", tA );
        EXPECT_EQ( tA, 1u );
        tFile.close_active_group();

        tFile.select_group( "grp2" );
        belfem::uint tB = 0;
        tFile.load_data( "value", tB );
        EXPECT_EQ( tB, 2u );
    }
}

TEST_F( HDF5Test, CloseActiveGroupNavigatesBack )
{
    this->reset_path();

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        const hid_t tOuter = tFile.create_group( "outer" );
        tFile.save_data( "outer_val", 10u );
        tFile.create_group( "inner", tOuter );
        tFile.save_data( "inner_val", 20u );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        tFile.select_group( "outer" );
        tFile.select_group( "inner" );
        EXPECT_EQ( tFile.tree(), "/outer/inner" );

        tFile.close_active_group();
        EXPECT_EQ( tFile.tree(), "/outer" );

        // verify I/O resolves against the parent group, not just the path string
        belfem::uint tOuterVal = 0;
        tFile.load_data( "outer_val", tOuterVal );
        EXPECT_EQ( tOuterVal, 10u );
    }
}

TEST_F( HDF5Test, CloseTreeReturnsToRoot )
{
    this->reset_path();

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.save_data( "root_val", 99u );
        const hid_t tOuter = tFile.create_group( "outer" );
        tFile.create_group( "inner", tOuter );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        tFile.select_group( "outer" );
        tFile.select_group( "inner" );
        EXPECT_EQ( tFile.tree(), "/outer/inner" );

        tFile.close_tree();
        EXPECT_EQ( tFile.tree(), "/" );

        // verify I/O resolves against root, not just the path string
        belfem::uint tRootVal = 0;
        tFile.load_data( "root_val", tRootVal );
        EXPECT_EQ( tRootVal, 99u );
    }
}

TEST_F( HDF5Test, TreePathTrackingCorrect )
{
    this->reset_path();

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        const hid_t tOuter = tFile.create_group( "outer" );
        tFile.create_group( "inner", tOuter );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );

        tFile.select_group( "outer" );
        EXPECT_EQ( tFile.tree(), "/outer" );

        tFile.select_group( "inner" );
        EXPECT_EQ( tFile.tree(), "/outer/inner" );

        tFile.close_active_group();
        EXPECT_EQ( tFile.tree(), "/outer" );

        tFile.close_tree();
        EXPECT_EQ( tFile.tree(), "/" );
    }
}

//-----------------------------------------------------------------------------
// 8. Error paths
//-----------------------------------------------------------------------------

TEST_F( HDF5Test, EmptyPathThrows )
{
    EXPECT_THROW( belfem::HDF5 tFile( "", belfem::FileMode::NEW ), std::runtime_error );
}

TEST_F( HDF5Test, SaveDuplicateLabelThrows )
{
    this->reset_path();

    belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
    tFile.save_data( "x", 1u );
    EXPECT_THROW( tFile.save_data( "x", 2u ), std::runtime_error );
}

TEST_F( HDF5Test, LoadMissingLabelThrows )
{
    this->reset_path();

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        belfem::real tValue = 0.0;
        EXPECT_THROW( tFile.load_data( "missing", tValue ), std::runtime_error );
    }
}

TEST_F( HDF5Test, OpenNonexistentFileRdonlyThrows )
{
    this->reset_path();
    remove_file_if_exists( mPath );

    EXPECT_THROW( belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY ), std::runtime_error );
}

TEST_F( HDF5Test, OpenNonexistentFileRdwrThrows )
{
    this->reset_path();
    remove_file_if_exists( mPath );

    EXPECT_THROW( belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDWR ), std::runtime_error );
}

TEST_F( HDF5Test, CreateDuplicateGroupThrows )
{
    this->reset_path();

    belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
    tFile.create_group( "grp" );
    tFile.close_active_group();
    EXPECT_THROW( tFile.create_group( "grp" ), std::runtime_error );
}

TEST_F( HDF5Test, SelectMissingGroupThrows )
{
    this->reset_path();

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        EXPECT_THROW( tFile.select_group( "nonexistent" ), std::runtime_error );
    }
}

//-----------------------------------------------------------------------------
// 9. Mixed payload integration tests
//-----------------------------------------------------------------------------

TEST_F( HDF5Test, MixedTypesInRootGroup )
{
    this->reset_path();

    const belfem::string tText = "root";
    const belfem::Vector< belfem::real > tVec = { 1.0, 2.0, 3.0 };
    belfem::Matrix< belfem::real > tMat( 2, 2 );
    tMat( 0, 0 ) = 1.5;
    tMat( 0, 1 ) = 2.5;
    tMat( 1, 0 ) = 3.5;
    tMat( 1, 1 ) = 4.5;
    const belfem::Cell< belfem::string > tCell = { "a", "b" };

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.save_data( "scalar_int", 9u );
        tFile.save_data( "scalar_real", 1.25 );
        tFile.save_data( "text", tText );
        tFile.save_data( "vec", tVec );
        tFile.save_data( "mat", tMat );
        tFile.save_data( "cell", tCell );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );

        belfem::uint tInt = 0;
        belfem::real tReal = 0.0;
        belfem::string tLoadedText;
        belfem::Vector< belfem::real > tLoadedVec;
        belfem::Matrix< belfem::real > tLoadedMat;
        belfem::Cell< belfem::string > tLoadedCell;

        tFile.load_data( "scalar_int", tInt );
        tFile.load_data( "scalar_real", tReal );
        tFile.load_data( "text", tLoadedText );
        tFile.load_data( "vec", tLoadedVec );
        tFile.load_data( "mat", tLoadedMat );
        tFile.load_data( "cell", tLoadedCell );

        EXPECT_EQ( tInt, 9u );
        EXPECT_NEAR( tReal, 1.25, tEps );
        EXPECT_EQ( tLoadedText, tText );
        expect_vector_near( tLoadedVec, tVec );
        expect_matrix_near( tLoadedMat, tMat );
        expect_cell_string_eq( tLoadedCell, tCell );
    }
}

TEST_F( HDF5Test, MixedTypesInSubGroup )
{
    this->reset_path();

    const belfem::Vector< belfem::sint > tVec = { -2, -1, 0, 1, 2 };
    belfem::Matrix< belfem::uint > tMat( 2, 2 );
    tMat( 0, 0 ) = 10u;
    tMat( 0, 1 ) = 20u;
    tMat( 1, 0 ) = 30u;
    tMat( 1, 1 ) = 40u;

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::NEW );
        tFile.create_group( "payload" );
        tFile.save_data( "flag", true );
        tFile.save_data( "name", belfem::string( "payload" ) );
        tFile.save_data( "vec", tVec );
        tFile.save_data( "mat", tMat );
    }

    {
        belfem::HDF5 tFile( mPath, belfem::FileMode::OPEN_RDONLY );
        tFile.select_group( "payload" );

        bool tFlag = false;
        belfem::string tName;
        belfem::Vector< belfem::sint > tLoadedVec;
        belfem::Matrix< belfem::uint > tLoadedMat;

        tFile.load_data( "flag", tFlag );
        tFile.load_data( "name", tName );
        tFile.load_data( "vec", tLoadedVec );
        tFile.load_data( "mat", tLoadedMat );

        EXPECT_TRUE( tFlag );
        EXPECT_EQ( tName, "payload" );
        expect_vector_eq( tLoadedVec, tVec );
        expect_matrix_eq( tLoadedMat, tMat );
    }
}

#endif // BELFEM_HDF5
