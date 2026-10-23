/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_HDF5_TOOLS_HPP
#define BELFEM_HDF5_TOOLS_HPP

#include <cstring>

#include "stringtools.hpp"
#include "assert.hpp"
#include "hdf5_types.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"

namespace belfem
{
    namespace hdf5
    {
//------------------------------------------------------------------------------

#ifdef BELFEM_HDF5
        /**
         * Scope owner of the three identifiers a dataset read or write opens.
         * Each member is assigned at the call that opens it; the destructor
         * closes every member that is a valid identifier, so a refusal that
         * throws between open and scope exit leaves no handle behind.
         * Close order: dataspace, datatype, dataset.
         */
        struct DatasetHandles
        {
            hid_t mDataSet   = H5I_INVALID_HID ;
            hid_t mDataType  = H5I_INVALID_HID ;
            hid_t mDataSpace = H5I_INVALID_HID ;

            DatasetHandles() = default ;

            // a copy or a move would close the same identifiers twice
            DatasetHandles( const DatasetHandles & ) = delete ;
            DatasetHandles & operator=( const DatasetHandles & ) = delete ;
            DatasetHandles( DatasetHandles && ) = delete ;
            DatasetHandles & operator=( DatasetHandles && ) = delete ;

            ~DatasetHandles()
            {
                if( mDataSpace >= 0 )
                {
                    H5Sclose( mDataSpace );
                }
                if( mDataType >= 0 )
                {
                    H5Tclose( mDataType );
                }
                if( mDataSet >= 0 )
                {
                    H5Dclose( mDataSet );
                }
            }
        };
#endif

//------------------------------------------------------------------------------

        /**
         * test if a dataset exists ( tests for a link of that name;
         * the link type is not checked )
         */
        inline bool
        dataset_exists( hid_t aLoc, const std::string & aLabel )
        {
#ifdef BELFEM_HDF5
            hid_t tDataSet = 0;
            return H5Lexists( aLoc, aLabel.c_str(), tDataSet );
#else
            return false;
#endif
        }

//------------------------------------------------------------------------------

        /**
         * test if a group exists ( tests for a link of that name;
         * the link type is not checked )
         */
        inline bool
        group_exists( hid_t aLoc, const std::string & aLabel )
        {
#ifdef BELFEM_HDF5
            hid_t tGroup = 0;
            return H5Lexists( aLoc, aLabel.c_str(), tGroup );
#else
            return false;
#endif
        }

//------------------------------------------------------------------------------

// The whole function is guarded, not just its body: its SIGNATURE names
// H5L_info2_t, which hdf5_types.hpp does not stub. With the body alone guarded
// this header stopped compiling under -DUSE_HDF5=OFF, and since cl_HDF5.hpp
// includes it and src/io always compiles cl_HDF5.cpp, that broke the whole
// tree rather than only the HDF5 callers
#ifdef BELFEM_HDF5

// collect_names() and get_groups() below are written against the 1.12 link and
// object interfaces ( H5L_info2_t, H5O_info2_t, H5Literate2,
// H5Oget_info_by_name3 ). belfem_find_package carries no version test, so the
// requirement is stated here, where it is created -- and this header is
// included across the whole tree, so on an older HDF5 the failure would
// otherwise be a wall of template errors in every module rather than one line
#if ! H5_VERSION_GE( 1, 12, 0 )
#error "BELFEM requires HDF5 1.12 or newer ( hdf5_tools.hpp uses H5Literate2 and H5Oget_info_by_name3 )."
#endif

        inline herr_t
        collect_names(
        hid_t                aLoc,
        const char         * aLabel,
        const H5L_info2_t  * aInfo,
        void               * aData )
        {
            auto * tLabels = static_cast<std::vector<std::string> *>( aData );

            // GROUPS only. H5Literate2 visits every link at the location,
            // datasets and named datatypes included, and a caller that opens
            // each returned name as a group fails on the first one that is not.
            // A soft or external link is skipped as well rather than followed:
            // resolving it can leave the location, which is not what listing a
            // location means
            if( aInfo != nullptr && aInfo->type != H5L_TYPE_HARD )
            {
                return 0;
            }

            H5O_info2_t tInfo ;

            // a failed query is an ERROR, not an entry to skip. Returning 0
            // here would let a permission failure or a corrupt object produce a
            // short list that the caller cannot tell from a complete one. The
            // callback contract is: negative aborts the iteration and becomes
            // H5Literate2's return value, which get_groups() checks
            if( H5Oget_info_by_name3( aLoc, aLabel, & tInfo,
                                      H5O_INFO_BASIC, H5P_DEFAULT ) < 0 )
            {
                return -1;
            }

            // not a group is not an error -- it is the thing this filter exists
            // to drop, so the iteration continues
            if( tInfo.type != H5O_TYPE_GROUP )
            {
                return 0;
            }

            tLabels->emplace_back( aLabel );

            return 0;
        }

#endif

        /**
         * list groups that exist in location
         */
        inline Cell< std::string >
        get_groups( hid_t aLoc )
        {
#ifdef BELFEM_HDF5
            Cell< std::string > tGroups ;

            hsize_t tIndex = 0 ;

            // H5_ITER_INC, not H5_ITER_NATIVE : native order is documented as
            // whatever the library finds fastest, and callers that treat the
            // first entry as a reference ( db2exo picks the grid of the first
            // table it sees ) would then depend on it
            const herr_t tStatus = H5Literate2(
                        aLoc,
                        H5_INDEX_NAME,
                        H5_ITER_INC,
                        & tIndex,
                        collect_names,
                        & tGroups.vector_data() );

            // runs once per file, so the always-active tier is affordable and
            // the message survives a release build. Silently returning a short
            // list here would be indistinguishable from a location that really
            // holds fewer groups
            BELFEM_ERROR( tStatus >= 0,
                "H5Literate2 failed while listing the groups of an HDF5 location." );

            return tGroups ;

#else
            return {};
#endif
        }

//------------------------------------------------------------------------------

        /**
         * this function returns true of both the HDF5 datatype
         * and the passed datatype have the same size
         */
        template < typename T >
        bool
        check_datatye_size()
        {
#ifdef BELFEM_HDF5
            return H5Tget_size( datatype<T>() ) == sizeof( T );
#else
            return false;
#endif
        }

//------------------------------------------------------------------------------

        /**
         * Guard for every read: the FILE datatype must be convertible into T.
         *
         * The loaders pass datatype<T>() as the MEMORY type, so HDF5 converts
         * between widths and a file written with another index width
         * ( index_t / int_t follow BELFEM_INT64, as in SpMatrix::save and the
         * .bfm meshes ) stays readable. Passing the file type instead would
         * mean "no conversion": HDF5 would copy at the file's element size
         * and overrun or half-fill the object, silently.
         *
         * This function refuses what conversion cannot fix: a different
         * datatype CLASS ( a float dataset read into an integer, say ), where
         * HDF5 would produce numerically meaningless values. Always-active
         * tier, because the mismatch shows up in release builds.
         */
        template < typename T >
        void
        check_read_datatype( const hid_t aFileType, const std::string & aLabel )
        {
#ifdef BELFEM_HDF5
            BELFEM_ERROR( H5Tget_class( aFileType ) == H5Tget_class( datatype<T>() ),
                "datatype class mismatch reading '%s': the file holds class %i, but %s is class %i.\n"
                "HDF5 can convert between widths, not between classes - this file does not hold what the reader expects.",
                aLabel.c_str(),
                ( int ) H5Tget_class( aFileType ),
                datatype_string<T>().c_str(),
                ( int ) H5Tget_class( datatype<T>() ) );
#endif
        }

//------------------------------------------------------------------------------

        /**
        * saves a scalar value to a file
                * file must be open
        *
        * @param[ in ]    aLoc     location ( file or group ) to write into
        * @param[ in ]    aLabel   dataset name
        * @param[ in ]    aValue   value that is to be stored
        * @param[ out ]   aStatus  HDF5 return code of the write
        */
        template < typename T >
        void
        save_scalar_to_file(
                hid_t               & aLoc,
                const std::string   & aLabel,
                const T             & aValue,
                herr_t              & aStatus
        )
        {
#ifdef BELFEM_HDF5
            BELFEM_ERROR( ! hdf5::dataset_exists( aLoc, aLabel ),
                    "Dataset %s of type %s does already exist.",
                    aLabel.c_str(),
                    datatype_string<T>().c_str() );

            BELFEM_ASSERT( hdf5::check_datatye_size<T>(),
                    "Error in datatype size of type %s.",
                        datatype_string<T>().c_str() );

            // FILE type: fixed width and little endian, so the on-disk
            // layout does not depend on this build's native type widths
            hid_t tDataType = H5Tcopy( filetype<T>() );

            hsize_t tDims[ 1 ] = { 1 };

            hid_t tDataSpace
                    = H5Screate_simple( 1, tDims, nullptr );

            hid_t tDataSet = H5Dcreate(
                    aLoc,
                    aLabel.c_str(),
                    tDataType,
                    tDataSpace,
                    H5P_DEFAULT,
                    H5P_DEFAULT,
                    H5P_DEFAULT );

            aStatus = H5Dwrite(
                    tDataSet,
                    datatype< T >(),   // MEMORY type: native, not the file type
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    & aValue );

            H5Sclose( tDataSpace );
            H5Tclose( tDataType );
            H5Dclose( tDataSet );

            BELFEM_ERROR( aStatus >= 0,
                    "Something went wrong while trying to store scalar %s of type %s",
                       aLabel.c_str(),
                       datatype_string<T>().c_str());
#endif
        }
//------------------------------------------------------------------------------

        template < typename T >
        void
        load_scalar_from_file(
                hid_t               & aLoc,
                const std::string   & aLabel,
                T                   & aValue,
                herr_t              & aStatus
        )
        {
#ifdef BELFEM_HDF5
            BELFEM_ERROR( hdf5::dataset_exists( aLoc, aLabel ),
                       "Dataset %s of type %s does not exist.",
                       aLabel.c_str(),
                       datatype_string<T>().c_str() );

            BELFEM_ASSERT( hdf5::check_datatye_size<T>(),
                        "Error in datatype size of type %s.",
                        datatype_string<T>().c_str() );

            DatasetHandles tHandles ;

            tHandles.mDataSet = H5Dopen1( aLoc, aLabel.c_str() );

            tHandles.mDataType = H5Dget_type( tHandles.mDataSet );

            // the file type must be convertible into T; the read
            // below then uses datatype<T>() so HDF5 does the conversion
            check_read_datatype< T >( tHandles.mDataType, aLabel );

            tHandles.mDataSpace = H5Dget_space( tHandles.mDataSet );

            aStatus = H5Dread(
                    tHandles.mDataSet,
                    datatype< T >(),   // MEMORY type: let HDF5 convert
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    &aValue );

            BELFEM_ERROR( aStatus >= 0,
                       "Something went wrong while trying to load scalar %s of type %s",
                       aLabel.c_str(),
                       datatype_string<T>().c_str() );

#endif
        }

//------------------------------------------------------------------------------

        void
        save_bool_to_file(
                hid_t               & aLoc,
                const std::string   & aLabel,
                const bool          & aValue,
                herr_t              & aStatus );

//------------------------------------------------------------------------------

        void
        load_bool_from_file(
                            hid_t   & aLoc,
                const std::string   & aLabel,
                             bool   & aValue,
                           herr_t   & aStatus );

//------------------------------------------------------------------------------

        void
        save_string_to_file(
                hid_t               & aLoc,
                const std::string   & aLabel,
                const std::string   & aValue,
                herr_t              & aStatus );

//------------------------------------------------------------------------------

        void
        load_string_from_file(
                      hid_t         & aLoc,
                const std::string   & aLabel,
                      std::string   & aValue,
                      herr_t        & aStatus );

//------------------------------------------------------------------------------
        template < typename T >
        void
        save_array_to_file(
                      hid_t      & aLoc,
                const std::string & aLabel,
                const T           * aData,
                const hsize_t     & aLength,
                      herr_t      & aStatus )
        {
#ifdef BELFEM_HDF5
            BELFEM_ERROR( ! hdf5::dataset_exists( aLoc, aLabel ),
                       "Dataset %s of type Vector<%s> does already exist.",
                       aLabel.c_str(),
                       datatype_string<T>().c_str() );

            BELFEM_ASSERT( hdf5::check_datatye_size<T>(),
                        "Error in datatype size of type %s.",
                        datatype_string<T>().c_str() );

            // FILE type: fixed width and little endian, so the on-disk
            // layout does not depend on this build's native type widths
            hid_t tDataType = H5Tcopy( filetype<T>() );

            hsize_t tDims[ 1 ];
            tDims[ 0 ] = aLength;

            hid_t  tDataSpace
                    = H5Screate_simple( 1, tDims, nullptr );

            hid_t tDataSet = H5Dcreate(
                    aLoc,
                    aLabel.c_str(),
                    tDataType,
                    tDataSpace,
                    H5P_DEFAULT,
                    H5P_DEFAULT,
                    H5P_DEFAULT );

            // an empty object skips the write, so its status cannot report a
            // failed create; the open handles are closed before the error
            if ( tDataSet < 0 )
            {
                H5Sclose( tDataSpace );
                H5Tclose( tDataType );
                BELFEM_ERROR( false, "Could not create dataset %s", aLabel.c_str() );
            }

            aStatus = 0 ;

            if( aLength > 0 )
            {
                aStatus = H5Dwrite(
                        tDataSet,
                        datatype< T >(),   // MEMORY type: native, not the file type
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,
                        &aData[ 0 ] );
            }

            H5Sclose( tDataSpace );
            H5Tclose( tDataType );
            H5Dclose( tDataSet );

            BELFEM_ERROR( aStatus >= 0,
                        "Something went wrong while trying to store vector %s of type %s",
                        aLabel.c_str(),
                        datatype_string<T>().c_str() );
#endif
        }

//------------------------------------------------------------------------------

        hsize_t
        get_array_size(
        hid_t       & aLoc,
        const std::string & aLabel  );

//------------------------------------------------------------------------------

        template < typename T >
        void
        load_array_from_file(
                      hid_t       & aLoc,
                const std::string & aLabel,
                      T           * aData,
                const hsize_t       aLength,
                      herr_t      & aStatus )
        {
#ifdef BELFEM_HDF5

            BELFEM_ERROR( hdf5::dataset_exists( aLoc, aLabel ),
                       "Dataset %s of type Vector<%s> does not exist.",
                       aLabel.c_str(),
                       datatype_string<T>().c_str() );

            BELFEM_ASSERT( hdf5::check_datatye_size<T>(),
                        "Error in datatype size of type %s.",
                        datatype_string<T>().c_str() );

            DatasetHandles tHandles ;

            tHandles.mDataSet = H5Dopen1( aLoc, aLabel.c_str() );

            tHandles.mDataType = H5Dget_type( tHandles.mDataSet );

            // the file type must be convertible into T; the read
            // below then uses datatype<T>() so HDF5 does the conversion
            check_read_datatype< T >( tHandles.mDataType, aLabel );

            tHandles.mDataSpace = H5Dget_space( tHandles.mDataSet );

            hsize_t tDims[ 1 ];

            // the extent buffer holds one entry; a dataset of another rank is
            // refused before the extents are read into it
            const int tRank = H5Sget_simple_extent_ndims( tHandles.mDataSpace );

            BELFEM_ERROR( tRank == 1, "Dataset %s is not a 1-D array (rank %i)",
                aLabel.c_str(), tRank );

            aStatus  = H5Sget_simple_extent_dims( tHandles.mDataSpace, tDims, nullptr );

            BELFEM_ERROR( tDims[ 0 ] == aLength, "Lengths do not match: is %lu, expect %lu.",
                         ( long unsigned int ) tDims[ 0 ], ( long unsigned int ) aLength );

            if( aLength > 0 )
            {
                aStatus = H5Dread(
                        tHandles.mDataSet,
                        datatype< T >(),   // MEMORY type: let HDF5 convert
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,
                        &aData[ 0 ] );
            }
            else
            {
                aStatus = 0;
            }

            BELFEM_ERROR( aStatus >= 0,
                        "Something went wrong while trying to load vector %s of type %s",
                        aLabel.c_str(),
                        datatype_string<T>().c_str() );
#endif
        }
//------------------------------------------------------------------------------

        template < typename T >
        void
        save_vector_to_file(
                hid_t         & aLoc,
                const std::string   & aLabel,
                const Vector< T >   & aVector,
                herr_t        & aStatus )
        {
#ifdef BELFEM_HDF5

            BELFEM_ERROR( ! hdf5::dataset_exists( aLoc, aLabel ),
                       "Dataset %s of type Vector<%s> does already exist.",
                       aLabel.c_str(),
                       datatype_string<T>().c_str() );

            BELFEM_ASSERT( hdf5::check_datatye_size<T>(),
                        "Error in datatype size of type %s.",
                        datatype_string<T>().c_str() );

            // FILE type: fixed width and little endian, so the on-disk
            // layout does not depend on this build's native type widths
            hid_t tDataType = H5Tcopy( filetype<T>() );
            hid_t tDataSet = 0;

            hsize_t tLength  = aVector.length();
            hsize_t tDims[ 1 ];
            tDims[ 0 ] = tLength;

            hid_t  tDataSpace
                    = H5Screate_simple( 1, tDims, nullptr );

            tDataSet = H5Dcreate(
                    aLoc,
                    aLabel.c_str(),
                    tDataType,
                    tDataSpace,
                    H5P_DEFAULT,
                    H5P_DEFAULT,
                    H5P_DEFAULT );

            // an empty object skips the write, so its status cannot report a
            // failed create; the open handles are closed before the error
            if ( tDataSet < 0 )
            {
                H5Sclose( tDataSpace );
                H5Tclose( tDataType );
                BELFEM_ERROR( false, "Could not create dataset %s", aLabel.c_str() );
            }

            aStatus = 0 ;

            if( tLength > 0 )
            {
                T* tData = ( T* ) malloc( tLength * sizeof( T ) );

                std::memcpy(
                        tData,
                        aVector.data(),
                        tLength * sizeof( T ) );

                aStatus = H5Dwrite(
                        tDataSet,
                        datatype< T >(),   // MEMORY type: native, not the file type
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,
                        & tData[ 0 ]);

                free( tData );
            }

            H5Sclose( tDataSpace );
            H5Tclose( tDataType );
            H5Dclose( tDataSet );

            BELFEM_ERROR( aStatus >= 0,
                        "Something went wrong while trying to store vector %s of type %s",
                        aLabel.c_str(),
                        datatype_string<T>().c_str() );
#endif
        }

//------------------------------------------------------------------------------

        template < typename T >
        void
        load_vector_from_file(
                hid_t         & aLoc,
                const std::string   & aLabel,
                Vector< T >   & aVector,
                herr_t        & aStatus )
        {
#ifdef BELFEM_HDF5
            BELFEM_ERROR( hdf5::dataset_exists( aLoc, aLabel ),
                       "Dataset %s of type Vector<%s> does not exist.",
                       aLabel.c_str(),
                       datatype_string<T>().c_str() );

            BELFEM_ASSERT( hdf5::check_datatye_size<T>(),
                        "Error in datatype size of type %s.",
                        datatype_string<T>().c_str() );

            DatasetHandles tHandles ;

            tHandles.mDataSet = H5Dopen1( aLoc, aLabel.c_str() );

            tHandles.mDataType = H5Dget_type( tHandles.mDataSet );

            // the file type must be convertible into T; the read
            // below then uses datatype<T>() so HDF5 does the conversion
            check_read_datatype< T >( tHandles.mDataType, aLabel );

            tHandles.mDataSpace = H5Dget_space( tHandles.mDataSet );

            hsize_t tDims[ 1 ];

            // the extent buffer holds one entry; a dataset of another rank is
            // refused before the extents are read into it
            const int tRank = H5Sget_simple_extent_ndims( tHandles.mDataSpace );

            BELFEM_ERROR( tRank == 1, "Dataset %s is not a 1-D array (rank %i)",
                aLabel.c_str(), tRank );

            aStatus  = H5Sget_simple_extent_dims( tHandles.mDataSpace, tDims, nullptr );

            hsize_t tLength = tDims[ 0 ];

            aVector.set_size( tLength );

            if( tLength > 0 )
            {
                T* tData = ( T* ) malloc( tLength * sizeof( T ) );

                aStatus = H5Dread(
                        tHandles.mDataSet,
                        datatype< T >(),   // MEMORY type: let HDF5 convert
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,
                        &tData[ 0 ] );

                std::memcpy(
                        aVector.data(),
                        tData,
                        tLength * sizeof( T ) );

                free( tData );
            }
            else
            {
                aStatus = 0;
            }

            BELFEM_ERROR( aStatus >= 0,
                        "Something went wrong while trying to load vector %s of type %s",
                        aLabel.c_str(),
                        datatype_string<T>().c_str() );
#endif
        }
//------------------------------------------------------------------------------

        template < typename T >
        void
        save_matrix_to_file(
                hid_t               & aLoc,
                const std::string   & aLabel,
                const Matrix< T >   & aMatrix,
                herr_t              & aStatus,
                const bool            aTranspose = false )
        {
#ifdef BELFEM_HDF5
            BELFEM_ERROR( ! hdf5::dataset_exists( aLoc, aLabel ),
                       "Dataset %s of type Matrix<%s> does already exist.",
                       aLabel.c_str(),
                       datatype_string<T>().c_str() );

            BELFEM_ASSERT( hdf5::check_datatye_size<T>(),
                        "Error in datatype size of type %s.",
                       datatype_string<T>().c_str() );

            // FILE type: fixed width and little endian, so the on-disk
            // layout does not depend on this build's native type widths
            hid_t tDataType = H5Tcopy( filetype<T>() );
            hid_t tDataSet = 0;

            hsize_t tDims[ 2 ];

            if ( aTranspose )
            {
                tDims[ 0 ] = aMatrix.n_cols();
                tDims[ 1 ] = aMatrix.n_rows();

            }
            else
            {
                tDims[ 0 ] = aMatrix.n_rows();
                tDims[ 1 ] = aMatrix.n_cols();
            }

            hid_t  tDataSpace
                    = H5Screate_simple( 2, tDims, nullptr );

            tDataSet = H5Dcreate(
                    aLoc,
                    aLabel.c_str(),
                    tDataType,
                    tDataSpace,
                    H5P_DEFAULT,
                    H5P_DEFAULT,
                    H5P_DEFAULT );

            // an empty object skips the write, so its status cannot report a
            // failed create; the open handles are closed before the error
            if ( tDataSet < 0 )
            {
                H5Sclose( tDataSpace );
                H5Tclose( tDataType );
                BELFEM_ERROR( false, "Could not create dataset %s", aLabel.c_str() );
            }

            aStatus = 0 ;

            if( tDims[ 0 ]*tDims[ 1 ] > 0 )
            {
                T** tData = ( T** ) malloc( tDims[ 0 ] * sizeof( T * ) );
                tData[ 0 ] = ( T* ) malloc( tDims[ 0 ]*tDims[ 1 ] * sizeof( T ) );

                for( hsize_t i=0; i<tDims[ 0 ]; ++i )
                {
                    tData[ i ] = tData[ 0 ]+ i*tDims[ 1 ];
                }

                if ( aTranspose )
                {
                    for ( hsize_t i = 0; i < tDims[ 0 ]; ++i )
                    {
                        for ( hsize_t j = 0; j < tDims[ 1 ]; ++j )
                        {
                            tData[ i ][ j ] = aMatrix( j, i );
                        }
                    }
                }
                else
                {
                    for ( hsize_t i = 0; i < tDims[ 0 ]; ++i )
                    {
                        for ( hsize_t j = 0; j < tDims[ 1 ]; ++j )
                        {
                            tData[ i ][ j ] = aMatrix( i, j );
                        }
                    }
                }
                aStatus = H5Dwrite(
                        tDataSet,
                        datatype< T >(),   // MEMORY type: native, not the file type
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,
                        & tData[ 0 ][ 0 ]);

                free( tData[ 0 ] );
                free( tData );
            }

            H5Sclose( tDataSpace );
            H5Tclose( tDataType );
            H5Dclose( tDataSet );

            BELFEM_ERROR( aStatus >= 0,
                       "Something went wrong while trying to store matrix %s of type %s",
                       aLabel.c_str(),
                       datatype_string<T>().c_str() );
#endif
        }

//------------------------------------------------------------------------------

        template < typename T >
        void
        load_matrix_from_file(
                hid_t               & aLoc,
                const std::string   & aLabel,
                      Matrix< T >   & aMatrix,
                herr_t              & aStatus,
                const bool            aTranspose = false )
        {
#ifdef BELFEM_HDF5

            BELFEM_ERROR( hdf5::dataset_exists( aLoc, aLabel ),
                       "Dataset %s of type Matrix<%s> does not exist.",
                       aLabel.c_str(),
                       datatype_string<T>().c_str() );

            BELFEM_ASSERT( hdf5::check_datatye_size<T>(),
                        "Error in datatype size of type %s.",
                        datatype_string<T>().c_str() );

            DatasetHandles tHandles ;

            tHandles.mDataSet = H5Dopen1( aLoc, aLabel.c_str() );

            tHandles.mDataType = H5Dget_type( tHandles.mDataSet );

            // the file type must be convertible into T; the read
            // below then uses datatype<T>() so HDF5 does the conversion
            check_read_datatype< T >( tHandles.mDataType, aLabel );

            tHandles.mDataSpace = H5Dget_space( tHandles.mDataSet );

            hsize_t tDims[ 2 ];

            aStatus  = H5Sget_simple_extent_dims( tHandles.mDataSpace, tDims, nullptr );

            // allocate memory for output (transpose swaps the on-disk dimensions)
            if ( aTranspose )
            {
                aMatrix.set_size( tDims[ 1 ], tDims[ 0 ] );
            }
            else
            {
                aMatrix.set_size( tDims[ 0 ], tDims[ 1 ] );
            }

            if( tDims[ 0 ]*tDims[ 1 ] > 0 )
            {
                T** tData = ( T** )
                        malloc( tDims[ 0 ]*sizeof( T* ) );

                tData[ 0 ] = ( T* )
                        malloc( tDims[ 0 ]*  tDims[ 1 ] * sizeof( T ) );

                for( hsize_t i=1; i<tDims[ 0 ]; ++i )
                {
                    tData[ i ] = tData[ 0 ]+ i*tDims[ 1 ];
                }

                aStatus = H5Dread(
                        tHandles.mDataSet,
                        datatype< T >(),   // MEMORY type: let HDF5 convert
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,
                        &tData[ 0 ][ 0 ] );

                if ( aTranspose )
                {
                    // on-disk (tDims[0] x tDims[1]) -> memory (tDims[1] x tDims[0])
                    for ( hsize_t i = 0; i < tDims[ 0 ]; ++i )
                    {
                        for ( hsize_t j = 0; j < tDims[ 1 ]; ++j )
                        {
                            aMatrix( j, i ) = tData[ i ][ j ];
                        }
                    }
                }
                else
                {
                    for ( hsize_t j = 0; j < tDims[ 1 ]; ++j )
                    {
                        for ( hsize_t i = 0; i < tDims[ 0 ]; ++i )
                        {
                            aMatrix( i, j ) = tData[ i ][ j ];
                        }
                    }
                }
                free( tData[ 0 ] );
                free( tData );
            }
            else if( aStatus == 2 )
            {
                aStatus = 0;
            }

            BELFEM_ERROR( aStatus >= 0,
                        "Something went wrong while trying to load matrix %s of type %s",
                        aLabel.c_str(),
                        datatype_string<T>().c_str() );
#endif
        }

        inline void
        save_strings_to_file( hid_t & aLoc,
                const std::string   & aLabel,
                const Cell< string >  & aCell,
                herr_t              & aStatus )
        {
#ifdef BELFEM_HDF5

            DatasetHandles tHandles ;

            hsize_t tNumStrings = aCell.size();
            hsize_t tDim[ 1 ] = { tNumStrings };
            tHandles.mDataSpace = H5Screate_simple( 1, tDim, NULL );

            tHandles.mDataType = H5Tcopy( H5T_C_S1 );
            H5Tset_size( tHandles.mDataType, H5T_VARIABLE );

            tHandles.mDataSet = H5Dcreate2(
                    aLoc,
                    aLabel.c_str(),
                    tHandles.mDataType,
                    tHandles.mDataSpace,
                    H5P_DEFAULT,
                    H5P_DEFAULT,
                    H5P_DEFAULT );

            BELFEM_ERROR( tHandles.mDataSet > -1, "Error creating dataset" );

            Cell< const char * > tBuffer( tNumStrings, nullptr );
            for ( hsize_t i=0; i<tNumStrings; ++i )
            {
                tBuffer( i ) = aCell( i ).c_str();
            }

            aStatus = H5Dwrite(
                    tHandles.mDataSet,
                    tHandles.mDataType,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    tBuffer.data() );

            BELFEM_ERROR( aStatus ==0, "Error writing data" );
#endif
        }

        inline void
        load_strings_from_file( hid_t & aLoc,
            const std::string   & aLabel,
            Cell< string >  & aCell,
            herr_t              & aStatus )
        {
#ifdef BELFEM_HDF5
            aCell.clear();

            BELFEM_ERROR( hdf5::dataset_exists( aLoc, aLabel ),
                "Dataset %s does not exist at expected location.", aLabel.c_str() );

            DatasetHandles tHandles ;

            tHandles.mDataSet = H5Dopen1( aLoc, aLabel.c_str() );

            tHandles.mDataSpace = H5Dget_space( tHandles.mDataSet );

            hsize_t tDim[ 1 ] ;

            // the extent buffer holds one entry; a dataset of another rank is
            // refused before the extents are read into it
            const int tRank = H5Sget_simple_extent_ndims( tHandles.mDataSpace );

            BELFEM_ERROR( tRank == 1, "Dataset %s is not a 1-D array (rank %i)",
                aLabel.c_str(), tRank );

            BELFEM_ERROR( H5Sget_simple_extent_dims( tHandles.mDataSpace, tDim, nullptr ) > -1 ,
                "Error getting dataspace dimensions" );

            hsize_t tNumStrings = tDim[ 0 ];

            // This must match the datatype used in writing -- INCLUDING the
            // character set: H5T_C_S1 defaults to ASCII, but h5py writes
            // vlen strings as UTF-8, and HDF5 refuses the UTF-8 -> ASCII
            // conversion path, failing the H5Dread below with no better
            // diagnostic than a nonzero status. Copy the cset from the
            // dataset itself so BELFEM-written ASCII files behave as before
            // and Python-written UTF-8 files become readable.
            tHandles.mDataType = H5Tcopy( H5T_C_S1 );
            H5Tset_size( tHandles.mDataType, H5T_VARIABLE );

            hid_t tFileType = H5Dget_type( tHandles.mDataSet );
            H5Tset_cset( tHandles.mDataType, H5Tget_cset( tFileType ) );
            H5Tclose( tFileType );

            // HDF5 allocates the strings during the read. The null slots let
            // the reclaim below free exactly what a failed read left allocated,
            // so it runs on both paths, before the refusal and before the owner
            // closes the type and the dataspace it needs. An empty Cell may
            // hand out a null data(), which H5Dread accepts at zero elements
            // and H5Dvlen_reclaim refuses, hence the size guard.
            Cell< char * > tBuffer( tNumStrings, nullptr );

            const herr_t tReadStatus = H5Dread( tHandles.mDataSet, tHandles.mDataType,
                    H5S_ALL, H5S_ALL, H5P_DEFAULT, tBuffer.data() );

            if( tReadStatus == 0 )
            {
                aCell.set_size( tNumStrings, "" );

                for ( hsize_t i=0; i<tNumStrings; ++i )
                {
                    aCell( i ) = std::string( tBuffer( i ) );
                }
            }

            if( tNumStrings > 0 )
            {
                H5Dvlen_reclaim( tHandles.mDataType, tHandles.mDataSpace, H5P_DEFAULT, tBuffer.data() );
            }

            BELFEM_ERROR( tReadStatus == 0, "Failed to read strings from file." );
#endif
        }
//------------------------------------------------------------------------------

        inline hid_t
        open_group( const std::string & aLabel,  hid_t aParent )
        {
#ifdef BELFEM_HDF5
            if ( hdf5::group_exists( aParent, aLabel ) )
            {

                return H5Gopen2(
                        aParent,
                        aLabel.c_str(),
                        H5P_DEFAULT );
            }
            else
            {
                BELFEM_ERROR( false,
                              "Group %s does not exist.",
                              aLabel.c_str() );

                return -1;
            }
#else
            return 0 ;
#endif
        }

//------------------------------------------------------------------------------

        inline herr_t
        close_group( hid_t aGroup )
        {
#ifdef BELFEM_HDF5
            return H5Gclose( aGroup );
#else
            return 0 ;
#endif
        }

//------------------------------------------------------------------------------

        string
        create_tree( const Cell< string > & aTree, const string & aLabel );

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
    } /* namespace hdf5 */
} /* namespace belfem */

#endif //BELFEM_HDF5_TOOLS_HPP
