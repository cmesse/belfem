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

#include "hdf5_tools.hpp"

namespace belfem
{
    namespace hdf5
    {
//------------------------------------------------------------------------------

        void
        save_bool_to_file(
                      hid_t       & aFileID,
                const std::string & aLabel,
                const bool        & aValue,
                      herr_t      & aStatus )
        {
#ifdef BELFEM_HDF5
            BELFEM_ERROR( ! hdf5::dataset_exists( aFileID, aLabel ),
                      "Dataset %s of type %s does already exist.",
                      aLabel.c_str(),
                      datatype_string<bool>().c_str());

            // FILE type: fixed width and little endian, so the on-disk
            // layout does not depend on this build's native type widths
            hid_t tDataType = H5Tcopy( filetype< bool >() );

            hsize_t tDims[1] = {1};

            hid_t tDataSpace
                    = H5Screate_simple(1, tDims, nullptr);

            hid_t tDataSet = H5Dcreate(
                    aFileID,
                    aLabel.c_str(),
                    tDataType,
                    tDataSpace,
                    H5P_DEFAULT,
                    H5P_DEFAULT,
                    H5P_DEFAULT);

            hbool_t tValue = ( hbool_t ) aValue;

            aStatus = H5Dwrite(
                    tDataSet,
                    H5T_NATIVE_HBOOL,   // MEMORY type: native, not the file type
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    & tValue );

            H5Sclose(tDataSpace);
            H5Tclose(tDataType);
            H5Dclose(tDataSet);

            BELFEM_ASSERT( aStatus == 0,
                      "Something went wrong while trying to store variable %s of type %s",
                      aLabel.c_str(),
                      datatype_string<bool>().c_str());
#endif
        }

//------------------------------------------------------------------------------

        void
        load_bool_from_file(
                  hid_t       & aFileID,
            const std::string & aLabel,
                  bool        & aValue,
                  herr_t      & aStatus )
        {
#ifdef BELFEM_HDF5
            BELFEM_ERROR( hdf5::dataset_exists( aFileID, aLabel ),
                      "Dataset %s of type %s does not exist.",
                      aLabel.c_str(),
                      datatype_string<bool>().c_str());

            hid_t tDataSet = H5Dopen1( aFileID, aLabel.c_str() );

            hid_t tDataType = H5Dget_type( tDataSet );

            // the file type must be convertible into a bool; the read
            // below then uses datatype<bool>() so HDF5 does the conversion
            check_read_datatype< bool >( tDataType, aLabel );

            hid_t tDataSpace = H5Dget_space( tDataSet );

            hbool_t tValue;

            aStatus = H5Dread(
                    tDataSet,
                    datatype< bool >(),   // MEMORY type: let HDF5 convert
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    &tValue );

            aValue = ( bool ) tValue;

            H5Tclose( tDataType );
            H5Dclose( tDataSet );
            H5Sclose( tDataSpace );

            BELFEM_ASSERT( aStatus == 0,
                       "Something went wrong while trying to load variable %s of type %s",
                       aLabel.c_str(),
                       datatype_string<bool>().c_str());
#endif
        }

//------------------------------------------------------------------------------

        void
        save_string_to_file(
                hid_t               & aFileID,
                const std::string   & aLabel,
                const std::string   & aValue,
                herr_t              & aStatus )
        {
#ifdef BELFEM_HDF5
            BELFEM_ERROR(! hdf5::dataset_exists( aFileID, aLabel ),
                      "Dataset %s of type %s does already exist.",
                      aLabel.c_str(),
                      datatype_string<std::string>().c_str());

            // HDF5 fixed-size strings require length > 0
            BELFEM_ERROR( ! aValue.empty(),
                      "Cannot save empty string as dataset %s. "
                      "HDF5 does not support zero-length fixed-size strings.",
                      aLabel.c_str() );

            hid_t tDataType = H5Tcopy( H5T_C_S1 );

            aStatus  = H5Tset_size( tDataType, aValue.length() );

            hsize_t tDims[ 1 ] = { 1 };

            hid_t tDataSpace
                    = H5Screate_simple( 1, tDims, nullptr );

            hid_t tDataSet = H5Dcreate(
                    aFileID,
                    aLabel.c_str(),
                    tDataType,
                    tDataSpace,
                    H5P_DEFAULT,
                    H5P_DEFAULT,
                    H5P_DEFAULT );

            aStatus = H5Dwrite(
                    tDataSet,
                    tDataType,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    aValue.c_str() );

            H5Sclose( tDataSpace );
            H5Tclose( tDataType );
            H5Dclose( tDataSet );

            BELFEM_ASSERT(aStatus == 0,
                      "Something went wrong while trying to store string %s.",
                      aLabel.c_str() );
#endif
        }

//------------------------------------------------------------------------------

        void
        load_string_from_file(
                      hid_t         & aFileID,
                const std::string   & aLabel,
                      std::string   & aValue,
                      herr_t        & aStatus )
        {
#ifdef BELFEM_HDF5
            BELFEM_ERROR( hdf5::dataset_exists( aFileID, aLabel ),
                      "Dataset %s of type %s does not exist.",
                      aLabel.c_str(),
                      datatype_string<std::string>().c_str());

            hid_t tDataSet = H5Dopen1( aFileID, aLabel.c_str() );

            hid_t tDataSpace = H5Dget_space( tDataSet );

            hid_t tDataType = H5Dget_type( tDataSet );

            hsize_t tSize = H5Dget_storage_size( tDataSet );

            char * tBuffer = ( char * ) malloc( tSize * sizeof( char ) );

            aStatus = H5Dread(
                    tDataSet,
                    tDataType,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    tBuffer );

            aValue.assign( tBuffer, tSize );

            free( tBuffer );

            H5Sclose( tDataSpace );
            H5Tclose( tDataType );
            H5Dclose( tDataSet );

            BELFEM_ASSERT( aStatus == 0,
                      "Something went wrong while trying to load string %s.",
                      aLabel.c_str() );
#endif
        }

//------------------------------------------------------------------------------

        hsize_t
        get_array_size(
                hid_t       & aFileID,
                const std::string & aLabel  )
        {
#ifdef BELFEM_HDF5
            hid_t tDataSet = H5Dopen1( aFileID, aLabel.c_str() );

            hid_t tDataSpace = H5Dget_space( tDataSet );

            hsize_t tDims[ 1 ];

            // this helper only sizes 1-D arrays. Verify the rank BEFORE reading
            // the extents into the single-element tDims buffer: H5Sget_simple_
            // extent_dims writes one entry per rank, so a rank>1 dataset would
            // overflow tDims. H5Sget_simple_extent_ndims returns the rank (>=0)
            // on success and a negative value on failure.
            const int tRank = H5Sget_simple_extent_ndims( tDataSpace );

            BELFEM_ERROR( tRank == 1, "Dataset %s is not a 1-D array (rank %i)",
                aLabel.c_str(), tRank );

            // now safe: exactly one dimension is written into tDims
            H5Sget_simple_extent_dims( tDataSpace, tDims, nullptr );

            H5Dclose( tDataSet );
            H5Sclose( tDataSpace );

            return tDims[ 0 ];
#else
            return 0 ;
#endif
        }

//------------------------------------------------------------------------------

        string
        create_tree( const Cell< string > & aTree, const string & aLabel )
        {
            string aOut = "/" ;

            uint n = aTree.size() ;

            for ( uint k=0; k<n; ++k )
            {
                aOut += aTree( k ) + "/" ;
            }
            aOut += aLabel ;
            return aOut ;
        }

//------------------------------------------------------------------------------
    } /* namespace hdf5 */
} /* namespace belfem */