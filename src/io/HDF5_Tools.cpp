//
// Created by Christian Messe on 2018-12-25.
//
#include "HDF5_Tools.hpp"

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
            // test if dataset exists
            BELFEM_ERROR( ! hdf5::dataset_exists( aFileID, aLabel ),
                      "Dataset %s of type %s does already exist.",
                      aLabel.c_str(),
                      get_datatype_string( aValue ).c_str());

            // select data type for matrix to save
            hid_t tDataType = H5Tcopy(H5T_NATIVE_HBOOL);

            // matrix dimensions
            hsize_t tDims[1] = {1};

            // create data space
            hid_t tDataSpace
                    = H5Screate_simple(1, tDims, nullptr);

            // create new dataset
            hid_t tDataSet = H5Dcreate(
                    aFileID,
                    aLabel.c_str(),
                    tDataType,
                    tDataSpace,
                    H5P_DEFAULT,
                    H5P_DEFAULT,
                    H5P_DEFAULT);

            // value to cast bool to
            hbool_t tValue = ( hbool_t ) aValue;

            // write data into dataset
            aStatus = H5Dwrite(
                    tDataSet,
                    tDataType,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    & tValue );

            // close open hids
            H5Sclose(tDataSpace);
            H5Tclose(tDataType);
            H5Dclose(tDataSet);

            // check for error
            BELFEM_ASSERT( aStatus == 0,
                      "Something went wrong while trying to store variable %s of type %s",
                      aLabel.c_str(),
                      get_datatype_string(aValue).c_str());
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
            // test if dataset exists
            BELFEM_ERROR( hdf5::dataset_exists( aFileID, aLabel ),
                      "Dataset %s of type %s does not exist.",
                      aLabel.c_str(),
                      get_datatype_string( aValue ).c_str());

            // open the data set
            hid_t tDataSet = H5Dopen1( aFileID, aLabel.c_str() );

            // get the data type of the set
            hid_t tDataType = H5Dget_type( tDataSet );

            // get handler to dataspace
            hid_t tDataSpace = H5Dget_space( tDataSet );

            // value to cast bool to
            hbool_t tValue;

            // read data from file
            aStatus = H5Dread(
                    tDataSet,
                    tDataType,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    &tValue );

            // cast output value
            aValue = ( bool ) tValue;

            // Close/release resources
            H5Tclose( tDataType );
            H5Dclose( tDataSet );
            H5Sclose( tDataSpace );

            // check for error
            BELFEM_ASSERT( aStatus == 0,
                       "Something went wrong while trying to load variable %s of type %s",
                       aLabel.c_str(),
                       get_datatype_string(aValue).c_str());
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
            // test if dataset exists
            BELFEM_ERROR(! hdf5::dataset_exists( aFileID, aLabel ),
                      "Dataset %s of type %s does already exist.",
                      aLabel.c_str(),
                      get_datatype_string( aValue ).c_str());

            // select data type for string
            hid_t tDataType = H5Tcopy( H5T_C_S1 );

            // set size of output type
            aStatus  = H5Tset_size( tDataType, aValue.length() );

            // matrix dimensions
            hsize_t tDims[ 1 ] = { 1 };

            // create data space
            hid_t tDataSpace
                    = H5Screate_simple( 1, tDims, nullptr );

            // create new dataset
            hid_t tDataSet = H5Dcreate(
                    aFileID,
                    aLabel.c_str(),
                    tDataType,
                    tDataSpace,
                    H5P_DEFAULT,
                    H5P_DEFAULT,
                    H5P_DEFAULT );

            // write data into dataset
            aStatus = H5Dwrite(
                    tDataSet,
                    tDataType,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    aValue.c_str() );

            // close open hids
            H5Sclose( tDataSpace );
            H5Tclose( tDataType );
            H5Dclose( tDataSet );

            // check for error
            BELFEM_ASSERT(aStatus == 0,
                      "Something went wrong while trying to store string %s-",
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
            // test if dataset exists
            BELFEM_ERROR( hdf5::dataset_exists( aFileID, aLabel ),
                      "Dataset %s of type %s does not exist.",
                      aLabel.c_str(),
                      get_datatype_string( aValue ).c_str());

            // open the data set
            hid_t tDataSet = H5Dopen1( aFileID, aLabel.c_str() );

            // get handler to dataspace
            hid_t tDataSpace = H5Dget_space( tDataSet );

            // get the data type of the set
            hid_t tDataType = H5Dget_type( tDataSet );

            // get length of string
            hsize_t tSize = H5Dget_storage_size( tDataSet );

            // allocate buffer for string
            char * tBuffer = ( char * ) malloc( tSize * sizeof( char ) );

            // load string into buffer
            aStatus = H5Dread(
                    tDataSet,
                    tDataType,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    tBuffer );

            // create string from buffer
            aValue.assign( tBuffer, tSize );

            // delete buffer
            free( tBuffer );

            // close open hids
            H5Sclose( tDataSpace );
            H5Tclose( tDataType );
            H5Dclose( tDataSet );

            // check for error
            BELFEM_ASSERT( aStatus == 0,
                      "Something went wrong while trying to store string %s.",
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
            // open the data set
            hid_t tDataSet = H5Dopen1( aFileID, aLabel.c_str() );

            // get handler to dataspace
            hid_t tDataSpace = H5Dget_space( tDataSet );

            // matrix dimensions
            hsize_t tDims[ 1 ];

            // ask hdf for dimensions
            herr_t tStatus  = H5Sget_simple_extent_dims( tDataSpace, tDims, nullptr );

            BELFEM_ERROR( tStatus == 0 , "HDF has thrown an error %i", (int) tStatus );

            // close the dataset
            H5Dclose( tDataSet );
            H5Sclose( tDataSpace );

            // return the length
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