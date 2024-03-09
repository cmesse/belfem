//
// Created by Christian Messe on 2018-12-23.
//

#ifndef BELFEM_HDF5_TOOLS_HPP
#define BELFEM_HDF5_TOOLS_HPP

#ifdef BELFEM_HDF5
#include <cstring>

#include <hdf5.h>
#else
namespace belfem
{
    typedef int hid_t;
    typedef int herr_t;
    typedef int hsize_t;
}
#endif

#include "stringtools.hpp"
#include "assert.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"

namespace belfem
{
    namespace hdf5
    {
//------------------------------------------------------------------------------

        /**
         * test if a dataset exists
         */
        inline bool
        dataset_exists( hid_t aFileID, const std::string & aLabel )
        {
#ifdef BELFEM_HDF5
            hid_t tDataSet = 0;
            return H5Lexists( aFileID, aLabel.c_str(), tDataSet );
#else
            return false;
#endif
        }

//------------------------------------------------------------------------------

        /**
         * test if a group exists
         */
        inline bool
        group_exists( hid_t aFileID, const std::string & aLabel )
        {
#ifdef BELFEM_HDF5
            std::string tLabel = "/" + aLabel;
            hid_t tGroup = 0;
            return H5Lexists( aFileID, tLabel.c_str(), tGroup );
#else
            return false;
#endif
        }


//------------------------------------------------------------------------------

        /**
         *
         * @brief                 returns a HDF5 enum defining the
         *                        data type that is to be communicated.
         *
         * @param[in] aSample     primitive data type with arbitrary value
         *
         * see also https://support.hdfgroup.org/HDF5/doc/H5.user/Datatypes.html
         */
        template < typename T >
        inline hid_t
        get_datatype( const T & aSample )
        {
            BELFEM_ERROR( false , "get_hdf5_datatype: unknown data type.");
            return 0;
        }

#ifdef BELFEM_HDF5
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template <>
        inline hid_t
        get_datatype( const int & aSample )
        {
            return H5T_NATIVE_INT;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template <>
        inline hid_t
        get_datatype( const long int & aSample )
        {
            return H5T_NATIVE_LONG;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template <>
        inline hid_t
        get_datatype( const unsigned int & aSample )
        {
            return H5T_NATIVE_UINT;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template <>
        inline hid_t
        get_datatype( const long unsigned int & aSample )
        {
            return H5T_NATIVE_ULONG;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template <>
        inline hid_t
        get_datatype( const double & aSample )
        {
            return H5T_NATIVE_DOUBLE;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template <>
        inline hid_t
        get_datatype( const long double & aSample )
        {
            return H5T_NATIVE_LDOUBLE;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<>
        inline hid_t
        get_datatype( const bool & aSample )
        {
            return H5T_NATIVE_HBOOL;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif /* BELFEM_HDF5 */

//------------------------------------------------------------------------------

        /**
         * this function returns true of both the HDF5 datatype
         * and the passed datatype have the same size
         */
        template < typename T >
        bool
        test_size_of_datatype( const T & aSample )
        {
#ifdef BELFEM_HDF5
            return     H5Tget_size( hdf5::get_datatype( aSample ) )
                       == sizeof( T );
#else
            return false;
#endif
        }

//------------------------------------------------------------------------------

        /**
        * saves a scalar value to a file
                * file must be open
        *
        * @param[ inout ] aFileID  handler to hdf5 file
        * @param[ in ]    aLabel   label of matrix to save
        * @param[ in ]    aValue   value that is to be stored
        * @param[ in ]    aStatus  error handler
        */
        template < typename T >
        void
        save_scalar_to_file(
                hid_t               & aFileID,
                const std::string   & aLabel,
                const T             & aValue,
                herr_t              & aStatus
        )
        {
#ifdef BELFEM_HDF5
            // test if dataset exists
            BELFEM_ERROR( ! hdf5::dataset_exists( aFileID, aLabel ),
                    "Dataset %s of type %s does already exist.",
                    aLabel.c_str(),
                    get_datatype_string( aValue ).c_str() );

            // check datatype
            BELFEM_ASSERT( hdf5::test_size_of_datatype( aValue ),
                    "Error in datatype size of type %s.",
                        get_datatype_string( aValue ).c_str() );

            // select datatype for data to save
            hid_t tDataType = H5Tcopy( hdf5::get_datatype( aValue ) );

            // set data type to little endian
            aStatus = H5Tset_order( tDataType, H5T_ORDER_LE );

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
                    & aValue );

            // close open hids
            H5Sclose( tDataSpace );
            H5Tclose( tDataType );
            H5Dclose( tDataSet );

            // check for error
            BELFEM_ASSERT( aStatus == 0,
                    "Something went wrong while trying to store scalar %s of type %s",
                       aLabel.c_str(),
                       get_datatype_string( aValue ).c_str() );
#endif
        }
//------------------------------------------------------------------------------

        template < typename T >
        void
        load_scalar_from_file(
                hid_t               & aFileID,
                const std::string   & aLabel,
                T                   & aValue,
                herr_t              & aStatus
        )
        {
#ifdef BELFEM_HDF5
            // test if dataset exists
            BELFEM_ERROR( hdf5::dataset_exists( aFileID, aLabel ),
                       "Dataset %s of type %s does not exist.",
                       aLabel.c_str(),
                       get_datatype_string( aValue ).c_str() );

            // check datatype
            BELFEM_ASSERT( hdf5::test_size_of_datatype( aValue ),
                        "Error in datatype size of type %s.",
                        get_datatype_string( aValue ).c_str() );

            // open the data set
            hid_t tDataSet = H5Dopen1( aFileID, aLabel.c_str() );

            // get the data type of the set
            hid_t tDataType = H5Dget_type( tDataSet );

            // get handler to dataspace
            hid_t tDataSpace = H5Dget_space( tDataSet );

            // read data from file
            aStatus = H5Dread(
                    tDataSet,
                    tDataType,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    &aValue );

            // Close/release resources
            H5Tclose( tDataType );
            H5Dclose( tDataSet );
            H5Sclose( tDataSpace );

            // check for error
            BELFEM_ASSERT( aStatus == 0,
                       "Something went wrong while trying to load scalar %s of type %s",
                       aLabel.c_str(),
                       get_datatype_string( aValue ).c_str() );

#endif
        }

//------------------------------------------------------------------------------

        void
        save_bool_to_file(
                hid_t               & aFileID,
                const std::string   & aLabel,
                const bool          & aValue,
                herr_t              & aStatus );

//------------------------------------------------------------------------------

        void
        load_bool_from_file(
                            hid_t   & aFileID,
                const std::string   & aLabel,
                             bool   & aValue,
                           herr_t   & aStatus );

//------------------------------------------------------------------------------

        void
        save_string_to_file(
                hid_t               & aFileID,
                const std::string   & aLabel,
                const std::string   & aValue,
                herr_t              & aStatus );

//------------------------------------------------------------------------------

        void
        load_string_from_file(
                      hid_t         & aFileID,
                const std::string   & aLabel,
                      std::string   & aValue,
                      herr_t        & aStatus );

//------------------------------------------------------------------------------
        template < typename T >
        void
        save_array_to_file(
                       hid_t      & aFileID,
                const std::string & aLabel,
                const T           * aData,
                const hsize_t     & aLength,
                      herr_t      & aStatus )
        {
#ifdef BELFEM_HDF5
            // create a sample
            T tSample = 0;

            // test if dataset exists
            BELFEM_ERROR( ! hdf5::dataset_exists( aFileID, aLabel ),
                       "Dataset %s of type Vector<%s> does already exist.",
                       aLabel.c_str(),
                       get_datatype_string( tSample ).c_str() );

            // check datatype
            BELFEM_ASSERT( hdf5::test_size_of_datatype( tSample ),
                        "Error in datatype size of type %s.",
                        get_datatype_string( tSample ).c_str() );

            // select datatype for data to save
            hid_t tDataType = H5Tcopy( hdf5::get_datatype( tSample ) );

            // set data type to little endian
            aStatus = H5Tset_order( tDataType, H5T_ORDER_LE );

            // matrix dimensions
            hsize_t tDims[ 1 ];
            tDims[ 0 ] = aLength;

            // create data space
            hid_t  tDataSpace
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

            // test if vector is not empty
            if( aLength > 0 )
            {
                // write data into dataset
                aStatus = H5Dwrite(
                        tDataSet,
                        tDataType,
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,
                        &aData[ 0 ] );
            }

            // close open hids
            H5Sclose( tDataSpace );
            H5Tclose( tDataType );
            H5Dclose( tDataSet );

            // check for error
            BELFEM_ASSERT( aStatus == 0,
                        "Something went wrong while trying to store vector %s of type %s",
                        aLabel.c_str(),
                        get_datatype_string( tSample ).c_str() );
#endif
        }

//------------------------------------------------------------------------------

        hsize_t
        get_array_size(
        hid_t       & aFileID,
        const std::string & aLabel  );

//------------------------------------------------------------------------------

        template < typename T >
        void
        load_array_from_file(
                      hid_t       & aFileID,
                const std::string & aLabel,
                      T           * aData,
                const hsize_t       aLength,
                      herr_t      & aStatus )
        {
#ifdef BELFEM_HDF5
            // create a sample
            T tSample = 0;

            // test if dataset exists
            BELFEM_ERROR( hdf5::dataset_exists( aFileID, aLabel ),
                       "Dataset %s of type Vector<%s> does not exist.",
                       aLabel.c_str(),
                       get_datatype_string( tSample ).c_str() );

            // check datatype
            BELFEM_ASSERT( hdf5::test_size_of_datatype( tSample ),
                        "Error in datatype size of type %s.",
                        get_datatype_string( tSample ).c_str() );

            // open the data set
            hid_t tDataSet = H5Dopen1( aFileID, aLabel.c_str() );

            // get the data type of the set
            hid_t tDataType = H5Dget_type( tDataSet );

            // get handler to dataspace
            hid_t tDataSpace = H5Dget_space( tDataSet );

            // matrix dimensions
            hsize_t tDims[ 1 ];

            // ask hdf for dimensions
            aStatus  = H5Sget_simple_extent_dims( tDataSpace, tDims, nullptr );

            // get length
            BELFEM_ERROR( tDims[ 0 ] == aLength, "Lengths do not match: is %lu, expect %lu.",
                         ( long unsigned int ) tDims[ 0 ], ( long unsigned int ) aLength );

            // test if vector is empty
            if( aLength > 0 )
            {
                // read data from file
                aStatus = H5Dread(
                        tDataSet,
                        tDataType,
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,
                        &aData[ 0 ] );
            }
            else
            {
                // all good. reset status
                aStatus = 0;
            }

            // Close/release resources
            H5Tclose( tDataType );
            H5Dclose( tDataSet );
            H5Sclose( tDataSpace );

            // check for error
            BELFEM_ASSERT( aStatus == 0,
                        "Something went wrong while trying to load vector %s of type %s",
                        aLabel.c_str(),
                        get_datatype_string( tSample ).c_str() );
#endif
        }
//------------------------------------------------------------------------------

        template < typename T >
        void
        save_vector_to_file(
                hid_t         & aFileID,
                const std::string   & aLabel,
                const Vector< T >   & aVector,
                herr_t        & aStatus )
        {
#ifdef BELFEM_HDF5
            // create a sample
            T tSample = 0;

            // test if dataset exists
            BELFEM_ERROR( ! hdf5::dataset_exists( aFileID, aLabel ),
                       "Dataset %s of type Vector<%s> does already exist.",
                       aLabel.c_str(),
                       get_datatype_string( tSample ).c_str() );

            // check datatype
            BELFEM_ASSERT( hdf5::test_size_of_datatype( tSample ),
                        "Error in datatype size of type %s.",
                        get_datatype_string( tSample ).c_str() );

            // select datatype for data to save
            hid_t tDataType = H5Tcopy( hdf5::get_datatype( tSample ) );

            // set data type to little endian
            aStatus = H5Tset_order( tDataType, H5T_ORDER_LE );
            hid_t tDataSet = 0;

            // matrix dimensions
            hsize_t tLength  = aVector.length();
            hsize_t tDims[ 1 ];
            tDims[ 0 ] = tLength;

            // create data space
            hid_t  tDataSpace
                    = H5Screate_simple( 1, tDims, nullptr );

            // create new dataset
            tDataSet = H5Dcreate(
                    aFileID,
                    aLabel.c_str(),
                    tDataType,
                    tDataSpace,
                    H5P_DEFAULT,
                    H5P_DEFAULT,
                    H5P_DEFAULT );

            // test if vector is not empty
            if( tLength > 0 )
            {
                // allocate memory for data
                T* tData = ( T* ) malloc( tLength * sizeof( T ) );

                // copy vector to data
                std::memcpy(
                        tData,
                        aVector.data(),
                        tLength * sizeof( T ) );

                // write data into dataset
                aStatus = H5Dwrite(
                        tDataSet,
                        tDataType,
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,
                        & tData[ 0 ]);

                // tidy up memory
                free( tData );
            }

            // close open hids
            H5Sclose( tDataSpace );
            H5Tclose( tDataType );
            H5Dclose( tDataSet );

            // check for error
            BELFEM_ASSERT( aStatus == 0,
                        "Something went wrong while trying to store vector %s of type %s",
                        aLabel.c_str(),
                        get_datatype_string( tSample ).c_str() );
#endif
        }

//------------------------------------------------------------------------------

        template < typename T >
        void
        load_vector_from_file(
                hid_t         & aFileID,
                const std::string   & aLabel,
                Vector< T >   & aVector,
                herr_t        & aStatus )
        {
#ifdef BELFEM_HDF5
            // create a sample
            T tSample = 0;

            // test if dataset exists
            BELFEM_ERROR( hdf5::dataset_exists( aFileID, aLabel ),
                       "Dataset %s of type Vector<%s> does not exist.",
                       aLabel.c_str(),
                       get_datatype_string( tSample ).c_str() );

            // check datatype
            BELFEM_ASSERT( hdf5::test_size_of_datatype( tSample ),
                        "Error in datatype size of type %s.",
                        get_datatype_string( tSample ).c_str() );

            // open the data set
            hid_t tDataSet = H5Dopen1( aFileID, aLabel.c_str() );

            // get the data type of the set
            hid_t tDataType = H5Dget_type( tDataSet );

            // get handler to dataspace
            hid_t tDataSpace = H5Dget_space( tDataSet );

            // matrix dimensions
            hsize_t tDims[ 1 ];

            // ask hdf for dimensions
            aStatus  = H5Sget_simple_extent_dims( tDataSpace, tDims, nullptr );

            // get length
            hsize_t tLength = tDims[ 0 ];

            // allocate length of output vector
            aVector.set_size( tLength );

            // test if vector is empty
            if( tLength > 0 )
            {
                // allocate buffer
                T* tData = ( T* ) malloc( tLength * sizeof( T ) );

                // read data from file
                aStatus = H5Dread(
                        tDataSet,
                        tDataType,
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,
                        &tData[ 0 ] );

                // copy data to vector
                std::memcpy(
                        aVector.data(),
                        tData,
                        tLength * sizeof( T ) );

                // tidy up memory
                free( tData );
            }
            else
            {
                // all good. reset status
                aStatus = 0;
            }

            // Close/release resources
            H5Tclose( tDataType );
            H5Dclose( tDataSet );
            H5Sclose( tDataSpace );

            // check for error
            BELFEM_ASSERT( aStatus == 0,
                        "Something went wrong while trying to load vector %s of type %s",
                        aLabel.c_str(),
                        get_datatype_string( tSample ).c_str() );
#endif
        }
//------------------------------------------------------------------------------

        template < typename T >
        void
        save_matrix_to_file(
                hid_t               & aFileID,
                const std::string   & aLabel,
                const Matrix< T >   & aMatrix,
                herr_t              & aStatus )
        {
#ifdef BELFEM_HDF5
            // create a sample
            T tSample = 0;

            // test if dataset exists
            BELFEM_ERROR( ! hdf5::dataset_exists( aFileID, aLabel ),
                       "Dataset %s of type Matrix<%s> does already exist.",
                       aLabel.c_str(),
                       get_datatype_string( tSample ).c_str() );

            // check datatype
            BELFEM_ASSERT( hdf5::test_size_of_datatype( tSample ),
                        "Error in datatype size of type %s.",
                        get_datatype_string( tSample ).c_str() );

            // select datatype for data to save
            hid_t tDataType = H5Tcopy( hdf5::get_datatype( tSample ) );

            // set data type to little endian
            aStatus = H5Tset_order( tDataType, H5T_ORDER_LE );
            hid_t tDataSet = 0;

            // matrix dimensions
            hsize_t tDims[ 2 ];
            tDims[ 0 ] = aMatrix.n_rows();
            tDims[ 1 ] = aMatrix.n_cols();

            // create data space
            hid_t  tDataSpace
                    = H5Screate_simple( 2, tDims, nullptr );

            // create new dataset
            tDataSet = H5Dcreate(
                    aFileID,
                    aLabel.c_str(),
                    tDataType,
                    tDataSpace,
                    H5P_DEFAULT,
                    H5P_DEFAULT,
                    H5P_DEFAULT );

            // test if vector is not empty
            if( tDims[ 0 ]*tDims[ 1 ] > 0 )
            {
                // allocate memory for data
                T** tData = ( T** ) malloc( tDims[ 0 ] * sizeof( T * ) );
                tData[ 0 ] = ( T* ) malloc( tDims[ 0 ]*tDims[ 1 ] * sizeof( T ) );

                // loop over all rows and allocate colums
                for( hsize_t i=0; i<tDims[ 0 ]; ++i )
                {
                    tData[ i ] = tData[ 0 ]+ i*tDims[ 1 ];
                }

                // convert matrix to array
                for ( hsize_t i = 0; i < tDims[ 0 ]; ++i )
                {
                    for ( hsize_t j = 0; j < tDims[ 1 ]; ++j )
                    {
                        tData[ i ][ j ] = aMatrix( i, j );
                    }
                }

                // write data into dataset
                aStatus = H5Dwrite(
                        tDataSet,
                        tDataType,
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,
                        & tData[ 0 ][ 0 ]);

                // tidy up memory
                free( tData[ 0 ] );
                free( tData );
            }

            // close open hids
            H5Sclose( tDataSpace );
            H5Tclose( tDataType );
            H5Dclose( tDataSet );

            // check for error
            BELFEM_ASSERT( aStatus == 0,
                       "Something went wrong while trying to store matrix %s of type %s",
                       aLabel.c_str(),
                       get_datatype_string( tSample ).c_str() );
#endif
        }

//------------------------------------------------------------------------------

        template < typename T >
        void
        load_matrix_from_file(
                hid_t               & aFileID,
                const std::string   & aLabel,
                      Matrix< T >   & aMatrix,
                herr_t              & aStatus )
        {
#ifdef BELFEM_HDF5
            // create a sample
            T tSample = 0;

            // test if dataset exists
            BELFEM_ERROR( hdf5::dataset_exists( aFileID, aLabel ),
                       "Dataset %s of type Matrix<%s> does not exist.",
                       aLabel.c_str(),
                       get_datatype_string( tSample ).c_str() );

            // check datatype
            BELFEM_ASSERT( hdf5::test_size_of_datatype( tSample ),
                        "Error in datatype size of type %s.",
                        get_datatype_string( tSample ).c_str() );


            // open the data set
            hid_t tDataSet = H5Dopen1( aFileID, aLabel.c_str() );

            // get the data type of the set
            hid_t tDataType = H5Dget_type( tDataSet );

            // test datatype compatibility
            BELFEM_ASSERT(   H5Tget_class( tDataType )
                          ==  H5Tget_class( hdf5::get_datatype( tSample ) ),
                        "Dataset %s does not seem to be %s.",
                        aLabel.c_str(),
                        get_datatype_string( tSample ).c_str() );

            // get handler to dataspace
            hid_t tDataSpace = H5Dget_space( tDataSet );

            // matrix dimensions
            hsize_t tDims[ 2 ];

            // ask hdf for dimensions
            aStatus  = H5Sget_simple_extent_dims( tDataSpace, tDims, nullptr );

            // allocate memory for output
            aMatrix.set_size( tDims[ 0 ], tDims[ 1 ] );

            // test if matrix is not empty
            if( tDims[ 0 ]*tDims[ 1 ] > 0 )
            {
                // allocate top level array which contains rows
                T** tData = ( T** )
                        malloc( tDims[ 0 ]*sizeof( T* ) );

                // allocate memory for data
                tData[ 0 ] = ( T* )
                        malloc( tDims[ 0 ]*  tDims[ 1 ] * sizeof( T ) );

                // loop over all rows and allocate colums
                for( hsize_t i=1; i<tDims[ 0 ]; ++i )
                {
                    tData[ i ] = tData[ 0 ]+ i*tDims[ 1 ];
                }


                // read data from file
                aStatus = H5Dread(
                        tDataSet,
                        tDataType,
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,
                        &tData[ 0 ][ 0 ] );

                // write values into matrix
                for ( hsize_t j = 0; j < tDims[ 1 ]; ++j )
                {
                    for ( hsize_t i = 0; i < tDims[ 0 ]; ++i )
                    {
                        aMatrix( i, j ) = tData[ i ][ j ];
                    }
                }

                // tidy up memory
                free( tData[ 0 ] );
                free( tData );
            }
            else if( aStatus == 2 )
            {
                // all good, reset status
                aStatus = 0;
            }

            // Close/release resources
            H5Tclose( tDataType );
            H5Dclose( tDataSet );
            H5Sclose( tDataSpace );

            // check for error
            BELFEM_ASSERT( aStatus == 0,
                        "Something went wrong while trying to load matrix %s of type %s",
                        aLabel.c_str(),
                        get_datatype_string( tSample ).c_str() );
#endif
        }

//------------------------------------------------------------------------------

        inline herr_t
        open_group( const std::string & aLabel,  hid_t aParent )
        {
#ifdef BELFEM_HDF5
            if ( hdf5::group_exists( aParent, aLabel ) )
            {
                // add backslash to label
                std::string tLabel = "/" + aLabel;

                return H5Gopen2(
                        aParent,
                        tLabel.c_str(),
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
    } /* namespace hdf5 */
} /* namespace belfem */

#endif //BELFEM_HDF5_TOOLS_HPP
