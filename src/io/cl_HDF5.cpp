//
// Created by Christian Messe on 2018-12-25.
//

#include "cl_HDF5.hpp"
#include "stringtools.hpp"
#include "assert.hpp"
#include "cl_Communicator.hpp"
#include "commtools.hpp"
//------------------------------------------------------------------------------
namespace belfem
{
//------------------------------------------------------------------------------

    HDF5::HDF5(
            const string & aPath,
            const enum FileMode aMode,
            const bool aParallelMode )
    {
#ifdef BELFEM_HDF5
        // make sure that file path is given
        BELFEM_ERROR( ! aPath.empty(), "No file path given." );

        if( aParallelMode )
        {
            // make path parallel
            mPath =  make_path_parallel( aPath );
        }
        else
        {
            mPath = aPath;
        }
        if( aParallelMode || gComm.rank() == 0 )
        {
            switch( aMode )
            {
                case( FileMode::NEW ) :
                {
                    // call HDF5 API
                    mFile = H5Fcreate(
                            mPath.c_str(),
                            H5F_ACC_TRUNC, // If file exists, erasing all existing data.
                            H5P_DEFAULT,   // File creation property list identifier
                            H5P_DEFAULT);  // Access property list identifier.

                    BELFEM_ERROR( mFile > 0,
                        "Something went wrong while trying to create file %s\nIs it in use?",
                         mPath.c_str() );

                    break;
                }
                case( FileMode::OPEN_RDONLY ) :
                {
                    // make sure that file exists
                    BELFEM_ERROR( file_exists( mPath ),
                               "File %s does not exist.",
                               mPath.c_str() );


                    // call HDF5 API
                    mFile = H5Fopen(
                            mPath.c_str(),
                            H5F_ACC_RDONLY,   // File creation property list identifier
                            H5P_DEFAULT);  // Access property list identifier.

                    BELFEM_ERROR( mFile > 0,
                                       "Something went wrong while trying to open file %s.\nIs it in use?",
                               mPath.c_str() );

                    break;
                }
                case( FileMode::OPEN_RDWR ) :
                {
                    // make sure that file exists
                    BELFEM_ERROR( file_exists( mPath ),
                               "File %s does not exist.",
                               mPath.c_str() );


                    // call HDF5 API
                    mFile = H5Fopen(
                            mPath.c_str(),
                            H5F_ACC_RDWR,   // File creation property list identifier
                            H5P_DEFAULT);  // Access property list identifier.

                    BELFEM_ERROR( mFile > 0,
                               "Something went wrong while trying to open file %s.\nIs it in use?",
                               mPath.c_str() );

                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "unknown filemode passed." );
                    break;
                }
            }

            // set flag for file is open
            mFileIsOpen = true;

            // copy my file id into active group
            mActiveGroup = mFile;

            mActiveGroupLabel = "";
        }
#endif
    }
//------------------------------------------------------------------------------

    /**
     * destructor
     */
    HDF5::~HDF5()
    {
        if( mFileIsOpen )
        {
            this->close();
        }
    }

//------------------------------------------------------------------------------

    void
    HDF5::close()
    {
#ifdef BELFEM_HDF5
        this->close_active_group();
        H5Fclose( mFile );
        mFileIsOpen = false;
#endif
    }

//------------------------------------------------------------------------------

    herr_t &
    HDF5::status()
    {
        return mStatus;
    }

//------------------------------------------------------------------------------

    /**
     * create a group
     */
     hid_t
     HDF5::create_group( const string & aLabel, const hid_t aParent )
     {
#ifdef BELFEM_HDF5

        hid_t tTarget = ( aParent == -1 ) ? ( mFile ) : ( aParent ) ;

        if( ! hdf5::group_exists(  tTarget, aLabel ) )
        {
            // add backslash to label
            std::string tLabel = "/" + aLabel;

            // check if any group is open
            if( mActiveGroup != mFile )
            {
                // close active group
                H5Gclose( mActiveGroup );
            }

            // create group and close it
            mActiveGroup =
                    H5Gcreate2(
                            tTarget,
                            tLabel.c_str(),
                            H5P_DEFAULT,
                            H5P_DEFAULT,
                            H5P_DEFAULT );

            // remember label
            mActiveGroupLabel = aLabel;

            return mActiveGroup;
        }
        else
        {
            BELFEM_ERROR( false,
                    "Group %s already exists in file %s",
                    aLabel.c_str(),
                    mPath.c_str() );

            return -1;
        }
#else
        return -1;
#endif
    }

//------------------------------------------------------------------------------

    hid_t
    HDF5::select_group(  const string & aLabel, const hid_t aParent )
    {
#ifdef BELFEM_HDF5

        hid_t tTarget = ( aParent == -1 ) ? ( mFile ) : ( aParent ) ;

        if ( hdf5::group_exists( tTarget, aLabel ) )
        {
            // add backslash to label
            std::string tLabel = "/" + aLabel;

            // check if any group is open
            if( mActiveGroup != mFile )
            {
                // close active group
                H5Gclose( mActiveGroup );
            }

            mActiveGroup = H5Gopen2(
                    tTarget,
                    tLabel.c_str(),
                    H5P_DEFAULT );

            // remember label
            mActiveGroupLabel = aLabel;

            return mActiveGroup;
        }
        else
        {
            BELFEM_ERROR( false,
                       "Group %s does not exist in file %s",
                       aLabel.c_str(),
                       mPath.c_str() );

            return -1;
        }
#else
        return -1;
#endif
    }


//------------------------------------------------------------------------------

    void
    HDF5::close_active_group()
    {
#ifdef BELFEM_HDF5
        if( mActiveGroup != mFile )
        {
            // close active group
            H5Gclose( mActiveGroup );
        }
#endif
        mActiveGroup = mFile;
        mActiveGroupLabel = "";
    }

//------------------------------------------------------------------------------
// Save Strings
//------------------------------------------------------------------------------

    void
    HDF5::save_data(const string &aLabel, const string & aValue)
    {
        // call interface
        hdf5::save_string_to_file(
                mActiveGroup,
                aLabel,
                aValue,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::save_data( const string &aLabel, const char* aValue )
    {
        // call interface
        hdf5::save_string_to_file(
                mActiveGroup,
                aLabel,
                std::string( aValue ),
                mStatus );
    }

//------------------------------------------------------------------------------
// Load Strings
//------------------------------------------------------------------------------

    void
    HDF5::load_data( const string & aLabel,
                           string & aValue)
    {
        // call interface
        hdf5::load_string_from_file(
            mActiveGroup,
            aLabel,
            aValue,
            mStatus );
    }

//------------------------------------------------------------------------------
// Save Scalars
//------------------------------------------------------------------------------

    void
    HDF5::save_data(const string & aLabel, const sint & aValue)
    {
        // call interface
        hdf5::save_scalar_to_file(
                mActiveGroup,
                aLabel,
                aValue,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::save_data(const string & aLabel, const uint & aValue)
    {
        // call interface
        hdf5::save_scalar_to_file(
                mActiveGroup,
                aLabel,
                aValue,
                mStatus );;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::save_data(const string & aLabel, const luint & aValue)
    {
        // call interface
        hdf5::save_scalar_to_file(
                mActiveGroup,
                aLabel,
                aValue,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::save_data(const std::string & aLabel, const real & aValue)
    {
        // call interface
        hdf5::save_scalar_to_file(
                mActiveGroup,
                aLabel,
                aValue,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::save_data(const std::string & aLabel,
                           const bool & aValue)
    {
        // call interface
        hdf5::save_bool_to_file(
                mActiveGroup,
                aLabel,
                aValue,
                mStatus );
    }

//------------------------------------------------------------------------------
// Load Scalars
//------------------------------------------------------------------------------

    void
    HDF5::load_data( const string & aLabel,
                           bool   & aValue )
    {
        hdf5::load_bool_from_file(
                mActiveGroup,
                aLabel,
                aValue,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::load_data( const string & aLabel,
                             sint & aValue )
    {
        // call interface
        hdf5::load_scalar_from_file(
                mActiveGroup,
                aLabel,
                aValue,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::load_data( const string & aLabel,
                             uint & aValue )
    {
        // call interface
        hdf5::load_scalar_from_file(
                mActiveGroup,
                aLabel,
                aValue,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::load_data( const string & aLabel,
                           luint  & aValue )
    {
        // call interface
        hdf5::load_scalar_from_file(
                mActiveGroup,
                aLabel,
                aValue,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::load_data( const string & aLabel,
                            real  & aValue )
    {
        // call interface
        hdf5::load_scalar_from_file(
                mActiveGroup,
                aLabel,
                aValue,
                mStatus );
    }
//------------------------------------------------------------------------------
// Save Vectors
//------------------------------------------------------------------------------

    void
    HDF5::save_data(   const string         & aLabel,
                       const Vector< sint > & aVector )
    {
        // call interface
        hdf5::save_vector_to_file(
                mActiveGroup,
                aLabel,
                aVector,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::save_data( const string         & aLabel,
                       const Vector< uint > & aVector )
    {
        // call interface
        hdf5::save_vector_to_file(
                mActiveGroup,
                aLabel,
                aVector,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::save_data( const string          & aLabel,
                       const Vector< luint > & aVector )
    {
        // call interface
        hdf5::save_vector_to_file(
                mActiveGroup,
                aLabel,
                aVector,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::save_data( const string         & aLabel,
                     const Vector< real > & aVector )
    {
        // call interface
        hdf5::save_vector_to_file(
                mActiveGroup,
                aLabel,
                aVector,
                mStatus );
    }
//------------------------------------------------------------------------------
// Load Vectors
//------------------------------------------------------------------------------

    void
    HDF5::load_data(   const string         & aLabel,
                             Vector< sint > & aVector )
    {
        // call interface
        hdf5::load_vector_from_file(
                mActiveGroup,
                aLabel,
                aVector,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::load_data(   const string   & aLabel,
                       Vector< uint > & aVector )
    {
        // call interface
        hdf5::load_vector_from_file(
                mActiveGroup,
                aLabel,
                aVector,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::load_data(   const string         & aLabel,
                       Vector< luint > & aVector )
    {
        // call interface
        hdf5::load_vector_from_file(
                mActiveGroup,
                aLabel,
                aVector,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::load_data(   const string   & aLabel,
                       Vector< real > & aVector )
    {
        // call interface
        hdf5::load_vector_from_file(
                mActiveGroup,
                aLabel,
                aVector,
                mStatus );
    }

//------------------------------------------------------------------------------
// Save Matrices
//------------------------------------------------------------------------------

    void
    HDF5::save_data( const string           & aLabel,
                       const Matrix< sint > & aMatrix )
    {
        // call interface
        hdf5::save_matrix_to_file(
                mActiveGroup,
                aLabel,
                aMatrix,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::save_data( const string         & aLabel,
                       const Matrix< uint > & aMatrix )
    {
        // call interface
        hdf5::save_matrix_to_file(
                mActiveGroup,
                aLabel,
                aMatrix,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::save_data( const string          & aLabel,
                       const Matrix< luint > & aMatrix )
    {
        // call interface
        hdf5::save_matrix_to_file(
                mActiveGroup,
                aLabel,
                aMatrix,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::save_data( const string          & aLabel,
                       const Matrix< real >  & aMatrix )
    {
        // call interface
        hdf5::save_matrix_to_file(
                mActiveGroup,
                aLabel,
                aMatrix,
                mStatus );
    }

//------------------------------------------------------------------------------
// Load Matrices
//------------------------------------------------------------------------------

    void
    HDF5::load_data( const string         & aLabel,
                           Matrix< sint > & aMatrix )
    {
        // call interface
        hdf5::load_matrix_from_file(
                mActiveGroup,
                aLabel,
                aMatrix,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::load_data( const string   & aLabel,
                     Matrix< uint > & aMatrix )
    {
        // call interface
        hdf5::load_matrix_from_file(
                mActiveGroup,
                aLabel,
                aMatrix,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::load_data( const string    & aLabel,
                     Matrix< luint > & aMatrix )
    {
        // call interface
        hdf5::load_matrix_from_file(
                mActiveGroup,
                aLabel,
                aMatrix,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::load_data( const string   & aLabel,
                     Matrix< real > & aMatrix )
    {
        // call interface
        hdf5::load_matrix_from_file(
                mActiveGroup,
                aLabel,
                aMatrix,
                mStatus );
    }

//------------------------------------------------------------------------------
}