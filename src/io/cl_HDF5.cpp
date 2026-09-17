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
                    BELFEM_ERROR( file_exists( mPath ),
                               "File %s does not exist.",
                               mPath.c_str() );

                    mFile = H5Fopen(
                            mPath.c_str(),
                            H5F_ACC_RDONLY,   // access flags
                            H5P_DEFAULT);  // Access property list identifier.

                    BELFEM_ERROR( mFile > 0,
                                       "Something went wrong while trying to open file %s.\nIs it in use?",
                               mPath.c_str() );

                    break;
                }
                case( FileMode::OPEN_RDWR ) :
                {
                    BELFEM_ERROR( file_exists( mPath ),
                               "File %s does not exist.",
                               mPath.c_str() );

                    mFile = H5Fopen(
                            mPath.c_str(),
                            H5F_ACC_RDWR,   // access flags
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

            mFileIsOpen = true;

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
        if ( mFileIsOpen )
        {
            this->close_active_group();
            H5Fclose( mFile );
            mFileIsOpen = false;
        }
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
     HDF5::create_group( const string & aLabel )
     {
#ifdef BELFEM_HDF5

        string tLabel = hdf5::create_tree( mTreeLabels, aLabel );

        if( ! hdf5::group_exists(  mActiveGroup, aLabel ) )
        {
            // create the HDF5 group under the target
            mActiveGroup =
                    H5Gcreate2(
                            mActiveGroup,
                            aLabel.c_str(),
                            H5P_DEFAULT,
                            H5P_DEFAULT,
                            H5P_DEFAULT );

            BELFEM_ERROR( mActiveGroup > 0,
                           "Error trying to create group %s in %s", tLabel.c_str(), mPath.c_str() );

            // update tree navigation state
            mActiveGroupLabel = tLabel;
            mTreeLabels.push( aLabel );
            mTree.push( mActiveGroup );

            return mActiveGroup;
        }
        else
        {
            BELFEM_ERROR( false,
                    "Group %s already exists in file %s",
                    tLabel.c_str(),
                    mPath.c_str() );

            return -1;
        }
#else
        return -1;
#endif
    }

//------------------------------------------------------------------------------

    hid_t
    HDF5::create_group( const string & aLabel, const hid_t aParent )
    {
#ifdef BELFEM_HDF5
        // the navigation state is a stack: a group can only be created
        // under the group that is currently active
        BELFEM_ERROR( aParent == mActiveGroup,
            "create_group( %s, parent ): the parent handle must be the active group ( tree is %s )",
            aLabel.c_str(),
            mActiveGroupLabel.c_str() );

        return this->create_group( aLabel );
#else
        return -1;
#endif
    }

//------------------------------------------------------------------------------

    hid_t
    HDF5::select_group(  const string & aLabel )
    {
#ifdef BELFEM_HDF5

        string tLabel = hdf5::create_tree( mTreeLabels, aLabel );

        if ( hdf5::group_exists( mActiveGroup, aLabel ) )
        {

            mActiveGroup = H5Gopen2(
                    mActiveGroup,
                    aLabel.c_str(),
                    H5P_DEFAULT );

            BELFEM_ERROR( mActiveGroup > 0,
                "Error trying to open group %s in %s", tLabel.c_str(), mPath.c_str() );

            mTreeLabels.push( aLabel );
            mTree.push( mActiveGroup );

            mActiveGroupLabel = tLabel;

            return mActiveGroup;
        }
        else
        {
            BELFEM_ERROR( false,
                       "Group %s does not exist in file %s",
                       tLabel.c_str(),
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
        if( mActiveGroup == mFile )
        {
            return ;
        }

        if( mTree.size() > 1 )
        {
            H5Gclose(  mTree.pop() );

            mTreeLabels.pop() ;

            mActiveGroup = mTree.last();

            string tLastLabel = mTreeLabels.pop() ;
            mActiveGroupLabel = hdf5::create_tree( mTreeLabels, tLastLabel );
            mTreeLabels.push( tLastLabel );
        }
        else
        {
            H5Gclose( mActiveGroup );
            mTree.clear() ;
            mTreeLabels.clear() ;
            mActiveGroupLabel = "/";
            mActiveGroup = mFile ;
        }

#endif
    }

//------------------------------------------------------------------------------

    void
    HDF5::close_tree()
    {
#ifdef BELFEM_HDF5
        if( mActiveGroup == mFile )
        {
            return ;
        }

        uint tNumGroups = mTree.size() ;
        for( uint g=0; g<tNumGroups; ++g )
        {
            H5Gclose(  mTree.pop() );
        }

        mTree.clear() ;
        mTreeLabels.clear() ;
        mActiveGroupLabel = "/";
        mActiveGroup = mFile ;

#endif
    }
//------------------------------------------------------------------------------
// Save Strings
//------------------------------------------------------------------------------

    void
    HDF5::save_data(const string &aLabel, const string & aValue)
    {
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
        hdf5::save_scalar_to_file(
                mActiveGroup,
                aLabel,
                aValue,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::save_data(const string & aLabel, const lluint & aValue)
    {
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
        hdf5::load_scalar_from_file(
                mActiveGroup,
                aLabel,
                aValue,
                mStatus );
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    HDF5::load_data( const string & aLabel,
                           lluint & aValue )
    {
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
        hdf5::load_scalar_from_file(
                mActiveGroup,
                aLabel,
                aValue,
                mStatus );
    }

//------------------------------------------------------------------------------

    hid_t
    HDF5::active_group() const
    {
        return mActiveGroup ;
    }

//------------------------------------------------------------------------------

    const string &
    HDF5::tree() const
    {
        return mActiveGroupLabel;
    }

//------------------------------------------------------------------------------

}
