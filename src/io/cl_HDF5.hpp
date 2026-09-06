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

#ifndef BELFEM_CL_HDF5_HPP
#define BELFEM_CL_HDF5_HPP


#include "typedefs.hpp"
#include "filetools.hpp"

#include "hdf5_tools.hpp"

#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"

//------------------------------------------------------------------------------

namespace belfem
{

//------------------------------------------------------------------------------

    /**
     * @brief Hierarchical data format I/O.
     *
     * @ingroup grp_io
     * @see @ref io_io_usage_guide
     */
    class HDF5
    {
        // path to data file
        string mPath;

        // pointer to file
        hid_t  mFile;

        // pointer to active group
        hid_t  mActiveGroup = -1;

        string mActiveGroupLabel = "";

        // error status of file
        herr_t mStatus = 0;

        // flag telling if file is open
        bool   mFileIsOpen = false;

        Cell< string > mTreeLabels ;
        Cell< hid_t  > mTree ;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        /**
        * constructor
        */
        HDF5(   const string & aPath,
             const enum FileMode aMode=FileMode::NEW,
             const bool aParallelMode=false );

//------------------------------------------------------------------------------

        /**
         * destructor
         */
        ~HDF5();

//------------------------------------------------------------------------------

        /**
         * close file
         */
        void
        close();

//------------------------------------------------------------------------------

        herr_t &
        status();

//------------------------------------------------------------------------------
// GROUP FUNCTIONS
//------------------------------------------------------------------------------

        /**
         * create a group and make it the active group
         */
        hid_t
        create_group( const string & aLabel );

//------------------------------------------------------------------------------

        /**
         * create a group as a child of a given parent group. The navigation
         * state is a stack, so the parent must be the active group — the
         * handle argument makes the intent explicit and is verified
         */
        hid_t
        create_group( const string & aLabel,
                      const hid_t    aParent );

//------------------------------------------------------------------------------

        /**
         * open an existing group and make it the active group
         */
        hid_t
        select_group( const string & aLabel );

//------------------------------------------------------------------------------

        void
        close_active_group();

//------------------------------------------------------------------------------

        void
        close_tree();

//------------------------------------------------------------------------------

        hid_t
        active_group() const;

//------------------------------------------------------------------------------

        const string &
        tree() const ;

//------------------------------------------------------------------------------
// Save Strings
//------------------------------------------------------------------------------

        /**
         * save a string into the active group
         */
        void
        save_data( const string & aLabel,
                   const string & aValue );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        save_data( const string & aLabel,
                   const char *   aValue );

//------------------------------------------------------------------------------
// Load Strings
//------------------------------------------------------------------------------

        void
        load_data( const string & aLabel,
                         string & aValue );

//------------------------------------------------------------------------------
// Save Scalars
//------------------------------------------------------------------------------

        /**
        * save an integer into the active group
        */
        void
        save_data( const string & aLabel,
                   const sint   & aValue );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * save an unsigned integer into the active group
         */
        void
        save_data( const string & aLabel,
                   const uint   & aValue );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * save a long unsigned integer into the active group
         */
        void
        save_data( const string & aLabel,
                   const luint  & aValue );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * save a long long unsigned integer into the active group.
         * Distinct from the luint overload on every platform: uint64_t is
         * `unsigned long` on linux but `unsigned long long` on macOS, and
         * without this overload the macOS spelling matches nothing.
         */
        void
        save_data( const string & aLabel,
                   const lluint & aValue );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * save a real into the active group
         */
        void
        save_data( const string & aLabel,
                   const real   & aValue );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        /**
         * save a bool into the active group
         */
        void
        save_data( const string & aLabel,
                   const bool   & aValue );

//------------------------------------------------------------------------------
// Load Scalars
//------------------------------------------------------------------------------

        void
        load_data( const string & aLabel,
                         bool   & aValue );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        load_data( const string & aLabel,
                         sint   & aValue );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        load_data( const string & aLabel,
                         uint   & aValue );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        load_data( const string & aLabel,
                         luint   & aValue );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        load_data( const string & aLabel,
                         lluint  & aValue );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        load_data( const string & aLabel,
                         real   & aValue );

//------------------------------------------------------------------------------
// Cells
//------------------------------------------------------------------------------

        template< typename T >
        void
        save_data( const string    & aLabel,
                   const Cell< T > & aCell )
        {
            hdf5::save_array_to_file(
                mActiveGroup,
                aLabel,
                aCell.data(),
                aCell.size(),
                mStatus );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template< typename T >
        void
        load_data( const string    & aLabel,
                         Cell< T > & aCell )
        {
            hsize_t tSize = hdf5::get_array_size( mActiveGroup, aLabel );

            aCell.set_size( tSize );

            hdf5::load_array_from_file(
                mActiveGroup,
                aLabel,
                aCell.data(),
                aCell.size(),
                mStatus );
        }

//------------------------------------------------------------------------------
// Vectors
//------------------------------------------------------------------------------

        template< typename T >
        void
        save_data( const string      & aLabel,
                   const Vector< T > & aVector )
        {
            // call interface
            hdf5::save_vector_to_file(
                    mActiveGroup,
                    aLabel,
                    aVector,
                    mStatus );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template< typename T >
        void
        load_data( const string      & aLabel,
                         Vector< T > & aVector )
        {
            hdf5::load_vector_from_file(
                    mActiveGroup,
                    aLabel,
                    aVector,
                    mStatus );
        }

//------------------------------------------------------------------------------
//  Matrices
//------------------------------------------------------------------------------

        template< typename T >
        void
        save_data( const string      & aLabel,
                   const Matrix< T > & aMatrix,
                   const bool            aTranspose = false )
        {
            // call interface
            hdf5::save_matrix_to_file(
                    mActiveGroup,
                    aLabel,
                    aMatrix,
                    mStatus,
                    aTranspose );
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template< typename T >
        void
        load_data( const string      & aLabel,
                         Matrix< T > & aMatrix,
                         const bool    aTranspose = false )
        {
            // call interface
            hdf5::load_matrix_from_file(
                    mActiveGroup,
                    aLabel,
                    aMatrix,
                    mStatus,
                    aTranspose );
        }

//------------------------------------------------------------------------------
//  Raw Arrays
//------------------------------------------------------------------------------

        template< typename T >
        void
        save_data(
            const string  & aLabel,
            const T       * aArray,
            const index_t   aMemorySize )
        {
            hdf5::save_array_to_file(
                mActiveGroup,
                aLabel,
                aArray,
                aMemorySize,
                mStatus );
        }

        template< typename T >
        void
        load_data(
            const string  & aLabel,
            T             * aArray,
            const index_t   aMemorySize ) // <-- must be known a priori
        {
            hdf5::load_array_from_file(
                mActiveGroup,
                aLabel,
                aArray,
                aMemorySize,
                mStatus );
        }

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
// Cells of strings
//------------------------------------------------------------------------------

    template<>
    inline void
    HDF5::save_data( const string         & aLabel,
                     const Cell< string > & aCell )
    {
#ifdef BELFEM_HDF5
        hdf5::save_strings_to_file( mActiveGroup, aLabel, aCell, mStatus );
#endif
    }

    template<>
    inline void
    HDF5::load_data(
        const string         & aLabel,
        Cell< string > & aCell )
    {
#ifdef BELFEM_HDF5
      hdf5::load_strings_from_file( mActiveGroup, aLabel, aCell, mStatus );
#endif
    }

//------------------------------------------------------------------------------
} /* namespace belfem */

//------------------------------------------------------------------------------
#endif //BELFEM_CL_HDF5_HPP
