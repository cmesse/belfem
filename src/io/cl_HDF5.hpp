//
// Created by Christian Messe on 2018-12-23.
//

#ifndef BELFEM_CL_HDF5_HPP
#define BELFEM_CL_HDF5_HPP


#include "typedefs.hpp"
#include "filetools.hpp"

#include "HDF5_Tools.hpp"

#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"

//------------------------------------------------------------------------------

namespace belfem
{

//------------------------------------------------------------------------------

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
        bool   mFileIsOpen;

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
         * create a dataset
         */
         hid_t
         create_group( const string & aLabel, const hid_t aParent=-1 );

//------------------------------------------------------------------------------

        /**
         * select a dataset
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
                         real   & aValue );

//------------------------------------------------------------------------------
// Save Vectors
//------------------------------------------------------------------------------

        void
        save_data( const string         & aLabel,
                   const Vector< sint > & aVector );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        save_data( const string         & aLabel,
                   const Vector< uint > & aVector );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        save_data( const string          & aLabel,
                   const Vector< luint > & aVector );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        save_data( const string         & aLabel,
                   const Vector< real > & aVector );

//------------------------------------------------------------------------------
// Load Vectors
//------------------------------------------------------------------------------

        void
        load_data( const string         & aLabel,
                         Vector< sint > & aVector );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        load_data( const string         & aLabel,
                         Vector< uint > & aVector );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        load_data( const string          & aLabel,
                         Vector< luint > & aVector );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        load_data( const string         & aLabel,
                         Vector< real > & aVector );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        load_data(  const string         & aLabel,
                    uint                 * aArray,
                    const index_t          aMemorySize );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        load_data(  const string         & aLabel,
                    real                 * aArray,
                    const index_t          aMemorySize );

//------------------------------------------------------------------------------
// Save Matrices
//------------------------------------------------------------------------------

        void
        save_data( const string         & aLabel,
                   const Matrix< sint > & aMatrix );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        save_data( const string         & aLabel,
                   const Matrix< uint > & aMatrix );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        save_data( const string          & aLabel,
                   const Matrix< luint > & aMatrix );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        save_data( const string         & aLabel,
                   const Matrix< real > & aMatrix );

//------------------------------------------------------------------------------
// Load Matrices
//------------------------------------------------------------------------------

        void
        load_data( const string         & aLabel,
                         Matrix< sint > & aMatrix );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        load_data( const string         & aLabel,
                         Matrix< uint > & aMatrix );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        load_data( const string          & aLabel,
                         Matrix< luint > & aMatrix );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        void
        load_data( const string         & aLabel,
                         Matrix< real > & aMatrix );

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
} /* namespace belfem */

//------------------------------------------------------------------------------
#endif //BELFEM_CL_HDF5_HPP
