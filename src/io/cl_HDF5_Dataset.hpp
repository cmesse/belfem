/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_CL_HDF5_DATASET_HPP
#define BELFEM_CL_HDF5_DATASET_HPP

#include "assert.hpp"
#include "hdf5_tools.hpp"

namespace belfem
{
    namespace hdf5
    {

#ifdef BELFEM_HDF5
        template< typename T >
        class Dataset
        {
            const hid_t mType ;
            const bool mWriteMode ;
            const string mLabel ;
            hsize_t mSize;

            hid_t mSpace ;
            hid_t mDataset ;
            hvl_t * mData = nullptr;
        public:

            // aForceWrite selects the write branch even when aSize == 0, so an
            // empty collection ( zero faces/facets ) can still be written. It is
            // additive: the default ( false ) preserves the size-based inference,
            // so Dataset( file, label ) still reads and Dataset( file, label, N )
            // still writes exactly as before.
            Dataset( hid_t aFile, const string & aLabel, const hsize_t aSize=0,
                     const bool aForceWrite=false ) :

                mType( H5Tvlen_create( hdf5::datatype<T>() ) ),
                mWriteMode( aForceWrite || aSize > 0 ),
                mLabel( aLabel ),
                mSize( aSize )
            {
                if ( ! mWriteMode )
                {
                    // read mode: open the existing dataset and read it.
                    // Every error path below closes the handles already opened
                    // before raising: a throwing constructor does not run
                    // ~Dataset(), so partially-built state must be torn down by
                    // hand (debug-only — a release build aborts instead of throwing, see assert.hpp).
                    if ( ! hdf5::dataset_exists( aFile, aLabel ) )
                    {
                        H5Tclose( mType );
                        BELFEM_ERROR( false, "Dataset %s does not exist", aLabel.c_str() );
                    }

                    mDataset = H5Dopen1( aFile, aLabel.c_str() );

                    // guard: the on-disk type must be a vlen of base type T. This
                    // is an I/O-integrity check ( a wrong/corrupt on-disk type
                    // would otherwise be read as silent garbage ), so it stays
                    // active in release as a BELFEM_ERROR. The cost is O(1) type
                    // metadata per open, negligible against the H5Dread itself.
                    //
                    // HDF5 stores standardized file types ( e.g. STD_I32LE ), so an
                    // H5Tequal against the native type would yield false negatives
                    // on every valid read. Compare the vlen base by class + size
                    // ( + sign for integers ) instead, which catches the real
                    // hazards ( 4- vs 8-byte, int-vs-float, signed-vs-unsigned )
                    // without tripping on representation.
                    {
                        hid_t tFileType = H5Dget_type( mDataset );
                        const H5T_class_t tOuterClass = H5Tget_class( tFileType );

                        hid_t tWantBase = hdf5::datatype< T >();
                        const H5T_class_t tWantClass = H5Tget_class( tWantBase );
                        const size_t      tWantBytes = H5Tget_size( tWantBase );
                        const H5T_sign_t  tWantSign  = ( tWantClass == H5T_INTEGER )
                                                       ? H5Tget_sign( tWantBase ) : H5T_SGN_NONE ;

                        // base comparands are only meaningful for a vlen type; if the
                        // on-disk type is not vlen, leave them so the outer check ( which
                        // runs first ) is what reports the mismatch
                        H5T_class_t tBaseClass = H5T_NO_CLASS ;
                        size_t      tBaseBytes = 0 ;
                        H5T_sign_t  tBaseSign  = H5T_SGN_NONE ;

                        if ( tOuterClass == H5T_VLEN )
                        {
                            hid_t tFileBase = H5Tget_super( tFileType );   // vlen element type
                            tBaseClass = H5Tget_class( tFileBase );
                            tBaseBytes = H5Tget_size( tFileBase );
                            tBaseSign  = ( tWantClass == H5T_INTEGER )
                                         ? H5Tget_sign( tFileBase ) : H5T_SGN_NONE ;
                            H5Tclose( tFileBase );
                        }

                        // close before checking so the throw path ( debug ) leaks nothing
                        H5Tclose( tFileType );

                        // on a type mismatch, close the dataset + vlen type before raising
                        auto tFail = [ & ] () { H5Dclose( mDataset ); H5Tclose( mType ); };

                        if ( tOuterClass != H5T_VLEN )
                        {
                            tFail();
                            BELFEM_ERROR( false,
                                "Dataset %s is not a variable-length array on disk", aLabel.c_str() );
                        }
                        if ( tBaseClass != tWantClass )
                        {
                            tFail();
                            BELFEM_ERROR( false,
                                "On-disk element type of dataset %s does not match the requested type", aLabel.c_str() );
                        }
                        if ( tBaseBytes != tWantBytes )
                        {
                            tFail();
                            BELFEM_ERROR( false,
                                "On-disk element size (%lu) of dataset %s does not match requested size (%lu)",
                                ( long unsigned int ) tBaseBytes, aLabel.c_str(), ( long unsigned int ) tWantBytes );
                        }
                        if ( tBaseSign != tWantSign )
                        {
                            tFail();
                            BELFEM_ERROR( false,
                                "On-disk element signedness of dataset %s does not match the requested type", aLabel.c_str() );
                        }
                    }

                    mSpace = H5Dget_space( mDataset );

                    // number of rows = number of points in the dataspace
                    hssize_t tNumPoints = H5Sget_simple_extent_npoints( mSpace );
                    if ( tNumPoints < 0 )
                    {
                        H5Sclose( mSpace );
                        H5Dclose( mDataset );
                        H5Tclose( mType );
                        BELFEM_ERROR( false, "Failed to read extent of dataset %s", aLabel.c_str() );
                    }
                    mSize = static_cast< hsize_t >( tNumPoints );

                    // value-initialize so each row starts as { len=0, p=null }
                    mData = new hvl_t[ mSize ]() ;

                    herr_t tStatus = H5Dread( mDataset, mType, H5S_ALL, H5S_ALL, H5P_DEFAULT, mData );
                    if ( tStatus != 0 )
                    {
                        // failed read: drop the row array. A partial read leaves
                        // indeterminate per-row pointers, so do not vlen-reclaim.
                        delete[] mData;
                        mData = nullptr ;
                        H5Sclose( mSpace );
                        H5Dclose( mDataset );
                        H5Tclose( mType );
                        BELFEM_ERROR( false, "Failed to read dataset %s", aLabel.c_str() );
                    }
                }
                else
                {
                    // write mode: create a new dataset with aSize rows
                    mSpace = H5Screate_simple( 1, &aSize, nullptr  );
                    mDataset = H5Dcreate2(
                        aFile,
                        aLabel.c_str(),
                        mType, mSpace,
                        H5P_DEFAULT,
                        H5P_DEFAULT,
                        H5P_DEFAULT );

                    // value-initialize so unallocated rows are { len=0, p=null }
                    mData = new hvl_t[ aSize ]() ;
                }
            }

            // owns raw buffers + HDF5 handles — non-copyable / non-movable
            Dataset( const Dataset & ) = delete ;
            Dataset & operator=( const Dataset & ) = delete ;
            Dataset( Dataset && ) = delete ;
            Dataset & operator=( Dataset && ) = delete ;

            T *
            set_size( const hsize_t aIndex, const hsize_t aMemory )
            {
                BELFEM_ASSERT( mWriteMode, "set_size() called on read-mode dataset %s", mLabel.c_str() );

                BELFEM_ASSERT( aIndex < mSize , "Index %lu out of range (must be < %lu)",
                    ( long unsigned int ) aIndex,
                    ( long unsigned int ) mSize );

                BELFEM_ASSERT( mData[ aIndex ].len == 0 , "Index %lu of dataset %s already allocated",
                    ( long unsigned int ) aIndex, mLabel.c_str() );

                mData[ aIndex ].p = new T[ aMemory ];
                mData[ aIndex ].len = aMemory;
                return static_cast< T * >( mData[ aIndex ].p );
            }
            T *
            operator[]( const hsize_t aIndex )
            {
                BELFEM_ASSERT( aIndex < mSize , "Index %lu out of range (must be < %lu)",
                    ( long unsigned int ) aIndex,
                    ( long unsigned int ) mSize );
                return static_cast< T * >( mData[ aIndex ].p );
            }

            const T *
            operator[]( const hsize_t aIndex ) const
            {
                BELFEM_ASSERT( aIndex < mSize , "Index %lu out of range (must be < %lu)",
                    ( long unsigned int ) aIndex,
                    ( long unsigned int ) mSize );
                return static_cast< const T * >( mData[ aIndex ].p );
            }


            // number of T entries stored in row aIndex
            hsize_t
            length( const hsize_t aIndex ) const
            {
                BELFEM_ASSERT( aIndex < mSize , "Index %lu out of range (must be < %lu)",
                    ( long unsigned int ) aIndex,
                    ( long unsigned int ) mSize );
                return mData[ aIndex ].len ;
            }

            void
            close()
            {
                if ( mData == nullptr ) return ;

                if ( mWriteMode )
                {
                    // we own the per-row buffers ( new T[] in set_size() )
                    for ( hsize_t k=0; k<mSize; ++k )
                    {
                        delete[] static_cast< T * >( mData[ k ].p );
                    }
                }
                else
                {
                    // per-row buffers were allocated by HDF5 during H5Dread
                    H5Dvlen_reclaim( mType, mSpace, H5P_DEFAULT, mData );
                }

                delete[] mData;
                mData = nullptr ;
                mSize = 0 ;

                H5Dclose( mDataset );
                H5Tclose( mType );
                H5Sclose( mSpace );
            }
            hsize_t
            size() const
            {
                return mSize;
            }

            ~Dataset()
            {
                this->close();
            }

            void
            save()
            {
                BELFEM_ASSERT( mWriteMode, "dataset %s not in write mode", mLabel.c_str() );

                herr_t tStatus = H5Dwrite( mDataset, mType, H5S_ALL, H5S_ALL, H5P_DEFAULT, mData );

                BELFEM_ERROR( tStatus == 0 , "Failed to write dataset %s", mLabel.c_str() );

                this->close();
            }
        };
#endif // BELFEM_HDF5
    }
}
#endif //BELFEM_CL_HDF5_DATASET_HPP
