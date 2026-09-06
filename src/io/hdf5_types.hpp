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

#ifndef BELFEM_HDF5_TYPES_HPP
#define BELFEM_HDF5_TYPES_HPP

#ifdef BELFEM_HDF5
#include <hdf5.h>
#else
namespace belfem
{
    typedef int hid_t;
    typedef int herr_t;
    typedef int hsize_t;
}
#endif

#include <cstddef>
#include <type_traits>

#include "assert.hpp"

namespace belfem
{
    namespace hdf5
    {
        template<typename T>
        hid_t datatype()
        {
#ifdef BELFEM_HDF5
            BELFEM_ERROR( false, "Unsupported datatype");
            return 0;
#else
            BELFEM_ERROR( false, "We are not linked against HDF5");
            return 0 ;
#endif
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// datatype<T>() above is the MEMORY type: the native C type, with whatever
// width and byte order this build happens to use. filetype<T>() below is the
// ON-DISK type, pinned to an explicit width and to little endian, so that the
// layout of a dataset is a property of the file format and not of the build
// configuration that produced it.
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<typename T>
        hid_t filetype()
        {
#ifdef BELFEM_HDF5
            BELFEM_ERROR( false, "Unsupported datatype");
            return 0;
#else
            BELFEM_ERROR( false, "We are not linked against HDF5");
            return 0 ;
#endif
        }
#ifdef BELFEM_HDF5
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// chars
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        template<> inline
        hid_t datatype<char>()
        {
            return H5T_NATIVE_CHAR;
        }

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<> inline
        hid_t datatype<signed char>()
        {
            return H5T_NATIVE_SCHAR;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<> inline
        hid_t datatype<unsigned char>()
        {
            return H5T_NATIVE_UCHAR;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// signed integers
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<> inline
        hid_t datatype<short>()
        {
            return H5T_NATIVE_SHORT;
        }

        template<> inline
        hid_t datatype<int>()
        {
            return H5T_NATIVE_INT;
        }

        template<> inline
        hid_t datatype<long>()
        {
            return H5T_NATIVE_LONG;
        }

        template<> inline
        hid_t datatype<long long>()
        {
            return H5T_NATIVE_LLONG;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// unsigned integers
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<> inline
        hid_t datatype<unsigned short>()
        {
            return H5T_NATIVE_USHORT;
        }

        template<> inline
        hid_t datatype<unsigned int>()
        {
            return H5T_NATIVE_UINT;
        }

        template<> inline
        hid_t datatype<unsigned long>()
        {
            return H5T_NATIVE_ULONG;
        }

        template<> inline
        hid_t datatype<unsigned long long>()
        {
            return H5T_NATIVE_ULLONG;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// floating point
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<> inline
        hid_t datatype<float>()
        {
            return H5T_NATIVE_FLOAT;
        }

        template<> inline
        hid_t datatype<double>()
        {
            return H5T_NATIVE_DOUBLE;
        }

        template<> inline
        hid_t datatype<long double>()
        {
            return H5T_NATIVE_LDOUBLE;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// bool
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<> inline
        hid_t datatype<bool>()
        {
            return H5T_NATIVE_HBOOL;
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// fixed-width file types
//
// The integer cases dispatch on size and signedness rather than on the C type
// name, because the name-to-width mapping is platform dependent: `long` is
// 8 bytes under LP64 and 4 under LLP64, and plain `char` is signed on x86 but
// unsigned on ARM. Dispatching on the property reproduces the native type
// exactly on every platform; dispatching on the name would not.
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        inline hid_t
        int_filetype( const std::size_t aSize, const bool aSigned )
        {
            switch ( aSize )
            {
                case 1 : return aSigned ? H5T_STD_I8LE  : H5T_STD_U8LE  ;
                case 2 : return aSigned ? H5T_STD_I16LE : H5T_STD_U16LE ;
                case 4 : return aSigned ? H5T_STD_I32LE : H5T_STD_U32LE ;
                case 8 : return aSigned ? H5T_STD_I64LE : H5T_STD_U64LE ;
                default :
                {
                    BELFEM_ERROR( false,
                            "No fixed-width HDF5 file type for a %d-byte integer",
                            ( int ) aSize );
                    return 0 ;
                }
            }
        }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// chars
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<> inline hid_t filetype<char>()
        { return int_filetype( sizeof( char ), std::is_signed< char >::value ); }

        template<> inline hid_t filetype<signed char>()
        { return int_filetype( sizeof( signed char ), true ); }

        template<> inline hid_t filetype<unsigned char>()
        { return int_filetype( sizeof( unsigned char ), false ); }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// signed integers
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<> inline hid_t filetype<short>()
        { return int_filetype( sizeof( short ), true ); }

        template<> inline hid_t filetype<int>()
        { return int_filetype( sizeof( int ), true ); }

        template<> inline hid_t filetype<long>()
        { return int_filetype( sizeof( long ), true ); }

        template<> inline hid_t filetype<long long>()
        { return int_filetype( sizeof( long long ), true ); }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// unsigned integers
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<> inline hid_t filetype<unsigned short>()
        { return int_filetype( sizeof( unsigned short ), false ); }

        template<> inline hid_t filetype<unsigned int>()
        { return int_filetype( sizeof( unsigned int ), false ); }

        template<> inline hid_t filetype<unsigned long>()
        { return int_filetype( sizeof( unsigned long ), false ); }

        template<> inline hid_t filetype<unsigned long long>()
        { return int_filetype( sizeof( unsigned long long ), false ); }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// floating point
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        template<> inline hid_t filetype<float>()
        { return H5T_IEEE_F32LE; }

        template<> inline hid_t filetype<double>()
        { return H5T_IEEE_F64LE; }

        // long double is the one type with no fixed-width IEEE equivalent: on
        // x86 it is 16 bytes carrying 80 bits of precision, and H5T_IEEE_F64LE
        // would silently truncate it. It therefore keeps the native layout.
        // No writer instantiates it today ( `real` is `double` ).
        template<> inline hid_t filetype<long double>()
        { return H5T_NATIVE_LDOUBLE; }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// bool
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // written through an hbool_t, so the file width follows that type
        template<> inline hid_t filetype<bool>()
        { return int_filetype( sizeof( hbool_t ), false ); }

#endif
    }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}
#endif // BELFEM_HDF5_TYPES_HPP
