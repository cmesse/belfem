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

#ifndef COMMTYPES_HPP
#define COMMTYPES_HPP

#include <complex>

#ifdef BELFEM_MPI
#include <mpi.h>
#endif

namespace belfem
{
#ifdef BELFEM_MPI
    typedef MPI_Datatype comm_t ;
#else
    // define fake type if there is no MPI
    typedef int comm_t;
#endif
    typedef int proc_t ;

//------------------------------------------------------------------------------

    /**
     * \brief returns the MPI datatype handle for T
     */
    template<typename T> comm_t
    comm_type()
    {
#ifdef BELFEM_MPI
        BELFEM_ERROR( false, "Unknown datatype");
        return MPI_DATATYPE_NULL;
#else
        return 0 ;
#endif
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef BELFEM_MPI
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline comm_t
    comm_type< char > ()
    {
        return MPI_CHAR;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline comm_t
    comm_type< int > ()
    {
        return MPI_INT;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline comm_t
    comm_type< long int > ()
    {
        return MPI_LONG;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline comm_t
    comm_type< short unsigned int > ()
    {
        return MPI_UNSIGNED_SHORT;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline comm_t
    comm_type< unsigned int > ()
    {
        return MPI_UNSIGNED;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline comm_t
    comm_type< long unsigned int > ()
    {
        return MPI_UNSIGNED_LONG;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline comm_t
    comm_type< double > ()
    {
        return MPI_DOUBLE ;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline comm_t
    comm_type< long double > ()
    {
        return MPI_LONG_DOUBLE ;
    }


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline comm_t
    comm_type< signed char > ()
    {
        return MPI_SIGNED_CHAR;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline comm_t
    comm_type< unsigned char > ()
    {
        return MPI_UNSIGNED_CHAR;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline comm_t
    comm_type< short int > ()
    {
        return MPI_SHORT;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline comm_t
    comm_type< long long int > ()
    {
        return MPI_LONG_LONG;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline comm_t
    comm_type< unsigned long long int > ()
    {
        return MPI_UNSIGNED_LONG_LONG;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline comm_t
    comm_type< float > ()
    {
        return MPI_FLOAT;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline comm_t
    comm_type< bool > ()
    {
        return MPI_CXX_BOOL;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef MPI_CXX_FLOAT_COMPLEX
    template<> inline comm_t
    comm_type< std::complex< float > > ()
    {
        return MPI_CXX_FLOAT_COMPLEX;
    }
#endif

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef MPI_CXX_DOUBLE_COMPLEX
    template<> inline comm_t
    comm_type< std::complex< double > > ()
    {
        return MPI_CXX_DOUBLE_COMPLEX ;
    }
#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef MPI_CXX_LONG_DOUBLE_COMPLEX
    template<> inline comm_t
    comm_type< std::complex< long double > > ()
    {
        return MPI_CXX_LONG_DOUBLE_COMPLEX ;
    }
#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif // BELFEM_MPI
//------------------------------------------------------------------------------
}

#endif //COMMTYPES_HPP
