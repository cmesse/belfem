//
// Created by Christian Messe on 2019-01-18.
//

#ifndef BELFEM_COMM_TYPEDEFS_HPP
#define BELFEM_COMM_TYPEDEFS_HPP

#ifdef BELFEM_MPI
#include <mpi.h>
typedef MPI_Datatype comm_data_t ;
#else
// define fake type if there is no MPI
typedef int comm_data_t;
#endif

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     *
     * @brief                 returns an MPI enum defining the
     *                        data type that is to be communicated.
     *
     * @param[in] aSample     primitive data type with arbitrary value
     *
     * see also http://mpitutorial.com/tutorials/mpi-send-and-receive/
     */
    template < typename T >
    inline comm_data_t
    get_comm_datatype( const T & aSample )
    {
        BELFEM_ASSERT( false, "unknown datatype" );
#ifdef BELFEM_MPI
        return MPI_DATATYPE_NULL;
#else
        return 0;
#endif
    }
#ifdef BELFEM_MPI

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline comm_data_t
    get_comm_datatype( const int & aSample )
    {
        return MPI_INT;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline comm_data_t
    get_comm_datatype( const long int & aSample )
    {
        return MPI_LONG;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline comm_data_t
    get_comm_datatype( const unsigned int & aSample )
    {
        return MPI_UNSIGNED;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline comm_data_t
    get_comm_datatype( const long unsigned int & aSample )
    {
        return MPI_UNSIGNED_LONG;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline comm_data_t
    get_comm_datatype( const double & aSample )
    {
        return MPI_DOUBLE;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline comm_data_t
    get_comm_datatype( const long double & aSample )
    {
        return MPI_LONG_DOUBLE;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// MPI_CXX_DOUBLE_COMPLEX is supported since MPI-3.0
#ifdef MPI_CXX_DOUBLE_COMPLEX
    template <>
    inline comm_data_t
    get_comm_datatype( const std::complex<double> & aSample )
    {
        return MPI_CXX_DOUBLE_COMPLEX;
    }
#endif

// MPI_CXX_LONG_DOUBLE_COMPLEX is supported since MPI-3.0
#ifdef MPI_CXX_LONG_DOUBLE_COMPLEX
    template <>
    inline comm_data_t
    get_comm_datatype( const std::complex<long double> & aSample )
    {
        return MPI_CXX_LONG_DOUBLE_COMPLEX;
    }
#endif
#endif // end BELFEM_MPI
//------------------------------------------------------------------------------
} /* namespace belfem */
#endif //BELFEM_COMM_TYPEDEFS_HPP
