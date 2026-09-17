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

#ifndef BELFEM_CL_CELL_HPP
#define BELFEM_CL_CELL_HPP

#include <vector>
#include <initializer_list>
#include <algorithm> // for unique and reverse
#include <iterator>  // for std::make_move_iterator
#include <sstream>   // for std::ostringstream (generic print)

#include "typedefs.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace cell
    {
        void
        error_out_of_bounds( const index_t aI, const index_t aN );
    }
//------------------------------------------------------------------------------

    /**
     * A bounds-checked wrapper around std::vector. A Cell< T > owns its
     * elements. A Cell< T* > stores pointers but does not own the pointees.
     * Destruction, clear() and set_size() do not delete them. Copying a
     * Cell< T* > copies the pointers. sort(), unique(), reverse() and
     * append() are free functions.
     *
     * @ingroup grp_containers
     * @see @ref containers_container_usage_guide
     */
    template< typename T >
    class Cell
    {
//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

        std::vector< T > mCell;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        Cell() = default ;

//------------------------------------------------------------------------------

        /** Creates an empty Cell with capacity for aReserve elements. size()
         *  is 0 and operator() is out of bounds until push() or set_size().
         *  The two-argument constructor creates elements. */
        Cell( const std::size_t aReserve )
        {
            mCell.reserve( aReserve );
        }

        Cell( const std::size_t aSize, const T aInitializationValue )
        {
            mCell.resize( aSize, aInitializationValue );
        }

//------------------------------------------------------------------------------

        Cell( const std::initializer_list< T > & aInitList ) :
                mCell( aInitList )
        {
        }

//------------------------------------------------------------------------------

        // Copy constructor - explicit implementation
        Cell( const Cell< T > & other ) : mCell( other.mCell )
        {
        }

//------------------------------------------------------------------------------

        // Move constructor - explicit implementation with noexcept
        Cell( Cell< T > && other ) noexcept : mCell( std::move( other.mCell ) )
        {
            // other.mCell is now in valid but unspecified state (empty for vector)
        }

//------------------------------------------------------------------------------

        // Copy assignment operator - explicit implementation
        Cell< T >& operator=( const Cell< T > & other )
        {
            if ( this != &other )
            {
                mCell = other.mCell;
            }
            return *this;
        }

//------------------------------------------------------------------------------

        // Move assignment operator - explicit implementation with noexcept
        Cell< T >& operator=( Cell< T > && other ) noexcept
        {
            if ( this != &other )
            {
                mCell = std::move( other.mCell );
                // other.mCell is now in valid but unspecified state
            }
            return *this;
        }

//------------------------------------------------------------------------------

        ~Cell() = default;

//------------------------------------------------------------------------------

        /** Returns a borrowed pointer into contiguous storage, possibly nullptr
         *  when the Cell is empty ( test size(), never data() ). Any call that
         *  may reallocate invalidates it. MPI accepts it at count 0. */
        T * data()
        {
            return mCell.data();
        }

//------------------------------------------------------------------------------

        const T * data() const
        {
            return mCell.data();
        }

//------------------------------------------------------------------------------

        /** The underlying vector. Resizing it here bypasses the bounds checks
         *  of operator(). */
        std::vector< T > &
        vector_data()
        {
            return mCell;
        }

//------------------------------------------------------------------------------

        const std::vector< T > &
        vector_data() const
        {
            return mCell;
        }
        
//------------------------------------------------------------------------------

        auto
        operator()( size_t aIndex )
        -> decltype( mCell[ aIndex ] )
        {
            BELFEM_ASSERT( aIndex < mCell.size() ,
                "Cell index out of bounds: %lu (expect < %lu)",
                (long unsigned int)aIndex,
                (long unsigned int)mCell.size() );
            return mCell[ aIndex ];
        }


//------------------------------------------------------------------------------

        auto
        operator()( size_t aIndex ) const
        -> decltype( mCell[ aIndex ] )
        {
            BELFEM_ASSERT( aIndex < mCell.size() ,
            "Cell index out of bounds: %lu (expect < %lu)",
            (long unsigned int)aIndex,
            (long unsigned int)mCell.size() );
                    return mCell[ aIndex ];
        }

//------------------------------------------------------------------------------

        size_t
        size() const
        {
            return mCell.size();
        }

//------------------------------------------------------------------------------

        /** New elements are value-initialized ( zero for arithmetic T,
         *  nullptr for pointers ). The two-argument overload sets a value. */
        void
        set_size( const size_t aSize )
        {
            mCell.resize( aSize );
        }

//------------------------------------------------------------------------------

        void
        set_size( const size_t aSize, const T & aValue )
        {
            mCell.resize( aSize, aValue );
        }

//------------------------------------------------------------------------------

        auto
        begin() -> decltype( mCell.begin() )
        {
            return mCell.begin();
        }

//------------------------------------------------------------------------------

        auto
        end() -> decltype( mCell.end() )
        {
            return mCell.end();
        }

//------------------------------------------------------------------------------

        auto
        begin() const -> decltype( mCell.begin() )
        {
            return mCell.begin();
        }

//------------------------------------------------------------------------------

        auto
        end() const -> decltype( mCell.end() )
        {
            return mCell.end();
        }

//------------------------------------------------------------------------------

        /** Empties the Cell. It retains its capacity; shrink_to_fit() asks the
         *  vector to release it. Pointees of a Cell< T* > are not deleted. */
        void
        clear()
        {
            mCell.clear();
        }

//------------------------------------------------------------------------------

        void
        reserve( const size_t aSize )
        {
            mCell.reserve( aSize );
        }

//------------------------------------------------------------------------------

        void
        push( const T & aValue )
        {
            mCell.push_back( aValue );
        }

//------------------------------------------------------------------------------

        void
        push( T && aValue )
        {
            mCell.push_back( std::move( aValue ) );
        }

//------------------------------------------------------------------------------

        template< typename... Args >
        void
        emplace( Args&&... args )
        {
            mCell.emplace_back( std::forward< Args >( args )... );
        }

//------------------------------------------------------------------------------

        T
        pop()
        {
            T aPop = std::move( mCell.back() );  // Move instead of copy
            mCell.pop_back();
            return aPop;
        }


//------------------------------------------------------------------------------

        void
        shrink_to_fit()
        {
            mCell.shrink_to_fit();
        }

//------------------------------------------------------------------------------

        T &
        first()
        {
#if !defined( NDEBUG ) || defined( DEBUG )
            return mCell.at( 0 );
#else
            return mCell[ 0 ];
#endif
        }

//------------------------------------------------------------------------------

        const T &
        first() const
        {
#if !defined( NDEBUG ) || defined( DEBUG )
            return mCell.at( 0 );
#else
            return mCell[ 0 ];
#endif
        }

//------------------------------------------------------------------------------

        T &
        last()
        {
#if !defined( NDEBUG ) || defined( DEBUG )
            return mCell.at( mCell.size()-1 );
#else
            return mCell[ mCell.size()-1 ];
#endif
        }

//------------------------------------------------------------------------------

        const T &
        last() const
        {
#if !defined( NDEBUG ) || defined( DEBUG )
            return mCell.at( mCell.size()-1 );
#else
            return mCell[ mCell.size()-1 ];
#endif
        }

//------------------------------------------------------------------------------

        /**
         * Swap contents with another Cell
         */
        void
        swap( Cell< T > & other ) noexcept
        {
            mCell.swap( other.mCell );
        }

//------------------------------------------------------------------------------

        bool
        empty() const
        {
            return mCell.empty();
        }

//------------------------------------------------------------------------------

        size_t
        capacity() const
        {
            return mCell.capacity();
        }

//------------------------------------------------------------------------------

        auto
        insert( typename std::vector< T >::const_iterator pos, const T & value )
            -> decltype( mCell.insert( pos, value ) )
        {
            return mCell.insert( pos, value );
        }

//------------------------------------------------------------------------------

        auto
        insert( typename std::vector< T >::const_iterator pos, T && value )
            -> decltype( mCell.insert( pos, std::move( value ) ) )
        {
            return mCell.insert( pos, std::move( value ) );
        }

//------------------------------------------------------------------------------

        auto
        erase( typename std::vector< T >::const_iterator pos )
            -> decltype( mCell.erase( pos ) )
        {
            return mCell.erase( pos );
        }

//------------------------------------------------------------------------------

        /**
         * Erase range of elements
         */
        auto
        erase( typename std::vector< T >::const_iterator first,
               typename std::vector< T >::const_iterator last )
            -> decltype( mCell.erase( first, last ) )
        {
            return mCell.erase( first, last );
        }
        
//------------------------------------------------------------------------------

        void
        print( const std::string aLabel="" ) const ;
        
//------------------------------------------------------------------------------
    };

    template< typename T >
    void
    sort( Cell< T > & aCell )
    {
        std::vector< T > & tVec = aCell.vector_data();

        std::sort( tVec.begin(), tVec.end() );
    }

//------------------------------------------------------------------------------

    template< typename T, class C >
    void
    sort( Cell< T > & aCell, C & aComp, const size_t aNumberOfItems=0 )
    {
        std::vector< T > & tVec = aCell.vector_data();

        if( aNumberOfItems==0 )
        {
            std::sort( tVec.begin(), tVec.end(), aComp );
        }
        else
        {
            std::sort( tVec.begin(), tVec.begin()+aNumberOfItems, aComp );
        }
    }

//------------------------------------------------------------------------------

    /** Sorts aCell in place and removes the duplicates. It does not preserve
     *  the order. */
    template< typename T >
    void
    unique( Cell< T > & aCell )
    {
        std::vector< T > & tVec = aCell.vector_data();

        std::sort( tVec.begin(), tVec.end() );

        tVec.erase( std::unique( tVec.begin(), tVec.end() ), tVec.end() );
    }

//------------------------------------------------------------------------------

    /**
     * returns the position of a member inside a cell.
     * The cell must be unique and sorted for this function to work.
     */
    template< typename T >
    index_t
    find_index_in_unique_cell( const Cell< T > & aCell,  const T & aMember )
    {
        auto tIterator = std::lower_bound( aCell.vector_data().begin(),
            aCell.vector_data().end(), aMember );

        BELFEM_ASSERT( tIterator != aCell.vector_data().end()
            && ! ( aMember < *tIterator ),
            "member not found in cell (is the cell sorted and unique?)" );

        return tIterator - aCell.vector_data().begin();
    }

//------------------------------------------------------------------------------

    template< typename T >
    void
    reverse( Cell< T > & aCell )
    {
        std::vector< T > & tVec = aCell.vector_data();

        std::reverse( tVec.begin(), tVec.end() );
    }

//------------------------------------------------------------------------------

    template< typename T >
    void
    append( Cell< T > & aA, const Cell< T > & aB )
    {
        aA.vector_data().reserve(aA.vector_data().size() + aB.vector_data().size());

        aA.vector_data().insert(
            aA.vector_data().end(),
            aB.vector_data().begin(),
            aB.vector_data().end()
            );
    }

//------------------------------------------------------------------------------

    template< typename T >
    void
    append_move( Cell< T > & aTarget, Cell< T > & aSource )
    {
        aTarget.vector_data().reserve(aTarget.vector_data().size() + aSource.vector_data().size());

        aTarget.vector_data().insert(
            aTarget.vector_data().end(),
            std::make_move_iterator(aSource.vector_data().begin()),
            std::make_move_iterator(aSource.vector_data().end()));

        aSource.clear();
    }

//------------------------------------------------------------------------------

    /**
     * Swap specialization for Cell
     */
    template< typename T >
    void
    swap( Cell< T > & aA, Cell< T > & aB ) noexcept
    {
        aA.swap( aB );
    }
    
//------------------------------------------------------------------------------
    
// turn off annoying waning
#ifdef BELFEM_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wformat-security"
#elif BELFEM_GCC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat"
#elif BELFEM_INTEL
#pragma warning push
#pragma warning disable 1595
#endif

    template < typename T > void
    Cell< T >::print( const std::string aLabel ) const
    {
        FILE * tOutFile = stdout;

        fprintf( tOutFile, "\n%s = [ ... \n", aLabel.c_str() );

        const Cell< T > & tThis = *this;

        uint tLength = tThis.size();

        for( uint i=0; i< tLength; ++i )
        {
            // use std::ostringstream for generic types (including std::string)
            std::ostringstream tStream;
            tStream << tThis( i );
            fprintf( tOutFile, "%s; ", tStream.str().c_str() );

            if( i < tLength-1 )
            {
                fprintf( tOutFile, "...\n" );
            }
            else
            {
                fprintf( tOutFile, "];\n" );
            }
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template <>
    inline void
    Cell< real >::print( const std::string aLabel ) const
    {
        FILE * tOutFile = stdout;

        fprintf( tOutFile, "%s = [ ... \n", aLabel.c_str() );

        const Cell< real > & tThis = *this;

        uint tLength = tThis.size();

        for( uint i=0; i< tLength; ++i )
        {

            fprintf( tOutFile, "%+.15e; ", ( double ) tThis( i ) );


            if( i < tLength-1 )
            {
                fprintf( tOutFile, "...\n" );
            }
            else
            {
                fprintf( tOutFile, "];\n" );
            }
        }
    }

#ifdef BELFEM_CLANG
#pragma clang diagnostic pop
#elif BELFEM_GCC
#pragma GCC diagnostic pop
#elif BELFEM_INTEL
#pragma warning pop
#endif

    
//------------------------------------------------------------------------------
} /* namespace belfem */

#endif //BELFEM_CL_CELL_HPP