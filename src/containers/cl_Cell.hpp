//
// Created by Christian Messe on 17.11.18.
//

#ifndef BELFEM_CL_CELL_HPP
#define BELFEM_CL_CELL_HPP

#include <vector>
#include <initializer_list>
#include <algorithm> // for unique and reverse

#include "typedefs.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * Cell is a wrapper around the standard vector.
     * It also supports sorting and unique making of entries.
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

        Cell()
        {
        }

//------------------------------------------------------------------------------

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

        ~Cell() = default;

//------------------------------------------------------------------------------

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
#if !defined( NDEBUG ) || defined( DEBUG )
        -> decltype( mCell.at( aIndex ) )
        {
            return mCell.at( aIndex );
        }
#else
        -> decltype( mCell[ aIndex ] )
        {
            return mCell[ aIndex ];
        }
#endif

//------------------------------------------------------------------------------

        auto
        operator()( size_t aIndex ) const
#if !defined( NDEBUG ) || defined( DEBUG )
        -> decltype( mCell.at( aIndex ) )
        {
            return mCell.at( aIndex );
        }
#else
        -> decltype( mCell[ aIndex ] )
        {
            return mCell[ aIndex ];
        }
#endif

//------------------------------------------------------------------------------

        /**
         * return the size of the Cell
         */
        size_t
        size() const
        {
            return mCell.size();
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
        begin() const -> decltype( mCell.begin() ) const
        {
            return mCell.begin();
        }

//------------------------------------------------------------------------------

        auto
        end() const -> decltype( mCell.end() ) const
        {
            return mCell.end();
        }

//------------------------------------------------------------------------------

        /**
         * clear the memory
         */
        void
        clear()
        {
            mCell.clear();
        }

//------------------------------------------------------------------------------

        /**
         * push an entry to the end of the cell
         */
        void
        push( const T & aValue )
        {
            mCell.push_back( aValue );
        }

//------------------------------------------------------------------------------

        /**
         * pop an entry from the Cell
         */
        T
        pop()
        {
            T aPop = mCell.back();
            mCell.pop_back();
            return aPop;
        }

//------------------------------------------------------------------------------

        /**
         * free unused memory
         */
        void
        shrink_to_fit()
        {
            mCell.shrink_to_fit();
        }

//------------------------------------------------------------------------------

    };

    template< typename T >
    void
    sort( Cell< T > & aCell )
    {
        // get ref to data
        std::vector< T > & tVec = aCell.vector_data();

        // sort data
        std::sort( tVec.begin(), tVec.end() );
    }

//------------------------------------------------------------------------------

    template< typename T, class C >
    void
    sort( Cell< T > & aCell, C & aComp, const size_t aNumberOfItems=0 )
    {
        // get ref to data
        std::vector< T > & tVec = aCell.vector_data();

        // sort data
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

    template< typename T >
    void
    unique( Cell< T > & aCell )
    {
        // get ref to data
        std::vector< T > & tVec = aCell.vector_data();

        // sort data
        std::sort( tVec.begin(), tVec.end() );

        // trim vector
        tVec.erase( std::unique( tVec.begin(), tVec.end() ), tVec.end() );
    }

//------------------------------------------------------------------------------

    template< typename T >
    void
    reverse( Cell< T > & aCell )
    {
        // get ref to data
        std::vector< T > & tVec = aCell.vector_data();

        // reverse data
        std::reverse( tVec.begin(), tVec.end() );
    }

//------------------------------------------------------------------------------
} /* namespace belfem */

#endif //BELFEM_CL_CELL_HPP