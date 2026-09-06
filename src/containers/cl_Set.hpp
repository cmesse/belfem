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


#ifndef BELFEM_CL_SET_HPP
#define BELFEM_CL_SET_HPP

#include <unordered_set>
#include <string>

#include "assert.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * @brief Hash set with set operations.
     *
     * @ingroup grp_containers
     * @see @ref containers_container_usage_guide
     */
    template< typename Key >
    class Set
    {
        std::unordered_set< Key > mSet;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        using iterator = typename std::unordered_set< Key >::iterator;
        using const_iterator = typename std::unordered_set< Key >::const_iterator;

        /**
         * Empty constructor
         */
        Set() = default;

//------------------------------------------------------------------------------

        /**
         * Constructor with initializer list
         */
        Set( std::initializer_list< Key > aInit ) : mSet( aInit )
        {
        }

//------------------------------------------------------------------------------

        /**
         * Constructor from iterators
         */
        template< typename InputIt >
        Set( InputIt aFirst, InputIt aLast ) : mSet( aFirst, aLast )
        {
        }

//------------------------------------------------------------------------------

        /**
         * Copy constructor
         */
        Set( const Set< Key > & aSet ) = default;

//------------------------------------------------------------------------------

        /**
         * Move constructor
         */
        Set( Set< Key > && aSet ) noexcept = default;

//------------------------------------------------------------------------------

        /**
         * Copy assignment operator
         */
        Set< Key > & operator=( const Set< Key > & aSet ) = default;

//------------------------------------------------------------------------------

        /**
         * Move assignment operator
         */
        Set< Key > & operator=( Set< Key > && aSet ) noexcept = default;

//------------------------------------------------------------------------------

        /**
         * Destructor
         */
        ~Set() = default;

//------------------------------------------------------------------------------

        /**
         * Clear the set
         */
        void
        clear()
        {
            mSet.clear();
        }

//------------------------------------------------------------------------------

        /**
         * Returns the size of the set
         */
        size_t
        size() const
        {
            return mSet.size();
        }

//------------------------------------------------------------------------------

        /**
         * Check if set is empty
         */
        bool
        empty() const
        {
            return mSet.empty();
        }

//------------------------------------------------------------------------------

        /**
         * Insert a key into the set
         * @return pair of iterator and bool (true if inserted, false if already existed)
         */
        std::pair< iterator, bool >
        insert( const Key & aKey )
        {
            return mSet.insert( aKey );
        }

//------------------------------------------------------------------------------

        /**
         * Insert a key using move semantics
         */
        std::pair< iterator, bool >
        insert( Key && aKey )
        {
            return mSet.insert( std::move( aKey ) );
        }

//------------------------------------------------------------------------------

        /**
         * Insert a range of elements
         */
        template< typename InputIt >
        void
        insert( InputIt aFirst, InputIt aLast )
        {
            mSet.insert( aFirst, aLast );
        }

//------------------------------------------------------------------------------

        /**
         * Emplace a key (construct in-place)
         */
        template< typename... Args >
        std::pair< iterator, bool >
        emplace( Args&&... args )
        {
            return mSet.emplace( std::forward< Args >( args )... );
        }

//------------------------------------------------------------------------------

        /**
         * Check if a key exists in the set
         */
        bool
        key_exists( const Key & aKey ) const
        {
            return mSet.find( aKey ) != mSet.end();
        }

//------------------------------------------------------------------------------

        /**
         * Alternative name for key_exists (more set-like)
         */
        bool
        contains( const Key & aKey ) const
        {
            return this->key_exists( aKey );
        }

//------------------------------------------------------------------------------

        /**
         * Count occurrences of key (0 or 1 for set)
         */
        size_t
        count( const Key & aKey ) const
        {
            return mSet.count( aKey );
        }

//------------------------------------------------------------------------------

        /**
         * Find an element
         */
        iterator
        find( const Key & aKey )
        {
            return mSet.find( aKey );
        }

//------------------------------------------------------------------------------

        /**
         * Find an element (const version)
         */
        const_iterator
        find( const Key & aKey ) const
        {
            return mSet.find( aKey );
        }

//------------------------------------------------------------------------------

        /**
         * Erase a key from the set
         * @return number of elements removed (0 or 1)
         */
        size_t
        erase( const Key & aKey )
        {
            return mSet.erase( aKey );
        }

//------------------------------------------------------------------------------

        /**
         * Erase an element by iterator
         */
        iterator
        erase( const_iterator aPos )
        {
            return mSet.erase( aPos );
        }

//------------------------------------------------------------------------------

        /**
         * Erase a range of elements
         */
        iterator
        erase( const_iterator aFirst, const_iterator aLast )
        {
            return mSet.erase( aFirst, aLast );
        }

//------------------------------------------------------------------------------

        /**
         * Reserve space for at least n elements
         */
        void
        reserve( size_t n )
        {
            mSet.reserve( n );
        }

//------------------------------------------------------------------------------

        /**
         * Get begin iterator
         */
        iterator
        begin()
        {
            return mSet.begin();
        }

//------------------------------------------------------------------------------

        /**
         * Get begin iterator (const)
         */
        const_iterator
        begin() const
        {
            return mSet.begin();
        }

//------------------------------------------------------------------------------

        /**
         * Get end iterator
         */
        iterator
        end()
        {
            return mSet.end();
        }

//------------------------------------------------------------------------------

        /**
         * Get end iterator (const)
         */
        const_iterator
        end() const
        {
            return mSet.end();
        }

//------------------------------------------------------------------------------

        /**
         * Expose the underlying container
         */
        std::unordered_set< Key > &
        set_data()
        {
            return mSet;
        }

//------------------------------------------------------------------------------

        /**
         * Expose the underlying container (const)
         */
        const std::unordered_set< Key > &
        set_data() const
        {
            return mSet;
        }

//------------------------------------------------------------------------------

        /**
         * Swap contents with another set
         */
        void
        swap( Set< Key > & aOther )
        {
            mSet.swap( aOther.mSet );
        }

//------------------------------------------------------------------------------

        /**
         * Set union operation
         * @return new set containing elements from both sets
         */
        Set< Key >
        operator|( const Set< Key > & aOther ) const
        {
            Set< Key > tResult( *this );
            tResult.insert( aOther.begin(), aOther.end() );
            return tResult;
        }

//------------------------------------------------------------------------------

        /**
         * Set intersection operation
         * @return new set containing only common elements
         */
        Set< Key >
        operator&( const Set< Key > & aOther ) const
        {
            Set< Key > tResult;
            for ( const auto & tKey : mSet )
            {
                if ( aOther.contains( tKey ) )
                {
                    tResult.insert( tKey );
                }
            }
            return tResult;
        }

//------------------------------------------------------------------------------

        /**
         * Set difference operation
         * @return new set containing elements in this but not in other
         */
        Set< Key >
        operator-( const Set< Key > & aOther ) const
        {
            Set< Key > tResult;
            for ( const auto & tKey : mSet )
            {
                if ( !aOther.contains( tKey ) )
                {
                    tResult.insert( tKey );
                }
            }
            return tResult;
        }

//------------------------------------------------------------------------------

        /**
         * Set symmetric difference operation
         * @return new set containing elements in either set but not both
         */
        Set< Key >
        operator^( const Set< Key > & aOther ) const
        {
            Set< Key > tResult;

            // Add elements from this set not in other
            for ( const auto & tKey : mSet )
            {
                if ( !aOther.contains( tKey ) )
                {
                    tResult.insert( tKey );
                }
            }

            // Add elements from other set not in this
            for ( const auto & tKey : aOther )
            {
                if ( !this->contains( tKey ) )
                {
                    tResult.insert( tKey );
                }
            }

            return tResult;
        }

//------------------------------------------------------------------------------

        /**
         * Check if this is a subset of another set
         */
        bool
        is_subset_of( const Set< Key > & aOther ) const
        {
            if ( this->size() > aOther.size() )
            {
                return false;
            }

            for ( const auto & tKey : mSet )
            {
                if ( !aOther.contains( tKey ) )
                {
                    return false;
                }
            }
            return true;
        }

//------------------------------------------------------------------------------

        /**
         * Check if this is a superset of another set
         */
        bool
        is_superset_of( const Set< Key > & aOther ) const
        {
            return aOther.is_subset_of( *this );
        }

//------------------------------------------------------------------------------

        /**
         * Equality operator
         */
        bool
        operator==( const Set< Key > & aOther ) const
        {
            if ( this->size() != aOther.size() )
            {
                return false;
            }

            for ( const auto & tKey : mSet )
            {
                if ( !aOther.contains( tKey ) )
                {
                    return false;
                }
            }
            return true;
        }

//------------------------------------------------------------------------------

        /**
         * Inequality operator
         */
        bool
        operator!=( const Set< Key > & aOther ) const
        {
            return !( *this == aOther );
        }

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
} // namespace belfem

#endif // BELFEM_CL_SET_HPP