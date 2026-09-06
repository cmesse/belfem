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

#ifndef BELFEM_CL_MAP_HPP
#define BELFEM_CL_MAP_HPP

#include <unordered_map>


#include <string>

#include "typedefs.hpp"
#include "assert.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    namespace map
    {
        template < typename T >
        std::string
        KeyToString( const T & aKey )
        {
            return "unknown";
        }

        template <>
        inline std::string
        KeyToString( const std::string & aKey )
        {
            return aKey;
        }

        template <>
        inline std::string
        KeyToString( const int & aKey )
        {
            return std::to_string( aKey );
        }

        template <>
        inline std::string
        KeyToString( const unsigned int & aKey )
        {
            return std::to_string( aKey );
        }

        template <>
        inline std::string
        KeyToString( const long unsigned int & aKey )
        {
            return std::to_string( aKey );
        }
    }

//------------------------------------------------------------------------------

    /**
     * @brief Hash map (unordered key-value).
     *
     * @ingroup grp_containers
     * @see @ref containers_container_usage_guide
     */
    template< typename Key, typename Value >
    class Map
    {

        std::unordered_map< Key, Value > mMap;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        using const_iterator = typename std::unordered_map< Key, Value >::const_iterator ;

        /**
         * empty constructor
         */
        explicit Map() = default;

//------------------------------------------------------------------------------

        /**
         * copy constructor
         */
        Map( const Map< Key, Value > & aMap ) = default;

//------------------------------------------------------------------------------

        /**
         * move constructor
         */
        Map(Map<Key, Value>&& other) noexcept = default;

//------------------------------------------------------------------------------

        /**
        * Copy assignment operator
        * */
        Map<Key, Value>& operator=(const Map<Key, Value>& other) = default;

//------------------------------------------------------------------------------

        /**
        * Move assignment operator constructor
        */
        Map<Key, Value>& operator=(Map<Key, Value>&& other) noexcept = default;

//------------------------------------------------------------------------------

        /**
         * default destructor
         */
         ~Map() = default;

//------------------------------------------------------------------------------

        /**
         * clear map
         */
         void
         clear()
         {
             mMap.clear();
         }

//------------------------------------------------------------------------------

        /**
         * returns the size of the map
         */
        size_t
        size() const
        {
            return mMap.size();
        }

//------------------------------------------------------------------------------

        /**
         * insert operator
         */
         Value &
         operator[]( const Key & aKey )
         {
            return mMap[ aKey ];
         }

//------------------------------------------------------------------------------

        bool
        key_exists(  const Key & aKey  ) const
        {
             return mMap.find( aKey ) != mMap.end();
        }

//------------------------------------------------------------------------------

        auto
        begin() const -> decltype( mMap.begin() )
        {
            return mMap.begin();
        }

        auto
        end() const -> decltype( mMap.end() )
        {
            return mMap.end();
        }

        auto
        find ( const Key & tKey ) const -> decltype( mMap.find(tKey) )
        {
            return mMap.find( tKey );
        }

        auto
        empty() const -> decltype( mMap.empty() )
        {
             return mMap.empty();
        }

        auto
        get_entry( const index_t aIndex ) const -> decltype( * mMap.begin() )
         {
             BELFEM_ASSERT( aIndex < mMap.size(),
                           "Map::get_entry() index %lu out of bounds ( must be < %lu )",
                           ( long unsigned int ) aIndex,
                           ( long unsigned int ) mMap.size() );

             auto tEntry = mMap.begin();
             std::advance( tEntry, aIndex );

             const auto & aEntry = * tEntry;
             return aEntry;
         }

//------------------------------------------------------------------------------

        void
        erase_key(  const Key & aKey  )
        {
             // remove key from map
             mMap.erase( aKey );
        }
//------------------------------------------------------------------------------

        // expose the data container
        auto
        map_data() -> decltype( mMap ) &
        {
            return mMap ;
        }

//------------------------------------------------------------------------------

        /**
         * find operator
         */
        Value &
        operator()( const Key & aKey )
        {
            // check if key exists
            auto tIterator = mMap.find( aKey );

#if !defined( NDEBUG ) || defined( DEBUG )
            BELFEM_ASSERT( tIterator != mMap.end(),
                        "Key %s not found in map.",
                        map::KeyToString( aKey ).c_str() );

#else
            BELFEM_ERROR( tIterator != mMap.end(),
                       "Key not found in map." );
#endif

            return tIterator->second;
        }

//------------------------------------------------------------------------------

        /**
         * find operator ( const variant )
         */
        const Value &
        operator()( const Key & aKey ) const
        {
// check if key exists
            auto tIterator = mMap.find( aKey );

#if !defined( NDEBUG ) || defined( DEBUG )

            BELFEM_ASSERT( tIterator != mMap.end(),
                        "Key %s not found in map.",
                        map::KeyToString( aKey ).c_str() );

#else
            BELFEM_ERROR( tIterator != mMap.end(),
                       "Key not found in map." );
#endif

            return tIterator->second;
        }
    };
//------------------------------------------------------------------------------
} /* namespace belfem */
#endif //BELFEM_CL_MAP_HPP
