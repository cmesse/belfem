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

#ifndef BELFEM_CL_ORDEREDMAP_HPP
#define BELFEM_CL_ORDEREDMAP_HPP

#include <map>


#include "cl_Map.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * @brief Sorted map (ordered key-value).
     *
     * @ingroup grp_containers
     * @see @ref containers_container_usage_guide
     */
    template< typename Key, typename Value >
    class OrderedMap
    {

        std::map< Key, Value > mMap;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------


        using const_iterator = typename std::map< Key, Value >::const_iterator ;

        /**
         * empty constructor
         */
        explicit OrderedMap() = default;

//------------------------------------------------------------------------------

        /**
         * copy constructor
         */
        OrderedMap( const OrderedMap< Key, Value > & aMap ) = default;

//------------------------------------------------------------------------------

        /**
         * move constructor
         */
        OrderedMap(OrderedMap<Key, Value>&& other) noexcept = default;

//------------------------------------------------------------------------------

        /**
        * Copy assignment operator
        * */
        OrderedMap<Key, Value>& operator=(const OrderedMap<Key, Value>& other) = default;

//------------------------------------------------------------------------------

        /**
        * Move assignment operator constructor
        */
        OrderedMap<Key, Value>& operator=( OrderedMap<Key, Value>&& other) noexcept = default;

//------------------------------------------------------------------------------

        /**
         * default destructor
         */
         ~OrderedMap() = default;

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
#endif //BELFEM_CL_ORDEREDMAP_HPP
