//
// Created by Christian Messe on 17.11.18.
//

#ifndef BELFEM_CL_MAP_HPP
#define BELFEM_CL_MAP_HPP

#include <map>

#include "typedefs.hpp"
#include "assert.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    namespace map
    {
        template < typename T >
        inline std::string
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

    template< typename Key, typename Value >
    class Map
    {
        std::map< Key, Value > mMap;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

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

        typename std::map< Key, Value >::const_iterator
        begin() const
        {
            return mMap.begin();
        }

        typename std::map< Key, Value >::const_iterator
        end() const
        {
            return mMap.end();
        }
//------------------------------------------------------------------------------

        void
        erase_key(  const Key & aKey  )
        {
            // check if key exists
            if( this->key_exists( aKey ) )
            {
                // remove key from map
                mMap.erase( aKey );
            }
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
} /* namespace ssf */
#endif //SSF_CL_MAP_HPP
