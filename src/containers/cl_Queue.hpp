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

#ifndef CL_QUEUE_HPP
#define CL_QUEUE_HPP

#include <queue>
#include "cl_Cell.hpp"

namespace belfem
{
    /**
     * @brief FIFO queue.
     *
     * @ingroup grp_containers
     * @see @ref containers_container_usage_guide
     */
    template< typename T >
    class Queue
    {
        std::queue< T > mQueue;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

        Queue() = default ;

        // Copy constructor
        Queue( const Queue< T > & other ) = default;

        // Move constructor
        Queue( Queue< T > && other ) noexcept = default;

        // Copy assignment operator
        Queue< T >& operator=( const Queue< T > & other ) = default;

        // Move assignment operator
        Queue< T >& operator=( Queue< T > && other ) noexcept = default;

        // special operator for converting a cell
        Queue( Cell< T > & aCell )
        {
            for( auto tValue : aCell )
            {
                mQueue.push( tValue );
            }
        }

//------------------------------------------------------------------------------

        ~Queue() = default;

//------------------------------------------------------------------------------

        void
        push( const T & aValue )
        {
            mQueue.push( aValue );
        }

//------------------------------------------------------------------------------

        T
        pop()
        {
            T aPop = mQueue.front();
            mQueue.pop();
            return aPop;
        }

//------------------------------------------------------------------------------

        size_t
        size() const
        {
            return mQueue.size();
        }

//------------------------------------------------------------------------------

        auto
        empty()-> decltype( mQueue.empty() ) const
        {
            return mQueue.empty();
        }
        
//------------------------------------------------------------------------------
    };
}

#endif //CL_QUEUE_HPP
