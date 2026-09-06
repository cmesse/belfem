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


#ifndef BELFEM_CL_SHIFTREGISTER_HPP
#define BELFEM_CL_SHIFTREGISTER_HPP

#include <algorithm>
#include <initializer_list>
#include <cstdlib>
#include <new>
#include <type_traits>
#include "typedefs.hpp"
#include "assert.hpp"

namespace belfem
{
    // Forward declarations for shift-register-safe specializations
    template< typename T > class Vector;
    template< typename T > class Matrix;

//------------------------------------------------------------------------------

    /**
     * Type trait: may T be stored in a ShiftRegister?
     *
     * ShiftRegister keeps a malloc'd buffer of (capacity+1) slots. How those
     * slots are managed depends on T:
     *
     *  - Trivially copyable T (int, real, index_t, …): slots are used raw.
     *    Assignment just overwrites bytes and free() owes no destructor pass.
     *    This is the zero-overhead fast path.
     *
     *  - Owning T (Vector<T>, Matrix<T>, …): the register manages element
     *    lifetime explicitly.  Every slot is placement-new default-constructed
     *    when the buffer is reserved and destroyed before the buffer is freed,
     *    so all element operations (push / shift / copy / revert) act on valid,
     *    constructed objects — exactly what a non-trivial operator= expects.
     *    This requires T to be default-constructible.
     *
     * Trivially copyable types qualify automatically.  Owning types must be
     * whitelisted below.  Vector<T> and Matrix<T> qualify: both the Armadillo
     * and Blaze backends have a valid default-constructed empty state, a
     * well-behaved copy/move operator=, and a moved-from state that is a valid
     * empty object.  Do NOT whitelist a type that is not default-constructible.
     */
    template< typename T >
    struct is_shift_register_safe : std::is_trivially_copyable< T > {};

    template< typename T >
    struct is_shift_register_safe< Vector< T > > : std::true_type {};

    template< typename T >
    struct is_shift_register_safe< Matrix< T > > : std::true_type {};

//------------------------------------------------------------------------------

    /**
     * @brief Fixed-capacity FIFO with history, for time-stepping.
     *
     * @ingroup grp_containers
     * @see @ref containers_container_usage_guide
     */
    template < typename T >
    class ShiftRegister
    {
        static_assert( is_shift_register_safe< T >::value,
            "ShiftRegister<T> uses malloc/free — T must be trivially copyable, "
            "or specialize belfem::is_shift_register_safe<T> for an owning type "
            "that is default-constructible and has a well-behaved operator= "
            "(see cl_ShiftRegister.hpp for details)" );

        // owning (non-trivial) types are lifetime-managed via placement-new /
        // explicit destruction, which requires a default constructor
        static_assert( std::is_trivially_copyable< T >::value
                       || std::is_default_constructible< T >::value,
            "ShiftRegister<T>: an owning T must be default-constructible so its "
            "slots can be default-constructed on reserve()" );

        // fixed size data container
        // we keep one more so that we can do a revert
        // the backup value is not passed during copy ( a copy cannot revert );
        // a move steals the whole buffer, backup included
        T * mData = nullptr ;

        // size with actually filled values
        size_t mSize = 0;

        size_t mCapacity = 0 ;

        enum class RevertState : uint8_t { None, CanRevert, CanRevertFull };
        RevertState mCanRevert = RevertState::None;

//------------------------------------------------------------------------------
            private:
//------------------------------------------------------------------------------

        // Default-construct all (mCapacity+1) slots of a freshly malloc'd
        // buffer, so every slot is a valid object before any assignment.
        // Compiles to nothing for trivially copyable T (slots used raw).
        void
        construct_slots()
        {
            if constexpr ( ! std::is_trivially_copyable< T >::value )
            {
                for ( size_t k = 0; k < mCapacity + 1; ++k )
                {
                    ::new ( static_cast< void * >( mData + k ) ) T();
                }
            }
        }

//------------------------------------------------------------------------------

        // Destroy all (mCapacity+1) slots before the buffer is freed, so
        // owning elements release their storage. Must be called while mData
        // and mCapacity still describe the live buffer.
        // Compiles to nothing for trivially copyable T.
        void
        destroy_slots()
        {
            if constexpr ( ! std::is_trivially_copyable< T >::value )
            {
                for ( size_t k = 0; k < mCapacity + 1; ++k )
                {
                    ( mData + k )->~T();
                }
            }
        }

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

        // Constructor with capacity only ( register starts empty )
        ShiftRegister( const size_t aCapacity )
        {
            this->reserve( aCapacity );
        }

        // Constructor with fill value
        ShiftRegister( const size_t aCapacity, const T & aInitValue )
        {
            this->reserve( aCapacity );
            std::fill( mData, mData + aCapacity, aInitValue );
            mSize = aCapacity ;
        }


        // Constructor with initializer list (fills from newest to oldest)
        ShiftRegister( const std::initializer_list< T > & aInitList )
        {
            this->reserve( aInitList.size() );

            uint k = 0;
            for( const T & value : aInitList )
            {
                mData[ k++ ] = value;
            }

            // remember size
            mSize = aInitList.size();
        }

//------------------------------------------------------------------------------

        // Copy constructor
        ShiftRegister( const ShiftRegister & aOther )
        {
            // reserve memory
            this->reserve( aOther.capacity() );

            // Copy only the valid elements
            std::copy( aOther.begin(), aOther.end(), mData );

            mCanRevert = RevertState::None ;

            // remember size
            mSize = aOther.size();
        }

//------------------------------------------------------------------------------

        // Move constructor
        ShiftRegister( ShiftRegister && aOther ) noexcept
        {
            // a freshly constructed object owns no buffer, so just steal
            // Move the pointer
            mData = aOther.data();
            mSize = aOther.size();
            mCapacity = aOther.capacity();
            mCanRevert = aOther.mCanRevert;

            // reset the source object
            aOther.mData = nullptr;
            aOther.mCapacity = 0;
            aOther.mSize = 0;
            aOther.mCanRevert = RevertState::None;
        }

        // Copy assignment operator
        ShiftRegister &
        operator=( const ShiftRegister & aOther )
        {
            if( this != &aOther )
            {
                this->reserve( aOther.capacity() );

                mSize = aOther.size();

                std::copy( aOther.begin(), aOther.end(), mData );
                mCanRevert = RevertState::None ;

            }
            return *this;
        }

//------------------------------------------------------------------------------

        // Move assignment operator
        ShiftRegister &
        operator=( ShiftRegister && aOther ) noexcept
        {
            if( this != &aOther )
            {
                if ( mData != nullptr )
                {
                    this->destroy_slots();
                    std::free( mData );
                }

                // Move the pointer
                mData = aOther.data();
                mSize = aOther.size();
                mCapacity = aOther.capacity();
                mCanRevert = aOther.mCanRevert;

                // reset the source object
                aOther.mData = nullptr;
                aOther.mCapacity = 0;
                aOther.mSize = 0;
                aOther.mCanRevert = RevertState::None;
            }
            return *this;
        }

//------------------------------------------------------------------------------

        ~ShiftRegister()
        {
            if ( mData != nullptr )
            {
                this->destroy_slots();
                std::free( mData ) ;
            }
        }

//------------------------------------------------------------------------------

        void
        push( T & aValue )
        {
            if( mSize > 0 )
            {
                // Determine how many elements to move
                index_t tElementsToMove = std::min( mSize, mCapacity-1 );

                // if we have a backup, we remember this value too
                if ( mSize == mCapacity ) ++tElementsToMove;

                // Move elements to the right using std::move_backward
                // This moves from [0, tElementsToMove) to [1, tElementsToMove+1)
                // including the backup value
                std::move_backward( mData, mData + tElementsToMove, mData + tElementsToMove + 1 );
            }

            // Insert new value at position 0
            mData[ 0 ] = aValue;

            mCanRevert = mSize == mCapacity ? RevertState::CanRevertFull : RevertState::CanRevert ;

            // Update size (cap at N)
            mSize = std::min( mSize+1, mCapacity );
        }

//------------------------------------------------------------------------------

        void
        push( const T & aValue )
        {
            if( mSize > 0 )
            {
                // Determine how many elements to move
                index_t tElementsToMove = std::min( mSize, mCapacity-1 );

                // if we have a backup, we remember this value too
                if ( mSize == mCapacity ) ++tElementsToMove;

                // Move elements to the right using std::move_backward
                // This moves from [0, tElementsToMove) to [1, tElementsToMove+1)
                std::move_backward( mData, mData + tElementsToMove, mData + tElementsToMove + 1 );
            }

            // Insert new value at position 0
            mData[ 0 ] = aValue;

            mCanRevert = mSize == mCapacity ? RevertState::CanRevertFull : RevertState::CanRevert ;

            // Update size (cap at N)
            mSize = std::min( mSize+1, mCapacity );
        }

//------------------------------------------------------------------------------

        void
        revert()
        {
            BELFEM_ERROR( mCanRevert != RevertState::None, "Cannot revert" );

            std::move( mData + 1, mData + mSize + ( mCanRevert == RevertState::CanRevertFull ? 1 : 0 ), mData );

            if ( mCanRevert == RevertState::CanRevert ) --mSize;

            mCanRevert = RevertState::None ;
        }
//------------------------------------------------------------------------------

        /**
         * Access element by index (0 = newest, N-1 = oldest)
         */
        T &
        operator()( const uint aIndex )
        {
            BELFEM_ASSERT( aIndex < mSize,
                "Index %u out of range (expect < %u)",
                (unsigned int) aIndex, (unsigned int) mSize );

            return mData[ aIndex ];
        }

//------------------------------------------------------------------------------

        /**
         * Access element by index (const version)
         */
        const T &
        operator()( const uint aIndex ) const
        {
            BELFEM_ASSERT( aIndex < mSize,
                "Index %u out of range (expect < %u)",
                (unsigned int) aIndex, (unsigned int) mSize );

            return mData[ aIndex ];
        }

//------------------------------------------------------------------------------

        /**
         * Get current number of elements
         */
        auto
        size() const -> decltype( mSize )
        {
            return mSize;
        }

//------------------------------------------------------------------------------

        /**
         * Get maximum capacity
         */
        auto
        capacity() const -> decltype( mCapacity )
        {
            return mCapacity ;
        }

//------------------------------------------------------------------------------

        /**
         * Check if the register is empty
         */
        bool
        empty() const
        {
            return mSize == 0;
        }

//------------------------------------------------------------------------------

        /**
         * Check if the register is full
         */
        bool
        full() const
        {
            return mSize == mCapacity ;
        }

//------------------------------------------------------------------------------

        /**
         * Clear all elements
         */
        void
        clear()
        {
            mSize = 0;
            mCanRevert = RevertState::None ;
        }

//------------------------------------------------------------------------------

        /**
         * Get raw pointer to data array
         */
        T *
        data()
        {
            return mData;
        }

//------------------------------------------------------------------------------

        /**
         * Get raw pointer to data array (const version)
         */
        const T *
        data() const
        {
            return mData;
        }

//------------------------------------------------------------------------------

        T *
        begin()
        {
            return mData;
        }

//------------------------------------------------------------------------------

        T *
        end()
        {
            return mData + mSize;
        }

//------------------------------------------------------------------------------

        const T *
        begin() const
        {
            return mData;
        }

//------------------------------------------------------------------------------

        const T*
        end() const
        {
            return mData + mSize;
        }

        void
        fill( const T aValue )
        {
            std::fill( mData, mData + mSize, aValue );
            mCanRevert = RevertState::None ;
        }

        void
        reserve( const index_t aCapacity )
        {
            // a zero-capacity register never allocates the backup slot, so a
            // later push() would write through nullptr — disallow it outright
            BELFEM_ERROR( aCapacity > 0,
                "ShiftRegister: capacity must be greater than zero" );

            if ( mCapacity == aCapacity ) return ;

            // tear down the existing buffer (destroy owning slots first)
            if ( mData != nullptr )
            {
                this->destroy_slots();
                std::free( mData );
            }

            // reserve one extra for reversion backup
            mCapacity = aCapacity ;
            mData = ( T * ) std::malloc( ( aCapacity + 1 ) * sizeof( T ) );

            // construct_slots() dereferences mData immediately, so a failed
            // allocation must be caught here rather than becoming UB
            BELFEM_ERROR( mData != nullptr,
                "ShiftRegister: failed to allocate %lu slots",
                ( unsigned long ) ( aCapacity + 1 ) );

            // bring every slot to a valid, constructed state
            this->construct_slots();

            mSize = 0 ;
            mCanRevert = RevertState::None ;
        }

        void
        free()
        {
            if ( mData != nullptr )
            {
                this->destroy_slots();
                std::free( mData );
                mData = nullptr ;
                mSize = 0 ;
                mCapacity = 0 ;
                mCanRevert = RevertState::None ;
            }
        }

//------------------------------------------------------------------------------
    };

}
#endif // BELFEM_CL_SHIFTREGISTER_HPP
