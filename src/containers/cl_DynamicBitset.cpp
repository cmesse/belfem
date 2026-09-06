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

#include <algorithm>
#include <cstring>
#include <string>
#include <sstream>
#include <iomanip>
#include <cctype> // for std::toupper
#if __cplusplus >= 202002L
#include <bit>   // for std::popcount (C++20)
#endif

#include "cl_DynamicBitset.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    namespace
    {
        //! position of the least significant set bit. The block must not be
        //! zero - every caller tests that in its loop condition.
        inline int
        ctz64( const uint64_t aBlock )
        {
#if __cplusplus >= 202002L
            return std::countr_zero( aBlock );
#elif defined(__GNUC__) || defined(__clang__)
            return __builtin_ctzll( aBlock );
#elif defined(_MSC_VER)
            unsigned long tIndex;
            _BitScanForward64( &tIndex, aBlock );
            return static_cast< int >( tIndex );
#else
#   error "DynamicBitset requires a count-trailing-zeros intrinsic"
#endif
        }

        //! allocates a zeroed block of aNumWords 64-bit words, or nullptr when
        //! aNumWords is zero ( malloc(0) may legally return null )
        inline uint64_t *
        alloc_words( const index_t aNumWords )
        {
            if ( aNumWords == 0 )
            {
                return nullptr ;
            }

            uint64_t * aBlock = ( uint64_t * ) malloc( aNumWords * sizeof( uint64_t ) );

            BELFEM_ERROR( aBlock != nullptr,
                "Failed to allocate %lu bytes for bitset",
                ( long unsigned int ) ( aNumWords * sizeof( uint64_t ) ) );

            std::memset( aBlock, 0, aNumWords * sizeof( uint64_t ) );

            return aBlock ;
        }
    }

//------------------------------------------------------------------------------

    DynamicBitset::DynamicBitset( const index_t aNumberOfBits ) :
         mNumberOfBits( aNumberOfBits ),
         mMemorySize( ( aNumberOfBits + 63) / 64 )
    {
          mSummary1Size = ( mMemorySize   + 63 ) / 64 ;
          mSummary2Size = ( mSummary1Size + 63 ) / 64 ;

          // allocate zeroed; note that reset() cannot be used here because it
          // walks the summaries, which do not exist yet
          mData     = alloc_words( mMemorySize );
          mSummary1 = alloc_words( mSummary1Size );
          mSummary2 = alloc_words( mSummary2Size );

          this->select_to_int_function();
    }

//------------------------------------------------------------------------------

    DynamicBitset::DynamicBitset( const DynamicBitset & aBitset ) :
            mNumberOfBits( aBitset.mNumberOfBits ),
            mMemorySize( aBitset.mMemorySize ),
            mSummary1Size( aBitset.mSummary1Size ),
            mSummary2Size( aBitset.mSummary2Size ),
            mIndex( aBitset.mIndex ),
            mHash( aBitset.mHash )
    {
        // zero-size bitsets are valid and skip the allocation
        mData     = alloc_words( mMemorySize );
        mSummary1 = alloc_words( mSummary1Size );
        mSummary2 = alloc_words( mSummary2Size );

        // the source summaries are tight, so copying them keeps them tight
        if ( mMemorySize > 0 )
        {
            std::memcpy( mData, aBitset.mData, mMemorySize * sizeof( uint64_t ) );
            std::memcpy( mSummary1, aBitset.mSummary1, mSummary1Size * sizeof( uint64_t ) );
            std::memcpy( mSummary2, aBitset.mSummary2, mSummary2Size * sizeof( uint64_t ) );
        }
        this->select_to_int_function();
    }

//------------------------------------------------------------------------------

    DynamicBitset::DynamicBitset( DynamicBitset && aBitset ) noexcept :
            mNumberOfBits( aBitset.mNumberOfBits ),
            mMemorySize( aBitset.mMemorySize ),
            mData( aBitset.mData ),
            mSummary1( aBitset.mSummary1 ),
            mSummary1Size( aBitset.mSummary1Size ),
            mSummary2( aBitset.mSummary2 ),
            mSummary2Size( aBitset.mSummary2Size ),
            mIndex( aBitset.mIndex ),
            mHash( aBitset.mHash ),
            mFunToInt( aBitset.mFunToInt )
    {
        aBitset.mNumberOfBits = 0;
        aBitset.mMemorySize = 0;
        aBitset.mData = nullptr;
        aBitset.mSummary1 = nullptr;
        aBitset.mSummary1Size = 0;
        aBitset.mSummary2 = nullptr;
        aBitset.mSummary2Size = 0;
        aBitset.mIndex = gNoIndex;
        aBitset.mHash = 0;
        aBitset.mFunToInt = nullptr;
    }

//------------------------------------------------------------------------------

    DynamicBitset::~DynamicBitset()
    {
          free( mData );
          free( mSummary1 );
          free( mSummary2 );
    }

//------------------------------------------------------------------------------

    void
    DynamicBitset::select_to_int_function()
    {
        index_t tNumMaxBits = std::min( sizeof( uint64_t ) * 8, sizeof( index_t ) * 8 );

        if ( mMemorySize == 0 )
        {
            // to_int_partial would read mData[ 0 ], which is null here
            mFunToInt = & DynamicBitset::to_int_zero ;
        }
        else if ( mNumberOfBits <  tNumMaxBits )
        {
            mFunToInt = & DynamicBitset::to_int_partial ;
        }
        else if ( mNumberOfBits == tNumMaxBits )
        {
            mFunToInt = & DynamicBitset::to_int_full ;
        }
        else
        {
            mFunToInt = & DynamicBitset::to_int_fail ;
        }
    }


//------------------------------------------------------------------------------

    index_t
    DynamicBitset::to_int_fail() const
    {
        BELFEM_ERROR( false, "Bitset is too big to be cast to index_t" );
        return gNoIndex ;
    }

//------------------------------------------------------------------------------

    index_t
    DynamicBitset::count() const
    {
        index_t aSum = 0;

        for (index_t i = 0; i < mMemorySize; ++i)
        {
            uint64_t tBlock = mData[i];

            // For the last block, mask out bits beyond mNumberOfBits
            if (i == mMemorySize - 1)
            {
                index_t bitsInLastBlock = mNumberOfBits % 64;
                if (bitsInLastBlock != 0)
                {
                    uint64_t mask = (uint64_t(1) << bitsInLastBlock) - 1;
                    tBlock &= mask;
                }
                // If bitsInLastBlock == 0, the last block is fully used; no masking needed
            }

#if __cplusplus >= 202002L
            // Use C++20 std::popcount if available
            aSum += std::popcount(tBlock);
#elif defined(__GNUC__) || defined(__clang__)
            // Use compiler built-in function for GCC/Clang
            aSum += __builtin_popcountll(tBlock);
#elif defined(_MSC_VER)
            // Use compiler intrinsic for MSVC
            aSum += __popcnt64(tBlock);
#else
            // Portable method using Brian Kernighan's algorithm
            while (tBlock)
            {
                tBlock &= (tBlock - 1);
                ++aSum;
            }
#endif
        }
        return aSum;
    }

//------------------------------------------------------------------------------

    string
    DynamicBitset::to_string() const
    {
        string aString ;
        aString.reserve( mNumberOfBits );
        for( index_t k=1; k<=mNumberOfBits; ++k )
        {
            if( this->test( mNumberOfBits - k ) )
            {
                aString.push_back( '1' );
            }
            else
            {
                aString.push_back( '0' );
            }
        }
        return aString ;
    }


//------------------------------------------------------------------------------

    std::string DynamicBitset::to_hex() const
    {
        if (mNumberOfBits == 0)
        {
            return "";
        }

        std::ostringstream tStream;

        // Calculate the total number of hex digits
        index_t tNumHexDigits = (mNumberOfBits + 3) / 4;

        // Process each 64-bit block in reverse order
        for (index_t b = 1; b <= mMemorySize; ++b)
        {
            uint64_t tBlock = mData[mMemorySize - b];

            // Determine how many bits are in this block
            index_t tNumBitsInBlock = 64;
            if (b == 1)
            {
                // Adjust for the last block, which may not be fully used
                tNumBitsInBlock = mNumberOfBits % 64;
                if (tNumBitsInBlock == 0)
                {
                    tNumBitsInBlock = 64;
                }

                // Mask the block to only include valid bits
                if (tNumBitsInBlock < 64)
                {
                    uint64_t mask = (uint64_t(1) << tNumBitsInBlock) - 1;
                    tBlock &= mask;
                }
            }

            // Calculate the number of hex digits in this block
            index_t tNumHexDigitsInBlock = (tNumBitsInBlock + 3) / 4;

            // Format the block as hexadecimal
            tStream << std::hex << std::setw(tNumHexDigitsInBlock) << std::setfill('0') << tBlock;
        }

        // Get the hexadecimal string
        std::string aHexStr = tStream.str();

        // If necessary, trim excess digits
        if (aHexStr.length() > tNumHexDigits)
        {
            aHexStr = aHexStr.substr(aHexStr.length() - tNumHexDigits);
        }

        // Convert to uppercase
        for (char& c : aHexStr)
        {
            c = std::toupper(static_cast<unsigned char>(c));
        }

        return aHexStr;
    }

//------------------------------------------------------------------------------

    std::string
    DynamicBitset::to_raw_string() const
    {
        // Create a string with enough capacity
        std::string aResult;
        aResult.reserve(mMemorySize * sizeof(uint64_t));

        for (index_t i = 0; i < mMemorySize; ++i)
        {
            // Directly append the raw bytes of mData[i]
            aResult.append(reinterpret_cast<const char*>(&mData[i]), sizeof(uint64_t));
        }

        return aResult ;
    }

//------------------------------------------------------------------------------

    void
    DynamicBitset::set_from_hex( const string & aString )
    {
        // Clear the current bitset
        this->reset();

        // Calculate the number of hex digits available in the input string
        size_t hexLength = aString.length();

        // Process each hex digit
        for (size_t i = 0; i < hexLength; ++i)
        {
            // Get the current hex character and convert to an integer value
            char hexChar = aString[hexLength - 1 - i]; // Process from the last character
            unsigned int hexValue;

            // Convert hex character to numeric value
            if (hexChar >= '0' && hexChar <= '9')
            {
                hexValue = hexChar - '0';
            }
            else if (hexChar >= 'A' && hexChar <= 'F')
            {
                hexValue = hexChar - 'A' + 10;
            }
            else if (hexChar >= 'a' && hexChar <= 'f')
            {
                hexValue = hexChar - 'a' + 10;
            }
            else
            {
                BELFEM_ERROR( false, "Invalid hexadecimal string" );
                hexValue = 0;
            }

            // Set the corresponding bits for the current hex value (4 bits for each hex digit)
            for (int j = 0; j < 4; ++j)
            {
                if (hexValue & (1 << j))
                {
                    size_t pos = (i * 4) + j;
                    if (pos < mNumberOfBits)
                    {
                        this->set(pos);
                    }
                }
            }
        }
    }

    void
    DynamicBitset::reset()
    {
        this->unlock() ;

        // walk the summaries and zero only the words that carry bits, then
        // zero the summaries behind us
        for ( index_t b = 0; b < mSummary2Size; ++b )
        {
            uint64_t tLevel2 = mSummary2[ b ];

            while ( tLevel2 != 0 )
            {
                const index_t tS1 = b * 64 + ctz64( tLevel2 );

                uint64_t tLevel1 = mSummary1[ tS1 ];

                while ( tLevel1 != 0 )
                {
                    mData[ tS1 * 64 + ctz64( tLevel1 ) ] = 0 ;
                    tLevel1 &= tLevel1 - 1 ;
                }

                mSummary1[ tS1 ] = 0 ;
                tLevel2 &= tLevel2 - 1 ;
            }

            mSummary2[ b ] = 0 ;
        }
    }

//------------------------------------------------------------------------------

    void
    DynamicBitset::recompute_summaries()
    {
        if ( mSummary1Size > 0 )
        {
            std::memset( mSummary1, 0, mSummary1Size * sizeof( uint64_t ) );
        }
        if ( mSummary2Size > 0 )
        {
            std::memset( mSummary2, 0, mSummary2Size * sizeof( uint64_t ) );
        }

        for ( index_t w = 0; w < mMemorySize; ++w )
        {
            if ( mData[ w ] != 0 )
            {
                mSummary1[ w / 64 ] |= uint64_t( 1 ) << ( w % 64 );
            }
        }

        for ( index_t s = 0; s < mSummary1Size; ++s )
        {
            if ( mSummary1[ s ] != 0 )
            {
                mSummary2[ s / 64 ] |= uint64_t( 1 ) << ( s % 64 );
            }
        }
    }

//------------------------------------------------------------------------------

    bool
    DynamicBitset::summaries_are_tight() const
    {
        for ( index_t w = 0; w < mMemorySize; ++w )
        {
            const bool tBitSet = ( mSummary1[ w / 64 ]
                                   & ( uint64_t( 1 ) << ( w % 64 ) ) ) != 0 ;

            if ( tBitSet != ( mData[ w ] != 0 ) )
            {
                return false ;
            }
        }

        for ( index_t s = 0; s < mSummary1Size; ++s )
        {
            const bool tBitSet = ( mSummary2[ s / 64 ]
                                   & ( uint64_t( 1 ) << ( s % 64 ) ) ) != 0 ;

            if ( tBitSet != ( mSummary1[ s ] != 0 ) )
            {
                return false ;
            }
        }

        // the safety rule: no summary bit may stand beyond its child array.
        // Both loops above already cover every in-range index, so only the
        // tail of the last summary word can still be dirty.
        for ( index_t w = mMemorySize; w < mSummary1Size * 64; ++w )
        {
            if ( ( mSummary1[ w / 64 ] & ( uint64_t( 1 ) << ( w % 64 ) ) ) != 0 )
            {
                return false ;
            }
        }

        for ( index_t s = mSummary1Size; s < mSummary2Size * 64; ++s )
        {
            if ( ( mSummary2[ s / 64 ] & ( uint64_t( 1 ) << ( s % 64 ) ) ) != 0 )
            {
                return false ;
            }
        }

        return true ;
    }

//------------------------------------------------------------------------------

    void
    DynamicBitset::where_dense( Cell< index_t > & aBits ) const
    {
        // the summary walk never visits a zero word, so the old
        // count-then-fill variant has nothing left to win over the sparse one
        this->where_sparse( aBits );
    }

//------------------------------------------------------------------------------

    void
    DynamicBitset::where_sparse( Cell< index_t > & aBits ) const
    {
        aBits.clear(); // Reset the container ( retains capacity )

        // ascending level-2 -> ascending level-1 -> ascending word ->
        // ascending bit, so the output is sorted by construction
        for ( index_t b = 0; b < mSummary2Size; ++b )
        {
            uint64_t tLevel2 = mSummary2[ b ];

            while ( tLevel2 != 0 )
            {
                const index_t tS1 = b * 64 + ctz64( tLevel2 );

                uint64_t tLevel1 = mSummary1[ tS1 ];

                while ( tLevel1 != 0 )
                {
                    const index_t tWord = tS1 * 64 + ctz64( tLevel1 );

                    uint64_t tBlock = mData[ tWord ];

                    while ( tBlock != 0 )
                    {
                        const index_t tIndex = tWord * 64 + ctz64( tBlock );

                        // bits past mNumberOfBits are always zero: the
                        // constructor zeroes the array, set() and flip( pos )
                        // assert their bounds, flip() masks the last block and
                        // set_from_hex() guards the position - so this cannot
                        // fire for a caller that respects the contract
                        BELFEM_ASSERT( tIndex < mNumberOfBits,
                            "Bit %lu is set beyond the end of the bitset (size %lu)",
                            ( long unsigned int ) tIndex,
                            ( long unsigned int ) mNumberOfBits );

                        aBits.push( tIndex );

                        tBlock &= tBlock - 1 ;
                    }

                    tLevel1 &= tLevel1 - 1 ;
                }

                tLevel2 &= tLevel2 - 1 ;
            }
        }
    }

//------------------------------------------------------------------------------

    void
    DynamicBitset::flip()
    {
        BELFEM_ASSERT(mHash == 0, "Can't flip bits on a locked bitset");

        // index_t is unsigned, so mMemorySize - 1 wraps for an empty bitset
        // and the loop below would march through a null mData
        if ( mMemorySize == 0 )
        {
            return;
        }

        // Flip all complete blocks
        for (index_t i = 0; i < mMemorySize - 1; ++i)
        {
            mData[i] = ~mData[i];
        }

        // Handle the last block specially to avoid flipping unused bits.
        // The empty case already returned above, so no size guard is needed.
        index_t tBitsInLastBlock = mNumberOfBits % 64;
        if (tBitsInLastBlock == 0)
        {
            // Last block is fully used
            mData[mMemorySize - 1] = ~mData[mMemorySize - 1];
        }
        else
        {
            // Create a mask for valid bits in the last block
            uint64_t mask = (uint64_t(1) << tBitsInLastBlock) - 1;
            mData[mMemorySize - 1] = (~mData[mMemorySize - 1]) & mask;
        }

        // a flip inverts nearly every word, so deriving the summaries from
        // their old state buys nothing - rebuild them. This is O( mMemorySize )
        // on top of an operation that is already O( mMemorySize ).
        this->recompute_summaries();
    }
//------------------------------------------------------------------------------

        //! Assignment operator
        DynamicBitset& DynamicBitset::operator=(const DynamicBitset& aRhs)
        {
            if ( this != &aRhs )
            {
                if ( this->memory() != aRhs.memory() )
                {
                    // Deallocate existing memory
                    free( mData );
                    free( mSummary1 );
                    free( mSummary2 );

                    // Update memory size and number of bits
                    mMemorySize = aRhs.memory();
                    mNumberOfBits = aRhs.size();
                    mSummary1Size = aRhs.mSummary1Size;
                    mSummary2Size = aRhs.mSummary2Size;

                    // Allocate new memory; zero-size bitsets skip the allocation
                    mData     = alloc_words( mMemorySize );
                    mSummary1 = alloc_words( mSummary1Size );
                    mSummary2 = alloc_words( mSummary2Size );
                }
                else
                {
                    // Update number of bits even if memory size is the same
                    mNumberOfBits = aRhs.size();
                }

                // Copy data from the right-hand side bitset. Its summaries are
                // tight, so the copies are tight too.
                if ( mMemorySize > 0 )
                {
                    std::memcpy(mData, aRhs.data(), mMemorySize * sizeof(uint64_t));
                    std::memcpy(mSummary1, aRhs.mSummary1, mSummary1Size * sizeof(uint64_t));
                    std::memcpy(mSummary2, aRhs.mSummary2, mSummary2Size * sizeof(uint64_t));
                }

                // Copy index and hash ( same members the copy constructor copies )
                mIndex = aRhs.mIndex;
                mHash = aRhs.is_locked() ? aRhs.hash() : 0;

                // update int function selection
                this->select_to_int_function();
            }

            return *this;
        }

//------------------------------------------------------------------------------

        //! Move assignment operator
        DynamicBitset& DynamicBitset::operator=(DynamicBitset&& aRhs) noexcept
        {
            if ( this != &aRhs )
            {
                // Deallocate existing memory
                free( mData );
                free( mSummary1 );
                free( mSummary2 );

                // Move members
                mNumberOfBits = aRhs.mNumberOfBits;
                mMemorySize = aRhs.mMemorySize;
                mData = aRhs.mData;
                mSummary1 = aRhs.mSummary1;
                mSummary1Size = aRhs.mSummary1Size;
                mSummary2 = aRhs.mSummary2;
                mSummary2Size = aRhs.mSummary2Size;
                mIndex = aRhs.mIndex;
                mHash = aRhs.mHash;
                mFunToInt = aRhs.mFunToInt;

                // Reset the source object to a valid state
                aRhs.mNumberOfBits = 0;
                aRhs.mMemorySize = 0;
                aRhs.mData = nullptr;
                aRhs.mSummary1 = nullptr;
                aRhs.mSummary1Size = 0;
                aRhs.mSummary2 = nullptr;
                aRhs.mSummary2Size = 0;
                aRhs.mIndex = gNoIndex;
                aRhs.mHash = 0;
                aRhs.mFunToInt = nullptr;
            }

            return *this;
        }
//------------------------------------------------------------------------------
}
