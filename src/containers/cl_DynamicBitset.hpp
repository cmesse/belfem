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

#ifndef CL_DYNAMICBITSET_HPP
#define CL_DYNAMICBITSET_HPP

#include <cstdint>  // For uint64_t
#include <cstring>  // For std::memcpy

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Cell.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    /**
     * @brief Runtime-sized bitset; one bit per flag, packed into 64-bit words.
     *
     * @ingroup grp_containers
     * @see @ref containers_container_usage_guide
     */
    class DynamicBitset
    {
        //! Number of bits in the bitset
        index_t mNumberOfBits;

        //! Number of 64-bit blocks required to store the bits
        index_t mMemorySize;

        //! Pointer to the array of 64-bit blocks
        uint64_t* mData;

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // two-level summary bitmap
        //
        // mSummary1 holds one bit per data word, mSummary2 one bit per
        // level-1 word. where() and reset() walk the summaries instead of
        // scanning mData, which turns both from O( mNumberOfBits / 64 ) into
        // O( mSummary2Size + number of touched words ). For a 2.1e6-bit
        // workspace that is 9 + O( k ) words instead of 32813.
        //
        // Correctness invariant ( one-directional ): every nonzero data word
        // has its level-1 bit set, and every nonzero level-1 word has its
        // level-2 bit set. A summary bit standing over a zero word is skipped
        // harmlessly by the extraction loops.
        //
        // Tightness invariant ( maintained by every mutator, asserted by
        // summaries_are_tight() ): a summary bit is set if and only if the
        // word below it is nonzero. Tightness is not needed for correctness,
        // but losing it silently degrades where() back towards a full scan,
        // so the tests pin it.
        //
        // Safety rule: no level-1 bit may stand for a word index >=
        // mMemorySize, and no level-2 bit for a level-1 index >=
        // mSummary1Size, or the walks read past the arrays.
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        //! level-1 summary: one bit per data word
        uint64_t* mSummary1 = nullptr ;

        //! number of 64-bit words in the level-1 summary
        index_t mSummary1Size = 0 ;

        //! level-2 summary: one bit per level-1 word
        uint64_t* mSummary2 = nullptr ;

        //! number of 64-bit words in the level-2 summary
        index_t mSummary2Size = 0 ;

        //! index of bitset, can also be set when bitset is locked
        index_t mIndex = gNoIndex ;

        // if hash is 0, we make the bitset writable
        // once the hash is computed, we lock the bitset to prevent errors
        size_t mHash = 0 ;

        // Function pointer for optimized to_integer conversion
        index_t ( DynamicBitset::*mFunToInt )() const ;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        //! Constructor: Initializes the bitset with the given number of bits
        DynamicBitset(const index_t aNumberOfBits );

//------------------------------------------------------------------------------

        //! copy constructor
        DynamicBitset( const DynamicBitset & aBitset );

//------------------------------------------------------------------------------

        //! move constructor
        DynamicBitset( DynamicBitset && aBitset ) noexcept ;

//------------------------------------------------------------------------------

        //! Destructor: Releases allocated memory
        ~DynamicBitset();

//------------------------------------------------------------------------------

        //! Returns the number of bits in the bitset
        index_t
        size() const;

//------------------------------------------------------------------------------

        //! Returns the number of 64-bit blocks used (memory size)
        index_t
        memory() const;

//------------------------------------------------------------------------------

        //! Resets (clears) the bit at the given position
        void
        reset(const index_t aPos);

//------------------------------------------------------------------------------

        //! Resets (clears) all bits
        void
        reset() ;

//------------------------------------------------------------------------------
        //! Sets the bit at the given position to 1
        void
        set(const index_t aPos);

//------------------------------------------------------------------------------

        //! Sets the bit at the given position to the specified value
        void
        set(const index_t aPos, const bool aValue);

//------------------------------------------------------------------------------

        //! Flips (toggles) the bit at the given position
        void
        flip(const index_t aPos);

//------------------------------------------------------------------------------

        //! Flips all bits at once
        void
        flip();

//------------------------------------------------------------------------------

        //! Tests whether the bit at the given position is set
        bool
        test(const index_t aPos) const;

//------------------------------------------------------------------------------

        //! Counts the number of bits set to 1
        index_t
        count() const;

//------------------------------------------------------------------------------

        //! Returns a pointer to the data array (const)
        //!
        //! @note there is deliberately no mutable overload: a write through
        //!       one would bypass the summary bitmaps, after which where() and
        //!       reset() silently miss the affected words. Mutate through
        //!       set() / reset() / flip() instead.
        const uint64_t *
        data() const;

//------------------------------------------------------------------------------

        //! Diagnostic: checks that both summary levels are tight, i.e. that a
        //! summary bit is set exactly when the word below it is nonzero, and
        //! that no summary bit stands beyond its child array. Intended for
        //! tests and debugging; costs one pass over the data.
        bool
        summaries_are_tight() const ;

//------------------------------------------------------------------------------

        //! makes the bitset non-writable and computes the hash
        void
        lock() ;

//------------------------------------------------------------------------------

        //! makes the bitset writable and resets the hash
        void
        unlock() ;

//------------------------------------------------------------------------------

        //! checks if the bitset is writable
        bool
        is_locked() const ;

//------------------------------------------------------------------------------

        //! Returns a hash function for fast comparison
        size_t hash() const;

//------------------------------------------------------------------------------

        string
        to_string() const ;

//------------------------------------------------------------------------------

        string
        to_hex() const ;

//------------------------------------------------------------------------------

        index_t
        to_int() const ;

//------------------------------------------------------------------------------

        string
        to_raw_string() const ;

//------------------------------------------------------------------------------

        void
        set_from_hex( const string & aString );

//------------------------------------------------------------------------------

        void
        set_index( const index_t aIndex );

//------------------------------------------------------------------------------

        index_t
        index() const ;

//------------------------------------------------------------------------------

        // Returns the indices of the set bits, strictly ascending and free of
        // duplicates. The scan walks the summary bitmaps, so it never visits a
        // zero data word: cost is O( level-2 words + touched words + set bits ),
        // i.e. essentially the number of set bits plus a size/262144 term.
        //
        // aAssumeSparse is retained for source compatibility and no longer
        // selects a different algorithm - there is no longer a dense variant
        // worth having, because the walk never visits a zero word.
        void
        where( Cell< index_t > & aBits, const bool aAssumeSparse = true ) const ;

//------------------------------------------------------------------------------

        //! Comparison operator: Checks if two bitsets are equal.
        //! Both bitsets must be locked: a debug build asserts on an unlocked
        //! operand, a release build compares the words regardless.
        //! A size mismatch is an error in every build.
        bool operator==(const DynamicBitset & aRhs) const
        {
            BELFEM_ERROR(this->size() == aRhs.size(),
                "Bitsets don't have the same size (%lu vs. %lu)",
                (long unsigned int)this->size(),
                (long unsigned int) aRhs.size());

            if ( this->hash() != aRhs.hash() )
            {
                return false;  // Quick reject if hashes differ
            }

            // Full comparison to resolve collisions
            for (index_t k = 0; k < mMemorySize; ++k)
            {
                if (mData[k] != aRhs.data()[k])
                {
                    return false;
                }
            }
            return true;
        }

//------------------------------------------------------------------------------

        //! Comparison operator: Checks if two bitsets are not equal
        bool operator!=(const DynamicBitset& aRhs) const
        {
            return !(*this == aRhs);
        }

//------------------------------------------------------------------------------

        //! Assignment operator
        DynamicBitset& operator=(const DynamicBitset& aRhs ) ;

        //! Move assignment operator
        DynamicBitset& operator=(DynamicBitset&& aRhs) noexcept ;

//------------------------------------------------------------------------------
    private:

        // Note: there is no mutable data() overload. A write through one would
        // bypass the summary bitmaps, after which where() and reset() silently
        // miss the affected words. Members use mData directly.

        //! Rebuilds both summary levels from mData in one pass, producing
        //! tight summaries by construction. Called from the bulk mutators
        //! ( flip(), ^=, &= ), which are O( mMemorySize ) anyway.
        void
        recompute_summaries();

        void
        where_dense( Cell< index_t > & aBits ) const ;

        void
        where_sparse( Cell< index_t > & aBits ) const ;

        void
        select_to_int_function();

        index_t
        to_int_partial() const ;

        index_t
        to_int_full() const ;

        //! an empty bitset has no mData to read - dispatching here keeps the
        //! zero check off the hot path, as to_int_fail already does
        index_t
        to_int_zero() const ;

        index_t
        to_int_fail() const ;


//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        //! Bitwise OR operator
        DynamicBitset operator|(const DynamicBitset& aRhs) const
        {
            BELFEM_ERROR(this->size() == aRhs.size(),
                "Bitsets don't have the same size (%lu vs. %lu)",
                (long unsigned int)this->size(),
                (long unsigned int)aRhs.size());

            DynamicBitset aResult(mNumberOfBits);

            const uint64_t* tSrcA = mData;
            const uint64_t* tSrcB = aRhs.data();
            uint64_t* tDst = aResult.mData;
            const uint64_t* tEnd = tSrcA + mMemorySize;

            while(tSrcA != tEnd)
            {
                *tDst++ = *tSrcA++ | *tSrcB++;
            }

            // a|b is nonzero exactly where a or b is, so OR-ing tight
            // summaries yields a tight summary - no rebuild needed
            for( index_t i = 0; i < mSummary1Size; ++i )
            {
                aResult.mSummary1[ i ] = mSummary1[ i ] | aRhs.mSummary1[ i ];
            }
            for( index_t i = 0; i < mSummary2Size; ++i )
            {
                aResult.mSummary2[ i ] = mSummary2[ i ] | aRhs.mSummary2[ i ];
            }

            return aResult;
        }

        //! Assignment or operator
        DynamicBitset & operator|=( const DynamicBitset& aRhs )
        {
            BELFEM_ERROR(this->size() == aRhs.size(),
                "Bitsets don't have the same size (%lu vs. %lu)",
                (long unsigned int)this->size(),
                (long unsigned int) aRhs.size());

            BELFEM_ERROR( mHash == 0, "Can't modify a locked bitset" );

            const uint64_t* tData = aRhs.data();

            for(index_t i = 0; i < mMemorySize; ++i)
            {
                mData[i] |= tData[i];
            }

            // OR of two tight summaries is tight ( see operator| )
            for( index_t i = 0; i < mSummary1Size; ++i )
            {
                mSummary1[ i ] |= aRhs.mSummary1[ i ];
            }
            for( index_t i = 0; i < mSummary2Size; ++i )
            {
                mSummary2[ i ] |= aRhs.mSummary2[ i ];
            }

            return *this;
        }

        //! Bitwise XOR operator
        DynamicBitset operator^(const DynamicBitset& aRhs) const
        {
            BELFEM_ERROR(this->size() == aRhs.size(),
                "Bitsets don't have the same size (%lu vs. %lu)",
                (long unsigned int)this->size(),
                (long unsigned int)aRhs.size());

            DynamicBitset aResult(mNumberOfBits);

            const uint64_t* tData = aRhs.data();

            for(index_t i = 0; i < mMemorySize; ++i)
            {
                aResult.mData[i] = mData[i] ^ tData[i];
            }

            // XOR can zero a word that was nonzero in both operands, so the
            // summaries cannot be derived from the inputs - rebuild them
            aResult.recompute_summaries();

            return aResult;
        }

        //! Bitwise XOR assignment operator
        DynamicBitset& operator^=(const DynamicBitset& aRhs)
        {
            BELFEM_ERROR(this->size() == aRhs.size(),
                "Bitsets don't have the same size (%lu vs. %lu)",
                (long unsigned int)this->size(),
                (long unsigned int)aRhs.size());

            BELFEM_ERROR( mHash == 0, "Can't modify a locked bitset");

            for(index_t i = 0; i < mMemorySize; ++i)
            {
                mData[i] ^= aRhs.mData[i];
            }

            // XOR can zero words ( see operator^ )
            this->recompute_summaries();

            return *this;
        }

        //! Bitwise AND operator
        DynamicBitset operator&(const DynamicBitset& aRhs) const
        {
            BELFEM_ERROR(this->size() == aRhs.size(),
                "Bitsets don't have the same size (%lu vs. %lu)",
                (long unsigned int)this->size(),
                (long unsigned int)aRhs.size());

            DynamicBitset tResult(mNumberOfBits);

            const uint64_t* tData = aRhs.data();

            for(index_t i = 0; i < mMemorySize; ++i)
            {
                tResult.mData[i] = mData[i] & tData[i];
            }

            // AND can zero a word that was nonzero in both operands
            tResult.recompute_summaries();

            return tResult;
        }

        //! Bitwise AND assignment operator
        DynamicBitset& operator&=(const DynamicBitset& aRhs)
        {
            BELFEM_ERROR(this->size() == aRhs.size(),
                "Bitsets don't have the same size (%lu vs. %lu)",
                (long unsigned int)this->size(),
                (long unsigned int)aRhs.size());

            const uint64_t* tData = aRhs.data();

            BELFEM_ERROR( mHash == 0, "Can't modify a locked bitset");

            for(index_t i = 0; i < mMemorySize; ++i)
            {
                mData[i] &= tData[i];
            }

            // AND can zero words ( see operator& )
            this->recompute_summaries();

            return *this;
        }

    };



//------------------------------------------------------------------------------

    // Inline implementations of member functions
    inline index_t DynamicBitset::size() const
    {
        return mNumberOfBits;
    }

//------------------------------------------------------------------------------

    inline index_t DynamicBitset::memory() const
    {
        return mMemorySize;
    }

//------------------------------------------------------------------------------

    inline void DynamicBitset::set(const index_t aPos)
    {
        BELFEM_ASSERT(aPos < mNumberOfBits,
            "Index %lu out of range (expect < %lu)",
            (long unsigned int)aPos,
            (long unsigned int)mNumberOfBits);

        BELFEM_ASSERT( mHash == 0,
            "can't write on a locked bitset");

        // setting a bit always makes the word nonzero, so both summary levels
        // can be updated unconditionally - no branch on this hot path
        const index_t tWord = aPos / 64 ;
        const index_t tS1   = tWord / 64 ;

        mData    [ tWord ]      |= uint64_t( 1 ) << ( aPos  % 64 );
        mSummary1[ tS1 ]        |= uint64_t( 1 ) << ( tWord % 64 );
        mSummary2[ tS1 / 64 ]   |= uint64_t( 1 ) << ( tS1   % 64 );
    }

//------------------------------------------------------------------------------

    inline void DynamicBitset::reset(const index_t aPos)
    {
        BELFEM_ASSERT(aPos < mNumberOfBits,
            "Index %lu out of range (expect < %lu)",
            (long unsigned int)aPos,
            (long unsigned int)mNumberOfBits);

        BELFEM_ASSERT( mHash == 0,
            "can't reset an single digit of a locked bitset" );

        // clearing the last bit of a word must clear the summary bits above
        // it. Leaving them stale would still be correct, but a caller that
        // clears bits one by one instead of calling reset() would saturate
        // the summaries and lose the whole benefit of the bitmap - which is
        // exactly what the sparsity-pattern builders used to do.
        const index_t tWord = aPos / 64 ;

        mData[ tWord ] &= ~( uint64_t( 1 ) << ( aPos % 64 ) );

        if ( mData[ tWord ] == 0 )
        {
            const index_t tS1 = tWord / 64 ;

            mSummary1[ tS1 ] &= ~( uint64_t( 1 ) << ( tWord % 64 ) );

            if ( mSummary1[ tS1 ] == 0 )
            {
                mSummary2[ tS1 / 64 ] &= ~( uint64_t( 1 ) << ( tS1 % 64 ) );
            }
        }
    }


//------------------------------------------------------------------------------

    // reset() is implemented in the .cpp: it walks the summary bitmaps rather
    // than memsetting the whole array, and shares the ctz helper with the
    // other walkers.

//------------------------------------------------------------------------------

    inline void DynamicBitset::set(const index_t aPos, const bool aValue)
    {
        if (aValue)
        {
            this->set(aPos);
        }
        else
        {
            this->reset(aPos);
        }
    }

//------------------------------------------------------------------------------

    inline void DynamicBitset::flip(const index_t aPos)
    {
        BELFEM_ASSERT(aPos < mNumberOfBits,
            "Index %lu out of range (expect < %lu)",
            (long unsigned int)aPos,
            (long unsigned int)mNumberOfBits);

        BELFEM_ASSERT( mHash == 0,
            "Can't flip a bit on a locked bitset" );

        // a flip can turn a bit ON in a word whose summary bit is clear, so
        // the summaries must be maintained here for correctness, not merely
        // for tightness: without the OR branch, where() would miss the bit
        const index_t tWord = aPos / 64 ;
        const index_t tS1   = tWord / 64 ;

        mData[ tWord ] ^= uint64_t( 1 ) << ( aPos % 64 );

        if ( mData[ tWord ] != 0 )
        {
            mSummary1[ tS1 ]      |= uint64_t( 1 ) << ( tWord % 64 );
            mSummary2[ tS1 / 64 ] |= uint64_t( 1 ) << ( tS1   % 64 );
        }
        else
        {
            mSummary1[ tS1 ] &= ~( uint64_t( 1 ) << ( tWord % 64 ) );

            if ( mSummary1[ tS1 ] == 0 )
            {
                mSummary2[ tS1 / 64 ] &= ~( uint64_t( 1 ) << ( tS1 % 64 ) );
            }
        }
    }

//------------------------------------------------------------------------------

    inline bool DynamicBitset::test(const index_t aPos) const
    {
        BELFEM_ASSERT(aPos < mNumberOfBits,
            "Index %lu out of range (expect < %lu)",
            (long unsigned int)aPos,
            (long unsigned int)mNumberOfBits);

        return ( mData[aPos / 64] & (uint64_t(1) << (aPos % 64))) != 0;
    }

//------------------------------------------------------------------------------

    inline const uint64_t* DynamicBitset::data() const
    {
        return mData;
    }

//------------------------------------------------------------------------------

    inline void
    DynamicBitset::lock()
    {
        mHash = 14695981039346656037ULL;  // FNV-1a 64-bit offset basis
        mHash ^= std::hash<index_t>{}(mNumberOfBits);
        mHash *= 1099511628211ULL;       // FNV-64 multiplier
        for (index_t i = 0; i < mMemorySize; ++i)
        {
            mHash ^= std::hash<uint64_t>{}(mData[i]);
            mHash *= 1099511628211ULL;
        }

        // zero is the "unlocked" sentinel, so a hash that happens to land on
        // it would make a locked bitset report itself writable
        if ( mHash == 0 )
        {
            mHash = 1 ;
        }
    }

//------------------------------------------------------------------------------

    inline void
    DynamicBitset::unlock()
    {
        mHash = 0 ;
    }

//------------------------------------------------------------------------------

    inline bool DynamicBitset::is_locked() const
    {
        return mHash != 0 ;
    }

//------------------------------------------------------------------------------

    //! returns the hash value
    inline size_t
    DynamicBitset::hash() const
    {
        BELFEM_ASSERT( mHash != 0, "can't return the hash of a writable bitset" );
        return mHash ;
    }

//------------------------------------------------------------------------------

    inline void
    DynamicBitset::where( Cell< index_t > & aBits, const bool aAssumeSparse ) const
    {
        if ( mNumberOfBits == 0 )
        {
            aBits.clear(); // Reset the container
            return;
        }
        // aAssumeSparse no longer selects a different algorithm: since the
        // scan walks the summary bitmaps it never visits a zero word, so the
        // count-then-fill variant has nothing left to win. The parameter is
        // kept so existing call sites still compile.
        ( void ) aAssumeSparse ;

        this->where_sparse( aBits ) ;
    }

//------------------------------------------------------------------------------

    inline index_t
    DynamicBitset::to_int() const
    {
        return ( this->*mFunToInt )();
    }

//------------------------------------------------------------------------------

    inline index_t
    DynamicBitset::to_int_partial() const
    {
        uint64_t mask = (uint64_t(1) << mNumberOfBits) - 1;
        return static_cast< index_t >(mData[0] & mask);
    }

//------------------------------------------------------------------------------

    inline index_t
    DynamicBitset::to_int_full() const
    {
        return static_cast<index_t>(mData[0]);
    }

//------------------------------------------------------------------------------

    inline index_t
    DynamicBitset::to_int_zero() const
    {
        return 0;
    }

//------------------------------------------------------------------------------

    inline void
    DynamicBitset::set_index( const index_t aIndex )
    {
        mIndex = aIndex ;
    }

//------------------------------------------------------------------------------

    inline index_t
    DynamicBitset::index() const
    {
        return mIndex;
    }

//------------------------------------------------------------------------------
} // namespace belfem

#endif // CL_DYNAMICBITSET_HPP
