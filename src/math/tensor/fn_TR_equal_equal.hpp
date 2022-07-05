//
// Created by Christian Messe on 17.12.20.
//

#ifndef BELFEM_FN_TR_EQUAL_EQUAL_HPP
#define BELFEM_FN_TR_EQUAL_EQUAL_HPP

namespace belfem
{
    namespace tensor
    {
//----------------------------------------------------------------------------

        template < typename T >
        inline bool
        equal_equal( const T * A, const T * B, const index_t aCapacity )
        {
            for( index_t k=0; k<aCapacity; ++k )
            {
                if( A[ k ] != B[ k ] )
                {
                    return false ;
                }
            }

            return true ;
        }

//----------------------------------------------------------------------------
    } /* namespace tensor */
} /* namespace belfem */
#endif //BELFEM_FN_TR_EQUAL_EQUAL_HPP
