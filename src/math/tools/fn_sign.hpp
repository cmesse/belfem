//
// Created by Christian Messe on 13.09.19.
//

#ifndef BELFEM_FN_SIGN_HPP
#define BELFEM_FN_SIGN_HPP

namespace belfem
{
    template < typename T >
    inline T sign( const T aX )
    {
        return ( aX > 0) ? 1 : ( ( aX < 0) ? -1 : 0);
    }

}

#endif //BELFEM_FN_SIGN_HPP
