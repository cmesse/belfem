//
// Created by Christian Messe on 24.08.20.
//

#ifndef BELFEM_FN_GM_HELMHOLZ_DERIVTEST_HPP
#define BELFEM_FN_GM_HELMHOLZ_DERIVTEST_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_GM_Helmholtz.hpp"

namespace belfem
{
    namespace gasmodels
    {
        void
        deriv_test( Helmholtz & aGas, Vector< real > & aR2 );

    }
}
#endif //BELFEM_FN_GM_HELMHOLZ_DERIVTEST_HPP
