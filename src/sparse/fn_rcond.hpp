//
// Created by Christian Messe on 07.09.22.
//

#ifndef BELFEM_FN_RCOND_HPP
#define BELFEM_FN_RCOND_HPP

#include "cl_SpMatrix.hpp"

namespace belfem
{
    real
    rcond( SpMatrix & aMatrix );
}
#endif //BELFEM_FN_RCOND_HPP
