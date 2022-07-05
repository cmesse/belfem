//
// Created by Christian Messe on 2019-07-02.
//

#ifndef BELFEM_FN_LINSPACE_HPP
#define BELFEM_FN_LINSPACE_HPP

#ifdef BELFEM_ARMADILLO
#include "fn_AR_linspace.hpp"
#elif  BELFEM_BLAZE
#include "fn_BZ_linspace.hpp"
#endif

#endif //BELFEM_FN_LINSPACE_HPP
