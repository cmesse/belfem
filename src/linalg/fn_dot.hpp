//
// Created by Christian Messe on 2019-02-07.
//

#ifndef BELFEM_FN_DOT_HPP
#define BELFEM_FN_DOT_HPP

#ifdef BELFEM_ARMADILLO
#include "fn_AR_dot.hpp"
#elif  BELFEM_BLAZE
#include "fn_BZ_dot.hpp"
#endif

#endif //BELFEM_FN_DOT_HPP
