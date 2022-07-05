//
// Created by Christian Messe on 2019-04-07.
//

#ifndef BELFEM_FN_GESV_HPP
#define BELFEM_FN_GESV_HPP

#ifdef BELFEM_ARMADILLO
#include "fn_AR_gesv.hpp"
#elif BELFEM_BLAZE
#include "fn_BZ_gesv.hpp"
#endif


#endif //BELFEM_FN_GESV_HPP
