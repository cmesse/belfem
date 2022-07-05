//
// Created by Christian Messe on 09.12.19.
//

#ifndef BELFEM_FN_POLYFIT_HPP
#define BELFEM_FN_POLYFIT_HPP

#ifdef BELFEM_ARMADILLO
#include "fn_AR_polyfit.hpp"
#elif  BELFEM_BLAZE
#include "fn_BZ_polyfit.hpp"
#endif

#endif //BELFEM_FN_POLYFIT_HPP
