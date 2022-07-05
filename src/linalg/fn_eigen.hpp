//
// Created by christian on 8/2/21.
//

#ifndef BELFEM_FN_EIGEN_HPP
#define BELFEM_FN_EIGEN_HPP

#ifdef BELFEM_ARMADILLO
#include "fn_AR_eigen.hpp"
#elif  BELFEM_BLAZE
#include "fn_BZ_eigen.hpp"
#endif

#endif //BELFEM_FN_EIGEN_HPP
