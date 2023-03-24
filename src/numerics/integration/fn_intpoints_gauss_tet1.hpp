//
// Created by Christian Messe on 2019-01-17.
//

#ifndef BELFEM_FN_INTPOINTS_GAUSS_TET1_HPP
#define BELFEM_FN_INTPOINTS_GAUSS_TET1_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace integration
    {
// ----------------------------------------------------------------------------

        inline void
        gauss_tet1(
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {
            aPoints.set_size( 4, 1, 0.25 );
            aWeights.set_size( 1, 1.0/6.0 );
        }

// ----------------------------------------------------------------------------
    } /* namespace integration */
} /* namespace belfem */


#endif //BELFEM_FN_INTPOINTS_GAUSS_TET1_HPP
