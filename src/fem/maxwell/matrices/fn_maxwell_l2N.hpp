//
// Created by christian on 10/24/24.
//

#ifndef FN_MAXWELL_L2N_HPP
#define FN_MAXWELL_L2N_HPP

#include "typedefs.hpp"
#include "cl_Matrix.hpp"

#include "cl_FEM_Calculator.hpp"

namespace belfem
{
    namespace fem
    {
        namespace maxwell
        {
            inline const Matrix< real > &
            l2N( Calculator * aCalc, const uint aIndex )
            {
                const Vector< real > & n = aCalc->Nvec( aIndex ) ;
                Matrix< real > & N = aCalc->matrix("L2N");

                uint nnodes = n.length();
                uint ndim = N.n_rows() ;
                for( uint j=0; j<ndim; ++j )
                {
                    uint c = j ;
                    for( uint i=0; i<nnodes; ++i )
                    {
                        N( j, c ) = n( i );
                        c+=ndim;
                    }
                }

                return N;
            }
        }
    }
}
#endif //FN_MAXWELL_L2N_HPP
