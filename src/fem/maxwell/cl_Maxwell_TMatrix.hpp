/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef CL_MAXWELL_TMATRIX_HPP
#define CL_MAXWELL_TMATRIX_HPP
#include "cl_FEM_SideSet.hpp"
#include "cl_FEM_Calculator.hpp"

namespace belfem
{
    namespace fem
    {
        namespace maxwell
        {
            class TMatrix
            {
                SideSet    * mGroup ;
                Calculator * mCalc ;

                Matrix< real > mNabla ;
                Matrix< uint > mIndexLookup ;
                Matrix< uint > mFacetLookup ;

                Matrix< real > mResult10 ;
                Vector< real > mResult6 ;
                Vector< real > mE1 ;
                Vector< real > mE2 ;
           public:

                TMatrix( Mesh * aMesh );

                ~TMatrix();

                const Vector< real > &
                process( mesh::Facet * aFacet );

            private:

                void
                compute_nabla( const uint aIndex );

                void
                compute_nedelec_function( const uint aI, const uint aJ, const uint aK, const uint aIndex );

            };
        }
    }
}
#endif //CL_MAXWELL_TMATRIX_HPP
