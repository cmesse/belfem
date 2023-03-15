//
// Created by christian on 2/28/23.
//

#ifndef BELFEM_FEM_GEOMETRY_HPP
#define BELFEM_FEM_GEOMETRY_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Group.hpp"

namespace belfem
{
    namespace fem
    {
        namespace geometry
        {
//------------------------------------------------------------------------------

            const Vector <real> &
                    normal_tri( Group * aGroup,
            const uint aMasterIndex
            );

//------------------------------------------------------------------------------

            const Vector <real> &
                    normal_tri( Group * aGroup,
            const uint aMasterIndex,
            const uint aIndex
            );

//------------------------------------------------------------------------------

            const Vector <real> &
                    normal_quad( Group * aGroup,
            const uint aMasterIndex
            );

//------------------------------------------------------------------------------

            const Vector <real> &
                    normal_quad( Group * aGroup,
            const uint aMasterIndex,
            const uint aIndex
            );

//------------------------------------------------------------------------------

            const Vector <real> &
                    normal_tet( Group * aGroup,
            const uint aMasterIndex
            );

//------------------------------------------------------------------------------

            const Vector <real> &
                    normal_tet( Group * aGroup,
            const uint aMasterIndex,
            const uint aIndex
            );

//------------------------------------------------------------------------------

            const Vector <real> &
                    normal_hex( Group * aGroup,
            const uint aMasterIndex,
            const uint aIndex
            );

//------------------------------------------------------------------------------

            /*
             * help function, computes the inverse of the Jacobian matrix
             */
            const Matrix< real > &
                    inv_J_2d( Group * aGroup );

//------------------------------------------------------------------------------

            /*
             * help function, computes the inverse of the Jacobian matrix
             */
            const Matrix< real > &
                    inv_J_3d( Group * aGroup );

//------------------------------------------------------------------------------

            /*
             * help function, computes the inverse of the Jacobian matrix
             */
            const Matrix< real > &
                    inv_J_tri3( Group * aGroup );

//------------------------------------------------------------------------------

            /*
             * help function, computes the inverse of the jacobian matrix
             */
            const Matrix< real > &
                    inv_J_tet4( Group * aGroup );

//------------------------------------------------------------------------------

            /**
             * special function only needed for test
             */
            real
            test_pipette( Group * aGroup ) ;

//------------------------------------------------------------------------------

            /**
            * special function only needed for test
            */
            bool
            test_intpoints( Group * aGroup,
                            const Matrix< real > & aPoints,
                            const Vector< real > & aWeights );

//------------------------------------------------------------------------------
        }
    }
}

#endif //BELFEM_FEM_GEOMETRY_HPP
