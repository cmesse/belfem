//
// Created by Christian Messe on 05.04.21.
//

#ifndef BELFEM_CL_TENSORMESHFACTORY_HPP
#define BELFEM_CL_TENSORMESHFACTORY_HPP

#include "typedefs.hpp"
#include "cl_Mesh.hpp"

namespace belfem
{
    class TensorMeshFactory
    {
//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        TensorMeshFactory() = default ;

//------------------------------------------------------------------------------

        ~TensorMeshFactory() = default ;

//------------------------------------------------------------------------------

        Mesh *
        create_tensor_mesh(
                const Vector< uint > & aNumElems,
                const Vector< real > & aMinPoint,
                const Vector< real > & aMaxPoint,
                const uint aOrder = 1 );

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

        void
        create_2d_mesh_linear(
                const Vector< uint > & aNumElems,
                const Vector< real > & aMinPoint,
                const Vector< real > & aMaxPoint,
                      Mesh           * aMesh );

//------------------------------------------------------------------------------

        void
        create_2d_mesh_quadratic(
                const Vector< uint > & aNumElems,
                const Vector< real > & aMinPoint,
                const Vector< real > & aMaxPoint,
                Mesh           * aMesh );

//------------------------------------------------------------------------------

        void
        create_3d_mesh_linear(
                const Vector< uint > & aNumElems,
                const Vector< real > & aMinPoint,
                const Vector< real > & aMaxPoint,
                Mesh           * aMesh );

//------------------------------------------------------------------------------

        void
        create_3d_mesh_quadratic(
                const Vector< uint > & aNumElems,
                const Vector< real > & aMinPoint,
                const Vector< real > & aMaxPoint,
                Mesh           * aMesh );

//------------------------------------------------------------------------------
    };
}

#endif //BELFEM_CL_TENSORMESHFACTORY_HPP
