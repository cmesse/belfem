//
// Created by christian on 2/7/23.
//

#ifndef BELFEM_CL_MAXWELL_THERMALMESHEXTRACTOR_HPP
#define BELFEM_CL_MAXWELL_THERMALMESHEXTRACTOR_HPP

#include "cl_Mesh.hpp"

namespace belfem
{
    class Spline;

    namespace fem
    {
        class ThermalMeshExtractor
        {
            //! ref to the original mesh
            Mesh & mInputMesh ;

            // vector with IDs of blocks on original mesh
            const Vector< id_t > & mGhostBlockIDs ;

            // node id map
            Map< id_t, index_t > mNodeIndexMap ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            ThermalMeshExtractor( Mesh * aInputMesh, const Vector< id_t > & aGhostBlockIDs );

//-----------------------------------------------------------------------------

            ~ThermalMeshExtractor() = default ;

//-----------------------------------------------------------------------------

            Mesh *
            create_mesh();

//-----------------------------------------------------------------------------
        private :
//-----------------------------------------------------------------------------

            Mesh *
            create_2d_mesh();

//-----------------------------------------------------------------------------

            Mesh *
            create_3d_mesh();

//-----------------------------------------------------------------------------

            // creates the nodes on the output mesh
            void
            create_nodes( Mesh * aMesh );

//-----------------------------------------------------------------------------

            void
            create_resistors( Mesh * aMesh );

//-----------------------------------------------------------------------------

            // compute the volume of the elements on the input mesh
            // and write into output
            void
            compute_element_volumes( Mesh * aMesh );

//-----------------------------------------------------------------------------
        };
    }
}

#endif //BELFEM_CL_MAXWELL_THERMALMESHEXTRACTOR_HPP
