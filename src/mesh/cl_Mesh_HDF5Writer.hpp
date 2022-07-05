//
// Created by Christian Messe on 27.10.19.
//

#ifndef BELFEM_CL_MESH_HDF5WRITER_HPP
#define BELFEM_CL_MESH_HDF5WRITER_HPP

#include "cl_Mesh.hpp"
#include "cl_HDF5.hpp"

namespace belfem
{
    namespace mesh
    {
        class HDF5Writer
        {
            const string mFilePath;

            HDF5 mFile;

            Mesh * mMesh;



//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

            HDF5Writer( const string & aFilePath, Mesh * aMesh  );

//------------------------------------------------------------------------------

            ~HDF5Writer() = default;

//------------------------------------------------------------------------------
            private:
//------------------------------------------------------------------------------

            void
            save_nodes();


//------------------------------------------------------------------------------

            void
            save_edges();

//------------------------------------------------------------------------------

            void
            save_faces();

//------------------------------------------------------------------------------

            void
            save_elements();

//------------------------------------------------------------------------------

            void
            save_sidesets();

//------------------------------------------------------------------------------

            void
            save_vertices();

//------------------------------------------------------------------------------

            void
            save_globals();

//------------------------------------------------------------------------------

            void
            save_fields();

//------------------------------------------------------------------------------

        };
    }
}
#endif //BELFEM_CL_MESH_HDF5WRITER_HPP
