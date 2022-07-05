//
// Created by Christian Messe on 22.06.20.
//

#ifndef BELFEM_CL_MESH_VTKWRITER_HPP
#define BELFEM_CL_MESH_VTKWRITER_HPP

#include <fstream>
#include "typedefs.hpp"
#include "cl_Mesh.hpp"
#include "cl_Map.hpp"

namespace belfem
{
    namespace mesh
    {
        class VtkWriter
        {
            const string mFilePath;

            Mesh * mMesh;

            std::ofstream mFile;

            // map that connects node IDs to indices in that order as they
            // have been written into the file
            Map< id_t, int > mNodeMap;

            int mNumberOfNodes ;

            int mNumberOfElements ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            VtkWriter( const string & aFilePath, Mesh * aMesh );

//------------------------------------------------------------------------------

            ~VtkWriter() = default;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            write_header();

//------------------------------------------------------------------------------

            void
            write_time();

//------------------------------------------------------------------------------

            void
            write_nodes();

//------------------------------------------------------------------------------

            void
            write_elements();

//------------------------------------------------------------------------------

            void
            write_element_fields();

//------------------------------------------------------------------------------

            void
            write_node_fields();

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

    }
//------------------------------------------------------------------------------
}
#endif //BELFEM_CL_MESH_VTKWRITER_HPP
