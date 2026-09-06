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
