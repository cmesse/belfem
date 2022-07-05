//
// Created by Christian Messe on 27.10.19.
//

#ifndef BELFEM_CL_MESH_HDF5READER_HPP
#define BELFEM_CL_MESH_HDF5READER_HPP

#include "cl_Mesh.hpp"
#include "cl_HDF5.hpp"

namespace belfem
{
    namespace mesh
    {
        class HDF5Reader
        {
            const string mFilePath;

            HDF5 mFile;

            Mesh * mMesh;

            /**
             * flag telling if mesh is deleted on destruction
             * once the mesh is called using this->mesh()
             * or passed in the constructor, it is the
             *   responsibility of the user to delete it
             */
            bool mOwnMesh = false;

            // header
            uint mNumberOfDimensions;
            uint mNumberOfBlocks;
            uint mNumberOfSideSets;
            uint mNumberOfGlobals;
            uint mNumberOfFields;
            uint mTimeStep;
            real mTimeStamp;

            uint mNumberOfNodes;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            HDF5Reader( const string & aFilePath, Mesh * aMesh = nullptr );

//------------------------------------------------------------------------------

            ~HDF5Reader();
//------------------------------------------------------------------------------

            /**
             * get the pointer to the mesh
             */
             Mesh *
             mesh();

//------------------------------------------------------------------------------

        private:

//------------------------------------------------------------------------------

            void
            read_header();

//------------------------------------------------------------------------------

            void
            read_nodes();

//------------------------------------------------------------------------------

            void
            read_blocks();

//------------------------------------------------------------------------------

            void
            read_sidesets();

//------------------------------------------------------------------------------

            void
            read_vertices();

//------------------------------------------------------------------------------

            void
            read_globals();

//------------------------------------------------------------------------------

            void
            read_fields();

//------------------------------------------------------------------------------
        };
    }
}

#endif //BELFEM_CL_MESH_HDF5READER_HPP
