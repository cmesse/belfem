//
// Created by Christian Messe on 2019-07-28.
//

#ifndef BELFEM_CL_MESH_EXODUSWRITER_HPP
#define BELFEM_CL_MESH_EXODUSWRITER_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Mesh.hpp"

namespace belfem
{
    namespace mesh
    {
        class ExodusWriter
        {
#ifdef BELFEM_EXODUS
            // pointer to mesh
            Mesh * mMesh;

            // path to file
            string mPath;

            // handler to file
            int mHandle = -1;

            // error handle
            int mError = 0;

            int mCpuWordSize;

            int mIoWordSize;

            int64_t mNumDim = 0;
            int64_t mNumNodes = 0;
            int64_t mNumElements = 0;

            int64_t mNumBlocks = 0;
            int64_t mNumSideSets = 0;
            int64_t mNumNodeSets = 0;

            int mTimeStep = 1;
            double mTimeValue = 0.0;

#endif

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            ExodusWriter( Mesh * aMesh );

//------------------------------------------------------------------------------

            ~ExodusWriter() = default;

//------------------------------------------------------------------------------

            void
            save( const string & aPath );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            check( const string & aRoutine );

//------------------------------------------------------------------------------
            void
            close_file();

//------------------------------------------------------------------------------

            void
            populate_header();

//------------------------------------------------------------------------------

            void
            populate_node_coords();

//------------------------------------------------------------------------------

            void
            populate_element_ids();

//------------------------------------------------------------------------------

            void
            populate_blocks();

//------------------------------------------------------------------------------

            string
            element_string_from_type( const ElementType & aType );

//------------------------------------------------------------------------------

            void
            populate_element_connectivity();

//------------------------------------------------------------------------------

            void
            populate_sidesets();

//------------------------------------------------------------------------------

            void
            populate_time();

//------------------------------------------------------------------------------

            void
            populate_global_variables();

//------------------------------------------------------------------------------

            void
            populate_fields();

//------------------------------------------------------------------------------

            void
            populate_node_fields( Cell< mesh::Field * > & aFields );

//------------------------------------------------------------------------------

            void
            populate_element_fields( Cell< mesh::Field * > & aFields );

//------------------------------------------------------------------------------

            /**
             * special treatment for unsupported higher order elements
             */
            int
            fix_num_nodes( const ElementType & aType );

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_MESH_EXODUSWRITER_HPP
