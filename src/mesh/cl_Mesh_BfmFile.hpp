/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_CL_MESH_BFMFILE_HPP
#define BELFEM_CL_MESH_BFMFILE_HPP

#include "cl_HDF5.hpp"

#ifndef BELFEM_HDF5
    typedef void hvl_t;
#endif


namespace belfem
{
    class Mesh ;

    namespace mesh
    {
        class ProtoMesh ;


        class BfmFile
        {
            const string mFilePath;

            Mesh * mMesh = nullptr;
            bool mOwnMesh = true ;

            ProtoMesh * mProto = nullptr ;

            HDF5 * mFile = nullptr;

            Map< id_t, uint > mNumNodesPerElement ;



//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            BfmFile(
                const string & aFilePath,
                Mesh * aMesh = nullptr );

//------------------------------------------------------------------------------

            ~BfmFile()  ;

//------------------------------------------------------------------------------

            size_t
            checksum() const ;

            //! settings fingerprint of the stored mesh, read WITHOUT loading the
            //! payload -- the reuse decision happens before load(), so it cannot
            //! wait for load_meta_data(). Returns 0 / empty on files written
            //! before the tag existed, which the caller must treat as a mismatch
            uint64_t
            config_tag() const ;

            string
            config_text() const ;

            Mesh *
            get() ;

//------------------------------------------------------------------------------

            void
            save();

            void
            load();

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------


            void
            save_meta_data();

            void
            load_meta_data();

            void
            save_node_data();

            void
            load_node_data();

            void
            save_group_data( const GroupType aType );

            void
            load_group_data( const GroupType aType );

            void
            save_element_data();

            void
            load_element_data();

            void
            save_facet_data( const bool aSaveTopology = false );

            void
            load_facet_data();

            void
            save_edge_data();

            void
            load_edge_data();

            void
            save_face_data();

            void
            load_face_data();

            void
            save_control_point_data();

            void
            load_control_point_data();

            void
            save_node_duplicate_data();

            void
            load_node_duplicate_data();

            void
            save_hanging_entities();

            void
            load_hanging_entities();

            void
            save_hanging_nodes( const hsize_t aNumNodes );

            void
            load_hanging_nodes();

            void
            save_hanging_edges( const hsize_t aNumEdges );

            void
            load_hanging_edges();

            void
            save_hanging_faces( const hsize_t aNumFaces );

            void
            load_hanging_faces();

            void
            save_hanging_facets( const hsize_t aNumFacets );

            void
            load_hanging_facets();

            void
            save_hanging_control_points( const hsize_t aNumControlPoints );

            void
            load_hanging_control_points();

            void
            save_periodicity_data();

            void
            save_periodic_planes();

            void
            save_periodic_nodes();

            void
            save_periodic_edges();

            void
            save_periodic_faces();

            void
            save_periodic_facets();

            void
            load_periodicity_data();

            void
            load_periodic_planes();

            void
            load_periodic_nodes();

            void
            load_periodic_edges();

            void
            load_periodic_faces();

            void
            load_periodic_facets();

            void
            save_thinshell_data();

            void
            load_thinshell_data();

            void
            save_vertex_data();

            void
            load_vertex_data();

            void
            save_curve_data();

            void
            load_curve_data();
        };

        inline
        Mesh *
        BfmFile::get()
        {
            mOwnMesh = false ;
            return mMesh ;
        }

    }
}
#endif // BELFEM_CL_MESH_BFMFILE_HPP
