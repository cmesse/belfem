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

#ifndef BELFEM_CL_VTK_MESHVIEW_HPP
#define BELFEM_CL_VTK_MESHVIEW_HPP

#include "vtktypes.hpp"
#include "cl_Mesh.hpp"

namespace belfem
{
    namespace vtk
    {
        /**
         * @class BlockActor
         * @brief Owns the VTK actor for a single mesh Block.
         *
         * Composition, not inheritance: the wrapper HOLDS a vtkActor (lifetime
         * managed by the smart pointer) rather than being one. The block's
         * volume elements become an unstructured grid, surface-extracted by a
         * vtkDataSetMapper. The grid shares the parent MeshView's point set;
         * only cell connectivity is stored here.
         */
        class BlockActor
        {
            mesh::Block     * mBlock ;

            UnstructuredGrid  mGrid ;
            DataSetMapper     mMapper ;
            Actor             mActor ;

        public:
            BlockActor( Points                   aPoints,
                        Map< id_t, vtkIdType > & aNodeMap,
                        Map< id_t, vtkIdType > & aElementMap,
                        mesh::Block *            aBlock );

            ~BlockActor() = default;

            inline mesh::Block * block() { return mBlock; }

            //! The VTK actor; ownership stays with this wrapper.
            inline Actor actor() { return mActor; }
        };

        class SideSetActor
        {
            mesh::SideSet     * mSideSet ;

            UnstructuredGrid  mGrid ;
            DataSetMapper     mMapper ;
            Actor             mActor ;

        public:
            SideSetActor( Points                   aPoints,
                        Map< id_t, vtkIdType > & aNodeMap,
                        Map< id_t, vtkIdType > & aElementMap,
                        mesh::SideSet *            aSideSet );

            ~SideSetActor() = default;

            inline mesh::SideSet * sideset() { return mSideSet; }

            //! The VTK actor; ownership stays with this wrapper.
            inline Actor actor() { return mActor; }
        };

        /**
         * @class MeshView
         * @brief The renderable view of a belfem::Mesh: one shared point set
         *        plus per-block and per-sideset actor wrappers.
         *
         * Source-agnostic: it knows nothing about what the mesh represents.
         * A consumer collects the actors via actors() and, when the mesh is
         * placed in a larger scene (e.g. a spacecraft in an orbit viewer),
         * moves the whole view as a rigid body via set_user_transform().
         *
         * @ingroup grp_visualizer
         * @see @ref visualizer_index
         */
        class MeshView
        {

            belfem::Mesh    * mMesh ;
            Points    mPoints ;
            Map< id_t, vtkIdType > mNodeMap ;
            Map< id_t, vtkIdType > mElementMap ;

            Cell< BlockActor * >   mBlocks ;
            Cell< SideSetActor * > mSideSets ;

        public:

            MeshView( belfem::Mesh * aMesh ) ;

            ~MeshView();

            //! The per-block actor wrappers; ownership stays with this view.
            inline Cell< BlockActor * > & blocks() { return mBlocks; }

            //! The per-sideset actor wrappers; ownership stays with this view.
            inline Cell< SideSetActor * > & sidesets() { return mSideSets; }

            //! All VTK actors (blocks + sidesets), for handing to a renderer.
            Cell< Actor > actors();

            //! Apply one user transform to every block/sideset actor, so the
            //! whole mesh moves as a rigid body. The caller keeps ownership
            //! of the transform and may keep mutating it.
            void set_user_transform( const Transform & aTransform );

        };
    }
}
#endif //BELFEM_CL_VTK_MESHVIEW_HPP
