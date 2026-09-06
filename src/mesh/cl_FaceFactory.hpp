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

#ifndef BELFEM_CL_FACEFACTORY_HPP
#define BELFEM_CL_FACEFACTORY_HPP
#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    class Mesh ;

    namespace mesh
    {
        class Node ;
        class Face ;
        class Element ;

//-----------------------------------------------------------------------

        class FaceFactory
        {
            const proc_t mCommRank ;

            // ref to mesh we work on
            Mesh & mMesh ;

            const key128_t mNumberOfNodes ;

//-----------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------

            FaceFactory( Mesh & aMesh );

//-----------------------------------------------------------------------

            FaceFactory( Mesh * aMesh );

//-----------------------------------------------------------------------

            ~FaceFactory() = default ;

//----------------------------------------------------------------------

            void
            create_faces(  const Vector< id_t > aNedelecBlocks = Vector< id_t >(),
                           const Vector< id_t > aNedelecSideSets = Vector< id_t >() );

//----------------------------------------------------------------------

            void
            print();

//-----------------------------------------------------------------------
        private:
//-----------------------------------------------------------------------

            key128_t
            face_key_2d(
                    Element           * aElement,
                    Cell< Node * >    & aWork  );

//-----------------------------------------------------------------------

            key128_t
            face_key_3d(
                    Element           * aElement,
                    const uint          aFaceIndex,
                    Cell< Node * >    & aWork );

//-----------------------------------------------------------------------

            void
            get_all_block_ids( Vector< id_t > & aBlockIDs ) ;

//-----------------------------------------------------------------------

            index_t
            count_faces(
                    const Vector< id_t >  & aBlockIDs,
                    const Vector< id_t >  & aSideSetIDs,
                          Map< key128_t, index_t >& aFaceMap );

//-----------------------------------------------------------------------

            void
            find_face_owners(
                    const Vector< id_t >        & aBlockIDs,
                    const Vector< id_t >        & aSideSetIDs,
                    const index_t               & aNumFaces,
                    const Map< key128_t, index_t > & aFaceMap,
                    Vector< id_t >              & aMasterIDs,
                    Vector< index_t >           & aMasterIndex,
                    Vector< id_t >              & aSlaveIDs,
                    Vector< index_t >           & aSlaveIndex );

//----------------------------------------------------------------------

            void
            allocate_face_containers( Vector< id_t > & aBlockIDs );

//----------------------------------------------------------------------


            // block faces in 2d
            void
            create_faces_2d( const Vector< id_t > & aBlockIDs );

//----------------------------------------------------------------------

            void
            create_faces_3d(
                        Vector< id_t >              & aMasterOwner,
                        Vector< index_t >           & aMasterIndex,
                        Vector< id_t >              & aSlaveOwner,
                        Vector< index_t >           & aSlaveIndex );

//----------------------------------------------------------------------

            // fix the ids for selected sidesets
            void
            set_face_ids( const Vector< id_t > & aSideSets, const Map< key128_t, index_t > & aFaceMap );

//----------------------------------------------------------------------
        };

//-----------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_FACEFACTORY_HPP
