//
// Created by christian on 7/12/21.
//

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
            const proc_t mRank ;

            // ref to mesh we work on
            Mesh & mMesh ;

            const luint mNumberOfNodes ;

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

            luint
            face_key_2d(
                    Element           * aElement,
                    Cell< Node * >    & aTriNodes,
                    Cell< Node * >    & aQuadNodes );

//-----------------------------------------------------------------------

            luint
            face_key_tri(
                    Element           * aElement,
                    const uint          aFaceIndex,
                    Cell< Node * >    & aAllNodes,
                    Cell< Node * >    & aCornerNodes );

//-----------------------------------------------------------------------

            luint
            face_key_quad(
                    Element           * aElement,
                    const uint          aFaceIndex,
                    Cell< Node * >    & aAllNodes,
                    Cell< Node * >    & aCornerNodes );

//-----------------------------------------------------------------------

            void
            get_all_block_ids( Vector< id_t > & aBlockIDs ) ;

//-----------------------------------------------------------------------

            index_t
            count_faces(
                    const Vector< id_t >  & aBlockIDs,
                    const Vector< id_t >  & aSideSetIDs,
                          Map< luint, index_t >& aFaceMap );

//-----------------------------------------------------------------------

            void
            find_face_owners(
                    const Vector< id_t >        & aBlockIDs,
                    const Vector< id_t >        & aSideSetIDs,
                    const index_t               & aNumFaces,
                    const Map< luint, index_t > & aFaceMap,
                    Vector< id_t >              & aMasterOwner,
                    Vector< index_t >           & aMasterIndex,
                    Vector< id_t >              & aSlaveOwner,
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
            set_face_ids( const Vector< id_t > & aSideSets, const Map< luint, index_t > & aFaceMap );

//----------------------------------------------------------------------
        };

//-----------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_FACEFACTORY_HPP
