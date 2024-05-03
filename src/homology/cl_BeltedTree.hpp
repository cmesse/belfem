//
// Created by grgia@ge.polymtl.ca on 31/01/24.
//

#include "cl_SimplicialComplex.hpp"
#include "cl_Chain.hpp"
#include "cl_Cochain.hpp"
#include "cl_Edge.hpp"
#include "cl_Element.hpp"

#ifndef BELFEM_CL_BELTEDTREE_HPP
#define BELFEM_CL_BELTEDTREE_HPP

namespace belfem
{
    namespace mesh
    {
//-----------------------------------------------------------------------------
        class BeltedTree
        {
            Cell< Chain * > m1HomologyGenerators ;

            Cell< Cochain * > m1CohomologyGenerators ;

            Cell < Edge * > mBeltFasteners ;

            Cell < id_t > mTree ;

            SimplicialComplex * mSimplicialComplex ;

            Mesh* mMesh ;

            bool mCycleFound = false;

            id_t mNodeToFind = 0;

            //-----------------------------------------------------------------------------
        public:

            //-----------------------------------------------------------------------------

            BeltedTree(Mesh * aMesh,  SimplicialComplex * aSimplicialComplex, Cell< Chain * > a1HomologyGenerators) ;

            //-----------------------------------------------------------------------------

            ~BeltedTree();

            //-----------------------------------------------------------------------------

            void
            select_belt_fasteners() ;

            //-----------------------------------------------------------------------------

            void
            create_tree() ;

            //-----------------------------------------------------------------------------

            void
            dfs(uint v, uint e, Map <id_t, bool > & tfinished, Map <id_t, bool > & tvisited) ;

            //-----------------------------------------------------------------------------

            void
            compute_cohomology() ;

            //-----------------------------------------------------------------------------

            void
            create_TreeField(string tFieldName) ;

            //-----------------------------------------------------------------------------

            void
            create_cohomologyField( Mesh* tMeshEdge ) ;

            //-----------------------------------------------------------------------------
        };
    }
}

#endif //BELFEM_CL_BELTEDTREE_HPP
