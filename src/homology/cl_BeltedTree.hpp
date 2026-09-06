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

            Cell < index_t > mTree ;

            SimplicialComplex * mSimplicialComplex ;

            Mesh* mMesh ;

//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            BeltedTree(Mesh * aMesh,
                SimplicialComplex * aSimplicialComplex,
                Cell< Chain * > a1HomologyGenerators) ;

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
            compute_cohomology() ;

//-----------------------------------------------------------------------------

            void
            create_TreeField(Mesh* tEdgeMesh , string tFieldName) ;

//-----------------------------------------------------------------------------

            void
            create_cohomologyField( Mesh* tMeshEdge ) ;

//-----------------------------------------------------------------------------

            Cell< Cochain * > &
            get_cohomology();

//-----------------------------------------------------------------------------
        };
    }
}

#endif //BELFEM_CL_BELTEDTREE_HPP
