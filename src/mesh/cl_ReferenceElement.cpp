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

#include "commtools.hpp"
#include "cl_ReferenceElement.hpp"
#include "cl_Mesh.hpp"
#include "random.hpp"

namespace belfem
{
    namespace mesh
    {
        ReferenceElement::ReferenceElement(Element * aElement, Cell< Node * > aNodes) :
        mElement(aElement),
        mNodes(aNodes)
        {

        }

        ReferenceElement::~ReferenceElement()
        {
            if ( mOwnPointers )
            {
                delete mElement;

                for ( Node * tNode : mNodes )
                {
                    delete tNode;
                }
            }
        }

        Mesh *
        ReferenceElement::export_mesh_and_pass_ownerships()
        {
            Mesh * aMesh = new Mesh( dimension( mElement->type() ) );

            aMesh->nodes() = mNodes ;
            Block * tBlock = new Block( 1, 1 );
            tBlock->insert_element( mElement );
            aMesh->add_block( tBlock );
            aMesh->finalize();
            mOwnPointers = false;

            return aMesh;
        }

    }
}
