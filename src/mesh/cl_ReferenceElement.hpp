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

#ifndef BELFEM_CL_REFERENCEELEMENT_HPP
#define BELFEM_CL_REFERENCEELEMENT_HPP


#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Node.hpp"
#include "cl_Element.hpp"

namespace belfem
{
    class Mesh ;
    namespace mesh
    {
        class ReferenceElement
        {
            Element * mElement = nullptr ;
            Cell< Node * > mNodes ;
            bool mOwnPointers = true ;
        public:

            ReferenceElement(Element * aElement, Cell< Node * > aNodes );

            ~ReferenceElement() ;

            Mesh *
            export_mesh_and_pass_ownerships();

            Element * element() { return mElement ; }

            const Element * element() const { return mElement ; }

            Cell< Node * > & nodes() { return mNodes ; }

            Node *
            node( const uint aIndex ) { return mNodes( aIndex ) ; }

            const Node *
            node( const uint aIndex ) const { return mNodes( aIndex ) ; }
        };
    }
}
#endif //BELFEM_CL_REFERENCEELEMENT_HPP