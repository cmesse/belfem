//
// Created by Christian Messe on 29.04.20.
//

#include "cl_BS_Element.hpp"

namespace belfem
{
    namespace bspline
    {
//------------------------------------------------------------------------------

        Element::Element( mesh::Element * aElement ) :
            mElement( aElement )
        {
            mBasis.set_size( aElement->number_of_nodes(), nullptr );
        }

//------------------------------------------------------------------------------

        void
        Element::set_basis(  Basis * aBasis, const uint aIndex )
        {
            mBasis( aIndex ) = aBasis;
        }

//------------------------------------------------------------------------------

        void
        Element::unflag_basis()
        {
            for( Basis * tBasis : mBasis )
            {
                tBasis->unflag();
            }
        }

//------------------------------------------------------------------------------
    }
}
