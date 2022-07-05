//
// Created by Christian Messe on 29.04.20.
//

#ifndef BELFEM_CL_BS_ELEMENT_HPP
#define BELFEM_CL_BS_ELEMENT_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Graph_Vertex.hpp"
#include "Mesh_Enums.hpp"
#include "cl_Element.hpp"
#include "cl_BS_Basis.hpp"
namespace belfem
{
    namespace bspline
    {
//------------------------------------------------------------------------------

        class Element
        {
            // background element on mesh
            mesh::Element * mElement ;
            Cell< Basis * > mBasis ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Element( mesh::Element * aElement );

//------------------------------------------------------------------------------

            void
            set_basis( Basis * aBasis, const uint & aIndex );

//------------------------------------------------------------------------------

            // expose element on mesh
            inline mesh::Element *
            element()
            {
                return mElement;
            }

//-----------------------------------------------------------------------------

            // expose basis container
            inline Cell< Basis * > &
            basis()
            {
                return mBasis;
            }

//-----------------------------------------------------------------------------
            // expose basis
            inline Basis *
            basis( const uint aIndex )
            {
                return mBasis( aIndex );
            }

//------------------------------------------------------------------------------

            // reset the flag of the basis
            void
            unflag_basis();

//------------------------------------------------------------------------------

        };

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_BS_ELEMENT_HPP
