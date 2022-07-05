//
// Created by Christian Messe on 29.04.20.
//

#ifndef BELFEM_CL_BS_BASIS_HPP
#define BELFEM_CL_BS_BASIS_HPP

#include "typedefs.hpp"
#include "cl_Graph_Vertex.hpp"
#include "cl_Cell.hpp"
namespace belfem
{
    namespace bspline
    {
        // forward declaration of element
        class Element ;

        class Basis : public graph::Vertex
        {
            Cell< Element * > mElements ;

            uint mElementCounter = 0;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Basis( const id_t aID ) ;

//------------------------------------------------------------------------------

            ~Basis();

//------------------------------------------------------------------------------

            inline void
            increment_element_counter()
            {
                ++mElementCounter ;
            }

//------------------------------------------------------------------------------

            void
            init_element_container();

//------------------------------------------------------------------------------

            void
            insert_element( Element * aElement );

//------------------------------------------------------------------------------

            void
            link_basis();

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_BS_BASIS_HPP
