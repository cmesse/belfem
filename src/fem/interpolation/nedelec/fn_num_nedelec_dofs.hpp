//
// Created by christian on 1/15/25.
//

#ifndef FN_NUM_NEDELEC_DOFS_HPP
#define FN_NUM_NEDELEC_DOFS_HPP
#include "typedefs.hpp"
#include "assert.hpp"
#include "Mesh_Enums.hpp"


namespace belfem
{
    namespace fem
    {
        uint
        num_nedelec_dofs( const ElementType aElementType )
        {
           switch( aElementType )
           {
                case ElementType::LINE2 : return  1 ;
                case ElementType::LINE3 : return  2 ;
                case ElementType::TRI3  : return  3 ;
                case ElementType::TRI6  : return  8 ;
                case ElementType::TET4  : return  6 ;
                case ElementType::TET10 : return 20 ;

                case ElementType::QUAD4TS : return 2 ;
                case ElementType::QUAD9TS : return 3 ;
                case ElementType::PENTA6TS : return  6 ;
                case ElementType::PENTA18TS : return 16 ;

                case ElementType::HEX8 : return 12 ;
                case ElementType::HEX8TS : return 8 ;
                case ElementType::HEX8TB : return 4 ;

                default : return 0 ;
           }
        }
    }
}
#endif //FN_NUM_NEDELEC_DOFS_HPP
