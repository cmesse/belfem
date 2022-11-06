//
// Created by Christian Messe on 2019-04-01.
//

#ifndef BELFEM_MESH_ENUMS_HPP
#define BELFEM_MESH_ENUMS_HPP

#include <string>

//------------------------------------------------------------------------------
namespace belfem
{
    /**
     * Element types. ( taken from gmsh )
     */
    enum class ElementType
    {
        EMPTY   =  0,
        LINE2   =  1,
        TRI3    =  2,
        QUAD4   =  3,
        TET4    =  4,
        HEX8    =  5,
        PENTA6  =  6,
        LINE3   =  8,
        TRI6    =  9,
        QUAD9   = 10,
        TET10   = 11,
        HEX27   = 12,
        PENTA18 = 13,
        VERTEX  = 15,
        QUAD8   = 16,
        HEX20   = 17,
        PENTA15 = 18,
        TRI10   = 21,
        TRI15   = 23,
        TRI21   = 25,
        LINE4   = 26,
        LINE5   = 27,
        LINE6   = 28,
        TET20   = 29,
        TET35   = 30,
        QUAD16  = 32,
        HEX64   = 92,
        UNDEFINED = 127
    };

//------------------------------------------------------------------------------

    enum class GeometryType
    {
        VERTEX = 0,
        LINE   = 1,
        TRI    = 2,
        QUAD   = 3,
        TET    = 4,
        HEX    = 5,
        PENTA  = 6,
        PYRA   = 7,
        UNDEFINED = 8,
    };

//------------------------------------------------------------------------------

    enum class InterpolationOrder
    {
        CONSTANT,
        LINEAR,
        QUADRATIC,
        SERENDIPITY,
        CUBIC,
        QUARTIC,
        QUINTIC,
        UNDEFINED
    };

//------------------------------------------------------------------------------

    enum class InterpolationType
    {
        LAGRANGE,
        HERMITE,
        BERNSTEIN,
        UNEFINED
    };

//------------------------------------------------------------------------------

    enum class EntityType
    {
        NODE      = 0,
        EDGE      = 1, // only used for DOFS, field is not saved in exodus
        FACE      = 2, // only used for DOFS, field is not saved in exodus
        CELL      = 3, // only used for DOFS, same as Element, but field is not saved
        FACET     = 4, // used for lambda dofs, field is not saved
        ELEMENT   = 5,
        UNDEFINED = 6
    };

//------------------------------------------------------------------------------

    enum class GroupType
    {
        BLOCK,
        SIDESET,
        CUT,
        SHELL,
        UNDEFINED
    };

//------------------------------------------------------------------------------

    enum class FieldType
    {
        SCALAR,
        UNDEFINED
    };

//------------------------------------------------------------------------------

    inline std::string
    to_string( const ElementType aElementType )
    {
        switch( aElementType )
        {
            case( ElementType::EMPTY ) :
            {
                return "empty" ;
            }
            case( ElementType::LINE2 ) :
            {
                return "line2" ;
            }
            case( ElementType::TRI3 ) :
            {
                return "tri3" ;
            }
            case( ElementType::QUAD4 ) :
            {
                return "quad4" ;
            }
            case( ElementType::TET4 ) :
            {
                return "tet4" ;
            }
            case( ElementType::HEX8 ) :
            {
                return "hex8" ;
            }
            case( ElementType::PENTA6 ) :
            {
                return "penta6" ;
            }
            case( ElementType::LINE3 ) :
            {
                return "line3" ;
            }
            case( ElementType::TRI6 ) :
            {
                return "tri6" ;
            }
            case( ElementType::QUAD9 ) :
            {
                return "quad9" ;
            }
            case( ElementType::TET10 ) :
            {
                return "tet10" ;
            }
            case( ElementType::HEX27 ) :
            {
                return "hex27" ;
            }
            case( ElementType::PENTA18 ) :
            {
                return "penta18" ;
            }
            case( ElementType::VERTEX ) :
            {
                return "vertex" ;
            }
            case( ElementType::QUAD8 ) :
            {
                return "quad8" ;
            }
            case( ElementType::HEX20 ) :
            {
                return "hex20" ;
            }
            case( ElementType::PENTA15 ) :
            {
                return "penta15" ;
            }
            case( ElementType::TRI10 ) :
            {
                return "tri10" ;
            }
            case( ElementType::TRI15 ) :
            {
                return "tri15" ;
            }
            case( ElementType::TRI21 ) :
            {
                return "tri21" ;
            }
            case( ElementType::LINE4 ) :
            {
                return "line4" ;
            }
            case( ElementType::LINE5 ) :
            {
                return "line5" ;
            }
            case( ElementType::LINE6 ) :
            {
                return "line6" ;
            }
            case( ElementType::TET20 ) :
            {
                return "tet20" ;
            }
            case( ElementType::TET35 ) :
            {
                return "tet35" ;
            }
            case( ElementType::QUAD16 ) :
            {
                return "quad16" ;
            }
            case( ElementType:: HEX64 ) :
            {
                return "hex64" ;
            }
            default :
            {
                return "unknown" ;
            }
        }
    }
//------------------------------------------------------------------------------
}
#endif //BELFEM_MESH_ENUMS_HPP
