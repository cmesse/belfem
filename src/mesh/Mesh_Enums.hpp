//
// Created by Christian Messe on 2019-04-01.
//

#ifndef BELFEM_MESH_ENUMS_HPP
#define BELFEM_MESH_ENUMS_HPP
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
        UNEFINED
    };

//------------------------------------------------------------------------------

    enum class EntityType
    {
        NODE,
        FACET, // used for lambda dofs, field is not saved
        ELEMENT,
        EDGE, // only used for DOFS, field is not saved in exodus
        FACE, // only used for DOFS, field is not saved in exodus
        CELL, // only used for DOFS, same as Element, but field is not saved
        UNDEFINED
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
}
#endif //BELFEM_MESH_ENUMS_HPP
