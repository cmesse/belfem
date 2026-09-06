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

#ifndef BELFEM_MESH_ENUMS_HPP
#define BELFEM_MESH_ENUMS_HPP

#include <string>

//------------------------------------------------------------------------------
namespace belfem
{
    /**
     * Element types. Numeric values follow gmsh where a gmsh type exists
     * ( QUAD16 = 32 deviates from gmsh's 36; element_type_from_gmsh() in
     * meshtools.cpp bridges that; the *TS and HEX8TB values are BELFEM
     * extensions ). must not be bigger than 255
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
        PYRA5   =  7,
        LINE3   =  8,
        TRI6    =  9,
        QUAD9   = 10,
        TET10   = 11,
        HEX27   = 12,
        PENTA18 = 13,
        PYRA14  = 14,
        VERTEX  = 15,
        QUAD8   = 16,
        HEX20   = 17,
        PENTA15 = 18,
        PYRA13  = 19,
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
        QUAD4TS = 103,
        QUAD9TS = 110,
        HEX8TS    = 105, // thin shell variant (not gmsh standard)
        PENTA6TS  = 106, // thin shell variant (not gmsh standard)
        PENTA18TS = 118, // thin shell variant (not gmsh standard)
        HEX8TB    = 125, // thin beam variant  (not gmsh standard)
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
        BubbleEdge0,
        BubbleEdge1,
        BubbleEdge2,
        BubbleFace0,
        BubbleFace1,
        BubbleFace2,
        BubbleFace3,
        UNEFINED
    };

//------------------------------------------------------------------------------

    enum class EntityType
    {
        NODE         = 0,
        EDGE         = 1, // only used for DOFS, field is not saved in exodus
        FACE         = 2, // only used for DOFS, field is not saved in exodus
        CELL         = 3, // only used for DOFS, same as Element, but field is not saved
        FACET        = 4, // used for lambda dofs, field is not saved
        ELEMENT      = 5,
        CONTROLPOINT = 6, // used for B-Splines
        UNDEFINED    = 7
    };

//------------------------------------------------------------------------------

    enum class GroupType
    {
        BLOCK,
        SIDESET,
        UNDEFINED
    };

//------------------------------------------------------------------------------

    enum class FieldType
    {
        SCALAR,
        UNDEFINED
    };

    enum class Connectivity
    {
        Compute = 0,
        NodeToVertex = 1,
        NodeToNode = 2,
        NodeToEdge = 3,
        NodeToFace = 4,
        NodeToFacet = 5,
        NodeToElement = 6,
        EdgeToVertex = 7,
        EdgeToNode = 8,
        EdgeToEdge = 9,
        EdgeToFace = 10,
        EdgeToFacet = 11,
        EdgeToElement = 12,
        FaceToVertex = 13,
        FaceToNode = 14,
        FaceToEdge = 15,
        FaceToFace = 16,
        FaceToFacet = 17,
        FaceToElement = 18,
        FacetToVertex = 19,
        FacetToNode = 20,
        FacetToEdge = 21,
        FacetToFace = 22,
        FacetToFacet = 23,
        FacetToElement = 24,
        ElementToVertex = 25,
        ElementToNode = 26,
        ElementToEdge = 27,
        ElementToFace = 28,
        ElementToFacet = 29,
        ElementToElement = 30,
        TsElementToTsElement = 31, // special for thin shells
        ShellToShell = 32,         // this actually tests neighbors (shared edges)
        ControlPointToControlPoint = 33,
        ElementToControlPoint = 34,
        ControlPointToElement = 35,
        UNDEFINED = 36
    };

//------------------------------------------------------------------------------

    inline std::string
    to_string( const GeometryType aGeometryType )
    {
        switch( aGeometryType )
        {
            case ( GeometryType::VERTEX ) :
            {
                return "vertex";
            }
            case ( GeometryType::LINE ) :
            {
                      return "line";
            }
            case ( GeometryType::TRI ) :
            {
                return "tri";
            }
            case ( GeometryType::QUAD ) :
            {
                return "quad";
            }
            case ( GeometryType::TET ) :
            {
                return "tet";
            }
            case ( GeometryType::PYRA ) :
            {
                return "pyra";
            }
            case ( GeometryType::PENTA ) :
            {
                return "penta";
            }
            case ( GeometryType::HEX ) :
            {
                return "hex";
            }
            default :
            {
                return "unknown";
            }
        }
    }

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
            case( ElementType::QUAD4TS ) :
            {
                return "quad4" ;
            }
            case( ElementType::TET4 ) :
            {
                return "tet4" ;
            }
            case( ElementType::HEX8 ) :
            case( ElementType::HEX8TS ) :
            case( ElementType::HEX8TB ) :
            {
                return "hex8" ;
            }
            case( ElementType::PENTA6 ) :
            case ( ElementType::PENTA6TS ) :
            {
                return "penta6" ;
            }
            case( ElementType::PYRA5 ) :
            {
                return "pyra5" ;
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
            case( ElementType::PYRA14 ) :
            {
                return "pyra14" ;
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
            case( ElementType::PYRA13 ) :
            {
                return "pyra13" ;
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


    inline ElementType element_type( const std::string & aStr )
    {
        uint tNumTypes = static_cast<uint>( ElementType::UNDEFINED );
        for ( uint k=0; k<tNumTypes; ++k )
        {
            if ( aStr == to_string( static_cast<ElementType>( k ) ) )
            {
                return static_cast<ElementType>( k );
            }
        }
        return ElementType::UNDEFINED;
    }

//------------------------------------------------------------------------------
}
#endif //BELFEM_MESH_ENUMS_HPP
