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

#include "meshtools.hpp"
#include "stringtools.hpp"
#include "assert.hpp"
#include "cl_Element.hpp"
#include "../physics/materials/cl_Material.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        unsigned int
        number_of_nodes(const enum ElementType aElementType )
        {
            switch ( aElementType )
            {
                case ElementType::VERTEX :
                {
                    return 1;
                }
                case ElementType::LINE2 :
                {
                    return 2;
                }
                case ElementType::LINE3 :
                {
                    return 3;
                }
                case ElementType::LINE4 :
                {
                    return 4;
                }
                case ElementType::TRI3 :
                {
                    return 3;
                }
                case ElementType::TRI6 :
                {
                    return 6;
                }
                case ElementType::TRI10 :
                {
                    return 10;
                }
                case ElementType::QUAD4 :
                case ElementType::QUAD4TS :
                {
                    return 4;
                }
                case ElementType::QUAD8 :
                {
                    return 8;
                }
                case ElementType::QUAD9 :
                {
                    return 9;
                }
                case ElementType::QUAD16 :
                {
                    return 16;
                }
                case ElementType::TET4 :
                {
                    return 4;
                }
                case ElementType::TET10 :
                {
                    return 10;
                }
                case ElementType::TET20 :
                {
                    return 20;
                }
                case ElementType::PENTA6 :
                case ElementType::PENTA6TS :
                {
                    return 6;
                }
                case ElementType::PENTA15 :
                {
                    return 15;
                }
                case ElementType::PENTA18 :
                case ElementType::PENTA18TS :
                {
                    return 18;
                }
                case ElementType::PYRA5 :
                {
                    return  5 ;
                }
                case ElementType::PYRA13 :
                {
                    return 13 ;
                }
                case ElementType::PYRA14 :
                {
                    return 14 ;
                }
                case ElementType::HEX8 :
                case ElementType::HEX8TS :
                case ElementType::HEX8TB :
                {
                    return 8;
                }
                case ElementType::HEX20 :
                {
                    return 20;
                }
                case ElementType::HEX27 :
                {
                    return 27;
                }
                case ElementType::HEX64 :
                {
                    return 64;
                }
                default:
                {
                    return 0;
                }
            }
        }

//------------------------------------------------------------------------------

        unsigned int
        number_of_edges(const enum ElementType aElementType )
        {
            switch ( aElementType )
            {
                case ElementType::VERTEX :
                {
                    return 0;
                }
                case ElementType::LINE2 :
                case ElementType::LINE3 :
                case ElementType::LINE4 :
                {
                    return 1;
                }
                case ElementType::QUAD4TS :
                {
                    return 2;
                }
                case ElementType::TRI3 :
                case ElementType::TRI6 :
                case ElementType::TRI10 :
                {
                    return 3;
                }
                case ElementType::QUAD4 :
                case ElementType::QUAD8 :
                case ElementType::QUAD9 :
                case ElementType::QUAD16 :
                {
                    return 4;
                }
                case ElementType::TET4 :
                case ElementType::TET10 :
                case ElementType::TET20 :
                case ElementType::PENTA6TS :
                {
                    return 6;
                }
                case ElementType::HEX8TS :
                {
                    return 8;
                }
                case ElementType::HEX8TB :
                {
                    return 4 ;
                }
                case ElementType::PENTA6 :
                case ElementType::PENTA15 :
                case ElementType::PENTA18 :
                case ElementType::PENTA18TS :
                {
                    return 9;
                }
                case ElementType::HEX8 :
                case ElementType::HEX20 :
                case ElementType::HEX27 :
                case ElementType::HEX64 :
                {
                    return 12;
                }
                default:
                {
                    return 0;
                }
            }
        }

//------------------------------------------------------------------------------

        belfem::InterpolationOrder
        interpolation_order(const enum ElementType  aElementType )
        {
            switch ( aElementType )
            {
                case ElementType::VERTEX :
                {
                    return InterpolationOrder::CONSTANT;
                }
                case ElementType::LINE2 :
                case ElementType::TRI3 :
                case ElementType::QUAD4 :
                case ElementType::QUAD4TS :
                case ElementType::TET4 :
                case ElementType::PENTA6 :
                case ElementType::PENTA6TS :
                case ElementType::HEX8 :
                case ElementType::HEX8TS :
                case ElementType::HEX8TB :
                {
                    return InterpolationOrder::LINEAR;
                }
                case ElementType::QUAD8 :
                case ElementType::PENTA15 :
                case ElementType::HEX20 :
                {
                    return InterpolationOrder::SERENDIPITY;
                }
                case ElementType::LINE3 :
                case ElementType::TRI6 :
                case ElementType::QUAD9 :
                case ElementType::QUAD9TS :
                case ElementType::TET10 :
                case ElementType::PENTA18 :
                case ElementType::PENTA18TS :
                case ElementType::HEX27 :
                {
                    return  belfem::InterpolationOrder::QUADRATIC;
                }
                case ElementType::LINE4 :
                case ElementType::TRI10 :
                case ElementType::QUAD16 :
                case ElementType::TET20 :
                case ElementType::HEX64 :
                {
                    return belfem::InterpolationOrder::CUBIC;
                }
                default:
                {
                    return belfem::InterpolationOrder::UNDEFINED;
                }
            }
        }

//------------------------------------------------------------------------------

        unsigned int
        interpolation_order_numeric(const enum ElementType  aElementType )
        {
            switch ( aElementType )
            {
                case ElementType::EMPTY :
                case ElementType::VERTEX :
                {
                    return 0;
                }
                case ElementType::LINE2 :
                case ElementType::TRI3 :
                case ElementType::QUAD4 :
                case ElementType::TET4 :
                case ElementType::PENTA6 :
                case ElementType::PYRA5 :
                case ElementType::QUAD4TS :
                case ElementType::PENTA6TS :
                case ElementType::HEX8 :
                case ElementType::HEX8TS :
                case ElementType::HEX8TB :
                {
                    return 1;
                }
                case ElementType::QUAD8 :
                case ElementType::PENTA15 :
                case ElementType::HEX20 :
                case ElementType::PYRA13 :
                case ElementType::PYRA14 :
                {
                    return 2;
                }
                case ElementType::LINE3 :
                case ElementType::TRI6 :
                case ElementType::QUAD9 :
                case ElementType::TET10 :
                case ElementType::PENTA18 :
                case ElementType::HEX27 :
                case ElementType::QUAD9TS :
                case ElementType::PENTA18TS :
                {
                    return 2;
                }
                case ElementType::LINE4 :
                case ElementType::TRI10 :
                case ElementType::QUAD16 :
                case ElementType::TET20 :
                case ElementType::HEX64 :
                {
                    return 3;
                }
                case ElementType::TET35 :
                {
                    return 4;
                }
                default:
                {
                    BELFEM_ERROR( false, "Unknown element type");
                    return BELFEM_UINT_MAX;
                }
            }
        }

//------------------------------------------------------------------------------

        belfem::GeometryType
        geometry_type( const enum ElementType aElementType )
        {
            switch ( aElementType )
            {
                case ElementType::VERTEX :
                {
                    return belfem::GeometryType::VERTEX;
                }
                case ElementType::LINE2 :
                case ElementType::LINE3 :
                case ElementType::LINE4 :
                {
                    return belfem::GeometryType::LINE;
                }
                case ElementType::TRI3 :
                case ElementType::TRI6 :
                case ElementType::TRI10 :
                case ElementType::TRI15 :
                case ElementType::TRI21 :
                {
                    return belfem::GeometryType::TRI;
                }
                case ElementType::QUAD4 :
                case ElementType::QUAD8 :
                case ElementType::QUAD9 :
                case ElementType::QUAD16 :
                case ElementType::QUAD4TS :
                case ElementType::QUAD9TS :
                {
                    return belfem::GeometryType::QUAD;
                }
                case ElementType::TET4 :
                case ElementType::TET10 :
                case ElementType::TET20 :
                case ElementType::TET35 :
                {
                    return belfem::GeometryType::TET;
                }
                case ElementType::PENTA6 :
                case ElementType::PENTA15 :
                case ElementType::PENTA18 :
                case ElementType::PENTA6TS :
                case ElementType::PENTA18TS :
                {
                    return belfem::GeometryType::PENTA;
                }
                case ElementType::PYRA5 :
                case ElementType::PYRA13 :
                case ElementType::PYRA14 :
                {
                    return belfem::GeometryType::PYRA ;
                }
                case ElementType::HEX8 :
                case ElementType::HEX20 :
                case ElementType::HEX27 :
                case ElementType::HEX64 :
                case ElementType::HEX8TS :
                case ElementType::HEX8TB :
                {
                    return belfem::GeometryType::HEX;
                }
                default:
                {
                    return belfem::GeometryType::UNDEFINED;
                }
            }
        }
//------------------------------------------------------------------------------

        int
        dimension( const enum GeometryType aGeometryType )
        {
            switch ( aGeometryType )
            {
                case GeometryType::VERTEX :
                {
                    return 0;
                }
                case GeometryType::LINE :
                {
                    return 1;
                }
                case GeometryType::TRI :
                case GeometryType::QUAD :
                {
                    return 2;
                }
                case GeometryType::TET :
                case GeometryType::PENTA :
                case GeometryType::PYRA :
                case GeometryType::HEX :
                {
                    return 3;
                }
                default:
                {
                    return -1;
                }
            }
        }

//------------------------------------------------------------------------------

        int
        dimension( const enum ElementType  aElementType )
        {
            return belfem::mesh::dimension(
                    belfem::mesh::geometry_type( aElementType ) );
        }

//------------------------------------------------------------------------------

        ElementType
        element_type_from_gmsh( const int  aGmshNumber )
        {
            // gmsh MSH type numbers ( GmshDefines.h ). BELFEM's enumerators
            // coincide with them except QUAD16: gmsh 36 is the 16-node quad,
            // while gmsh 32 is a 22-node tet BELFEM does not have
            switch ( aGmshNumber )
            {
                case  1 : return ElementType::LINE2 ;
                case  2 : return ElementType::TRI3 ;
                case  3 : return ElementType::QUAD4 ;
                case  4 : return ElementType::TET4 ;
                case  5 : return ElementType::HEX8 ;
                case  6 : return ElementType::PENTA6 ;
                case  7 : return ElementType::PYRA5 ;
                case  8 : return ElementType::LINE3 ;
                case  9 : return ElementType::TRI6 ;
                case 10 : return ElementType::QUAD9 ;
                case 11 : return ElementType::TET10 ;
                case 12 : return ElementType::HEX27 ;
                case 13 : return ElementType::PENTA18 ;
                case 14 : return ElementType::PYRA14 ;
                case 15 : return ElementType::VERTEX ;
                case 16 : return ElementType::QUAD8 ;
                case 17 : return ElementType::HEX20 ;
                case 18 : return ElementType::PENTA15 ;
                case 19 : return ElementType::PYRA13 ;
                case 21 : return ElementType::TRI10 ;
                case 23 : return ElementType::TRI15 ;
                case 25 : return ElementType::TRI21 ;
                case 26 : return ElementType::LINE4 ;
                case 27 : return ElementType::LINE5 ;
                case 28 : return ElementType::LINE6 ;
                case 29 : return ElementType::TET20 ;
                case 30 : return ElementType::TET35 ;
                case 36 : return ElementType::QUAD16 ;
                case 92 : return ElementType::HEX64 ;
                default :
                {
                    BELFEM_ERROR( false,
                        "gmsh element type %d is not supported by BELFEM",
                        aGmshNumber );
                    return ElementType::UNDEFINED ;
                }
            }
        }

//------------------------------------------------------------------------------

        int
        gmsh_from_element_type( const ElementType  aElementType )
        {
            switch ( aElementType )
            {
                case ElementType::QUAD16 : return 36 ;
                case ElementType::QUAD4TS :
                case ElementType::QUAD9TS :
                case ElementType::HEX8TS :
                case ElementType::PENTA6TS :
                case ElementType::PENTA18TS :
                case ElementType::HEX8TB :
                case ElementType::UNDEFINED :
                case ElementType::EMPTY :
                {
                    BELFEM_ERROR( false,
                        "element type %d has no gmsh equivalent",
                        static_cast< int >( aElementType ) );
                    return -1 ;
                }
                default : return static_cast< int >( aElementType );
            }
        }

//------------------------------------------------------------------------------

        ElementType
        linear_element_type( const ElementType  aElementType )
        {
            switch ( aElementType )
            {
                case ElementType::VERTEX :
                {
                    return ElementType::VERTEX ;
                }
                case ElementType::LINE2 :
                case ElementType::LINE3 :
                case ElementType::LINE4 :
                case ElementType::LINE5 :
                case ElementType::LINE6 :
                {
                    return ElementType::LINE2 ;
                }
                case ElementType::TRI3  :
                case ElementType::TRI6  :
                case ElementType::TRI10 :
                case ElementType::TRI15 :
                case ElementType::TRI21 :
                {
                    return ElementType::TRI3 ;
                }
                case ElementType::QUAD4  :
                case ElementType::QUAD8  :
                case ElementType::QUAD9  :
                case ElementType::QUAD16 :
                case ElementType::QUAD4TS  :
                case ElementType::QUAD9TS  :
                {
                    return ElementType::QUAD4 ;
                }
                case ElementType::TET4  :
                case ElementType::TET10 :
                case ElementType::TET20 :
                case ElementType::TET35 :
                {
                    return ElementType::TET4 ;
                }
                case ElementType::PENTA6  :
                case ElementType::PENTA15 :
                case ElementType::PENTA18 :
                case ElementType::PENTA6TS :
                case ElementType::PENTA18TS :
                {
                    return ElementType::PENTA6 ;
                }
                case ElementType::PYRA5 :
                case ElementType::PYRA13 :
                case ElementType::PYRA14 :
                {
                    return ElementType::PYRA5 ;
                }
                case ElementType::HEX8  :
                case ElementType::HEX20 :
                case ElementType::HEX27 :
                case ElementType::HEX64 :
                {
                    return ElementType::HEX8 ;
                }
                case ElementType::HEX8TS :
                {
                    // HEX8TS is already linear; it's its own linear type.
                    return ElementType::HEX8TS ;
                }
                case ElementType::HEX8TB :
                {
                    // same for the thin beam variant
                    return ElementType::HEX8TB ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Unknown element type");
                    return ElementType::EMPTY ;
                }
            }
        }

//------------------------------------------------------------------------------

        unsigned int
        number_of_corner_nodes(  const enum ElementType  aElementType )
        {
            // get the geometry type
            belfem::GeometryType tGeoType = geometry_type( aElementType );

            switch( tGeoType )
            {
                case belfem::GeometryType::VERTEX :
                {
                    return 1;
                }
                case belfem::GeometryType::LINE :
                {
                    return 2;
                }
                case belfem::GeometryType::TRI :
                {
                    return 3;
                }
                case belfem::GeometryType::QUAD :
                {
                    return 4;
                }
                case belfem::GeometryType::TET :
                {
                    return 4;
                }
                case belfem::GeometryType::PYRA :
                {
                    return 5;
                }
                case belfem::GeometryType::PENTA :
                {
                    return 6;
                }
                case belfem::GeometryType::HEX :
                {
                    return 8;
                }
                default:
                {
                    return 0;
                }
            }
        }

//------------------------------------------------------------------------------

        unsigned int
        number_of_facets(  const enum ElementType  aElementType )
        {
            // get the geometry type
            switch(  geometry_type( aElementType ) )
            {
                case belfem::GeometryType::VERTEX :
                {
                    return 0;
                }
                case belfem::GeometryType::LINE :
                {
                    return 0;
                }
                case belfem::GeometryType::TRI :
                {
                    return 3;
                }
                case belfem::GeometryType::QUAD :
                {
                    return 4;
                }
                case belfem::GeometryType::TET :
                {
                    return 4;
                }
                case belfem::GeometryType::PYRA :
                {
                    return 5;
                }
                case belfem::GeometryType::PENTA :
                {
                    return 5;
                }
                case belfem::GeometryType::HEX :
                {
                    return 6;
                }
                default:
                {
                    return 0;
                }
            }
        }

//------------------------------------------------------------------------------

        unsigned int
        bottom_facet_index( const enum ElementType aElementType )
        {
            switch ( geometry_type( aElementType ) )
            {
                // 2D thin shells: bottom = facet 0 (CCW mirror of volume QUAD)
                case GeometryType::QUAD : return 0 ;

                // 3D PENTA thin shells: bottom triangle sits at the
                // first tri slot after the three side quads
                case GeometryType::PENTA : return 3 ;

                // 3D HEX thin shell: bottom quad sits at slot 4
                case GeometryType::HEX : return 4 ;

                default :
                {
                    BELFEM_ERROR( false,
                        "bottom_facet_index: element type %d is not a thin-shell element",
                        ( int ) aElementType );
                    return 0 ; // unreachable
                }
            }
        }

//------------------------------------------------------------------------------

        unsigned int
        top_facet_index( const enum ElementType aElementType )
        {
            switch ( geometry_type( aElementType ) )
            {
                // 2D thin shells: top = facet 2 (CCW mirror of volume QUAD)
                case GeometryType::QUAD : return 2 ;

                // 3D PENTA thin shells: top triangle sits at slot 4
                case GeometryType::PENTA : return 4 ;

                // 3D HEX thin shell: top quad sits at slot 5
                case GeometryType::HEX : return 5 ;

                default :
                {
                    BELFEM_ERROR( false,
                        "top_facet_index: element type %d is not a thin-shell element",
                        ( int ) aElementType );
                    return 0 ; // unreachable
                }
            }
        }

//------------------------------------------------------------------------------

        unsigned int
        number_of_faces ( const enum ElementType  aElementType )
        {
            // Thin-shell "face" count matches the element template's F
            // parameter (midsurface and, for quadratic shells, layer faces)
            // — NOT the volume element's boundary-face count. Keep these
            // special cases in sync with each element header's template.
            switch ( aElementType )
            {
                case ElementType::QUAD4TS :     // <4,4,2,4,1>
                case ElementType::QUAD9TS :     // <9,4,3,4,1>
                case ElementType::PENTA6TS :    // <6,6,6,5,1>
                case ElementType::HEX8TS :      // <8,8,8,6,1>
                case ElementType::HEX8TB :      // <8,8,4,6,1>
                    return 1 ;
                case ElementType::PENTA18TS :   // <18,6,9,5,3>
                    return 3 ;
                default :
                    break ;
            }

            // get the geometry type
            switch( geometry_type( aElementType ) )
            {
                case belfem::GeometryType::VERTEX :
                {
                    return 0;
                }
                case belfem::GeometryType::LINE :
                {
                    return 0;
                }
                case belfem::GeometryType::TRI :
                {
                    return 1;
                }
                case belfem::GeometryType::QUAD :
                {
                    return 1;
                }
                case belfem::GeometryType::TET :
                {
                    return 4;
                }
                case belfem::GeometryType::PYRA :
                {
                    return 5;
                }
                case belfem::GeometryType::PENTA :
                {
                    return 5;
                }
                case belfem::GeometryType::HEX :
                {
                    return 6;
                }
                default:
                {
                    return 0;
                }
            }
        }

//------------------------------------------------------------------------------

        ElementType
        element_type_of_facet(
                const ElementType  aElementType,
                const  unsigned int  aFacetIndex )
        {
            // get the geometry type
            switch ( aElementType )
            {
                case ElementType::TRI3 :
                case ElementType::QUAD4 :
                {
                    return ElementType::LINE2;
                }
                case ElementType::TRI6  :
                case ElementType::QUAD8 :
                case ElementType::QUAD9 :
                {
                    return ElementType::LINE3;
                }
                case ElementType::TRI10 :
                case ElementType::QUAD16 :
                {
                    return ElementType::LINE4;
                }
                case ElementType::TRI15 :
                {
                    return ElementType::LINE5;
                }
                case ElementType::TRI21 :
                {
                    return ElementType::LINE6;
                }
                case ElementType::TET4 :
                {
                    return ElementType::TRI3;
                }
                case ElementType::TET10 :
                {
                    return ElementType::TRI6;
                }
                case ElementType::TET20 :
                {
                    return ElementType::TRI10;
                }
                case ElementType::HEX8 :
                {
                    return ElementType::QUAD4;
                }
                case ElementType::HEX20 :
                {
                    return ElementType::QUAD8;
                }
                case ElementType::HEX27 :
                {
                    return ElementType::QUAD9;
                }
                case ElementType::HEX64 :
                {
                    return ElementType::QUAD16;
                }
                case ElementType::PENTA6 :
                {
                    if( aFacetIndex < 3 )
                    {
                        return ElementType::QUAD4;
                    }
                    else
                    {
                        return ElementType::TRI3;
                    }
                }
                case ElementType::PENTA15 :
                {
                    if( aFacetIndex < 3 )
                    {
                        return ElementType::QUAD8;
                    }
                    else
                    {
                        return ElementType::TRI6;
                    }
                }
                case ElementType::PENTA18 :
                {
                    if( aFacetIndex < 3 )
                    {
                        return ElementType::QUAD9;
                    }
                    else
                    {
                        return ElementType::TRI6;
                    }
                }
                case ElementType::PYRA5 :
                {
                    if( aFacetIndex < 4 )
                    {
                        return ElementType::TRI3;
                    }
                    else
                    {
                        return ElementType::QUAD4;
                    }
                }
                case ElementType::PYRA13 :
                {
                    if( aFacetIndex < 4 )
                    {
                        return ElementType::TRI6;
                    }
                    else
                    {
                        return ElementType::QUAD8;
                    }
                }
                case ElementType::PYRA14 :
                {
                    if( aFacetIndex < 4 )
                    {
                        return ElementType::TRI6;
                    }
                    else
                    {
                        return ElementType::QUAD9;
                    }
                }
                default:
                {
                    return ElementType::UNDEFINED;
                }
            }
        }

//------------------------------------------------------------------------------

        uint
        number_of_orientations( const ElementType aElementType,
                                const uint        aFacetNumber )
        {
            // PENTA*TS and HEX8TS now share their facet topology with the
            // corresponding volume element, so they fall through to the
            // geometry-type switch below. QUAD*TS is the only thin shell
            // that needs a dedicated entry: its LINE facets need orientation
            // tracking (2 per facet) where volume TRI/QUAD return 1.
            switch ( aElementType )
            {
                case ElementType::QUAD4TS :
                case ElementType::QUAD9TS :
                    return 2 ;
                default:
                    break ;
            }

            switch ( mesh::geometry_type( aElementType ) )
            {
                case GeometryType::TRI :
                case GeometryType::QUAD :
                {
                    return 1 ;
                }
                case GeometryType::TET :
                {
                    return 3 ;
                }
                case GeometryType::PENTA :
                {
                    if( aFacetNumber < 3 )
                    {
                        return  4 ;
                    }
                    else
                    {
                        return  3 ;
                    }
                }
                case GeometryType::PYRA :
                {
                    if( aFacetNumber < 4 )
                    {
                        return  3 ;
                    }
                    else
                    {
                        return  4 ;
                    }
                }
                case GeometryType::HEX :
                {
                    return 4 ;
                }
                default:
                {
                    BELFEM_ERROR( false, "unsupported geometry type");
                    return 0 ;
                }
            }
        }

//------------------------------------------------------------------------------

        uint
        number_of_nedelec_dofs( const ElementType aElementType )
        {
            switch ( aElementType )
            {
                case ElementType::LINE2 :
                {
                    return 1 ;
                }
                case ElementType::LINE3 :
                {
                    return 2 ;
                }
                case ElementType::TRI3 :
                {
                    return 3 ;
                }
                case ElementType::TRI6 :
                {
                    return 8 ;
                }
                case ElementType::TRI10 :
                {
                    return 15 ;
                }
                case ElementType::QUAD4 :
                {
                    return 4 ;
                }
                case ElementType::QUAD4TS :
                {
                    return 2 ;
                }
                case ElementType::QUAD8 :
                case ElementType::QUAD9 :
                {
                    return 12 ;
                }
                case ElementType::TET4 :
                {
                    return 6 ;
                }
                case ElementType::TET10 :
                {
                    return 20 ;
                }
                case ElementType::TET20 :
                {
                    return 45 ;
                }
                case ElementType::PENTA6 :
                {
                    return 9 ;
                }
                case ElementType::PENTA6TS :
                {
                    return 6 ;
                }
                case ElementType::PENTA15 :
                case ElementType::PENTA18 :
                {
                    return 36 ;
                }
                case ElementType::HEX8 :
                {
                    return 12 ;
                }
                case ElementType::HEX8TS :
                {
                    return 8 ;
                }
                case ElementType::HEX8TB :
                {
                    return 4 ;
                }
                case ElementType::HEX20 :
                case ElementType::HEX27 :
                {
                    return 54 ;
                }
                default:
                {
                    BELFEM_ERROR( false, "unsupported element type: %s", to_string( aElementType ).c_str() );
                    return 0 ;
                }
            }
        }

//------------------------------------------------------------------------------

        ElementType
        element_type_from_numnodes( const uint aDimension, const uint aNumNodes, const bool aAsThinshell )
        {
            // special check for thinshells
            if ( aAsThinshell )
            {
                switch( aNumNodes )
                {
                    case 4 :
                    {
                        return ElementType::QUAD4TS ;
                    }
                    case 9 :
                    {
                        return ElementType::QUAD9TS ;
                    }
                    case 6 :
                    {
                        return ElementType::PENTA6TS ;
                    }
                    case 18 :
                    {
                        return ElementType::PENTA18TS ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Invalid number of nodes '%u' for thin shell element of dimension %u",
                                      ( unsigned int ) aDimension, ( unsigned int ) aNumNodes );
                        return ElementType::UNDEFINED ;
                    }
                }
            }

            if( aDimension == 1 )
            {
                switch( aNumNodes )
                {
                    case 2 :
                    {
                        return ElementType::LINE2 ;
                    }
                    case 3 :
                    {
                        return ElementType::LINE3 ;
                    }
                    case 4 :
                    {
                        return ElementType::LINE4 ;
                    }
                    case 5 :
                    {
                        return ElementType::LINE5 ;
                    }
                    case 6 :
                    {
                        return ElementType::LINE6 ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Invalid number of nodes '%u' for dimension %u",
                                      ( unsigned int ) aDimension, ( unsigned int ) aNumNodes );
                        return ElementType::UNDEFINED ;
                    }
                }
            }
            else if ( aDimension == 2 )
            {
                switch( aNumNodes )
                {
                    case 3 :
                    {
                        return ElementType::TRI3 ;
                    }
                    case 4 :
                    {
                        return ElementType::QUAD4 ;
                    }
                    case 6 :
                    {
                        return ElementType::TRI6 ;
                    }
                    case 8 :
                    {
                        return ElementType::QUAD8 ;
                    }
                    case 9 :
                    {
                        return ElementType::QUAD9 ;
                    }
                    case 10 :
                    {
                        return ElementType::TRI10 ;
                    }
                    case 15 :
                    {
                        return ElementType::TRI15 ;
                    }
                    case 16 :
                    {
                        return ElementType::QUAD16 ;
                    }
                    case 21 :
                    {
                        return ElementType::TRI21 ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Invalid number of nodes '%u' for dimension %u",
                                      ( unsigned int ) aDimension, ( unsigned int ) aNumNodes );
                        return ElementType::UNDEFINED ;
                    }
                }
            }
            else if ( aDimension == 3 )
            {
                switch( aNumNodes )
                {
                    case 4 :
                    {
                        return ElementType::TET4 ;
                    }
                    case 5 :
                    {
                        return ElementType::PYRA5 ;
                    }
                    case 6 :
                    {
                        return ElementType::PENTA6 ;
                    }
                    case 8 :
                    {
                        return ElementType::HEX8 ;
                    }
                    case 10 :
                    {
                        return ElementType::TET10 ;
                    }
                    case 13 :
                    {
                        return ElementType::PYRA13 ;
                    }
                    case 14 :
                    {
                        return ElementType::PYRA14 ;
                    }
                    case 15 :
                    {
                        return ElementType::PENTA15 ;
                    }
                    case 18 :
                    {
                        return ElementType::PENTA18 ;
                    }
                    case 20 :
                    {
                        BELFEM_ERROR( false, "Ambiguous number of nodes '%u' for dimension %u, could be either TET20 or HEX20",
                                      ( unsigned int ) aDimension, ( unsigned int ) aNumNodes );
                        return ElementType::UNDEFINED ;
                    }
                    case 27 :
                    {
                        return ElementType::HEX27 ;
                    }
                    case 35 :
                    {
                        return ElementType::TET35 ;
                    }
                    case 64 :
                    {
                        return ElementType::HEX64 ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Invalid number of nodes '%u' for dimension %u",
                                      ( unsigned int ) aDimension, ( unsigned int ) aNumNodes );
                        return ElementType::UNDEFINED ;
                    }
                }
            }
            else
            {

                BELFEM_ERROR( false, "Invalid number of dimensions %u",
                              ( unsigned int ) aDimension );
                return ElementType::UNDEFINED ;

            }
        }

//------------------------------------------------------------------------------

        void
        get_bottom_nodes( Element * aElement, Cell< Node * > & aNodes )
        {
            switch ( aElement->type() )
            {
                case ElementType::QUAD4 :
                case ElementType::QUAD4TS :
                {
                    aNodes.set_size( 2, nullptr );
                    aNodes( 0 ) = aElement->node( 0 );
                    aNodes( 1 ) = aElement->node( 1 );
                    break;
                }
                case ElementType::QUAD9 :
                case ElementType::QUAD9TS :
                {
                    aNodes.set_size( 3, nullptr );
                    aNodes( 0 ) = aElement->node( 0 );
                    aNodes( 1 ) = aElement->node( 1 );
                    aNodes( 2 ) = aElement->node( 5 );

                    break;
                }
                case ElementType::PENTA6 :
                case ElementType::PENTA6TS :
                {
                    aNodes.set_size( 3, nullptr );
                    aNodes( 0 ) = aElement->node( 0 );
                    aNodes( 1 ) = aElement->node( 1 );
                    aNodes( 2 ) = aElement->node( 2 );
                    break ;
                }
                case ElementType::PENTA15 :
                case ElementType::PENTA18 :
                case ElementType::PENTA18TS :
                {
                    aNodes.set_size( 6, nullptr );
                    aNodes( 0 ) = aElement->node( 0 );
                    aNodes( 1 ) = aElement->node( 1 );
                    aNodes( 2 ) = aElement->node( 2 );
                    aNodes( 3 ) = aElement->node( 6 );
                    aNodes( 4 ) = aElement->node( 7 );
                    aNodes( 5 ) = aElement->node( 8 );
                    break ;
                }
                case ElementType::HEX8 :
                case ElementType::HEX8TS :
                case ElementType::HEX8TB :
                {
                    aNodes.set_size( 4, nullptr );
                    aNodes( 0 ) = aElement->node( 0 );
                    aNodes( 1 ) = aElement->node( 1 );
                    aNodes( 2 ) = aElement->node( 2 );
                    aNodes( 3 ) = aElement->node( 3 );
                    break;
                }
                case ElementType::HEX20 :
                {
                    aNodes.set_size( 9, nullptr );
                    aNodes( 0 ) = aElement->node( 0 );
                    aNodes( 1 ) = aElement->node( 1 );
                    aNodes( 2 ) = aElement->node( 2 );
                    aNodes( 3 ) = aElement->node( 3 );
                    aNodes( 4 ) = aElement->node( 8 );
                    aNodes( 5 ) = aElement->node( 9 );
                    aNodes( 6 ) = aElement->node( 10 );
                    aNodes( 7 ) = aElement->node( 11 );
                    aNodes( 8 ) = aElement->node( 21 );
                    break;
                }
                case ElementType::HEX27 :
                {
                    aNodes.set_size( 9, nullptr );
                    aNodes( 0 ) = aElement->node( 0 );
                    aNodes( 1 ) = aElement->node( 1 );
                    aNodes( 2 ) = aElement->node( 2 );
                    aNodes( 3 ) = aElement->node( 3 );
                    aNodes( 4 ) = aElement->node( 8 );
                    aNodes( 5 ) = aElement->node( 9 );
                    aNodes( 6 ) = aElement->node( 10 );
                    aNodes( 7 ) = aElement->node( 11 );
                    aNodes( 8 ) = aElement->node( 21 );
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Unsupported element type: %s", to_string( aElement->type() ).c_str() );
                }
            }
        }

        void
        get_top_nodes( Element * aElement, Cell< Node * > & aNodes )
        {
            switch ( geometry_type( aElement->type() ) )
            {
                case GeometryType::QUAD:
                {
                    // the facet-2 traversal {2,3} runs opposite to the bottom
                    // curve; the thin-shell hang functions need position-wise
                    // alignment with the master orientation (node 3 sits above
                    // node 0), so swap the corner nodes
                    aElement->get_nodes_of_facet( 2, aNodes );
                    Node * tSwap = aNodes( 0 );
                    aNodes( 0 ) = aNodes( 1 );
                    aNodes( 1 ) = tSwap ;
                    break;
                }
                case GeometryType::PENTA :
                {
                    aElement->get_nodes_of_facet( 4, aNodes );
                    break;
                }
                case GeometryType::HEX  :
                {
                    aElement->get_nodes_of_facet( 5, aNodes );
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Unsupported element type: %s", to_string( aElement->type() ).c_str() );
                }
            }
        }

        void
        get_bottom_edges( Element * aElement, Cell< Edge * > & aEdges )
        {
            switch ( geometry_type( aElement->type() ) )
            {
                case GeometryType::QUAD:
                {
                    // 2D thin-shell: the bottom facet is a single LINE edge
                    // (the bottom curve). QUAD*TS only stores two Nédélec-
                    // carrying edges total (bottom + top curve), so attempting
                    // to read edge(2) or edge(3) would run past the element's
                    // edge container. edge(0) is the bottom curve by
                    // construction for QUAD4TS/QUAD9TS.
                    aEdges.set_size( 1, nullptr );
                    aEdges( 0 ) = aElement->edge( 0 );
                    break;
                }
                case GeometryType::PENTA :
                {
                    aEdges.set_size( 3, nullptr );
                    aEdges( 0 ) = aElement->edge( 0 );
                    aEdges( 1 ) = aElement->edge( 1 );
                    aEdges( 2 ) = aElement->edge( 2 );
                    break;
                }
                case GeometryType::HEX  :
                {
                    aEdges.set_size( 4, nullptr );
                    aEdges( 0 ) = aElement->edge( 0 );
                    aEdges( 1 ) = aElement->edge( 1 );
                    aEdges( 2 ) = aElement->edge( 2 );
                    aEdges( 3 ) = aElement->edge( 3 );
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Unsupported element type: %s", to_string( aElement->type() ).c_str() );
                }
            }
        }

        void
        get_top_edges( Element * aElement, Cell< Edge * > & aEdges )
        {
            switch ( geometry_type( aElement->type() ) )
            {
                case GeometryType::QUAD:
                {
                    aElement->get_edges_of_facet( 2, aEdges );
                    break;
                }
                case GeometryType::PENTA :
                {
                    aElement->get_edges_of_facet( 4, aEdges );
                    break;
                }
                case GeometryType::HEX  :
                {
                    aElement->get_edges_of_facet( 5, aEdges );
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Unsupported element type: %s", to_string( aElement->type() ).c_str() );
                }
            }
        }

//------------------------------------------------------------------------------
    }
}
