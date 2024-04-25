//
// Created by Christian Messe on 2018-12-28.
//

#include "meshtools.hpp"
#include "stringtools.hpp"
#include "assert.hpp"

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
                case ( ElementType::VERTEX ) :
                {
                    return 1;
                }
                case ( ElementType::LINE2 ) :
                {
                    return 2;
                }
                case ( ElementType::LINE3 ) :
                {
                    return 3;
                }
                case ( ElementType::LINE4 ) :
                {
                    return 4;
                }
                case ( ElementType::TRI3 ) :
                {
                    return 3;
                }
                case ( ElementType::TRI6 ) :
                {
                    return 6;
                }
                case ( ElementType::TRI10 ) :
                {
                    return 10;
                }
                case ( ElementType::QUAD4 ) :
                {
                    return 4;
                }
                case ( ElementType::QUAD8 ) :
                {
                    return 8;
                }
                case ( ElementType::QUAD9 ) :
                {
                    return 9;
                }
                case ( ElementType::QUAD16 ) :
                {
                    return 16;
                }
                case ( ElementType::TET4 ) :
                {
                    return 4;
                }
                case ( ElementType::TET10 ) :
                {
                    return 10;
                }
                case ( ElementType::TET20 ) :
                {
                    return 20;
                }
                case ( ElementType::PENTA6 ) :
                {
                    return 6;
                }
                case ( ElementType::PENTA15 ) :
                {
                    return 15;
                }
                case ( ElementType::PENTA18 ) :
                {
                    return 18;
                }
                case ( ElementType::HEX8 ) :
                {
                    return 8;
                }
                case ( ElementType::HEX20 ) :
                {
                    return 20;
                }
                case ( ElementType::HEX27 ) :
                {
                    return 27;
                }
                case ( ElementType::HEX64 ) :
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
                case ( ElementType::VERTEX ) :
                {
                    return 0;
                }
                case ( ElementType::LINE2 ) :
                case ( ElementType::LINE3 ) :
                case ( ElementType::LINE4 ) :
                {
                    return 1;
                }
                case ( ElementType::TRI3 ) :
                case ( ElementType::TRI6 ) :
                case ( ElementType::TRI10 ) :
                {
                    return 3;
                }
                case ( ElementType::QUAD4 ) :
                case ( ElementType::QUAD8 ) :
                case ( ElementType::QUAD9 ) :
                case ( ElementType::QUAD16 ) :
                {
                    return 4;
                }
                case ( ElementType::TET4 ) :
                case ( ElementType::TET10 ) :
                case ( ElementType::TET20 ) :
                {
                    return 6;
                }
                case ( ElementType::PENTA6 ) :
                case ( ElementType::PENTA15 ) :
                case ( ElementType::PENTA18 ) :
                {
                    return 9;
                }
                case ( ElementType::HEX8 ) :
                case ( ElementType::HEX20 ) :
                case ( ElementType::HEX27 ) :
                case ( ElementType::HEX64 ) :
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
                case ( ElementType::VERTEX ) :
                {
                    return belfem::InterpolationOrder::CONSTANT;
                }
                case ( ElementType::LINE2 ) :
                case ( ElementType::TRI3 ) :
                case ( ElementType::QUAD4 ) :
                case ( ElementType::TET4 ) :
                case ( ElementType::PENTA6 ) :
                case ( ElementType::HEX8 ) :
                {
                    return belfem::InterpolationOrder::LINEAR;
                }
                case ( ElementType::QUAD8 ) :
                case ( ElementType::PENTA15 ) :
                case ( ElementType::HEX20 ) :
                {
                    return belfem::InterpolationOrder::SERENDIPITY;
                }
                case ( ElementType::LINE3 ) :
                case ( ElementType::TRI6 ) :
                case ( ElementType::QUAD9 ) :
                case ( ElementType::TET10 ) :
                case ( ElementType::PENTA18 ) :
                case ( ElementType::HEX27 ) :
                {
                    return  belfem::InterpolationOrder::QUADRATIC;
                }
                case ( ElementType::LINE4 ) :
                case ( ElementType::TRI10 ) :
                case ( ElementType::QUAD16 ):
                case ( ElementType::TET20 ) :
                case ( ElementType::HEX64 ) :
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
                case( ElementType::EMPTY ) :
                case ( ElementType::VERTEX ) :
                {
                    return 0;
                }
                case ( ElementType::LINE2 ) :
                case ( ElementType::TRI3 ) :
                case ( ElementType::QUAD4 ) :
                case ( ElementType::TET4 ) :
                case ( ElementType::PENTA6 ) :
                case ( ElementType::HEX8 ) :
                {
                    return 1;
                }
                case ( ElementType::QUAD8 ) :
                case ( ElementType::PENTA15 ) :
                case ( ElementType::HEX20 ) :
                {
                    return 2;
                }
                case ( ElementType::LINE3 ) :
                case ( ElementType::TRI6 ) :
                case ( ElementType::QUAD9 ) :
                case ( ElementType::TET10 ) :
                case ( ElementType::PENTA18 ) :
                case ( ElementType::HEX27 ) :
                {
                    return 2;
                }
                case ( ElementType::LINE4 ) :
                case ( ElementType::TRI10 ) :
                case ( ElementType::QUAD16 ):
                case ( ElementType::TET20 ) :
                case ( ElementType::HEX64 ) :
                {
                    return 3;
                }
                default:
                {
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
                case ( ElementType::VERTEX ) :
                {
                    return belfem::GeometryType::VERTEX;
                }
                case ( ElementType::LINE2 ) :
                case ( ElementType::LINE3 ) :
                case ( ElementType::LINE4 ) :
                {
                    return belfem::GeometryType::LINE;
                }
                case ( ElementType::TRI3 ) :
                case ( ElementType::TRI6 ) :
                case ( ElementType::TRI10 ) :
                {
                    return belfem::GeometryType::TRI;
                }
                case ( ElementType::QUAD4 ) :
                case ( ElementType::QUAD8 ) :
                case ( ElementType::QUAD9 ) :
                case ( ElementType::QUAD16 ) :
                {
                    return belfem::GeometryType::QUAD;
                }
                case ( ElementType::TET4 ) :
                case ( ElementType::TET10 ) :
                case ( ElementType::TET20 ) :
                {
                    return belfem::GeometryType::TET;
                }
                case ( ElementType::PENTA6 ) :
                case ( ElementType::PENTA15 ) :
                case ( ElementType::PENTA18 ) :
                {
                    return belfem::GeometryType::PENTA;
                }
                case ( ElementType::HEX8 ) :
                case ( ElementType::HEX20 ) :
                case ( ElementType::HEX27 ) :
                case ( ElementType::HEX64 ) :
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
                case ( GeometryType::VERTEX ) :
                {
                    return 0;
                }
                case ( GeometryType::LINE ) :
                {
                    return 1;
                }
                case ( GeometryType::TRI ) :
                case ( GeometryType::QUAD ) :
                {
                    return 2;
                }
                case ( GeometryType::TET ) :
                case ( GeometryType::PENTA ) :
                case ( GeometryType::HEX ) :
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
            return static_cast< ElementType >( aGmshNumber );
        }

//------------------------------------------------------------------------------

        int
        gmsh_from_element_type( const ElementType  aElementType )
        {
            return static_cast< int >( aElementType );
        }

//------------------------------------------------------------------------------

        ElementType
        linear_element_type( const ElementType  aElementType )
        {
            switch ( aElementType )
            {
                case ( ElementType::VERTEX ) :
                {
                    return ElementType::VERTEX ;
                }
                case ( ElementType::LINE2 ) :
                case ( ElementType::LINE3 ) :
                case ( ElementType::LINE4 ) :
                case ( ElementType::LINE5 ) :
                case ( ElementType::LINE6 ) :
                {
                    return ElementType::LINE2 ;
                }
                case ( ElementType::TRI3  ) :
                case ( ElementType::TRI6  ) :
                case ( ElementType::TRI10 ) :
                case ( ElementType::TRI15 ) :
                case ( ElementType::TRI21 ) :
                {
                    return ElementType::TRI3 ;
                }
                case ( ElementType::QUAD4  ) :
                case ( ElementType::QUAD8  ) :
                case ( ElementType::QUAD9  ) :
                case ( ElementType::QUAD16 ) :
                {
                    return ElementType::QUAD4 ;
                }
                case ( ElementType::TET4  ) :
                case ( ElementType::TET10 ) :
                case ( ElementType::TET20 ) :
                case ( ElementType::TET35 ) :
                {
                    return ElementType::TET4 ;
                }
                case ( ElementType::PENTA6  ) :
                case ( ElementType::PENTA15 ) :
                case ( ElementType::PENTA18 ) :
                {
                    return ElementType::PENTA6 ;
                }
                case ( ElementType::HEX8  ) :
                case ( ElementType::HEX20 ) :
                case ( ElementType::HEX27 ) :
                case ( ElementType::HEX64 ) :
                {
                    return ElementType::HEX8 ;
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
                case( belfem::GeometryType::VERTEX ) :
                {
                    return 1;
                }
                case( belfem::GeometryType::LINE ) :
                {
                    return 2;
                }
                case( belfem::GeometryType::TRI ) :
                {
                    return 3;
                }
                case( belfem::GeometryType::QUAD ) :
                {
                    return 4;
                }
                case( belfem::GeometryType::TET ) :
                {
                    return 4;
                }
                case( belfem::GeometryType::PYRA ) :
                {
                    return 5;
                }
                case( belfem::GeometryType::PENTA ) :
                {
                    return 6;
                }
                case( belfem::GeometryType::HEX ) :
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
                case( belfem::GeometryType::VERTEX ) :
                {
                    return 0;
                }
                case( belfem::GeometryType::LINE ) :
                {
                    return 0;
                }
                case( belfem::GeometryType::TRI ) :
                {
                    return 3;
                }
                case( belfem::GeometryType::QUAD ) :
                {
                    return 4;
                }
                case( belfem::GeometryType::TET ) :
                {
                    return 4;
                }
                case( belfem::GeometryType::PYRA ) :
                {
                    return 5;
                }
                case( belfem::GeometryType::PENTA ) :
                {
                    return 5;
                }
                case( belfem::GeometryType::HEX ) :
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
        number_of_faces ( const enum ElementType  aElementType )
        {
            // get the geometry type
            switch( geometry_type( aElementType ) )
            {
                case( belfem::GeometryType::VERTEX ) :
                {
                    return 0;
                }
                case( belfem::GeometryType::LINE ) :
                {
                    return 0;
                }
                case( belfem::GeometryType::TRI ) :
                {
                    return 1;
                }
                case( belfem::GeometryType::QUAD ) :
                {
                    return 1;
                }
                case( belfem::GeometryType::TET ) :
                {
                    return 4;
                }
                case( belfem::GeometryType::PYRA ) :
                {
                    return 5;
                }
                case( belfem::GeometryType::PENTA ) :
                {
                    return 5;
                }
                case( belfem::GeometryType::HEX ) :
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
                case( ElementType::TRI3 ) :
                case( ElementType::QUAD4 ) :
                {
                    return ElementType::LINE2;
                }
                case( ElementType::TRI6 )  :
                case( ElementType::QUAD8 ) :
                case( ElementType::QUAD9 ) :
                {
                    return ElementType::LINE3;
                }
                case( ElementType::TRI10 ):
                case( ElementType::QUAD16 ):
                {
                    return ElementType::LINE4;
                }
                case( ElementType::TRI15 ) :
                {
                    return ElementType::LINE5;
                }
                case( ElementType::TRI21 ) :
                {
                    return ElementType::LINE6;
                }
                case( ElementType::TET4 ) :
                {
                    return ElementType::TRI3;
                }
                case( ElementType::TET10 ) :
                {
                    return ElementType::TRI6;
                }
                case( ElementType::TET20 ) :
                {
                    return ElementType::TRI10;
                }
                case( ElementType::HEX8 ) :
                {
                    return ElementType::QUAD4;
                }
                case( ElementType::HEX20 ) :
                {
                    return ElementType::QUAD8;
                }
                case( ElementType::HEX27 ) :
                {
                    return ElementType::QUAD9;
                }
                case( ElementType::HEX64 ) :
                {
                    return ElementType::QUAD16;
                }
                case( ElementType::PENTA6 ) :
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
                case ( ElementType::PENTA15 ) :
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
                case ( ElementType::PENTA18 ) :
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
            switch ( mesh::geometry_type( aElementType ) )
            {
                case( GeometryType::TRI ):
                case( GeometryType::QUAD ):
                {
                    return 1 ;
                }
                case( GeometryType::TET ):
                {
                    return 3 ;
                }
                case( GeometryType::PENTA ):
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
                case( GeometryType::HEX ):
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
                case( ElementType::LINE2 ) :
                {
                    return 1 ;
                }
                case( ElementType::LINE3 ) :
                {
                    return 2 ;
                }
                case( ElementType::TRI3 ) :
                {
                    return 3 ;
                }
                case( ElementType::TRI6 ) :
                {
                    return 8 ;
                }
                case( ElementType::TRI10 ) :
                {
                    return 15 ;
                }
                case( ElementType::QUAD4 ) :
                {
                    return 4 ;
                }
                case( ElementType::QUAD8 ) :
                case( ElementType::QUAD9 ) :
                {
                    return 12 ;
                }
                case( ElementType::TET4 ) :
                {
                    return 6 ;
                }
                case( ElementType::TET10 ) :
                {
                    return 20 ;
                }
                case( ElementType::TET20 ) :
                {
                    return 45 ;
                }
                case( ElementType::PENTA6 ) :
                {
                    return 9 ;
                }
                case( ElementType::PENTA15 ) :
                case( ElementType::PENTA18 ) :
                {
                    return 36 ;
                }
                case( ElementType::HEX8 ) :
                {
                    return 12 ;
                }
                case( ElementType::HEX20 ) :
                case( ElementType::HEX27 ) :
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
        element_type_from_numnodes( const uint aDimension, const uint aNumNodes )
        {
            if( aDimension == 1 )
            {
                switch( aNumNodes )
                {
                    case( 2 ) :
                    {
                        return ElementType::LINE2 ;
                    }
                    case( 3 ) :
                    {
                        return ElementType::LINE3 ;
                    }
                    case( 4 ) :
                    {
                        return ElementType::LINE4 ;
                    }
                    case( 5 ) :
                    {
                        return ElementType::LINE5 ;
                    }
                    case( 6 ) :
                    {
                        return ElementType::LINE6 ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Invalid number of nodes '%u% for dimension %u",
                                      ( unsigned int ) aDimension, ( unsigned int ) aNumNodes );
                        return ElementType::UNDEFINED ;
                    }
                }
            }
            else if ( aDimension == 2 )
            {
                switch( aNumNodes )
                {
                    case( 3 ) :
                    {
                        return ElementType::TRI3 ;
                    }
                    case( 4 ) :
                    {
                        return ElementType::QUAD4 ;
                    }
                    case( 6 ) :
                    {
                        return ElementType::TRI6 ;
                    }
                    case( 8 ) :
                    {
                        return ElementType::QUAD8 ;
                    }
                    case( 9 ) :
                    {
                        return ElementType::QUAD9 ;
                    }
                    case( 10 ) :
                    {
                        return ElementType::TRI10 ;
                    }
                    case( 15 ) :
                    {
                        return ElementType::TRI15 ;
                    }
                    case( 16 ) :
                    {
                        return ElementType::QUAD16 ;
                    }
                    case( 21 ) :
                    {
                        return ElementType::TRI21 ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Invalid number of nodes '%u% for dimension %u",
                                      ( unsigned int ) aDimension, ( unsigned int ) aNumNodes );
                        return ElementType::UNDEFINED ;
                    }
                }
            }
            else if ( aDimension == 3 )
            {
                switch( aNumNodes )
                {
                    case( 4 ) :
                    {
                        return ElementType::TET4 ;
                    }
                    case( 6 ) :
                    {
                        return ElementType::PENTA6 ;
                    }
                    case( 8 ) :
                    {
                        return ElementType::HEX8 ;
                    }
                    case( 10 ) :
                    {
                        return ElementType::TET10 ;
                    }
                    case( 15 ) :
                    {
                        return ElementType::PENTA15 ;
                    }
                    case( 18 ) :
                    {
                        return ElementType::PENTA18 ;
                    }
                    case( 20 ) :
                    {
                        BELFEM_ERROR( false, "Ambiguous number of nodes '%u% for dimension %u, could be either TET20 or HEX20" );
                        return ElementType::UNDEFINED ;
                    }
                    case( 27 ) :
                    {
                        return ElementType::HEX27 ;
                    }
                    case( 35 ) :
                    {
                        return ElementType::TET35 ;
                    }
                    case( 64 ) :
                    {
                        return ElementType::HEX64 ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Invalid number of nodes '%u% for dimension %u",
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
    }
}