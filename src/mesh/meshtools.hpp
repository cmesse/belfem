//
// Created by Christian Messe on 2018-12-28.
//

#ifndef BELFEM_ElementTypeS_HPP
#define BELFEM_ElementTypeS_HPP
#include <string>
#include "Mesh_Enums.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        unsigned int
        number_of_nodes( const enum ElementType aElementType );

//------------------------------------------------------------------------------

        unsigned int
        number_of_edges( const enum ElementType aElementType );
//------------------------------------------------------------------------------

        unsigned int
        number_of_corner_nodes(  const enum ElementType aElementType );

//------------------------------------------------------------------------------

        unsigned int
        number_of_facets(  const enum ElementType aElementType );

//------------------------------------------------------------------------------

        unsigned int
        number_of_faces(  const enum ElementType aElementType );

//------------------------------------------------------------------------------

        InterpolationOrder
        interpolation_order( const enum ElementType aElementType );

//------------------------------------------------------------------------------

        unsigned int
        interpolation_order_numeric( const enum ElementType aElementType );

//------------------------------------------------------------------------------

        GeometryType
        geometry_type( const enum ElementType aElementType );

//------------------------------------------------------------------------------

        int
        dimension( const enum GeometryType  aGeometryType );

//------------------------------------------------------------------------------

        int
        dimension( const enum ElementType aElementType );

//------------------------------------------------------------------------------

        ElementType
        element_type_from_gmsh( const int  aGmshNumber );

//------------------------------------------------------------------------------

        int
        gmsh_from_element_type( const ElementType aElementType );

//------------------------------------------------------------------------------

        ElementType
        linear_element_type( const ElementType aElementType );

//------------------------------------------------------------------------------

        ElementType
        element_type_of_facet(
                const ElementType  aElementType,
                const  unsigned int  aFacetIndex );


//------------------------------------------------------------------------------

        double
        length_scale_factor( const std::string  aUnitLabel );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_ElementTypeS_HPP
