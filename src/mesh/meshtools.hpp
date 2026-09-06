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

#ifndef BELFEM_ElementTypeS_HPP
#define BELFEM_ElementTypeS_HPP
#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "Mesh_Enums.hpp"
#include "cl_Node.hpp"

namespace belfem
{
    namespace mesh
    {
        // forward declarations for the thin-shell layer helpers below.
        // Node is also forward-declared here (even though cl_Node.hpp is
        // included above) because cl_Node.hpp itself re-includes this
        // header, so the include guard can leave Node incomplete when this
        // header is read from inside cl_Node.hpp's include chain.
        class Node ;
        class Element ;
        class Edge ;

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

        /**
         * Facet index of the "lower" face of a thin-shell element
         * (the face whose outward normal points away from the shell on
         * the bottom-layer side). Dispatches on geometry type: volume
         * QUAD/PENTA/HEX types return the same slot as their thin-shell
         * twin; errors for geometries without a thin-shell variant
         * ( TRI, TET, PYRA, LINE, VERTEX ).
         */
        unsigned int
        bottom_facet_index( const enum ElementType aElementType );

//------------------------------------------------------------------------------

        /**
         * Facet index of the "upper" face of a thin-shell element
         * (the face whose outward normal points away from the shell on
         * the top-layer side). Dispatches on geometry type: volume
         * QUAD/PENTA/HEX types return the same slot as their thin-shell
         * twin; errors for geometries without a thin-shell variant
         * ( TRI, TET, PYRA, LINE, VERTEX ).
         */
        unsigned int
        top_facet_index( const enum ElementType aElementType );

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

        uint
        number_of_orientations(
                const ElementType aElementType,
                const uint        aFacetNumber );

//------------------------------------------------------------------------------

        uint
        number_of_nedelec_dofs( const ElementType aElementType );

//------------------------------------------------------------------------------

        ElementType
        element_type_from_numnodes( const uint aDimension, const uint aNumNodes, const bool aAsThinshell = false );

//------------------------------------------------------------------------------

        void
        get_bottom_nodes( Element * aElement, Cell< Node * > & aNodes );

        void
        get_top_nodes( Element * aElement, Cell< Node * > & aNodes );

        void
        get_bottom_edges( Element * aElement, Cell< Edge * > & aEdges );

        void
        get_top_edges( Element * aElement, Cell< Edge * > & aEdges );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_ElementTypeS_HPP
