//
// Created by Christian Messe on 22.06.20.
//

#include "cl_Vertex.hpp"
#include "cl_Node.hpp"
#include "cl_Element.hpp"

#include "vtktools.hpp"
#include "meshtools.hpp"
namespace belfem
{
    namespace vtk
    {
//------------------------------------------------------------------------------

        int
        vtk_type( const ElementType & aElementType )
        {
            // see vtkCellType.h in VTK source
            switch ( aElementType )
            {
                case( ElementType::VERTEX ) :
                {
                    return  1 ; // VTK_VERTEX
                }
                case( ElementType::LINE2 ) :
                {
                    return  3 ; // VTK_LINE
                }
                case( ElementType::LINE3 ) :
                {
                    return 21 ; // VTK_QUADRATIC_EDGE
                }
                case( ElementType::LINE4 ) :
                {
                    return 35 ; // VTK_CUBIC_LINE
                }
                case( ElementType::TRI3 ) :
                {
                    return  5 ;  // VTK_TRIANGLE
                }
                case( ElementType::TRI6 ) :
                {
                    return 22 ; // VTK_QUADRATIC_TRIANGLE
                }
                case( ElementType::QUAD4 ) :
                {
                    return  9 ; // VTK_QUAD
                }
                case( ElementType::QUAD8 ) :
                {
                    return 23 ; // VTK_QUADRATIC_QUAD ?
                }
                case( ElementType::QUAD9 ) :
                {
                    return 28 ; // VTK_BIQUADRATIC_QUAD
                }
                case( ElementType::TET4 ) :
                {
                    return 10 ; // VTK_TETRA
                }
                case( ElementType::TET10 ) :
                {
                    return 24 ; // VTK_QUADRATIC_TETRA
                }
                case( ElementType::PENTA6 ) :
                {
                    return 13 ; // VTK_WEDGE
                }
                case( ElementType::PENTA15 ) :
                {
                    return 26 ; // VTK_QUADRATIC_WEDGE ?
                }
                case( ElementType::PENTA18 ) :
                {
                    return 32 ; // VTK_BIQUADRATIC_QUADRATIC_WEDGE
                }
                case( ElementType::HEX8 ) :
                {
                    return 12 ; /// VTK_HEXAHEDRON
                }
                case( ElementType::HEX20 ) :
                {
                    return 25 ; // VTK_QUADRATIC_HEXAHEDRON
                }
                case( ElementType::HEX27 ) :
                {
                    return 29 ; // VTK_TRIQUADRATIC_HEXAHEDRON
                }
                default:
                {
                    BELFEM_ERROR( false, "Unknown Element Type" );
                    return 0 ;
                }
            }
        }
//------------------------------------------------------------------------------

        void
        get_node_ids( mesh::Element * aElement, Vector< id_t > & aNodeIDs )
        {
            // get element type
            ElementType tType = aElement->type() ;

            BELFEM_ASSERT( aNodeIDs.length() == mesh::number_of_nodes( tType ),
                "Lenth of aNodeIDs does not match ( is %u but should be %u )",
                          ( unsigned int ) aNodeIDs.length(),
                          ( unsigned int ) mesh::number_of_nodes( tType ) );

            // see vtkCellType.h in VTK source
            switch ( tType )
            {
                case( ElementType::VERTEX ) :
                {
                    aNodeIDs( 0 ) = aElement->node( 0 )->id() ;
                    break ;
                }
                case( ElementType::LINE2 ) :
                {
                    aNodeIDs( 0 ) = aElement->node( 0 )->id() ;
                    aNodeIDs( 1 ) = aElement->node( 1 )->id() ;
                    break ;
                }
                case( ElementType::LINE3 ) :
                {
                    aNodeIDs( 0 ) = aElement->node( 0 )->id() ;
                    aNodeIDs( 1 ) = aElement->node( 1 )->id() ;
                    aNodeIDs( 2 ) = aElement->node( 2 )->id() ;
                    break ;
                }
                case( ElementType::LINE4 ) :
                {
                    aNodeIDs( 0 ) = aElement->node( 0 )->id() ;
                    aNodeIDs( 1 ) = aElement->node( 1 )->id() ;
                    aNodeIDs( 2 ) = aElement->node( 2 )->id() ;
                    aNodeIDs( 3 ) = aElement->node( 3 )->id() ;
                    break ;
                }
                case( ElementType::TRI3 ) :
                {
                    aNodeIDs( 0 ) = aElement->node( 0 )->id() ;
                    aNodeIDs( 1 ) = aElement->node( 1 )->id() ;
                    aNodeIDs( 2 ) = aElement->node( 2 )->id() ;
                    break ;
                }
                case( ElementType::TRI6 ) :
                {
                    aNodeIDs( 0 ) = aElement->node( 0 )->id() ;
                    aNodeIDs( 1 ) = aElement->node( 1 )->id() ;
                    aNodeIDs( 2 ) = aElement->node( 2 )->id() ;
                    aNodeIDs( 3 ) = aElement->node( 3 )->id() ;
                    aNodeIDs( 4 ) = aElement->node( 4 )->id() ;
                    aNodeIDs( 5 ) = aElement->node( 5 )->id() ;
                    break ;
                }
                case( ElementType::QUAD4 ) :
                {
                    aNodeIDs( 0 ) = aElement->node( 0 )->id() ;
                    aNodeIDs( 1 ) = aElement->node( 1 )->id() ;
                    aNodeIDs( 2 ) = aElement->node( 2 )->id() ;
                    aNodeIDs( 3 ) = aElement->node( 3 )->id() ;
                    break ;
                }
                case( ElementType::QUAD8 ) :
                {
                    aNodeIDs( 0 ) = aElement->node( 0 )->id() ;
                    aNodeIDs( 1 ) = aElement->node( 1 )->id() ;
                    aNodeIDs( 2 ) = aElement->node( 2 )->id() ;
                    aNodeIDs( 3 ) = aElement->node( 3 )->id() ;
                    aNodeIDs( 4 ) = aElement->node( 4 )->id() ;
                    aNodeIDs( 5 ) = aElement->node( 5 )->id() ;
                    aNodeIDs( 6 ) = aElement->node( 6 )->id() ;
                    aNodeIDs( 7 ) = aElement->node( 7 )->id() ;
                    break ;
                }
                case( ElementType::QUAD9 ) :
                {
                    aNodeIDs( 0 ) = aElement->node( 0 )->id() ;
                    aNodeIDs( 1 ) = aElement->node( 1 )->id() ;
                    aNodeIDs( 2 ) = aElement->node( 2 )->id() ;
                    aNodeIDs( 3 ) = aElement->node( 3 )->id() ;
                    aNodeIDs( 4 ) = aElement->node( 4 )->id() ;
                    aNodeIDs( 5 ) = aElement->node( 5 )->id() ;
                    aNodeIDs( 6 ) = aElement->node( 6 )->id() ;
                    aNodeIDs( 7 ) = aElement->node( 7 )->id() ;
                    aNodeIDs( 8 ) = aElement->node( 8 )->id() ;
                    break ;
                }
                case( ElementType::TET4 ) :
                {
                    aNodeIDs( 0 ) = aElement->node( 0 )->id() ;
                    aNodeIDs( 1 ) = aElement->node( 1 )->id() ;
                    aNodeIDs( 2 ) = aElement->node( 2 )->id() ;
                    aNodeIDs( 3 ) = aElement->node( 3 )->id() ;
                    break ;
                }
                case( ElementType::TET10 ) :
                {
                    aNodeIDs( 0 ) = aElement->node( 0 )->id() ;
                    aNodeIDs( 1 ) = aElement->node( 1 )->id() ;
                    aNodeIDs( 2 ) = aElement->node( 2 )->id() ;
                    aNodeIDs( 3 ) = aElement->node( 3 )->id() ;
                    aNodeIDs( 4 ) = aElement->node( 4 )->id() ;
                    aNodeIDs( 5 ) = aElement->node( 5 )->id() ;
                    aNodeIDs( 6 ) = aElement->node( 6 )->id() ;
                    aNodeIDs( 7 ) = aElement->node( 7 )->id() ;
                    aNodeIDs( 8 ) = aElement->node( 8 )->id() ;
                    aNodeIDs( 9 ) = aElement->node( 9 )->id() ;
                    break ;
                }
                case( ElementType::PENTA6 ) :
                {
                    aNodeIDs( 0 ) = aElement->node( 0 )->id() ;
                    aNodeIDs( 1 ) = aElement->node( 1 )->id() ;
                    aNodeIDs( 2 ) = aElement->node( 2 )->id() ;
                    aNodeIDs( 3 ) = aElement->node( 3 )->id() ;
                    aNodeIDs( 4 ) = aElement->node( 4 )->id() ;
                    aNodeIDs( 5 ) = aElement->node( 5 )->id() ;
                    break ;
                }
                case( ElementType::PENTA15 ) :
                {
                    aNodeIDs(  0 ) = aElement->node(  0 )->id(); //  1
                    aNodeIDs(  1 ) = aElement->node(  2 )->id(); //  3
                    aNodeIDs(  2 ) = aElement->node(  1 )->id(); //  2
                    aNodeIDs(  3 ) = aElement->node(  3 )->id(); //  4
                    aNodeIDs(  4 ) = aElement->node(  5 )->id(); //  6
                    aNodeIDs(  5 ) = aElement->node(  4 )->id(); //  5
                    aNodeIDs(  6 ) = aElement->node(  8 )->id(); //  8
                    aNodeIDs(  7 ) = aElement->node(  7 )->id(); // 10
                    aNodeIDs(  8 ) = aElement->node(  6 )->id(); //  7
                    aNodeIDs(  9 ) = aElement->node( 14 )->id(); // 14
                    aNodeIDs( 10 ) = aElement->node( 13 )->id(); // 15
                    aNodeIDs( 11 ) = aElement->node( 12 )->id(); // 13
                    aNodeIDs( 12 ) = aElement->node(  9 )->id(); //  9
                    aNodeIDs( 13 ) = aElement->node( 11 )->id(); // 12
                    aNodeIDs( 14 ) = aElement->node( 10 )->id(); // 11
                    break ;
                }
                case( ElementType::PENTA18 ) :
                {
                    aNodeIDs(  0 ) = aElement->node(  0 )->id(); //  1
                    aNodeIDs(  1 ) = aElement->node(  2 )->id(); //  3
                    aNodeIDs(  2 ) = aElement->node(  1 )->id(); //  2
                    aNodeIDs(  3 ) = aElement->node(  3 )->id(); //  4
                    aNodeIDs(  4 ) = aElement->node(  5 )->id(); //  6
                    aNodeIDs(  5 ) = aElement->node(  4 )->id(); //  5
                    aNodeIDs(  6 ) = aElement->node(  8 )->id(); //  8
                    aNodeIDs(  7 ) = aElement->node(  7 )->id(); // 10
                    aNodeIDs(  8 ) = aElement->node(  6 )->id(); //  7
                    aNodeIDs(  9 ) = aElement->node( 14 )->id(); // 14
                    aNodeIDs( 10 ) = aElement->node( 13 )->id(); // 15
                    aNodeIDs( 11 ) = aElement->node( 12 )->id(); // 13
                    aNodeIDs( 12 ) = aElement->node(  9 )->id(); //  9
                    aNodeIDs( 13 ) = aElement->node( 11 )->id(); // 12
                    aNodeIDs( 14 ) = aElement->node( 10 )->id(); // 11
                    aNodeIDs( 15 ) = aElement->node( 17 )->id(); // 17
                    aNodeIDs( 16 ) = aElement->node( 16 )->id(); // 18
                    aNodeIDs( 17 ) = aElement->node( 15 )->id(); // 16
                    break ;
                }
                case( ElementType::HEX8 ) :
                {
                    aNodeIDs( 0 ) = aElement->node( 0 )->id() ;
                    aNodeIDs( 1 ) = aElement->node( 1 )->id() ;
                    aNodeIDs( 2 ) = aElement->node( 2 )->id() ;
                    aNodeIDs( 3 ) = aElement->node( 3 )->id() ;
                    aNodeIDs( 4 ) = aElement->node( 4 )->id() ;
                    aNodeIDs( 5 ) = aElement->node( 5 )->id() ;
                    aNodeIDs( 6 ) = aElement->node( 6 )->id() ;
                    aNodeIDs( 7 ) = aElement->node( 7 )->id() ;
                    break ;
                }
                case( ElementType::HEX20 ) :
                {
                    aNodeIDs(  0 ) =  aElement->node(  0 )->id();
                    aNodeIDs(  1 ) =  aElement->node(  1 )->id();
                    aNodeIDs(  2 ) =  aElement->node(  2 )->id();
                    aNodeIDs(  3 ) =  aElement->node(  3 )->id();
                    aNodeIDs(  4 ) =  aElement->node(  4 )->id();
                    aNodeIDs(  5 ) =  aElement->node(  5 )->id();
                    aNodeIDs(  6 ) =  aElement->node(  6 )->id();
                    aNodeIDs(  7 ) =  aElement->node(  7 )->id();
                    aNodeIDs(  8 ) =  aElement->node(  8 )->id();
                    aNodeIDs(  9 ) =  aElement->node(  9 )->id();
                    aNodeIDs( 10 ) =  aElement->node( 10 )->id();
                    aNodeIDs( 11 ) =  aElement->node( 11 )->id();
                    aNodeIDs( 12 ) =  aElement->node( 16 )->id();
                    aNodeIDs( 13 ) =  aElement->node( 17 )->id();
                    aNodeIDs( 14 ) =  aElement->node( 18 )->id();
                    aNodeIDs( 15 ) =  aElement->node( 19 )->id();
                    aNodeIDs( 16 ) =  aElement->node( 12 )->id();
                    aNodeIDs( 17 ) =  aElement->node( 13 )->id();
                    aNodeIDs( 18 ) =  aElement->node( 14 )->id();
                    aNodeIDs( 19 ) =  aElement->node( 15 )->id();
                    break ;
                }
                case( ElementType::HEX27 ) :
                {
                    aNodeIDs(  0 ) =  aElement->node(  0 )->id();
                    aNodeIDs(  1 ) =  aElement->node(  1 )->id();
                    aNodeIDs(  2 ) =  aElement->node(  2 )->id();
                    aNodeIDs(  3 ) =  aElement->node(  3 )->id();
                    aNodeIDs(  4 ) =  aElement->node(  4 )->id();
                    aNodeIDs(  5 ) =  aElement->node(  5 )->id();
                    aNodeIDs(  6 ) =  aElement->node(  6 )->id();
                    aNodeIDs(  7 ) =  aElement->node(  7 )->id();
                    aNodeIDs(  8 ) =  aElement->node(  8 )->id();
                    aNodeIDs(  9 ) =  aElement->node(  9 )->id();
                    aNodeIDs( 10 ) =  aElement->node( 10 )->id();
                    aNodeIDs( 11 ) =  aElement->node( 11 )->id();
                    aNodeIDs( 12 ) =  aElement->node( 16 )->id();
                    aNodeIDs( 13 ) =  aElement->node( 17 )->id();
                    aNodeIDs( 14 ) =  aElement->node( 18 )->id();
                    aNodeIDs( 15 ) =  aElement->node( 19 )->id();
                    aNodeIDs( 16 ) =  aElement->node( 12 )->id();
                    aNodeIDs( 17 ) =  aElement->node( 13 )->id();
                    aNodeIDs( 18 ) =  aElement->node( 14 )->id();
                    aNodeIDs( 19 ) =  aElement->node( 15 )->id();
                    aNodeIDs( 20 ) =  aElement->node( 23 )->id();
                    aNodeIDs( 21 ) =  aElement->node( 24 )->id();
                    aNodeIDs( 22 ) =  aElement->node( 25 )->id();
                    aNodeIDs( 23 ) =  aElement->node( 26 )->id();
                    aNodeIDs( 24 ) =  aElement->node( 21 )->id();
                    aNodeIDs( 25 ) =  aElement->node( 22 )->id();
                    aNodeIDs( 26 ) =  aElement->node( 20 )->id();
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Unknown Element Type" );
                    break ;
                }
            }
        }

//------------------------------------------------------------------------------
    }
}
