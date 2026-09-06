//
// Created by Christian Messe on 25.10.19.
//

#include "typedefs.hpp"
#include "assert.hpp"
#include "Mesh_Enums.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"
#include "cl_IF_InterpolationFunctionTemplate.hpp"
#include "meshtools.hpp"
#include "lagrange/cl_IF_LINE2.hpp"
#include "lagrange/cl_IF_LINE3.hpp"
#include "lagrange/cl_IF_LINE4.hpp"
#include "lagrange/cl_IF_LINE5.hpp"
#include "lagrange/cl_IF_TRI3.hpp"
#include "lagrange/cl_IF_TRI6.hpp"
#include "lagrange/cl_IF_TRI10.hpp"
#include "lagrange/cl_IF_TRI15.hpp"
#include "lagrange/cl_IF_QUAD4.hpp"
#include "lagrange/cl_IF_QUAD8.hpp"
#include "lagrange/cl_IF_QUAD9.hpp"
#include "lagrange/cl_IF_QUAD16.hpp"
#include "lagrange/cl_IF_TET4.hpp"
#include "lagrange/cl_IF_TET10.hpp"
#include "lagrange/cl_IF_TET20.hpp"
#include "lagrange/cl_IF_TET35.hpp"
#include "lagrange/cl_IF_HEX8.hpp"
#include "lagrange/cl_IF_HEX20.hpp"
#include "lagrange/cl_IF_HEX27.hpp"
#include "lagrange/cl_IF_HEX64.hpp"

#include "lagrange/cl_IF_PYRA5.hpp"
#include "lagrange/cl_IF_PYRA13.hpp"
#include "lagrange/cl_IF_PYRA14.hpp"

#include "lagrange/cl_IF_PENTA6.hpp"
#include "lagrange/cl_IF_PENTA15.hpp"
#include "lagrange/cl_IF_PENTA18.hpp"

#include "hermite/cl_IF_BEAM.hpp"
#include "hermite/cl_IF_PLATE.hpp"
#include "bernstein/cl_IFB_LINE3.hpp"
#include "bernstein/cl_IFB_TRI6.hpp"

#include "bubble/cl_IFG_TRI3A.hpp"
#include "bubble/cl_IFG_TRI3B.hpp"
#include "bubble/cl_IFG_TRI3C.hpp"

#include "bubble/cl_IFG_TRI6A.hpp"
#include "bubble/cl_IFG_TRI6B.hpp"
#include "bubble/cl_IFG_TRI6C.hpp"

#include "bubble/cl_IFG_TET4A.hpp"
#include "bubble/cl_IFG_TET4B.hpp"
#include "bubble/cl_IFG_TET4C.hpp"
#include "bubble/cl_IFG_TET4D.hpp"

#include "bubble/cl_IFG_TET10A.hpp"
#include "bubble/cl_IFG_TET10B.hpp"
#include "bubble/cl_IFG_TET10C.hpp"
#include "bubble/cl_IFG_TET10D.hpp"

namespace belfem
{
    namespace fem
    {

//------------------------------------------------------------------------------

        InterpolationFunction *
        InterpolationFunctionFactory::create_function(
                const ElementType aElementType,
                const InterpolationType aType )
        {
            switch ( aType )
            {
                case( InterpolationType::LAGRANGE ) :
                {
                    return this->create_lagrange_function( aElementType );
                }
                case( InterpolationType::BERNSTEIN ) :
                {
                    return this->create_bernstein_function( aElementType );
                }
                case( InterpolationType::HERMITE ) :
                {
                    return this->create_hermite_function( aElementType );
                }
                default:
                {
                    BELFEM_ERROR( false, "unknown interpolation funciton type");
                    return nullptr ;
                }
            }
        }

//------------------------------------------------------------------------------

        InterpolationFunction *
        InterpolationFunctionFactory::create_lagrange_function( const ElementType aElementType )
        {
            switch( aElementType )
            {
                case ( ElementType::LINE2 ) :
                {
                    return new InterpolationFunctionTemplate<GeometryType::LINE, InterpolationType::LAGRANGE, 1, 2 >();
                }
                case ( ElementType::LINE3 ) :
                {
                    return new InterpolationFunctionTemplate<GeometryType::LINE, InterpolationType::LAGRANGE, 1, 3 >();
                }
                case ( ElementType::LINE4 ) :
                {
                    return new InterpolationFunctionTemplate<GeometryType::LINE, InterpolationType::LAGRANGE, 1, 4 >();
                }
                case ( ElementType::LINE5 ) :
                {
                    return new InterpolationFunctionTemplate<GeometryType::LINE, InterpolationType::LAGRANGE, 1, 5 >();
                }
                case( ElementType::TRI3 ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::TRI, InterpolationType::LAGRANGE, 2, 3 >();
                }
                case( ElementType::TRI6 ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::TRI, InterpolationType::LAGRANGE, 2, 6 >();
                }
                case( ElementType::TRI10 ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::TRI, InterpolationType::LAGRANGE, 2, 10 >();
                }
                case( ElementType::TRI15 ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::TRI, InterpolationType::LAGRANGE, 2, 15 >();
                }
                case( ElementType::QUAD4 ):
                case( ElementType::QUAD4TS ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 4 >();
                }
                case( ElementType::QUAD8 ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 8 >();
                }
                case( ElementType::QUAD9 ):
                case( ElementType::QUAD9TS ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 9 >();
                }
                case( ElementType::QUAD16 ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 16 >();
                }
                case( ElementType::TET4 ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::TET, InterpolationType::LAGRANGE, 3, 4 >();
                }
                case( ElementType::TET10 ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::TET, InterpolationType::LAGRANGE, 3, 10 >();
                }
                case( ElementType::TET20 ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::TET, InterpolationType::LAGRANGE, 3, 20 >();
                }
                case( ElementType::TET35 ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::TET, InterpolationType::LAGRANGE, 3, 35 >();
                }
                case( ElementType::PENTA6 ):
                case ( ElementType::PENTA6TS ) :
                {
                    return new InterpolationFunctionTemplate<GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 6 >();
                }
                case( ElementType::PENTA15 ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 15 >();
                }
                case( ElementType::PENTA18 ):
                case ( ElementType::PENTA18TS ) :
                {
                    return new InterpolationFunctionTemplate<GeometryType::PENTA, InterpolationType::LAGRANGE, 3, 18 >();
                }
                case( ElementType::PYRA5 ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 5 >();
                }
                case( ElementType::PYRA13 ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 13 >();
                }
                case( ElementType::PYRA14 ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::PYRA, InterpolationType::LAGRANGE, 3, 14 >();
                }
                case( ElementType::HEX8 ):
                case( ElementType::HEX8TS ) :
                case( ElementType::HEX8TB ) :
                {
                    return new InterpolationFunctionTemplate<GeometryType::HEX, InterpolationType::LAGRANGE, 3, 8 >();
                }
                case( ElementType::HEX20 ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::HEX, InterpolationType::LAGRANGE, 3, 20 >();
                }
                case( ElementType::HEX27 ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::HEX, InterpolationType::LAGRANGE, 3, 27 >();
                }
                case( ElementType::HEX64 ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::HEX, InterpolationType::LAGRANGE, 3, 64 >();
                }
                default :
                {
                    BELFEM_ERROR( false,
                            "no lagrange function available for selected ElementType" );
                    return nullptr;
                }
            }
        }

        InterpolationFunction *
        InterpolationFunctionFactory::create_bernstein_function( const ElementType aElementType )
        {
            switch( aElementType )
            {
                case ( ElementType::LINE2 ) :
                {
                    return new InterpolationFunctionTemplate<GeometryType::LINE, InterpolationType::LAGRANGE, 1, 2 >();
                }
                case ( ElementType::LINE3 ) :
                {
                    return new InterpolationFunctionTemplate<GeometryType::LINE, InterpolationType::BERNSTEIN, 1, 3 >();
                }
                case( ElementType::TRI3 ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::TRI, InterpolationType::LAGRANGE, 2, 3 >();
                }
                case( ElementType::TRI6 ):
                {
                    return new InterpolationFunctionTemplate<GeometryType::TRI, InterpolationType::BERNSTEIN, 2, 6 >();
                }
                default :
                {
                    BELFEM_ERROR( false,
                                  "no bernstein function available for selected ElementType" );
                    return nullptr;
                }
            }
        }

//------------------------------------------------------------------------------

        InterpolationFunction *
        InterpolationFunctionFactory::create_hermite_function( const ElementType aElementType )
        {
            switch( aElementType )
            {
                case ( ElementType::LINE2 ) :
                {
                    return new InterpolationFunctionTemplate<GeometryType::LINE, InterpolationType::HERMITE, 1, 4 >();
                }
                case ( ElementType::QUAD4 ) :
                {
                    return new InterpolationFunctionTemplate<GeometryType::QUAD, InterpolationType::HERMITE, 2, 16 >();
                }
                default :
                {
                    BELFEM_ERROR( false,
                                 "no hermite function available for selected ElementType" );
                    return nullptr;
                }
            }
        }
        
//------------------------------------------------------------------------------
        
        InterpolationFunction *
        InterpolationFunctionFactory::create_bubble_function( const ElementType aElementType, const uint aFacet )
        {
            switch( aElementType )
            {
                case( ElementType::TRI3 ) :
                {
                    switch ( aFacet )
                    {
                        case( 0 ):
                        {
                            return new InterpolationFunctionTemplate<GeometryType::TRI, InterpolationType::BubbleEdge0, 2, 1 >();
                        }
                        case( 1 ):
                        {
                            return new InterpolationFunctionTemplate<GeometryType::TRI, InterpolationType::BubbleEdge1, 2, 1 >();
                        }
                        case( 2 ):
                        {
                            return new InterpolationFunctionTemplate<GeometryType::TRI, InterpolationType::BubbleEdge2, 2, 1 >();
                        }
                        default:
                        {
                            BELFEM_ERROR( false, "invalid facet index for triangle: %u", ( unsigned int ) aFacet );
                            return nullptr;
                        }
                    }
                } // end case tri3
                case( ElementType::TRI6 ) :
                {
                    switch ( aFacet )
                    {
                        case( 0 ):
                        {
                            return new InterpolationFunctionTemplate<GeometryType::TRI, InterpolationType::BubbleEdge0, 2, 2 >();
                        }
                        case( 1 ):
                        {
                            return new InterpolationFunctionTemplate<GeometryType::TRI, InterpolationType::BubbleEdge1, 2, 2 >();
                        }
                        case( 2 ):
                        {
                            return new InterpolationFunctionTemplate<GeometryType::TRI, InterpolationType::BubbleEdge2, 2, 2 >();
                        }
                        default:
                        {
                            BELFEM_ERROR( false, "invalid facet index for triangle: %u", ( unsigned int ) aFacet );
                            return nullptr;
                        }
                    }
                } // end case tri6
                case( ElementType::TET4 ) :
                {
                    switch ( aFacet )
                    {
                        case( 0 ):
                        {
                            return new InterpolationFunctionTemplate<GeometryType::TET, InterpolationType::BubbleFace0, 3, 1 >();
                        }
                        case( 1 ):
                        {
                            return new InterpolationFunctionTemplate<GeometryType::TET, InterpolationType::BubbleFace1, 3, 1 >();
                        }
                        case( 2 ):
                        {
                            return new InterpolationFunctionTemplate<GeometryType::TET, InterpolationType::BubbleFace2, 3, 1 >();
                        }
                        case( 3 ):
                        {
                            return new InterpolationFunctionTemplate<GeometryType::TET, InterpolationType::BubbleFace3, 3, 1 >();
                        }
                        default:
                        {
                            BELFEM_ERROR( false, "invalid facet index for tet: %u", ( unsigned int ) aFacet );
                            return nullptr;
                        }
                    }
                } // end case tet4
                case( ElementType::TET10 ) :
                {
                    switch ( aFacet )
                    {
                        case( 0 ):
                        {
                            return new InterpolationFunctionTemplate<GeometryType::TET, InterpolationType::BubbleFace0, 3, 3 >();
                        }
                        case( 1 ):
                        {
                            return new InterpolationFunctionTemplate<GeometryType::TET, InterpolationType::BubbleFace1, 3, 3 >();
                        }
                        case( 2 ):
                        {
                            return new InterpolationFunctionTemplate<GeometryType::TET, InterpolationType::BubbleFace2, 3, 3 >();
                        }
                        case( 3 ):
                        {
                            return new InterpolationFunctionTemplate<GeometryType::TET, InterpolationType::BubbleFace3, 3, 3 >();
                        }
                        default:
                        {
                            BELFEM_ERROR( false, "invalid facet index for tet: %u", ( unsigned int ) aFacet );
                            return nullptr;
                        }
                    }
                } // end case tet10
                default:
                {
                    BELFEM_ERROR( false, "no bubble function available for selected ElementType" );
                    return nullptr;
                }
            }
        }
        
//------------------------------------------------------------------------------
    }
}