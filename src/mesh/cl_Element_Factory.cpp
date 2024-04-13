//
// Created by Christian Messe on 2019-07-28.
//

#include "cl_Element_Factory.hpp"
#include "meshtools.hpp"

#include "cl_ElementTemplate.hpp"

#include "cl_Element_Vertex.hpp"
#include "cl_Element_LINE2.hpp"
#include "cl_Element_LINE3.hpp"
#include "cl_Element_LINE4.hpp"
#include "cl_Element_LINE5.hpp"
#include "cl_Element_TRI3.hpp"
#include "cl_Element_TRI6.hpp"
#include "cl_Element_TRI10.hpp"
#include "cl_Element_TRI15.hpp"
#include "cl_Element_QUAD4.hpp"
#include "cl_Element_QUAD8.hpp"
#include "cl_Element_QUAD9.hpp"
#include "cl_Element_QUAD16.hpp"
#include "cl_Element_TET4.hpp"
#include "cl_Element_TET10.hpp"
#include "cl_Element_PENTA6.hpp"
#include "cl_Element_PENTA15.hpp"
#include "cl_Element_PENTA18.hpp"
#include "cl_Element_HEX8.hpp"
#include "cl_Element_HEX20.hpp"
#include "cl_Element_HEX27.hpp"
#include "cl_Element_HEX64.hpp"
namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        Element *
        ElementFactory::create_element(
                const ElementType aType, const id_t aID ) const
        {
            switch ( aType )
            {
                case ( ElementType::VERTEX ) :
                {
                    return new ElementTemplate< 1, 1, 0, 0, 0 >( aID );
                }
                case ( ElementType::LINE2 ) :
                {
                    return new ElementTemplate< 2, 2, 1, 0, 0 >( aID );
                }
                case ( ElementType::LINE3 ) :
                {
                    return new ElementTemplate< 3, 2, 1, 0, 0 >( aID );
                }
                case ( ElementType::LINE4 ) :
                {
                    return new ElementTemplate< 4, 2, 1, 0, 0 >( aID );
                }
                case ( ElementType::LINE5 ) :
                {
                    return new ElementTemplate< 5, 2, 1, 0, 0 >( aID );
                }
                case ( ElementType::TRI3 ) :
                {
                    return new ElementTemplate< 3, 3, 3, 3, 1 >( aID );
                }
                case ( ElementType::TRI6 ) :
                {
                    return new ElementTemplate< 6, 3, 3, 3, 1 >( aID );
                }
                case ( ElementType::TRI10 ) :
                {
                    return new ElementTemplate< 10, 3, 3, 3, 1 >( aID );
                }
                case ( ElementType::TRI15 ) :
                {
                    return new ElementTemplate< 15, 3, 3, 3, 1 >( aID );
                }
                case ( ElementType::QUAD4 ) :
                {
                    return new ElementTemplate< 4, 4, 4, 4, 1 >( aID );
                }
                case ( ElementType::QUAD8 ) :
                {
                    return new ElementTemplate< 8, 4, 4, 4, 1 >( aID );
                }
                case ( ElementType::QUAD9 ) :
                {
                    return new ElementTemplate< 9, 4, 4, 4, 1 >( aID );
                }
                case ( ElementType::QUAD16 ) :
                {
                    return new ElementTemplate< 16, 4, 4, 4, 1 >( aID );
                }
                case ( ElementType::TET4 ) :
                {
                    return new ElementTemplate< 4, 4, 6, 4, 4 >( aID );
                }
                case ( ElementType::TET10 ) :
                {
                    return new ElementTemplate< 10, 4, 6, 4, 4 >( aID );
                }
                case ( ElementType::TET20 ) :
                {
                    return new ElementTemplate< 20, 4, 6, 4, 4 >( aID );
                }
                case ( ElementType::TET35 ) :
                {
                    return new ElementTemplate< 35, 4, 6, 4, 4 >( aID );
                }
                case ( ElementType::PENTA6 ) :
                {
                    return new ElementTemplate< 6, 6, 9, 5, 5 >( aID );
                }
                case ( ElementType::PENTA15 ) :
                {
                    return new ElementTemplate< 15, 6, 9, 5, 5 >( aID );
                }
                case ( ElementType::PENTA18 ) :
                {
                    return new ElementTemplate< 18, 6, 9, 5, 5 >( aID );
                }
                case ( ElementType::HEX8 ) :
                {
                    return new ElementTemplate< 8, 8, 12, 6, 6 >( aID );
                }
                case ( ElementType::HEX20 ) :
                {
                    return new ElementTemplate< 20, 8, 12, 6, 6 >( aID );
                }
                case ( ElementType::HEX27 ) :
                {
                    return new ElementTemplate< 27, 8, 12, 6, 6 >( aID );
                }
                case ( ElementType::HEX64 ) :
                {
                    return new ElementTemplate< 64, 8, 12, 6, 6 >( aID );
                }
                default:
                {
                    BELFEM_ERROR( false, "Unknown element type" );
                    return nullptr ;
                }
            }
        }

//------------------------------------------------------------------------------
    }
}
