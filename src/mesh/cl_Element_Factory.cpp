//
// Created by Christian Messe on 2019-07-28.
//

#include "cl_Element_Factory.hpp"
#include "meshtools.hpp"

#include "cl_LagrangeElement.hpp"

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
        ElementFactory::create_lagrange_element(
                const ElementType & aType, const id_t & aID ) const
        {
            // create an element pointer
            Element * aElement = nullptr;

            switch ( aType )
            {
                case ( ElementType::VERTEX ) :
                {
                    aElement = new LagrangeElement< 1, 1, 0, 0, 0 >( aID );
                    break;
                }
                case ( ElementType::LINE2 ) :
                {
                    aElement = new LagrangeElement< 2, 2, 1, 0, 0 >( aID );
                    break;
                }
                case ( ElementType::LINE3 ) :
                {
                    aElement = new LagrangeElement< 3, 2, 1, 0, 0 >( aID );
                    break;
                }
                case ( ElementType::LINE4 ) :
                {
                    aElement = new LagrangeElement< 4, 2, 1, 0, 0 >( aID );
                    break;
                }
                case ( ElementType::LINE5 ) :
                {
                    aElement = new LagrangeElement< 5, 2, 1, 0, 0 >( aID );
                    break;
                }
                case ( ElementType::TRI3 ) :
                {
                    aElement = new LagrangeElement< 3, 3, 3, 3, 1 >( aID );
                    break;
                }
                case ( ElementType::TRI6 ) :
                {
                    aElement = new LagrangeElement< 6, 3, 3, 3, 1 >( aID );
                    break;
                }
                case ( ElementType::TRI10 ) :
                {
                    aElement = new LagrangeElement< 10, 3, 3, 3, 1 >( aID );
                    break;
                }
                case ( ElementType::TRI15 ) :
                {
                    aElement = new LagrangeElement< 15, 3, 3, 3, 1 >( aID );
                    break;
                }
                case ( ElementType::QUAD4 ) :
                {
                    aElement = new LagrangeElement< 4, 4, 4, 4, 1 >( aID );
                    break;
                }
                case ( ElementType::QUAD8 ) :
                {
                    aElement = new LagrangeElement< 8, 4, 4, 4, 1 >( aID );
                    break;
                }
                case ( ElementType::QUAD9 ) :
                {
                    aElement = new LagrangeElement< 9, 4, 4, 4, 1 >( aID );
                    break;
                }
                case ( ElementType::QUAD16 ) :
                {
                    aElement = new LagrangeElement< 16, 4, 4, 4, 1 >( aID );
                    break;
                }
                case ( ElementType::TET4 ) :
                {
                    aElement = new LagrangeElement< 4, 4, 6, 4, 4 >( aID );
                    break;
                }
                case ( ElementType::TET10 ) :
                {
                    aElement = new LagrangeElement< 10, 4, 6, 4, 4 >( aID );
                    break;
                }
                case ( ElementType::TET20 ) :
                {
                    aElement = new LagrangeElement< 20, 4, 6, 4, 4 >( aID );
                    break;
                }
                case ( ElementType::TET35 ) :
                {
                    aElement = new LagrangeElement< 35, 4, 6, 4, 4 >( aID );
                    break;
                }
                case ( ElementType::PENTA6 ) :
                {
                    aElement = new LagrangeElement< 6, 6, 9, 5, 5 >( aID );
                    break;
                }
                case ( ElementType::PENTA15 ) :
                {
                    aElement = new LagrangeElement< 15, 6, 9, 5, 5 >( aID );
                    break;
                }
                case ( ElementType::PENTA18 ) :
                {
                    aElement = new LagrangeElement< 18, 6, 9, 5, 5 >( aID );
                    break;
                }
                case ( ElementType::HEX8 ) :
                {
                    aElement = new LagrangeElement< 8, 8, 12, 6, 6 >( aID );
                    break;
                }
                case ( ElementType::HEX20 ) :
                {
                    aElement = new LagrangeElement< 20, 8, 12, 6, 6 >( aID );
                    break;
                }
                case ( ElementType::HEX27 ) :
                {
                    aElement = new LagrangeElement< 27, 8, 12, 6, 6 >( aID );
                    break;
                }
                case ( ElementType::HEX64 ) :
                {
                    aElement = new LagrangeElement< 64, 8, 12, 6, 6 >( aID );
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Unknown element type" );
                }
            }

            return aElement;
        }

//------------------------------------------------------------------------------
    }
}
