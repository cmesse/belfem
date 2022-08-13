//
// Created by Christian Messe on 2019-07-25.
//

#include "cl_Element.hpp"
#include "cl_Facet.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        Element::Element( const id_t aID )
        {
            this->set_id( aID );
            mElements = nullptr ;
            mFacets   = nullptr ;
        }

//------------------------------------------------------------------------------

        Element::~Element()
        {
            if( mElements != nullptr )
            {
                free( mElements );
            }
            if( mFacets != nullptr )
            {
                free( mFacets );
            }
        }

//------------------------------------------------------------------------------

        inline void
        Element::increment_element_counter()
        {
            ++mNumberOfElements;
        }

//------------------------------------------------------------------------------

        void
        Element::allocate_element_container( const uint aSize )
        {
            if( aSize != 0 )
            {
                // allocate container
                mElements = ( Element ** ) malloc( aSize * sizeof( Element ) );

                for( uint e=0; e<aSize; ++e )
                {
                    mElements[ e ] = nullptr ;
                }

                // reset counter
                mNumberOfElements = 0;
            }
        }

//------------------------------------------------------------------------------

        void
        Element::allocate_facet_container( const uint aSize )
        {
            if( aSize != 0 )
            {
                // allocate container
                mFacets = ( Facet ** ) malloc( aSize * sizeof( Facet ) );

                for( uint f=0; f<aSize; ++f )
                {
                    mElements[ f ] = nullptr ;
                }
            }
        }


//------------------------------------------------------------------------------

        void
        Element::insert_element( Element * aElement )
        {
            mElements[ mNumberOfElements++ ] = aElement;
        }

//------------------------------------------------------------------------------

        void
        Element::insert_facet( Facet * aFacet, const uint aIndex )
        {
            mFacets[ aIndex ] = aFacet ;
        }

//------------------------------------------------------------------------------

        void
        Element::reset_element_container()
        {
            mNumberOfElements = 0 ;
            if( mElements != nullptr )
            {
                free( mElements );
                mElements = nullptr ;
            }
        }

//------------------------------------------------------------------------------

        void
        Element::print() const
        {
            BELFEM_ERROR( false,
                          "invalid call of base class function print() from element %lu",
                          ( long unsigned int ) this->id() );
        }

//------------------------------------------------------------------------------
    }
}
