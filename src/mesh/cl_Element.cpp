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
            if( mNeighbors != nullptr )
            {
                free( mNeighbors );
            }
            if ( mControlPoints != nullptr )
            {
                free( mControlPoints );
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
            if( mElements != nullptr )
            {
                free( mElements );
                mElements = nullptr ;
            }

            // always reset: a stale counter over a freed container makes
            // number_of_elements() lie about a null array
            mNumberOfElements = 0;

            // the count field is narrow; refuse before anything is allocated
            BELFEM_ERROR( aSize <= std::numeric_limits< decltype( mNumberOfElements ) >::max(),
                "element count of element %lu overflows its counter width ( %u requested )",
                ( long unsigned int ) this->id(), ( unsigned int ) aSize );

            if( aSize != 0 )
            {
                // allocate container
                mElements = ( Element ** ) malloc( aSize * sizeof( Element * ) );

                for( uint e=0; e<aSize; ++e )
                {
                    mElements[ e ] = nullptr ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Element::allocate_facet_container( const uint aSize )
        {
            if( mFacets != nullptr )
            {
                free( mFacets );
                mFacets = nullptr ;
            }

            if( aSize != 0 )
            {
                // memory() charges one slot per facet
                BELFEM_ERROR( aSize == this->number_of_facets(),
                    "facet container of element %lu must hold one slot per facet",
                    ( long unsigned int ) this->id() );

                // allocate container
                mFacets = ( Facet ** ) malloc( aSize * sizeof( Facet * ) );

                for( uint f=0; f<aSize; ++f )
                {
                    mFacets[ f ] = nullptr ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Element::allocate_neighbor_container()
        {
            if( this->number_of_facets() == 0 )
            {
                return;

            }
            if( mNeighbors != nullptr )
            {
                free( mNeighbors );
            }
            mNeighbors = ( Element ** ) malloc( this->number_of_facets() * sizeof( Element * ) );
            for( uint f=0; f<this->number_of_facets(); ++f )
            {
                mNeighbors[ f ] = nullptr ;
            }
        }

//------------------------------------------------------------------------------

        void
        Element::insert_neighbor( Element * aNeighbor, const uint aIndex )
        {
            mNeighbors[ aIndex ] = aNeighbor ;
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
