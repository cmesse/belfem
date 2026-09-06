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

#ifndef BELFEM_CL_SIDESET_HPP
#define BELFEM_CL_SIDESET_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Element.hpp"
#include "cl_Facet.hpp"
#include "cl_Map.hpp"
#include "en_DomainType.hpp"

namespace belfem
{
    namespace mesh
    {
        class GmshReader;

        class SideSet
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            id_t mID;
            index_t mIndex = gNoIndex ;

            index_t mFacetCounter = 0;

            DomainType mDomainType = DomainType::Default ;

            string mLabel;

            Cell< Facet * > mFacets;
            Cell< Node *  > mNodes ;

            friend GmshReader;

            // this flag hides the sideset from exodus
            bool mIsHidden = false ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            SideSet( const id_t aID, const index_t aNumElements );

//------------------------------------------------------------------------------

            ~SideSet();

//------------------------------------------------------------------------------

            id_t
            id() const;

//------------------------------------------------------------------------------

            index_t
            index() const ;

            void
            set_index( const index_t aIndex ) ;
//------------------------------------------------------------------------------

            string &
            label();

//------------------------------------------------------------------------------

            void
            insert_facet( Facet * aFacet );

//------------------------------------------------------------------------------

            Facet *
            facet_by_index( const index_t aIndex );

//------------------------------------------------------------------------------

            //Facet *
            //facet( const id_t aID );

//------------------------------------------------------------------------------

            Cell < Facet * > &
            facets();

//------------------------------------------------------------------------------

            index_t
            number_of_facets() const;

//------------------------------------------------------------------------------

            void
            unflag_all_nodes();

//------------------------------------------------------------------------------

            void
            flag_all_nodes();

//------------------------------------------------------------------------------

            void
            flag_corner_nodes();

//------------------------------------------------------------------------------

            void
            unflag_all_facets();

//------------------------------------------------------------------------------

            void
            flag_all_facets();

//------------------------------------------------------------------------------

            void
            flag_edges();

//------------------------------------------------------------------------------

            void
            unflag_edges();

//------------------------------------------------------------------------------

            void
            flag_faces();

//------------------------------------------------------------------------------

            void
            unflag_faces();

//------------------------------------------------------------------------------

            void
            collect_nodes() ;

//------------------------------------------------------------------------------

            /**
             * expose list with all nodes on sideset
             * @return
             */
            Cell< Node * > &
            nodes() ;

//------------------------------------------------------------------------------

            void
            reset_node_container();

//------------------------------------------------------------------------------

            void
            reset_facet_container();

//------------------------------------------------------------------------------

            ElementType
            element_type() const ;

//------------------------------------------------------------------------------

            void
            hide( const bool aSwitch = true );

//------------------------------------------------------------------------------

            bool
            is_hidden() const ;

//------------------------------------------------------------------------------

            DomainType
            domain_type() const ;

            void
            set_domain_type( const DomainType aType );

            void
            set_facet_counter( const index_t aCounter );

            size_t
            memory() const ;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline id_t
        SideSet::id() const
        {
            return mID;
        }

//------------------------------------------------------------------------------

        inline index_t
        SideSet::index() const
        {
            return mIndex;
        }

//------------------------------------------------------------------------------

        inline void
        SideSet::set_index( const index_t aIndex )
        {
            mIndex = aIndex ;
        }

//------------------------------------------------------------------------------

        inline Facet *
        SideSet::facet_by_index( const index_t aIndex )
        {
            BELFEM_ASSERT( aIndex < mFacets.size(),
                          "Invalid Facet index" );

            return mFacets( aIndex );
        }

//------------------------------------------------------------------------------

        inline index_t
        SideSet::number_of_facets() const
        {
            return mFacets.size();
        }

//------------------------------------------------------------------------------

        inline string &
        SideSet::label()
        {
            return mLabel;
        }

//------------------------------------------------------------------------------

        inline Cell < Facet * > &
        SideSet::facets()
        {
            return mFacets;
        }

//------------------------------------------------------------------------------

        inline Cell < Node * > &
        SideSet::nodes()
        {
            return mNodes;
        }

//------------------------------------------------------------------------------

        inline ElementType
        SideSet::element_type() const
        {
            if( mFacets.size() == 0 )
            {
                return ElementType::EMPTY ;
            }
            else
            {
                return mFacets( 0 )->element()->type() ;
            }
        }

//------------------------------------------------------------------------------

        inline void
        SideSet::hide( const bool aSwitch )
        {
            mIsHidden = aSwitch ;
        }

//------------------------------------------------------------------------------

        inline bool
        SideSet::is_hidden() const
        {
            return mIsHidden ;
        }

//------------------------------------------------------------------------------

        inline DomainType
        SideSet::domain_type() const
        {
            return mDomainType ;
        }

//------------------------------------------------------------------------------

        inline
        void SideSet::set_domain_type( const DomainType aType )
        {
            mDomainType = aType ;
        }

//------------------------------------------------------------------------------

        inline void SideSet::set_facet_counter( const index_t aCounter )
        {
            mFacetCounter = aCounter ;
        }

//------------------------------------------------------------------------------

        inline size_t SideSet::memory() const
        {
            return sizeof( SideSet )
                + mFacets.size() * sizeof( Facet * )
                + mNodes.size() * sizeof( Node * )
                + mLabel.capacity() * sizeof( char ) ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_SIDESET_HPP
