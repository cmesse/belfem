//
// Created by Christian Messe on 2019-07-28.
//

#ifndef BELFEM_CL_SIDESET_HPP
#define BELFEM_CL_SIDESET_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Element.hpp"
#include "cl_Facet.hpp"
#include "cl_Map.hpp"
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

            index_t mFacetCounter = 0;

            string mLabel;

            Cell< Facet * > mFacets;
            Cell< Node *  > mNodes ;

            Map< id_t, Facet * > mFacetMap ;

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

            string &
            label();

//------------------------------------------------------------------------------

            void
            insert_facet( Facet * aFacet );

//------------------------------------------------------------------------------

            Facet *
            facet_by_index( const index_t aIndex );

//------------------------------------------------------------------------------

            Facet *
            facet( const id_t aID );

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
            flag_faces();

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

            inline void
            hide( const bool aSwitch = true );

//------------------------------------------------------------------------------

            inline bool
            is_hidden() const ;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline id_t
        SideSet::id() const
        {
            return mID;
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

        inline Facet *
        SideSet::facet( const id_t aID )
        {
            return mFacetMap( aID );
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

        void
        SideSet::hide( const bool aSwitch )
        {
            mIsHidden = true ;
        }

//------------------------------------------------------------------------------

        bool
        SideSet::is_hidden() const
        {
            return mIsHidden ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_SIDESET_HPP
