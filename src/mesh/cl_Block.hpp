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

#ifndef BELFEM_CL_BLOCK_HPP
#define BELFEM_CL_BLOCK_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Element.hpp"
#include "en_DomainType.hpp"

namespace belfem
{
    namespace mesh
    {
        class GmshReader;

        class Block
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            id_t mID;
            index_t mIndex = gNoIndex ;

            index_t mElementCounter = 0;

            Cell< Element * > mElements;

            DomainType mDomainType = DomainType::Default ;

            string mLabel;

            //! flag telling if elements on this block have edges
            bool mHasEdges = false ;

            //!  flag telling if elements on this block have faces
            bool mHasFaces = false ;

            //! flag telling if block is visualized in VTK tool (not exodus output)
            bool mIsHidden = false ;

            //! thickness, only used for thin shells
            real mThickness = BELFEM_QUIET_NAN  ;

            friend GmshReader;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Block( const id_t aID, const index_t aNumElements );

//------------------------------------------------------------------------------

            ~Block() ;

//------------------------------------------------------------------------------

            inline id_t
            id() const
            {
                return mID;
            }

            inline index_t
            index() const
            {
                return mIndex;
            }

            inline void
            set_index( const index_t aIndex )
            {
                mIndex = aIndex;
            }
//------------------------------------------------------------------------------

            void
            insert_element( Element * aElement );

//------------------------------------------------------------------------------

            inline Element *
            element( const index_t aIndex )
            {
                BELFEM_ASSERT( aIndex < mElements.size(),
                              "Invalid Element index %lu in block %lu",
                              ( long unsigned int ) aIndex,
                              ( long unsigned int ) mID );

                return mElements( aIndex );
            }

//------------------------------------------------------------------------------

            inline Cell< Element * > &
            elements()
            {
                return mElements;
            }

//------------------------------------------------------------------------------

            inline index_t
            number_of_elements() const
            {
                return mElements.size();
            }

//------------------------------------------------------------------------------

            inline string &
            label()
            {
                return mLabel;
            }

//------------------------------------------------------------------------------

            inline ElementType
            element_type() const
            {
                BELFEM_ASSERT( mElements.size() > 0, "Block %u seems to be empty.",
                              ( unsigned int ) mID );

                return mElements( 0 )->type();
            }


//------------------------------------------------------------------------------

            /**
             * function telling if elements in this group have edges
             * ( default : false )
             */
            inline bool
            has_edges() const
            {
                return mHasEdges;
            }

//------------------------------------------------------------------------------

            /**
             * function telling if elements in this group have faces
             * ( default : false )
             */
            inline bool
            has_faces() const
            {
                return mHasFaces;
            }

//------------------------------------------------------------------------------

            inline bool
            is_hidden() const
            {
                return mIsHidden ;
            }

//------------------------------------------------------------------------------

            void
            set_edges_flag( const bool aFlag=true );

            void
            set_faces_flag( const bool aFlag=true );

//------------------------------------------------------------------------------

            void
            flag_elements( const uint8_t aFlagIndex = 0 );

//------------------------------------------------------------------------------

            void
            unflag_elements( const uint8_t aFlagIndex = 0 );

//------------------------------------------------------------------------------

            void
            flag_nodes();

//------------------------------------------------------------------------------

            void
            unflag_nodes();

//------------------------------------------------------------------------------

            void
            flag_corner_nodes();

//------------------------------------------------------------------------------

            void
            flag_edges();

//------------------------------------------------------------------------------

            void
            flag_faces();

//------------------------------------------------------------------------------

            DomainType
            domain_type() const ;

//------------------------------------------------------------------------------

            void
            set_domain_type( const DomainType aType );

            void
            set_thickness( const real aThickness );

            real
            thickness() const ;

            void
            set_element_counter( const index_t aCounter );

//------------------------------------------------------------------------------

            size_t
            memory() const ;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline DomainType
        Block::domain_type() const
        {
            return mDomainType ;
        }

//------------------------------------------------------------------------------

        inline
        void Block::set_domain_type( const DomainType aType )
        {
            mDomainType = aType ;
        }

//------------------------------------------------------------------------------

        inline void
        Block::set_thickness( const real aThickness )
        {
            mThickness = aThickness ;
        }

//------------------------------------------------------------------------------

        inline real
        Block::thickness() const
        {
            return mThickness ;
        }

//------------------------------------------------------------------------------

        inline size_t Block::memory() const
        {
            return sizeof( Block )
                + mElements.size() * sizeof( Element * )
                + mLabel.capacity() * sizeof( char ) ;
        }

//------------------------------------------------------------------------------

        inline void Block::set_element_counter( const index_t aCounter )
        {
            mElementCounter = aCounter ;
        }

//------------------------------------------------------------------------------

    }
}
#endif //BELFEM_CL_BLOCK_HPP
