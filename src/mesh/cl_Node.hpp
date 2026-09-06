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

#ifndef BELFEM_CL_NODE_HPP
#define BELFEM_CL_NODE_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "meshtools.hpp"
#include "Mesh_Enums.hpp"

#include "cl_Vector.hpp"
#include "cl_Vertex.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    namespace mesh
    {
        class Node : public Vertex
        {
            // coordinates of this vector
            real mCoords[ 3 ];

            // container for duplicates
            Node ** mDuplicates = nullptr ;
            int mNumberOfDuplicates = 0 ;

            Node * mPeriodic = nullptr ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Node( const id_t & aID,
                  const real aX=0.0,
                  const real aY=0.0,
                  const real aZ=0.0 );

//------------------------------------------------------------------------------

            ~Node() override;

//------------------------------------------------------------------------------

            EntityType
            entity_type() const override ;

//------------------------------------------------------------------------------

            real
            x() const;

//------------------------------------------------------------------------------

            real
            x( const uint aDimension ) const;

//------------------------------------------------------------------------------

            real
            y() const;

//------------------------------------------------------------------------------

            real
            z() const;

//------------------------------------------------------------------------------

            Vector<real>
            coords() const ;

//------------------------------------------------------------------------------

            /**
             * change the coordinates of this node
             */
            void
            set_coords( const Vector< real > & aCoords );

            void
            get_coords( Vector< real > & aCoords );

            void
            set_coords( const real aX, const real aY );

            void
            set_coords( const real aX, const real aY, const real & aZ  );

            void
            allocate_duplicate_container( const uint aSize );

            void
            reset_duplicate_container();

            void
            add_duplicate( Node * aNode );

//------------------------------------------------------------------------------

            /**
             * returns the number of duplicates of this node
             */
            uint
            number_of_duplicates() const ;

//------------------------------------------------------------------------------

            bool
            is_duplicate() const ;

//------------------------------------------------------------------------------

            /**
             * returns duplicates of this node
             */
             Node *
             duplicate( const uint aIndex );

//------------------------------------------------------------------------------

            /**
             * returns the original of this node
             */
            Node *
            original();

            const Node *
            original() const ;

//------------------------------------------------------------------------------

            void
            set_original( Node * aNode );

//------------------------------------------------------------------------------

            void
            unlink_from_original();

//------------------------------------------------------------------------------

            size_t
            memory() const ;

//------------------------------------------------------------------------------

            void
            set_periodic( Node * aNode ) ;

            Node *
            periodic() ;

            const Node *
            periodic() const ;

            bool
            is_periodic() const ;

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        inline EntityType
        Node::entity_type() const
        {
            return EntityType::NODE ;
        }

//------------------------------------------------------------------------------

        inline real
        Node::x() const
        {
            return mCoords[ 0 ];
        }

//------------------------------------------------------------------------------

        inline real
        Node::y() const
        {
            return mCoords[ 1 ];
        }

//------------------------------------------------------------------------------

        inline real
        Node::z() const
        {
            return mCoords[ 2 ];
        }

//------------------------------------------------------------------------------

        inline Vector<real>
        Node::coords() const
        {
            return Vector< real >( {mCoords[0],mCoords[1], mCoords[2]} );
        }

//------------------------------------------------------------------------------

        inline real
        Node::x( const uint aDimension ) const
        {
            return mCoords[ aDimension ];
        }

//------------------------------------------------------------------------------

        /**
         * returns the number of duplicates of this node
         */
        inline uint
        Node::number_of_duplicates() const
        {
            return mNumberOfDuplicates == -1 ? 0 : mNumberOfDuplicates ;
        }

//------------------------------------------------------------------------------

        inline bool
        Node::is_duplicate() const
        {
            return mNumberOfDuplicates == -1 ;
        }

//------------------------------------------------------------------------------

        /**
         * returns duplicates of this node
         */
        inline Node *
        Node::duplicate( const uint aIndex )
        {
            BELFEM_ASSERT( aIndex < std::abs( mNumberOfDuplicates ) ,
                           "Duplicate Index %d for node %lu out of bounds, which must be less than %d",
                           ( int ) aIndex,
                           ( long unsigned int ) this->id(),
                           ( int ) mNumberOfDuplicates );

            return mDuplicates[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline void
        Node::allocate_duplicate_container( const uint aSize )
        {
            BELFEM_ASSERT( ! this->is_duplicate(), "Duplicates of a duplicate node are forbidden (node %lu, original: %lu )",
                ( long unsigned int ) this->id(),( long unsigned int ) this->original()->id());


            BELFEM_ASSERT( mNumberOfDuplicates == 0, "duplicate container of node %lu already allocated", (luint) this->id() );

            if ( aSize > 0 )
            {
                mDuplicates = new Node * [ aSize ];
            }
            mNumberOfDuplicates = 0;
        }

//------------------------------------------------------------------------------

        inline void
        Node::reset_duplicate_container()
        {
            if( mDuplicates != nullptr )
            {
                delete[] mDuplicates;
                mDuplicates = nullptr;
                mNumberOfDuplicates = 0;
            }
        }
//------------------------------------------------------------------------------

        inline void
        Node::add_duplicate( Node * aNode )
        {
            mDuplicates[ mNumberOfDuplicates++ ] = aNode ;
        }

//------------------------------------------------------------------------------

        inline void
        Node::set_original( Node * aNode )
        {
            BELFEM_ASSERT( mNumberOfDuplicates == 0, "can't set original if this node has duplicates" );

            mNumberOfDuplicates = -1 ;
            mDuplicates = new Node * [ 1 ];
            mDuplicates[ 0 ] = aNode ;
        }

        inline void
        Node::unlink_from_original()
        {
            BELFEM_ASSERT( this->is_duplicate(), "can't unlink from original if this node is not a duplicate" );
            mNumberOfDuplicates = 0 ;
            delete[] mDuplicates;
            mDuplicates = nullptr ;
        }

//------------------------------------------------------------------------------

        inline Node *
        Node::original()
        {
            return mNumberOfDuplicates == -1 ? mDuplicates[ 0 ] : this ;
        }
//------------------------------------------------------------------------------

        inline const Node *
        Node::original() const
        {
            return mNumberOfDuplicates == -1 ? mDuplicates[ 0 ] : this ;
        }

//------------------------------------------------------------------------------

        inline size_t
        Node::memory() const
        {
            // mCoords[3] is a member array, already counted in sizeof(Node).
            // the duplicate array holds the duplicates of an original, or the
            // one slot a duplicate keeps for its original ( counter -1 );
            // its capacity is not stored, this is the fill count
            return sizeof( Node ) + this->array_memory()
                 + std::abs( mNumberOfDuplicates ) * sizeof( Node * );
        }

//------------------------------------------------------------------------------

        inline void
        Node::set_periodic( Node * aNode )
        {
            mPeriodic = aNode ;
        }

        inline Node *
        Node::periodic()
        {
            return mPeriodic ;
        }

        inline const Node *
        Node::periodic() const
        {
            return mPeriodic ;
        }

        inline bool
        Node::is_periodic() const
        {
            return mPeriodic != nullptr ;
        }

        inline void Node::get_coords( Vector< real > & aCoords )
        {
            aCoords.set_size( 3 );
            aCoords( 0 ) = mCoords[ 0 ];
            aCoords( 1 ) = mCoords[ 1 ];
            aCoords( 2 ) = mCoords[ 2 ];
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_NODE_HPP
