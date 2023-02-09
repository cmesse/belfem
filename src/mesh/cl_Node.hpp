//
// Created by Christian Messe on 2019-07-22.
//

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
            Node **mDuplicates = nullptr ;
            uint mNumberOfDuplicates = 0 ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Node( const id_t & aID,
                  const real aX=0.0,
                  const real aY=0.0,
                  const real aZ=0.0 );

//------------------------------------------------------------------------------

            ~Node();

//------------------------------------------------------------------------------

            EntityType
            entity_type() const ;

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
            set_coords( const real aX, const real aY );

            void
            set_coords( const real aX, const real aY, const real & aZ  );

            void
            allocate_duplicate_container( const uint aSize );

            void
            add_duplicate( Node * aNode );

//------------------------------------------------------------------------------

            /**
             * returns the number of duplicates of this node
             */
            uint
            number_of_duplicates() const ;

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
            return mNumberOfDuplicates ;
        }

//------------------------------------------------------------------------------

        /**
         * returns duplicates of this node
         */
        inline Node *
        Node::duplicate( const uint aIndex )
        {
            BELFEM_ASSERT( aIndex < mNumberOfDuplicates,
                           "Duplicate Index %d out of bounds, which must be less than %d",
                           ( int ) aIndex,
                           ( int ) mNumberOfDuplicates );

            return mDuplicates[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline void
        Node::allocate_duplicate_container( const uint aSize )
        {
            if ( aSize > 0 )
            {
                mDuplicates = new Node * [ aSize ];
            }
            mNumberOfDuplicates = 0;
        }

//------------------------------------------------------------------------------

        inline void
        Node::add_duplicate( Node * aNode )
        {
            mDuplicates[ mNumberOfDuplicates++ ] = aNode ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_NODE_HPP
