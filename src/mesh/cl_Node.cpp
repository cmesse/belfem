//
// Created by Christian Messe on 2019-07-25.
//

#include "cl_Node.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        Node::Node(
                const id_t & aID,
                const real aX,
                const real aY,
                const real aZ ) :
                Vertex()
        {
            // set the id
            this->set_id( aID );

            // copy node coordinates into coordinate vector
            mCoords[ 0 ] = aX;
            mCoords[ 1 ] = aY;
            mCoords[ 2 ] = aZ;
        }

//------------------------------------------------------------------------------

        Node::~Node()
        {
            if( mDuplicates != nullptr )
            {
                delete[] mDuplicates;
                mNumberOfDuplicates = 0;
            }
            this->delete_containers();
        }

//------------------------------------------------------------------------------

        void
        Node::set_coords( const Vector< real > & aCoords )
        {
            std::copy( aCoords.data(),
                       aCoords.data()
                        + aCoords.length(), mCoords );
        }

//------------------------------------------------------------------------------

        void
        Node::set_coords( const real aX, const real aY )
        {
            mCoords[ 0 ] = aX;
            mCoords[ 1 ] = aY;
        }

//------------------------------------------------------------------------------

        void
        Node::set_coords( const real aX, const real aY, const real & aZ )
        {
            mCoords[ 0 ] = aX;
            mCoords[ 1 ] = aY;
            mCoords[ 2 ] = aZ;
        }


//------------------------------------------------------------------------------

        Node *
        Node::original()
        {
            if( mNumberOfDuplicates == 0 )
            {
                return this ;
            }
            else
            {
                // this might look slow, but we don't call this
                // very often, and we save one pointer!
                id_t tID = this->id() ;
                uint tIndex = 0 ;
                for( uint d=0; d<mNumberOfDuplicates; ++d )
                {
                    if( mDuplicates[ d ]->id() < tID )
                    {
                        tID = mDuplicates[ d ]->id() ;
                        tIndex = d ;
                    }
                }
                return mDuplicates[ tIndex ];
            }
        }

//------------------------------------------------------------------------------
    }
}