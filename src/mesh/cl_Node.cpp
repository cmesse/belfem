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

            // set size of coordinate vector
            mCoords.set_size( 3 );

            // copy node coordinates into coordinate vector
            mCoords( 0 ) = aX;
            mCoords( 1 ) = aY;
            mCoords( 2 ) = aZ;
        }

//------------------------------------------------------------------------------

        Node::~Node()
        {
            this->delete_containers();
        }

//------------------------------------------------------------------------------

        void
        Node::set_coords( const Vector< real > & aCoords )
        {
            mCoords = aCoords;
        }

//------------------------------------------------------------------------------

        void
        Node::set_coords( const real aX, const real aY )
        {
            mCoords( 0 ) = aX;
            mCoords( 1 ) = aY;
        }

//------------------------------------------------------------------------------

        void
        Node::set_coords( const real aX, const real aY, const real & aZ )
        {
            mCoords( 0 ) = aX;
            mCoords( 1 ) = aY;
            mCoords( 2 ) = aZ;
        }

//------------------------------------------------------------------------------
    }
}