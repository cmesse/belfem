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

#include "cl_ControlPoint.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        ControlPoint::ControlPoint(
                const id_t aID,
                const real aX,
                const real aY,
                const real aZ ) :
                Basis()
        {
            // set the id
            this->set_id( aID );

            // copy coordinates into coordinate vector
            mCoords[ 0 ] = aX;
            mCoords[ 1 ] = aY;
            mCoords[ 2 ] = aZ;

            // initialize pointers
            mElements = nullptr;
            mControlPoints = nullptr;
        }

//------------------------------------------------------------------------------

        ControlPoint::~ControlPoint()
        {
            this->reset_element_container();
            this->reset_control_point_container();
        }

//------------------------------------------------------------------------------

        void
        ControlPoint::set_coords( const real aX, const real aY, const real aZ )
        {
            mCoords[ 0 ] = aX;
            mCoords[ 1 ] = aY;
            mCoords[ 2 ] = aZ;
        }

//------------------------------------------------------------------------------

        void
        ControlPoint::allocate_element_container( )
        {
            BELFEM_ASSERT( mElements == nullptr, "Elements are already allocated." );
            if( mElementCounter > 0 )
            {
                mElements = ( Element ** ) malloc( mElementCounter * sizeof( Element * ) );
                mElementCounter = 0 ;
            }
        }

//------------------------------------------------------------------------------

        void
        ControlPoint::allocate_control_point_container()
        {
            BELFEM_ASSERT( mControlPoints == nullptr, "Control points are already allocated." );

            if( mControlPointCounter > 0 )
            {
                mControlPoints = ( ControlPoint ** ) malloc( mControlPointCounter * sizeof( ControlPoint * ) );
                mControlPointCounter = 0 ;
            }
        }

//------------------------------------------------------------------------------

        void
        ControlPoint::reset_element_container()
        {
            if ( mElements != nullptr )
            {
                free( mElements );
                mElements = nullptr;
                mElementCounter = 0;
            }
        }

//------------------------------------------------------------------------------

        void
        ControlPoint::reset_control_point_container()
        {
            if ( mControlPoints != nullptr )
            {
                free( mControlPoints );
                mControlPoints = nullptr;
                mControlPointCounter = 0;
            }
        }

//------------------------------------------------------------------------------
    }
}
