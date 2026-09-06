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

#ifndef BELFEM_CL_CONTROLPOINT_HPP
#define BELFEM_CL_CONTROLPOINT_HPP

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Mesh_Basis.hpp"
#include "Mesh_Enums.hpp"

namespace belfem
{
    namespace mesh
    {
        class Element ;

        class ControlPoint : public Basis
        {
            // coordinates of this control point
            real mCoords[ 3 ];

            //! pointer to elements
            Element **mElements;

            //! number of elements connected to this vertex
            uint mElementCounter = 0 ;

            //! connected control points
            ControlPoint **mControlPoints;

            //! number of control points connected to this vertex
            uint mControlPointCounter = 0 ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            ControlPoint( const id_t aID, const real aX=0.0, const real aY=0.0, const real aZ=0.0 ) ;

//------------------------------------------------------------------------------

            ~ControlPoint() override ;

//------------------------------------------------------------------------------

            EntityType
            entity_type() const override;

//------------------------------------------------------------------------------

            real
            x() const;

//------------------------------------------------------------------------------

            real
            y() const;

//------------------------------------------------------------------------------

            real
            z() const;

//------------------------------------------------------------------------------

            real
            x( const uint aIndex ) const;

//------------------------------------------------------------------------------

            Vector<real>
            coords() const;

//------------------------------------------------------------------------------

            void
            set_coords( const real aX, const real aY, const real aZ );

//------------------------------------------------------------------------------

            void
            allocate_element_container();

//------------------------------------------------------------------------------

            void
            allocate_control_point_container();

//------------------------------------------------------------------------------

            void
            reset_element_container() override;

//------------------------------------------------------------------------------

            void
            reset_control_point_container();

//------------------------------------------------------------------------------

            void
            increment_element_counter();

//------------------------------------------------------------------------------

            void
            increment_control_point_counter();
//------------------------------------------------------------------------------

            void
            add_element( Element * aElement );

//------------------------------------------------------------------------------

            void
            add_control_point( ControlPoint * aControlPoint );

//------------------------------------------------------------------------------

            uint
            number_of_elements() const override;

//------------------------------------------------------------------------------

            uint
            number_of_control_points() const;

//------------------------------------------------------------------------------

            Element *
            element( const uint aIndex ) override;

//------------------------------------------------------------------------------

            const Element *
            element( const uint aIndex ) const;

//------------------------------------------------------------------------------

            ControlPoint *
            control_point( const uint aIndex );

//------------------------------------------------------------------------------

            const ControlPoint *
            control_point( const uint aIndex ) const;

//------------------------------------------------------------------------------

            size_t memory() const;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline EntityType
        ControlPoint::entity_type() const
        {
            return EntityType::CONTROLPOINT;
        }

//------------------------------------------------------------------------------

        inline real
        ControlPoint::x() const
        {
            return mCoords[ 0 ];
        }

//------------------------------------------------------------------------------

        inline real
        ControlPoint::y() const
        {
            return mCoords[ 1 ];
        }

//------------------------------------------------------------------------------

        inline real
        ControlPoint::z() const
        {
            return mCoords[ 2 ];
        }

//------------------------------------------------------------------------------

        inline real
        ControlPoint::x( const uint aIndex ) const
        {
            BELFEM_ASSERT( aIndex < 3, "Invalid index %u for control point %lu.", ( unsigned int ) aIndex, ( long unsigned int ) id() );
            return mCoords[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline Vector<real>
        ControlPoint::coords() const
        {
            return Vector<real>( { mCoords[0], mCoords[1], mCoords[2] } );
        }

//------------------------------------------------------------------------------

        inline uint
        ControlPoint::number_of_elements() const
        {
            return mElementCounter;
        }

//------------------------------------------------------------------------------

        inline uint
        ControlPoint::number_of_control_points() const
        {
            return mControlPointCounter;
        }

//------------------------------------------------------------------------------

        inline Element *
        ControlPoint::element( const uint aIndex )
        {
            BELFEM_ASSERT(  aIndex < mElementCounter,
                "Invalid control point index %u for element %lu. (expect < %u)",
                ( unsigned int ) aIndex, ( long unsigned int ) id(),
                ( unsigned int ) mElementCounter );
            return mElements[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline const Element *
        ControlPoint::element( const uint aIndex ) const
        {
            BELFEM_ASSERT(  aIndex < mElementCounter,
                "Invalid control point index %u for element %lu. (expect < %u)",
                ( unsigned int ) aIndex, ( long unsigned int ) id(), ( unsigned int )
                mElementCounter );
            return mElements[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline ControlPoint *
        ControlPoint::control_point( const uint aIndex )
        {
            BELFEM_ASSERT( aIndex < mControlPointCounter,
                "Invalid control point index %u for control point %lu. (expect < %u)",
                ( unsigned int ) aIndex, ( long unsigned int ) id(),
                ( unsigned int ) mControlPointCounter );
            return mControlPoints[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline const ControlPoint *
        ControlPoint::control_point( const uint aIndex ) const
        {
            BELFEM_ASSERT( aIndex < mControlPointCounter,
            "Invalid control point index %u for control point %lu. (expect < %u)",
            ( unsigned int ) aIndex, ( long unsigned int ) id(),
            ( unsigned int ) mControlPointCounter );
            return mControlPoints[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline void
        ControlPoint::increment_element_counter()
        {
            ++mElementCounter;
        }

//------------------------------------------------------------------------------

        inline void ControlPoint::increment_control_point_counter()
        {
            ++mControlPointCounter;
        }

//------------------------------------------------------------------------------

        inline void ControlPoint::add_element( Element *aElement )
        {
            mElements[ mElementCounter++ ] = aElement;
        }

//------------------------------------------------------------------------------

        inline void ControlPoint::add_control_point( ControlPoint *aControlPoint )
        {
            mControlPoints[ mControlPointCounter++ ] = aControlPoint;
        }

//------------------------------------------------------------------------------

        inline size_t
        ControlPoint::memory() const
        {
            return sizeof( ControlPoint )
                + ( mElementCounter ) * sizeof( Element * )
                + ( mControlPointCounter ) * sizeof( ControlPoint * );
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_CONTROLPOINT_HPP