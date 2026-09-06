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

#ifndef BELFEM_CL_PIPETTE_HPP
#define BELFEM_CL_PIPETTE_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Element.hpp"
#include "cl_Facet.hpp"
#include "cl_IF_IntegrationData.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        /**
         * the pipette is a class that measures the volume of an element
         */
        class Pipette
        {

            ElementType mType = ElementType::UNDEFINED ;
            uint mNumDim = 0 ;
            uint mNumNodes = 0 ;
            uint mNumCornerNodes = 0 ;

            uint mNumIntPoints = 0 ;
            uint mNumIntPointsLinear = 0 ;


            // help fector
            real * mW = nullptr ;

            // raw vector for coordinates
            real * mX = nullptr ;
            real * mY = nullptr ;
            real * mZ = nullptr ;

            Matrix< real > mNodeCoords ;
            Matrix< real > mNodeCoordsLinear ;

            // for thin shells
            Vector< real > mN ; // normal
            Matrix< real > mJ ; // jacobian

            // the default integration function
            fem::IntegrationData * mIntegrationData = nullptr ;

            // special integration function for non-curved elements
            // used for all elements except tri and tet
            fem::IntegrationData * mIntegrationDataLinear = nullptr ;

            real
            ( Pipette::*mVolumeFunction )( const Element * aElement );

            // function that we can use if the element is not curved
            real
            ( Pipette::*mVolumeFunctionLinear )( const Element * aElement );

            real
            ( Pipette::*mSurfaceFunction )( const Facet * aFacet );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Pipette();

//------------------------------------------------------------------------------

            ~Pipette();

//------------------------------------------------------------------------------

            void
            set_element_type( const ElementType aType );

//------------------------------------------------------------------------------

            void
            set_facet_type( const ElementType aType );

//------------------------------------------------------------------------------

            real
            measure( const Element * aElement );

//------------------------------------------------------------------------------

            real
            measure( const Facet * aFacet );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            reset_containers();

//------------------------------------------------------------------------------

            void
            collect_node_coords( const Element * aElement );

//------------------------------------------------------------------------------

            real
            measure_tri3( const Element * aElement );

//------------------------------------------------------------------------------

            real
            measure_quad4( const Element * aElement );

//------------------------------------------------------------------------------

            real
            measure_tet4( const Element * aElement );

//------------------------------------------------------------------------------

            real
            measure_linear( const Element * aElement );

//------------------------------------------------------------------------------

            real
            measure_quad4ts( const Element * aElement );

//------------------------------------------------------------------------------

            real
            measure_higher_order( const Element * aElement );

//------------------------------------------------------------------------------

            real
            measure_surface_tri3( const Facet * aFacet );

            real
            measure_surface_line2( const Facet * aFacet );

            real
            measure_surface_higher_order( const Facet * aFacet );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        inline real
        Pipette::measure( const Element * aElement )
        {
            return ( this->*mVolumeFunction )( aElement );
        }

//------------------------------------------------------------------------------

        inline real
        Pipette::measure( const Facet * aFacet )
        {
            return ( this->*mSurfaceFunction )( aFacet );
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_PIPETTE_HPP
