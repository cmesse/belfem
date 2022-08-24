//
// Created by christian on 8/23/22.
//

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
            uint mNumIntPoints = 0 ;
            bool mHaveW = false ;

            real * mW = nullptr ;
            real * mX = nullptr ;
            real * mY = nullptr ;
            real * mZ = nullptr ;

            Matrix< real > mNodeCoords ;
            Matrix< real > mNodeCoordsFacet ;

            belfem::fem::IntegrationData * mIntegrationData = nullptr ;

            real
            ( Pipette::*mVolumeFunction )( const Element * aElement );

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

            void
            collect_node_coords( const Facet * aFacet );

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
            measure_penta6( const Element * aElement );

//------------------------------------------------------------------------------

            real
            measure_hex8( const Element * aElement );

//------------------------------------------------------------------------------

            real
            measure_higher_order( const Element * aElement );

//------------------------------------------------------------------------------

            real
            measure_higher_order_quad( const Facet * aFacet );

//------------------------------------------------------------------------------

            real
            measure_linear_tri( const Facet * aFacet );

//------------------------------------------------------------------------------

            real
            measure_higher_order_tri( const Facet * aFacet );
//------------------------------------------------------------------------------


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
