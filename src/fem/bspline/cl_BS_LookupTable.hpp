//
// Created by Christian Messe on 30.04.20.
//

#ifndef BELFEM_CL_BS_LOOKUPTABLE_HPP
#define BELFEM_CL_BS_LOOKUPTABLE_HPP

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Mesh.hpp"
#include "cl_Element.hpp"
#include "cl_IF_InterpolationFunction.hpp"

namespace belfem
{
    namespace bspline
    {
//------------------------------------------------------------------------------

        class LookupTable
        {
            // mesh containing interpolation data
            Mesh * mMesh = nullptr;

            const int mNumberOfDimensions;

            // element type of first block
            const ElementType mElementType;

            const uint mNumberOfNodesPerElement ;

            // elements of first block
            Cell< mesh::Element * > & mElements;

            // interpolation function for table
            fem::InterpolationFunction * mInterpolationFunction;

            // map to field lables with indices
            Map< string, index_t > mFieldMap ;

            // first point of bounding box
            Vector< real > mPmin;

            // last point of bounding box
            Vector< real > mPmax;

            // how many elements exist
            Vector< index_t > mNumberOfElementsPerDirection;

            // size of one element
            Vector< real > mElementLength;

            // factors for derivatives
            real mScaleDX    = BELFEM_QUIET_NAN ;
            real mScaleDY    = BELFEM_QUIET_NAN ;
            real mScaleDZ    = BELFEM_QUIET_NAN ;
            real mScaleD2X2  = BELFEM_QUIET_NAN ;
            real mScaleD2Y2  = BELFEM_QUIET_NAN ;
            real mScaleD2Z2  = BELFEM_QUIET_NAN ;

            real mScaleDYDZ  = BELFEM_QUIET_NAN ;
            real mScaleDXDZ  = BELFEM_QUIET_NAN ;
            real mScaleDXDY  = BELFEM_QUIET_NAN ;

            // current values for x, y and z
            real mX = BELFEM_REAL_MAX ;
            real mY = BELFEM_REAL_MAX ;
            real mZ = BELFEM_REAL_MAX ;

            // current element
            mesh::Element * mElement = nullptr ;

            // parameter coordinates on element
            Vector< real > mXi ;

            // vector containing values for interpolation function
            Matrix< real > mN ;
            Matrix< real > mdNdxi ;
            Matrix< real > md2Ndxi2 ;

            // vector containing values for one element
            Vector< real > mValues ;

            // remmember index of stored field
            uint mFieldIndex = BELFEM_INT_MAX ;

            // remember id of selected element when field was colleted
            index_t mElementID = 0 ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            LookupTable( const string aPath );

//------------------------------------------------------------------------------

            ~LookupTable();

//------------------------------------------------------------------------------

            // get a field index from a lanbel
            inline index_t
            field_index( const string & aLabel );

//------------------------------------------------------------------------------
            real
            compute_value(
                    const index_t & aFieldIndex,
                    const real aX );

//------------------------------------------------------------------------------

            real
            compute_value(
                    const index_t & aFieldIndex,
                    const real aX,
                    const real aY );

//------------------------------------------------------------------------------

            real
            compute_value(
                    const index_t & aFieldIndex,
                    const real aX,
                    const real aY,
                    const real & aZ );

//------------------------------------------------------------------------------

            real
            compute_derivative(
                    const index_t aFieldIndex,
                    const real aX );

//------------------------------------------------------------------------------

            Vector< real >
            compute_derivative(
                    const index_t aFieldIndex,
                    const real aX,
                    const real aY );

//------------------------------------------------------------------------------

            Vector< real >
            compute_derivative(
                    const index_t aFieldIndex,
                    const real aX,
                    const real aY,
                    const real aZ );

//------------------------------------------------------------------------------

            real
            compute_second_derivative(
                    const index_t aFieldIndex,
                    const real aX );

//------------------------------------------------------------------------------

            Vector< real >
            compute_second_derivative(
                    const index_t aFieldIndex,
                    const real aX,
                    const real aY );

//------------------------------------------------------------------------------

            Vector< real >
            compute_second_derivative(
                    const index_t aFieldIndex,
                    const real aX,
                    const real aY,
                    const real aZ );

//------------------------------------------------------------------------------

            inline const real  &
            min( const uint aIndex ) const
            {
                return mPmin( aIndex );
            }

//------------------------------------------------------------------------------

            inline const real  &
            max( const uint aIndex ) const
            {
                return mPmax( aIndex );
            }

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            create_field_map();

//------------------------------------------------------------------------------

            void
            create_interpolation_function();

//------------------------------------------------------------------------------

            void
            polpulate_globals();

//------------------------------------------------------------------------------

            void
            compute_derivative_factors();

//------------------------------------------------------------------------------

            // element grabber for 1D problem
            inline mesh::Element *
            get_element( const index_t & aI );

//------------------------------------------------------------------------------

            // element grabber for 2D problem
            inline mesh::Element *
            get_element( const index_t & aI, const index_t & aJ );

//------------------------------------------------------------------------------

            // element grabber for 3D problem
            inline mesh::Element *
            get_element( const index_t & aI, const index_t & aJ, const index_t & aK );

//------------------------------------------------------------------------------

            // element finder for 1D problem
            inline void
            select_element( const real aX );

//------------------------------------------------------------------------------

            // element finder for 2D problem
            inline void
            select_element( const real aX, const real aY );

//------------------------------------------------------------------------------

            // element finder for 3D problem
            inline void
            select_element( const real aX, const real aY, const real & aZ );

//------------------------------------------------------------------------------

            void
            update_values( const index_t aFieldIndex );

//------------------------------------------------------------------------------

            inline void
            compute_xi( const uint & aI, const real aX );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
// Inline functions
//------------------------------------------------------------------------------

        // get a field index from a lanbel
        inline index_t
        LookupTable::field_index( const string & aLabel )
        {
            return mFieldMap( aLabel );
        }

//------------------------------------------------------------------------------
        // element grabber for 1D problem
        mesh::Element *
        LookupTable::get_element( const index_t & aI )
        {
            BELFEM_ASSERT( aI < mNumberOfElementsPerDirection( 0 ),
                          "index i out of bounds: %lu > %lu",
                          ( long unsigned int ) aI,
                          ( long unsigned int ) mNumberOfElementsPerDirection( 0 ) );

            return mElements( aI );

        }

//------------------------------------------------------------------------------

        // element grabber for 2D problem
        mesh::Element *
        LookupTable::get_element( const index_t & aI, const index_t & aJ )
        {
            BELFEM_ASSERT( aI < mNumberOfElementsPerDirection( 0 ),
                          "index i out of bounds: %lu > %lu",
                          ( long unsigned int ) aI,
                          ( long unsigned int ) mNumberOfElementsPerDirection( 0 ) );

            BELFEM_ASSERT( aJ < mNumberOfElementsPerDirection( 1 ),
                          "index j out of bounds: %lu > %lu",
                          ( long unsigned int ) aJ,
                          ( long unsigned int ) mNumberOfElementsPerDirection( 1 ) );

            return mElements( aJ*mNumberOfElementsPerDirection( 0 ) + aI );
        }


//------------------------------------------------------------------------------

        // element grabber for 1D problem
        mesh::Element *
        LookupTable::get_element( const index_t & aI, const index_t & aJ, const index_t & aK )
        {
            BELFEM_ASSERT( aI < mNumberOfElementsPerDirection( 0 ),
                          "index i out of bounds: %lu > %lu",
                          ( long unsigned int ) aI,
                          ( long unsigned int ) mNumberOfElementsPerDirection( 0 ) );

            BELFEM_ASSERT( aJ < mNumberOfElementsPerDirection( 1 ),
                          "index j out of bounds: %lu > %lu",
                          ( long unsigned int ) aJ,
                          ( long unsigned int ) mNumberOfElementsPerDirection( 1 ) );

            BELFEM_ASSERT( aK < mNumberOfElementsPerDirection( 2 ),
                          "index j out of bounds: %lu > %lu",
                          ( long unsigned int ) aK,
                          ( long unsigned int ) mNumberOfElementsPerDirection( 2 ) );

            return mElements( mNumberOfElementsPerDirection( 0 )
                *( mNumberOfElementsPerDirection( 1 ) * aK + aJ ) + aI );
        }

//------------------------------------------------------------------------------

        // element finder for 1D problem
        void
        LookupTable::select_element( const real aX )
        {
            mElement = this->get_element( std::floor( ( aX - mPmin( 0 ) ) / mElementLength( 0 ) ) );
        }

//------------------------------------------------------------------------------

        // element finder for 2D problem
        void
        LookupTable::select_element( const real aX, const real aY )
        {
            mElement = this->get_element(
                    std::floor( ( aX - mPmin( 0 ) ) / mElementLength( 0 ) ),
                    std::floor( ( aY - mPmin( 1 ) ) / mElementLength( 1 ) ) );
        }

//------------------------------------------------------------------------------

        // element finder for 3D problem
        void
        LookupTable::select_element( const real aX, const real aY, const real & aZ )
        {
            mElement = this->get_element(
                    std::floor( ( aX - mPmin( 0 ) ) / mElementLength( 0 ) ),
                    std::floor( ( aY - mPmin( 1 ) ) / mElementLength( 1 ) ),
                    std::floor( ( aZ - mPmin( 2 ) ) / mElementLength( 2 ) ) );
        }

//------------------------------------------------------------------------------

        void
        LookupTable::compute_xi( const uint & aI, const real aX )
        {
            mXi( aI ) = 2.0 * ( aX - mElement->node( 0 )->x( aI )) /
                        mElementLength( aI ) - 1.0;

            BELFEM_ASSERT( mXi( aI ) >= -1.0 && mXi( aI ) <= 1.0,
                    "Wrong element selected" );
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_BS_LOOKUPTABLE_HPP
