//
// Created by Christian Messe on 30.04.20.
//
#include "commtools.hpp"
#include "cl_BS_LookupTable.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"
#include "meshtools.hpp"

namespace belfem
{
    namespace bspline
    {
//------------------------------------------------------------------------------

        LookupTable::LookupTable( const string aPath ) :
            mMesh( new Mesh( aPath, comm_rank(), false ) ),
            mNumberOfDimensions( mMesh->number_of_dimensions() ),
            mElementType( mMesh->block( 1 )->element_type() ),
            mNumberOfNodesPerElement( mesh::number_of_nodes( mElementType ) ),
            mNumberOfFields( mMesh->number_of_fields() ),
            mElements( mMesh->block( 1 )->elements() )
        {
            this->create_interpolation_function() ;
            this->create_field_map();
            this->polpulate_globals() ;
            this->compute_derivative_factors() ;
        }

//------------------------------------------------------------------------------

        LookupTable::~LookupTable()
        {
            delete mInterpolationFunction ;
            delete mMesh ;
        }

//------------------------------------------------------------------------------

        void
        LookupTable::create_interpolation_function()
        {
            // create a factory
            fem::InterpolationFunctionFactory tFactory;

            // create the function
            mInterpolationFunction = tFactory.create_lagrange_function( mElementType );

            // container for parameter coordinates
            mXi.set_size( mNumberOfDimensions );

            // vector for shape funtion
            mN.set_size( 1, mNumberOfNodesPerElement );

            // matrix for first derivatives
            mdNdxi.set_size( mNumberOfDimensions, mNumberOfNodesPerElement );

            // matrix for second derivatives
            md2Ndxi2.set_size(
                    std::pow( 3, mNumberOfDimensions),
                    mNumberOfNodesPerElement );


            // allocate values vector
            mValues.set_size( mNumberOfNodesPerElement );
        }

//------------------------------------------------------------------------------

        void
        LookupTable::create_field_map()
        {
            mFieldMap.clear();

            uint tNumFields = mMesh->number_of_fields() ;

            for ( uint k=0; k<tNumFields; ++k )
            {
                mFieldMap[ mMesh->field( k )->label() ] = k ;
            }
        }

//------------------------------------------------------------------------------

        void
        LookupTable::polpulate_globals()
        {
            // allocate vectors
            mNumberOfElementsPerDirection.set_size( mNumberOfDimensions );

            mPmin.set_size( mNumberOfDimensions );
            mPmax.set_size( mNumberOfDimensions );
            mElementLength.set_size( mNumberOfDimensions );

            mNumberOfElementsPerDirection( 0 ) = ( index_t ) mMesh->global_variable_data("numElemsX" );
            mPmin( 0 ) = mMesh->global_variable_data("Xmin" );
            mPmax( 0 ) = mMesh->global_variable_data("Xmax" );
            mElementLength( 0 ) = mMesh->global_variable_data("DeltaX");

            if( mNumberOfDimensions > 1 )
            {
                mNumberOfElementsPerDirection( 1 ) = ( index_t ) mMesh->global_variable_data("numElemsY" );
                mPmin( 1 ) = mMesh->global_variable_data("Ymin" );
                mPmax( 1 ) = mMesh->global_variable_data("Ymax" );
                mElementLength( 1 ) = mMesh->global_variable_data("DeltaY");

                if( mNumberOfDimensions > 2 )
                {
                    mNumberOfElementsPerDirection( 2 ) = ( index_t ) mMesh->global_variable_data("numElemsZ" );
                    mPmin( 2 ) = mMesh->global_variable_data("Zmin" );
                    mPmax( 2 ) = mMesh->global_variable_data("Zmax" );
                    mElementLength( 2 ) = mMesh->global_variable_data("DeltaZ");
                }
            }
        }

//------------------------------------------------------------------------------

        void
        LookupTable::update_values( const index_t aFieldIndex )
        {
            // check if values are up to data
            if( aFieldIndex != mFieldIndex || mElementID != mElement->id() )
            {
                // get field from mesh
                const Vector< real > & tField = mMesh->field( aFieldIndex )->data();

                // populate data
                for ( uint k = 0; k < mNumberOfNodesPerElement; ++k )
                {
                    mValues( k ) = tField( mElement->node( k )->index());
                }

                // remember field index
                mFieldIndex = aFieldIndex;

                // remember ID of element
                mElementID = mElement->id();
            }
        }

//------------------------------------------------------------------------------

        void
        LookupTable::compute_derivative_factors()
        {
            mScaleDX    = 2.0 / mElementLength( 0 );
            mScaleD2X2  = 4.0 / ( mElementLength( 0 ) * mElementLength( 0 ) );

            if( mNumberOfDimensions > 1 )
            {
                mScaleDY    = 2.0 / mElementLength( 1 );
                mScaleD2Y2  = 4.0 / ( mElementLength( 1 ) * mElementLength( 1 ) );
                mScaleDXDY  = 4.0 / ( mElementLength( 0 ) * mElementLength( 1 ) );

                if( mNumberOfDimensions > 2 )
                {
                    mScaleDZ    = 2.0 / mElementLength( 2 );
                    mScaleD2Z2  = 4.0 / ( mElementLength( 2 ) * mElementLength( 2 ) );
                    mScaleDXDZ  = 4.0 / ( mElementLength( 0 ) * mElementLength( 2 ) );
                    mScaleDYDZ  = 4.0 / ( mElementLength( 1 ) * mElementLength( 2 ) );
                }
            }
        }

//------------------------------------------------------------------------------

        real
        LookupTable::compute_value( const index_t & aFieldIndex, const real aX )
        {
            BELFEM_ASSERT( mNumberOfDimensions == 1, "Table must be of dimension 1" );

            // check if element must be updated
            if( aX != mX )
            {
                // select the element from the mesh
                this->select_element( aX );

                // compute parameter coordinate
                this->compute_xi( 0, aX );

                // remember x value
                mX = aX ;
            }

            // update field if element or field was changed
            this->update_values( aFieldIndex );

            // compute interpolation function
            mInterpolationFunction->N( mXi, mN );

            // compute the value
            Vector< real > tVal( ( mN * mValues ) );

            // return the value
            return  tVal( 0 );
        }

//------------------------------------------------------------------------------

        real
        LookupTable::compute_value(
                const index_t & aFieldIndex,
                const real aX,
                const real aY )
        {
            BELFEM_ASSERT( mNumberOfDimensions == 2, "Table must be of dimension 2" );

            // check if element must be updated
            if( aX != mX || aY != mY )
            {
                // select the element from the mesh
                this->select_element( aX, aY );

                // compute parameter coordinate
                this->compute_xi( 0, aX );
                this->compute_xi( 1, aY );

                // remember x and y values
                mX = aX ;
                mY = aY ;
            }

            // update field if element or field was changed
            this->update_values( aFieldIndex );

            // compute interpolation function
            mInterpolationFunction->N( mXi, mN );

            // compute the value
            Vector< real > tVal( ( mN * mValues ) );

            // return the value
            return  tVal( 0 );
        }

//------------------------------------------------------------------------------

        real
        LookupTable::compute_value(
                const index_t & aFieldIndex,
                const real aX,
                const real aY,
                const real & aZ )
        {
            BELFEM_ASSERT( mNumberOfDimensions == 3, "Table must be of dimension 3" );

            // check if element must be updated
            if( aX != mX || aY != mY )
            {
                // select the element from the mesh
                this->select_element( aX, aY, aZ );

                // compute parameter coordinate
                this->compute_xi( 0, aX );
                this->compute_xi( 1, aY );
                this->compute_xi( 2, aZ );

                // remember x and y values
                mX = aX ;
                mY = aY ;
                mZ = aZ ;
            }

            // update field if element or field was changed
            this->update_values( aFieldIndex );

            // compute interpolation function
            mInterpolationFunction->N( mXi, mN );

            // compute the value
            Vector< real > tVal( ( mN * mValues ) );

            // return the value
            return  tVal( 0 );
        }

//------------------------------------------------------------------------------

        real
        LookupTable::compute_derivative(
                const index_t  aFieldIndex,
                const real     aX )
        {
            BELFEM_ASSERT( mNumberOfDimensions == 1, "Table must be of dimension 1" );

            // check if element must be updated
            if( aX != mX )
            {
                // select the element from the mesh
                this->select_element( aX );

                // compute parameter coordinate
                this->compute_xi( 0, aX );

                // remember x value
                mX = aX ;
            }

            // update field if element or field was changed
            this->update_values( aFieldIndex );

            // compute interpolation function
            mInterpolationFunction->dNdXi( mXi, mdNdxi );

            // compute the value
            Vector< real > aDerivative(mdNdxi * mValues * mScaleDX );
            return aDerivative( 0 );
        }

//------------------------------------------------------------------------------

        Vector< real >
        LookupTable::compute_derivative(
                const index_t aFieldIndex,
                const real aX,
                const real aY )
        {
            BELFEM_ASSERT( mNumberOfDimensions == 2, "Table must be of dimension 2" );

            // check if element must be updated
            if( aX != mX || aY != mY )
            {
                // select the element from the mesh
                this->select_element( aX, aY );

                // compute parameter coordinate
                this->compute_xi( 0, aX );
                this->compute_xi( 1, aY );

                // remember x and y values
                mX = aX ;
                mY = aY ;
            }

            // update field if element or field was changed
            this->update_values( aFieldIndex );

            // compute interpolation function
            mInterpolationFunction->dNdXi( mXi, mdNdxi );

            // compute the values
            Vector< real > aDerivative( mdNdxi * mValues );

            // transform derivatives
            aDerivative( 0 ) *= mScaleDX ;
            aDerivative( 1 ) *= mScaleDY ;

            return aDerivative ;
        }

//------------------------------------------------------------------------------

        Vector< real >
        LookupTable::compute_derivative(
                const index_t aFieldIndex,
                const real aX,
                const real aY,
                const real aZ )
        {
            BELFEM_ASSERT( mNumberOfDimensions == 3, "Table must be of dimension 3" );

            // check if element must be updated
            if ( aX != mX || aY != mY )
            {
                // select the element from the mesh
                this->select_element( aX, aY, aZ );

                // compute parameter coordinate
                this->compute_xi( 0, aX );
                this->compute_xi( 1, aY );
                this->compute_xi( 2, aZ );

                // remember x and y values
                mX = aX;
                mY = aY;
                mZ = aZ;
            }

            // update field if element or field was changed
            this->update_values( aFieldIndex );

            // compute interpolation function
            mInterpolationFunction->dNdXi( mXi, mdNdxi );

            // compute the values
            Vector< real > aDerivative( mdNdxi * mValues );


            // transform derivatives
            aDerivative( 0 ) *= mScaleDX ;
            aDerivative( 1 ) *= mScaleDY ;
            aDerivative( 2 ) *= mScaleDZ ;

            return aDerivative ;
        }
//------------------------------------------------------------------------------

        real
        LookupTable::compute_second_derivative(
                const index_t    aFieldIndex,
                const real       aX )
        {
            BELFEM_ASSERT( mNumberOfDimensions == 1, "Table must be of dimension 1" );

            // check if element must be updated
            if( aX != mX )
            {
                // select the element from the mesh
                this->select_element( aX );

                // compute parameter coordinate
                this->compute_xi( 0, aX );

                // remember x value
                mX = aX ;
            }

            // update field if element or field was changed
            this->update_values( aFieldIndex );

            // compute interpolation function
            mInterpolationFunction->d2NdXi2( mXi, md2Ndxi2 );

            // compute the value
            Vector< real > aDerivative( md2Ndxi2 * mValues * mScaleD2X2 );

            return aDerivative( 0 );
        }

//------------------------------------------------------------------------------

        Vector< real >
        LookupTable::compute_second_derivative(
                const index_t aFieldIndex,
                const real aX,
                const real aY )
        {
            BELFEM_ASSERT( mNumberOfDimensions == 2, "Table must be of dimension 2" );

            // check if element must be updated
            if( aX != mX || aY != mY )
            {
                // select the element from the mesh
                this->select_element( aX, aY );

                // compute parameter coordinate
                this->compute_xi( 0, aX );
                this->compute_xi( 1, aY );

                // remember x and y values
                mX = aX ;
                mY = aY ;
            }

            // update field if element or field was changed
            this->update_values( aFieldIndex );

            // compute interpolation function
            mInterpolationFunction->d2NdXi2( mXi, md2Ndxi2 );

            // compute the values
            Vector< real > aDerivative( md2Ndxi2 * mValues );

            // transform derivatives
            aDerivative( 0 ) *= mScaleD2X2 ;
            aDerivative( 1 ) *= mScaleD2Y2 ;
            aDerivative( 2 ) *= mScaleDXDY ;

            return aDerivative;
        }

//------------------------------------------------------------------------------

        Vector< real >
        LookupTable::compute_second_derivative(
                const index_t aFieldIndex,
                const real aX,
                const real aY,
                const real aZ )
        {
            BELFEM_ASSERT( mNumberOfDimensions == 3, "Table must be of dimension 3" );

            // check if element must be updated
            if ( aX != mX || aY != mY )
            {
                // select the element from the mesh
                this->select_element( aX, aY, aZ );

                // compute parameter coordinate
                this->compute_xi( 0, aX );
                this->compute_xi( 1, aY );
                this->compute_xi( 2, aZ );

                // remember x and y values
                mX = aX;
                mY = aY;
                mZ = aZ;
            }

            // update field if element or field was changed
            this->update_values( aFieldIndex );

            // compute interpolation function
            mInterpolationFunction->d2NdXi2( mXi, md2Ndxi2 );

            // compute the values
            Vector< real > aDerivative( md2Ndxi2 * mValues );


            // transform derivatives
            aDerivative( 0 ) *= mScaleD2X2 ;
            aDerivative( 1 ) *= mScaleD2Y2 ;
            aDerivative( 2 ) *= mScaleD2Z2 ;
            aDerivative( 3 ) *= mScaleDYDZ ;
            aDerivative( 4 ) *= mScaleDXDZ ;
            aDerivative( 5 ) *= mScaleDXDY ;

            return aDerivative ;
        }

//------------------------------------------------------------------------------
    }
}