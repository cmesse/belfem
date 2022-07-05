//
// Created by Christian Messe on 04.11.19.
//

#include "fn_FEM_initialize_shape_function.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        InterpolationFunction *
        initialize_shape_function(
                const ElementType      & aElementType,
                const Matrix< real >   & aXi,
                Cell< Matrix< real > > & aN,
                Cell< Matrix< real > > & adNdXi,
                Cell< Matrix< real > > & ad2NdXi2  )
        {
            // number of integration points
            uint tNumberOfPoints = aXi.n_cols();

            // the factory that creates the shape function
            InterpolationFunctionFactory tFactory;

            // create the shape function
            InterpolationFunction * aShapeFunction
                = tFactory.create_lagrange_function( aElementType );

            // allocate the matrices
            // ( even if we don't need all of them,
            //  initializing them does not hurt )

            // an empty matrix for the cell initialization
            Matrix< real > tEmpty;

            aN.set_size( tNumberOfPoints, tEmpty );
            adNdXi.set_size( tNumberOfPoints, tEmpty );
            ad2NdXi2.set_size( tNumberOfPoints, tEmpty );

            // allocate the matrices for each integration point
            for( uint k=0; k<tNumberOfPoints; ++k )
            {
                // get the parameter coordinates from the integration point
                Vector< real > tPoint ( aXi.col( k ) );

                // evaluate the shape function
                aShapeFunction->N( tPoint, aN( k ) );

                // evaluate the first derivative
                aShapeFunction->dNdXi( tPoint, adNdXi( k ) );

                // evaluate the second derivative
                aShapeFunction->d2NdXi2( tPoint, ad2NdXi2( k ) );
            }

            return aShapeFunction;
        }

//------------------------------------------------------------------------------
    }
}

