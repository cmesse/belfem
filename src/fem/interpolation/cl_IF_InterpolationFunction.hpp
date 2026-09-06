//
// Created by Christian Messe on 24.10.19.
//

#ifndef BELFEM_CL_IF_INTERPOLATIONFUNCTION_HPP
#define BELFEM_CL_IF_INTERPOLATIONFUNCTION_HPP

#include "typedefs.hpp"
#include "Mesh_Enums.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

namespace belfem
{
    namespace fem
    {
        /**
         * the shape function base class
          *
          * @ingroup grp_fem_interpolation
          * @see @ref fem_interpolation_interpolation_usage_guide
         */
         class InterpolationFunction
         {
//------------------------------------------------------------------------------
         public:
//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            InterpolationFunction() = default;

//------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            virtual ~InterpolationFunction() = default;

//------------------------------------------------------------------------------

            /**
              * returns a matrix containing the parameter coordinates of the nodes
              * < number of dimensions x number of nodes >
              */
            virtual void
            param_coords( Matrix< real > & aXiHat ) const = 0;

//------------------------------------------------------------------------------

            /**
             * evaluates the shape function at a given point
             *
             * @param[ in ]  aXi parameter coordinates
             *                   ( \<number of dimensions\>  x 1 )
             *
             * @param[ out ] aN  shape function as
             *                   ( 1 x \<number of bases\> )
             *
             */
            virtual void
            N( const Vector< real > & aXi,
                     Matrix< real > & aN  ) const = 0;

//------------------------------------------------------------------------------

            /**
             * calculates the first derivative of the shape function
             * in parameter space
             *
             * @param[ in ] aXi    point where function is evaluated
             *                     ( \<number of dimensions\>  x 1 )
             *
             * @param[ out ] adNdXi ( \<number of dimensions\> x \<number of bases\> )
             *
             *
             */
            virtual void
            dNdXi( const Vector< real > & aXi,
                         Matrix< real > & adNdXi  ) const = 0;

//------------------------------------------------------------------------------

            /**
              * calculates the second derivative of the shape function
              * in parameter space
              * @param[ in ] aXi    point where function is evaluated
              *                     ( \<number of dimensions\>  x 1 )
              *
              * @param[ out ] ad2NdXi2 ( \<number of second-derivative components\> x \<number of bases\> ):
              *              1 x n in 1D, 3 x n in 2D ( xi², eta², xi·eta ),
              *              6 x n in 3D ( xi², eta², zeta², eta·zeta, xi·zeta, xi·eta )
              *
              */
            virtual void
            d2NdXi2( const Vector< real > & aXi,
                           Matrix< real > & ad2NdXi2  ) const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the number of bases for this shape function
             */
            virtual uint
            number_of_bases() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the number of dimensions for this shape function
             */
            virtual uint
            number_of_dimensions() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the interpolation order
             */
            virtual InterpolationOrder
            interpolation_order() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the interpolation type
             */
             virtual InterpolationType
             interpolation_type() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the geometry type
             */
            virtual GeometryType
            geometry_type() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the element type
             */
            virtual ElementType
            element_type() const = 0;

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* end namespace fem */
}  /* end namespace belfem */
#endif //BELFEM_CL_IF_INTERPOLATIONFUNCTION_HPP
