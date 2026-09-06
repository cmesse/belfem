//
// Created by Christian Messe on 25.10.19.
//

#ifndef BELFEM_CL_IF_INTERPOLATIONFUNCTIONTEMPLATE_HPP
#define BELFEM_CL_IF_INTERPOLATIONFUNCTIONTEMPLATE_HPP

#include "cl_IF_InterpolationFunction.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------
        /**
         * shape function templated class
         * G : Geometry
         * T : Type
         * D : Dimension
         * B : Number of Basis
         */
        template< GeometryType G, InterpolationType T, uint D, uint B >
        class InterpolationFunctionTemplate : public InterpolationFunction
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            InterpolationFunctionTemplate() = default;

//------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
             ~InterpolationFunctionTemplate() override = default;

//------------------------------------------------------------------------------

            /**
              * returns a matrix containing the parameter coordinates of the nodes
              * < number of dimensions x number of nodes >
              */
            void
            param_coords( Matrix< real > & aXiHat ) const override
            {
                BELFEM_ERROR( false,
                "illegal call to template function InterpolationFunctionTemplate::param_coords()" );
            }

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
            void
            N( const Vector< real > & aXi,
                     Matrix< real > & aN  ) const override
            {
                BELFEM_ERROR( false,
                "illegal call to template function InterpolationFunctionTemplate::N()" );
            }

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
            void
            dNdXi( const Vector< real > & aXi,
                       Matrix< real > & adNdXi  ) const override
            {
                BELFEM_ERROR( false,
                "illegal call to template function InterpolationFunctionTemplate::adNdXi()" );
            }

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
            void
            d2NdXi2( const Vector< real > & aXi,
                           Matrix< real > & ad2NdXi2  ) const override
            {
                BELFEM_ERROR( false,
                 "illegal call to template function InterpolationFunctionTemplate::ad2NdXi2()" );
            }

//------------------------------------------------------------------------------

            /**
             * returns the number of bases for this shape function
             */
            uint
            number_of_bases() const override
            {
                return B;
            }

//------------------------------------------------------------------------------

            /**
             * returns the number of dimensions for this shape function
             */
            uint
            number_of_dimensions() const override
            {
                return D;
            }

//------------------------------------------------------------------------------

            /**
             * returns the interpolation order
             */
            InterpolationOrder
            interpolation_order() const override
            {
                BELFEM_ERROR( false,
                             "illegal call to template function InterpolationFunctionTemplate::interpolation_order()" );

                return InterpolationOrder::UNDEFINED;
            }

//------------------------------------------------------------------------------

            /**
             * returns the interpolation type
             */
            InterpolationType
            interpolation_type() const override
            {
                return T;
            }

//------------------------------------------------------------------------------

            /**
             * returns the geometry type
             */
            GeometryType
            geometry_type() const override
            {
                return G;
            }

//------------------------------------------------------------------------------

            /**
             * returns the element type
             */
            ElementType
            element_type() const override
            {
                BELFEM_ERROR( false,
                             "illegal call to template function InterpolationFunctionTemplate::element_type()" );

                return ElementType::UNDEFINED;
            }

//------------------------------------------------------------------------------

        };

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_IF_INTERPOLATIONFUNCTIONTEMPLATE_HPP
