//
// Created by Christian Messe on 25.10.19.
//

#ifndef BELFEM_CL_IF_QUAD16_HPP
#define BELFEM_CL_IF_QUAD16_HPP

#include "cl_IF_InterpolationFunctionTemplate.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        InterpolationOrder
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 16 >
        ::interpolation_order() const
        {
            return InterpolationOrder::CUBIC;
        }

//------------------------------------------------------------------------------

        template<>
        ElementType
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 16 >
        ::element_type() const
        {
            return ElementType::QUAD16;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 16 >
        ::param_coords( Matrix< real > & aXiHat )  const
        {
            aXiHat.set_size( 2, 16 );

            const real c = 1.0/3.0;
            aXiHat( 0,  0 ) = -1.000000;
            aXiHat( 1,  0 ) = -1.000000;

            aXiHat( 0,  1 ) =  1.000000;
            aXiHat( 1,  1 ) = -1.000000;

            aXiHat( 0,  2 ) =  1.000000;
            aXiHat( 1,  2 ) =  1.000000;

            aXiHat( 0,  3 ) = -1.000000;
            aXiHat( 1,  3 ) =  1.000000;

            aXiHat( 0,  4 ) = -c;
            aXiHat( 1,  4 ) = -1.000000;

            aXiHat( 0,  5 ) =  c;
            aXiHat( 1,  5 ) = -1.000000;

            aXiHat( 0,  6 ) =  1.000000;
            aXiHat( 1,  6 ) = -c;

            aXiHat( 0,  7 ) =  1.000000;
            aXiHat( 1,  7 ) =  c;

            aXiHat( 0,  8 ) =  c;
            aXiHat( 1,  8 ) =  1.000000;

            aXiHat( 0,  9 ) = -c;
            aXiHat( 1,  9 ) =  1.000000;

            aXiHat( 0, 10 ) = -1.000000;
            aXiHat( 1, 10 ) =  c;

            aXiHat( 0, 11 ) = -1.000000;
            aXiHat( 1, 11 ) = -c;

            aXiHat( 0, 12 ) = -c;
            aXiHat( 1, 12 ) = -c;

            aXiHat( 0, 13 ) =  c;
            aXiHat( 1, 13 ) = -c;

            aXiHat( 0, 14 ) =  c;
            aXiHat( 1, 14 ) =  c;

            aXiHat( 0, 15 ) = -c;
            aXiHat( 1, 15 ) =  c;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 16 >::N(
                const Vector< real > & aXi,
                Matrix< real > & aN  ) const
        {
            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );

            const real a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 )*0.0625;
            const real a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) )*0.0625;
            const real a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) )*0.0625;
            const real a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 )*0.0625;

            const real b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 )*0.0625;
            const real b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) )*0.0625;
            const real b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) )*0.0625;
            const real b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 )*0.0625;

            // populate matrix with values
            aN.set_size(1,16);
            aN(  0, 0 ) = a0*b0;
            aN(  0, 1 ) = a3*b0;
            aN(  0, 2 ) = a3*b3;
            aN(  0, 3 ) = a0*b3;
            aN(  0, 4 ) = a1*b0;
            aN(  0, 5 ) = a2*b0;
            aN(  0, 6 ) = a3*b1;
            aN(  0, 7 ) = a3*b2;
            aN(  0, 8 ) = a2*b3;
            aN(  0, 9 ) = a1*b3;
            aN( 0, 10 ) = a0*b2;
            aN( 0, 11 ) = a0*b1;
            aN( 0, 12 ) = a1*b1;
            aN( 0, 13 ) = a2*b1;
            aN( 0, 14 ) = a2*b2;
            aN( 0, 15 ) = a1*b2;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 16 >::dNdXi(
                const Vector< real > & aXi,
                Matrix< real > & adNdXi  ) const
        {
            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );

            // often used parameters
            const real a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 ) * 0.0625;
            const real a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) ) * 0.0625;
            const real a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) ) * 0.0625;
            const real a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 ) * 0.0625;

            const real b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 ) * 0.0625;
            const real b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) ) * 0.0625;
            const real b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) ) * 0.0625;
            const real b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 ) * 0.0625;

            const real da0 = (   1.0 + xi*( 18.0 - 27.0*xi )) * 0.0625;
            const real da1 = ( -27.0 - xi*( 18.0 - 81.0*xi )) * 0.0625;
            const real da2 = (  27.0 - xi*( 18.0 + 81.0*xi )) * 0.0625;
            const real da3 = (  -1.0 + xi*( 18.0 + 27.0*xi )) * 0.0625;

            const real db0 = (   1.0 + eta*( 18.0 - 27.0*eta )) * 0.0625;
            const real db1 = ( -27.0 - eta*( 18.0 - 81.0*eta )) * 0.0625;
            const real db2 = (  27.0 - eta*( 18.0 + 81.0*eta )) * 0.0625;
            const real db3 = (  -1.0 + eta*( 18.0 + 27.0*eta )) * 0.0625;

            // populate output matrix
            adNdXi.set_size( 2, 16 );

            adNdXi( 0,  0 ) = da0*b0;
            adNdXi( 1,  0 ) = a0*db0;

            adNdXi( 0,  1 ) = da3*b0;
            adNdXi( 1,  1 ) = a3*db0;

            adNdXi( 0,  2 ) = da3*b3;
            adNdXi( 1,  2 ) = a3*db3;

            adNdXi( 0,  3 ) = da0*b3;
            adNdXi( 1,  3 ) = a0*db3;

            adNdXi( 0,  4 ) = da1*b0;
            adNdXi( 1,  4 ) = a1*db0;

            adNdXi( 0,  5 ) = da2*b0;
            adNdXi( 1,  5 ) = a2*db0;

            adNdXi( 0,  6 ) = da3*b1;
            adNdXi( 1,  6 ) = a3*db1;

            adNdXi( 0,  7 ) = da3*b2;
            adNdXi( 1,  7 ) = a3*db2;

            adNdXi( 0,  8 ) = da2*b3;
            adNdXi( 1,  8 ) = a2*db3;

            adNdXi( 0,  9 ) = da1*b3;
            adNdXi( 1,  9 ) = a1*db3;

            adNdXi( 0, 10 ) = da0*b2;
            adNdXi( 1, 10 ) = a0*db2;

            adNdXi( 0, 11 ) = da0*b1;
            adNdXi( 1, 11 ) = a0*db1;

            adNdXi( 0, 12 ) = da1*b1;
            adNdXi( 1, 12 ) = a1*db1;

            adNdXi( 0, 13 ) = da2*b1;
            adNdXi( 1, 13 ) = a2*db1;

            adNdXi( 0, 14 ) = da2*b2;
            adNdXi( 1, 14 ) = a2*db2;

            adNdXi( 0, 15 ) = da1*b2;
            adNdXi( 1, 15 ) = a1*db2;
        }

//------------------------------------------------------------------------------

        template<>
        void
        InterpolationFunctionTemplate<
                GeometryType::QUAD, InterpolationType::LAGRANGE, 2, 16 >::d2NdXi2(
                const Vector< real > & aXi,
                Matrix< real > & ad2NdXi2  ) const
        {
            // unpack xi and eta from input vector
            const real  xi = aXi( 0 );
            const real eta = aXi( 1 );

            // often used parameters
            const real a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 ) * 0.0625;
            const real a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) ) * 0.0625;
            const real a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) ) * 0.0625;
            const real a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 ) * 0.0625;

            const real b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 ) * 0.0625;
            const real b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) ) * 0.0625;
            const real b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) ) * 0.0625;
            const real b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 ) * 0.0625;

            const real da0 = (   1.0 + xi*( 18.0 - 27.0*xi )) * 0.0625;
            const real da1 = ( -27.0 - xi*( 18.0 - 81.0*xi )) * 0.0625;
            const real da2 = (  27.0 - xi*( 18.0 + 81.0*xi )) * 0.0625;
            const real da3 = (  -1.0 + xi*( 18.0 + 27.0*xi )) * 0.0625;

            const real dda0 = ( 18.0 - 54.0*xi ) * 0.0625;
            const real dda1 = ( 162.0*xi - 18.0 ) * 0.0625;
            const real dda2 = ( - 162.0*xi - 18.0 ) * 0.0625;
            const real dda3 = ( 54.0*xi + 18.0 ) * 0.0625;

            const real db0 = (   1.0 + eta*( 18.0 - 27.0*eta )) * 0.0625;
            const real db1 = ( -27.0 - eta*( 18.0 - 81.0*eta )) * 0.0625;
            const real db2 = (  27.0 - eta*( 18.0 + 81.0*eta )) * 0.0625;
            const real db3 = (  -1.0 + eta*( 18.0 + 27.0*eta )) * 0.0625;

            const real ddb0 = ( 18.0 - 54.0*eta ) * 0.0625;
            const real ddb1 = ( 162.0*eta - 18.0 ) * 0.0625;
            const real ddb2 = ( - 162.0*eta - 18.0 ) * 0.0625;
            const real ddb3 = ( 54.0*eta + 18.0 ) * 0.0625;

            ad2NdXi2.set_size( 3, 16 );
            ad2NdXi2( 0,   0 ) = dda0*b0;
            ad2NdXi2( 1,   0 ) = a0*ddb0;
            ad2NdXi2( 2,   0 ) = da0*db0;

            ad2NdXi2( 0,   1 ) = dda3*b0;
            ad2NdXi2( 1,   1 ) = a3*ddb0;
            ad2NdXi2( 2,   1 ) = da3*db0;

            ad2NdXi2( 0,   2 ) = dda3*b3;
            ad2NdXi2( 1,   2 ) = a3*ddb3;
            ad2NdXi2( 2,   2 ) = da3*db3;

            ad2NdXi2( 0,   3 ) = dda0*b3;
            ad2NdXi2( 1,   3 ) = a0*ddb3;
            ad2NdXi2( 2,   3 ) = da0*db3;

            ad2NdXi2( 0,   4 ) = dda1*b0;
            ad2NdXi2( 1,   4 ) = a1*ddb0;
            ad2NdXi2( 2,   4 ) = da1*db0;

            ad2NdXi2( 0,   5 ) = dda2*b0;
            ad2NdXi2( 1,   5 ) = a2*ddb0;
            ad2NdXi2( 2,   5 ) = da2*db0;

            ad2NdXi2( 0,   6 ) = dda3*b1;
            ad2NdXi2( 1,   6 ) = a3*ddb1;
            ad2NdXi2( 2,   6 ) = da3*db1;

            ad2NdXi2( 0,   7 ) = dda3*b2;
            ad2NdXi2( 1,   7 ) = a3*ddb2;
            ad2NdXi2( 2,   7 ) = da3*db2;

            ad2NdXi2( 0,   8 ) = dda2*b3;
            ad2NdXi2( 1,   8 ) = a2*ddb3;
            ad2NdXi2( 2,   8 ) = da2*db3;

            ad2NdXi2( 0,   9 ) = dda1*b3;
            ad2NdXi2( 1,   9 ) = a1*ddb3;
            ad2NdXi2( 2,   9 ) = da1*db3;

            ad2NdXi2( 0,  10 ) = dda0*b2;
            ad2NdXi2( 1,  10 ) = a0*ddb2;
            ad2NdXi2( 2,  10 ) = da0*db2;

            ad2NdXi2( 0,  11 ) = dda0*b1;
            ad2NdXi2( 1,  11 ) = a0*ddb1;
            ad2NdXi2( 2,  11 ) = da0*db1;

            ad2NdXi2( 0,  12 ) = dda1*b1;
            ad2NdXi2( 1,  12 ) = a1*ddb1;
            ad2NdXi2( 2,  12 ) = da1*db1;

            ad2NdXi2( 0,  13 ) = dda2*b1;
            ad2NdXi2( 1,  13 ) = a2*ddb1;
            ad2NdXi2( 2,  13 ) = da2*db1;

            ad2NdXi2( 0,  14 ) = dda2*b2;
            ad2NdXi2( 1,  14 ) = a2*ddb2;
            ad2NdXi2( 2,  14 ) = da2*db2;

            ad2NdXi2( 0,  15 ) = dda1*b2;
            ad2NdXi2( 1,  15 ) = a1*ddb2;
            ad2NdXi2( 2,  15 ) = da1*db2;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_IF_QUAD16_HPP
