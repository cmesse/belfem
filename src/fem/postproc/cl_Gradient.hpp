//
// Created by christian on 2/10/25.
//

#ifndef CL_GRADIENT_HPP
#define CL_GRADIENT_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Cell.hpp"
#include "cl_Mesh.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_Postprocessor.hpp"
#include "fn_trans.hpp"

namespace belfem
{
    namespace fem
    {
        /**
         * @brief Gradient recovery on a solved mesh.
         *
         * @ingroup grp_fem_postproc
         * @see @ref fem_postproc_index
         */
        class Gradient : public Postprocessor
        {
            const int mNumDimensions ;

            // flag telling if we need the negative gradient
            bool mFlipSign = false ;

            // tells which dof manager is used, default: 0
            const uint mDofMaganerIndex ;

            Map< id_t, Block * > mBlocks ;
            Cell< mesh::Node * > mNodes ;
            Matrix< real > mX ; // node coordinates
            Matrix< real > mP ; // polynomial vector
            Matrix< real > mC ; // coefficient matrix
            Matrix< real > mV ; // vandermonde matrix
            Vector< int_t >  mPivot ; // pivot for lapack
            Vector< real > mG ; // gradient

            Matrix< real > mJ ;
            Matrix< real > mB ;
            Vector< real > mPhi ;
            Matrix< real > mElX ;

            int mN ; // number of coefficients

            uint mOrder = 0 ;

            string mScalarField ;
            Cell< string > mGradientFields ;

            void
            ( Gradient::*mFunComputePoly )( const Matrix< real > & aX );

            void
            ( Gradient::*mFunPoly2D )( const real x, const real y );

            void
            ( Gradient::*mFunPoly3D )( const real x, const real y, const real z );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

             Gradient(
                 Kernel * aKernel,
                 const Vector< id_t > & aBlocksIDs,
                 const bool aFlipSign = false,
                 const uint aDofManagerIndex = 0 );

             ~Gradient() override = default;

            void
            set_fields( const string & aScalarField, const string & aGradientField );

            virtual void
            process_node( mesh::Node * aNode );

            void
            set_order( const uint aOrder );

            void
            run() override;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            compute_poly( const Matrix< real > & aX );

            void
            compute_poly_2d( const Matrix< real > & aX );

            void
            compute_poly_3d( const Matrix< real > & aX );

            void
            poly1_2d( const real x, const real y );

            void
            poly2_2d( const real x, const real y );

            void
            poly3_2d( const real x, const real y );

            void
            poly4_2d( const real x, const real y );

            void
            poly1_3d( const real x, const real y, const real z  );

            void
            poly2_3d( const real x, const real y, const real z  );

            void
            poly3_3d( const real x, const real y, const real z  );

            void
            poly4_3d( const real x, const real y, const real z  );

            uint
            subprocess_node( mesh::Node * aNode );
        };

        inline void
        Gradient::compute_poly( const Matrix< real > & aX )
        {
            (this->*mFunComputePoly)( aX );
        }

        inline void
        Gradient::compute_poly_2d( const Matrix< real > & aX )
        {
            (this->*mFunPoly2D )( aX( 0, 0 ), aX( 0, 1 ) );
        }

        inline void
        Gradient::compute_poly_3d( const Matrix< real > & aX )
        {
            (this->*mFunPoly3D )( aX( 0, 0 ), aX( 0, 1 ), aX( 0, 2 ) );
        }

        inline void
        Gradient::poly1_2d( const real x, const real y )
        {
            mP( 0, 0 ) = 1.0 ;
            mP( 1, 0 ) = x ;
            mP( 2, 0 ) = y ;
        }

        inline void
        Gradient::poly2_2d( const real x, const real y )
        {
            mP( 0, 0 ) = 1.0 ;
            mP( 1, 0 ) = x ;
            mP( 2, 0 ) = y ;
            mP( 3, 0 ) = x*x ;
            mP( 4, 0 ) = x*y ;
            mP( 5, 0 ) = y*y ;
        }

        inline void
        Gradient::poly3_2d( const real x, const real y )
        {
            mP( 0, 0 ) = 1.0 ;
            mP( 1, 0 ) = x ;
            mP( 2, 0 ) = y ;
            mP( 3, 0 ) = x*x ;
            mP( 4, 0 ) = x*y ;
            mP( 5, 0 ) = y*y ;
            mP( 6, 0 ) = mP( 3, 0 )*x ;
            mP( 7, 0 ) = mP( 3, 0 )*y ;
            mP( 8, 0 ) = x * mP( 5, 0 );
            mP( 9, 0 ) = y * mP( 5, 0 );
        }

        inline void
        Gradient::poly4_2d( const real x, const real y )
        {
            mP(  0, 0 ) = 1.0 ;
            mP(  1, 0 ) = x ;
            mP(  2, 0 ) = y ;
            mP(  3, 0 ) = x*x ;
            mP(  4, 0 ) = x*y ;
            mP(  5, 0 ) = y*y ;
            mP(  6, 0 ) = mP( 3, 0 )*x ;  // x^3
            mP(  7, 0 ) = mP( 3, 0 )*y ;  // x^2 * y
            mP(  8, 0 ) = x * mP( 5, 0 ); // x * y^2
            mP(  9, 0 ) = y * mP( 5, 0 ); // y^3
            mP( 10, 0 ) = mP( 3, 0 ) * mP( 3, 0 ); // x^4
            mP( 11, 0 ) = mP( 6, 0 ) * y ; // x^3 * y
            mP( 12, 0 ) = mP( 3, 0 ) * mP( 5, 0 ) ; // x^2*y^2
            mP( 13, 0 ) = x * mP( 9, 0 ) ; // x * y^3
            mP( 14, 0 ) =  mP( 5, 0 ) * mP( 5, 0 ) ; // y^4
        }

        inline void
        Gradient::poly1_3d( const real x, const real y, const real z )
        {
            mP( 0, 0 ) = 1.0 ;
            mP( 1, 0 ) = x ;
            mP( 2, 0 ) = y ;
            mP( 3, 0 ) = z ;
        }

        inline void
        Gradient::poly2_3d( const real x, const real y, const real z )
        {
            mP(  0, 0 ) = 1.0 ;
            mP(  1, 0 ) = x ;
            mP(  2, 0 ) = y ;
            mP(  3, 0 ) = z ;
            mP(  4, 0 ) = x*x ;
            mP(  5, 0 ) = x*y ;
            mP(  6, 0 ) = y*y ;
            mP(  7, 0 ) = y*z ;
            mP(  8, 0 ) = z*z ;
            mP(  9, 0 ) = z*x ;
        }

        inline void
        Gradient::poly3_3d( const real x, const real y, const real z )
        {
            mP(  0, 0 ) = 1.0 ;
            mP(  1, 0 ) = x ;
            mP(  2, 0 ) = y ;
            mP(  3, 0 ) = z ;
            mP(  4, 0 ) = x*x ;
            mP(  5, 0 ) = x*y ;
            mP(  6, 0 ) = y*y ;
            mP(  7, 0 ) = y*z ;
            mP(  8, 0 ) = z*z ;
            mP(  9, 0 ) = z*x ;
            mP( 10, 0 ) = mP( 4, 0 ) * x ; // x^3
            mP( 11, 0 ) = mP( 4, 0 ) * y ; // x^2 * y
            mP( 12, 0 ) = x * mP( 6, 0 ) ; // x * y^2 ;
            mP( 13, 0 ) = y * mP( 6, 0 ) ; // y^3
            mP( 14, 0 ) = z * mP( 6, 0 ) ; // y^2 * z
            mP( 15, 0 ) = y * mP( 8, 0 ) ; // y * z^2
            mP( 16, 0 ) = z * mP( 8, 0 ) ; // z^3
            mP( 17, 0 ) = x * mP( 8, 0 ) ; // z^2 * x
            mP( 18, 0 ) = mP( 4, 0 ) * z ; // x^2 * z
            mP( 19, 0 ) = x * y * z ;
        }

        inline void
        Gradient::poly4_3d( const real x, const real y, const real z )
        {
            mP(  0, 0 ) = 1.0 ;
            mP(  1, 0 ) = x ;
            mP(  2, 0 ) = y ;
            mP(  3, 0 ) = z ;
            mP(  4, 0 ) = x*x ;
            mP(  5, 0 ) = x*y ;
            mP(  6, 0 ) = y*y ;
            mP(  7, 0 ) = y*z ;
            mP(  8, 0 ) = z*z ;
            mP(  9, 0 ) = z*x ;
            mP( 10, 0 ) = mP( 4, 0 ) * x ; // x^3
            mP( 11, 0 ) = mP( 4, 0 ) * y ; // x^2 * y
            mP( 12, 0 ) = x * mP( 6, 0 ) ; // x * y^2 ;
            mP( 13, 0 ) = y * mP( 6, 0 ) ; // y^3
            mP( 14, 0 ) = z * mP( 6, 0 ) ; // y^2 * z
            mP( 15, 0 ) = y * mP( 8, 0 ) ; // y * z^2
            mP( 16, 0 ) = z * mP( 8, 0 ) ; // z^3
            mP( 17, 0 ) = x * mP( 8, 0 ) ; // z^2 * x
            mP( 18, 0 ) = mP( 4, 0 ) * z ; // x^2 * z
            mP( 19, 0 ) = x * y * z ;
            mP( 20, 0 ) = mP( 4, 0 ) * mP( 4, 0 ) ; // x^4
            mP( 21, 0 ) = mP( 4, 0 ) * mP( 5, 0 ) ; // x^3 * y
            mP( 22, 0 ) = mP( 4, 0 ) * mP( 6, 0 ) ; // x^2 * y^2
            mP( 23, 0 ) = mP( 5, 0 ) * mP( 6, 0 ) ; // x * y^3
            mP( 24, 0 ) = mP( 6, 0 ) * mP( 6, 0 ) ; // y^4
            mP( 25, 0 ) = mP( 6, 0 ) * mP( 7, 0 ) ; // y^3 * z
            mP( 26, 0 ) = mP( 6, 0 ) * mP( 8, 0 ) ; // y^2 * z^2
            mP( 27, 0 ) = mP( 7, 0 ) * mP( 8, 0 ) ; // y * z^3
            mP( 28, 0 ) = mP( 8, 0 ) * mP( 8, 0 ) ; // z^4
            mP( 29, 0 ) = mP( 9, 0 ) * mP( 8, 0 ) ; // z^3 * x
            mP( 30, 0 ) = mP( 4, 0 ) * mP( 8, 0 ) ; // x^2 * z^2
            mP( 31, 0 ) = mP( 9, 0 ) * mP( 4, 0 ) ; // z * x^3
            mP( 32, 0 ) = mP(  4, 0 ) * mP( 7, 0  ) ; // x^2 * y * z
            mP( 33, 0 ) = mP( 6, 0 ) * mP( 9, 0 ) ; // x * y^2 * z
            mP( 34, 0 ) = mP( 5, 0 ) * mP( 8, 0 ) ; // x * y * z^2
        }
    }
}
#endif //CL_GRADIENT_HPP
