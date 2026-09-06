/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_CL_IF_ELEMENTMAPPER_HPP
#define BELFEM_CL_IF_ELEMENTMAPPER_HPP

#include "cl_Element.hpp"
#include "cl_Map.hpp"
#include "cl_IF_InterpolationFunction.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"

namespace belfem
{
    namespace fem
    {
        /**
         * Inverse element map: given physical coordinates aX inside an
         * element, solves x(xi) = aX for the natural coordinates xi.
         * Uses closed-form inverses for TRI3 / QUAD4 / TET4 and Newton
         * iteration (seeded by the linear / bilinear inverse where
         * available) for higher-order or curved elements.
         */
        class ElementMapper
        {
//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            ElementType mElementType = ElementType::UNDEFINED ;

            InterpolationFunctionFactory * mFactory ;
            Map< ElementType, InterpolationFunction * > mFunctions ;

            Matrix< real > mX ;     // node coordinates  ( <num nodes> x <spatial dim> )
            Matrix< real > mJ ;     // jacobian          ( <natural dim> x <natural dim> )
            Matrix< real > mInvJ ;  // jacobian inverse  ( <natural dim> x <natural dim> )
            Matrix< real > mA ;     // bilinear quad coefficients
                                    // rows = {x, y}; cols = {1, xi, eta, xi*eta}
            Matrix< real > mN ;     // shape values      ( 1 x <num nodes> )
            Matrix< real > mNxi ;   // shape derivatives ( <natural dim> x <num nodes> )

            Vector< real > mXi0 ;   // centroid seed for Newton
            Vector< real > mRHS ;   // residual / right-hand side
            bool mIsAffin = false ;

            uint mDim = 0 ;         // spatial dim of mX; 0 = auto-detect from element

            mesh::Element * mElement = nullptr ;
            InterpolationFunction * mFunction = nullptr ;

            // dispatch: closed-form solver, initial guess, in-domain predicate
            bool ( ElementMapper::* mFunEval  )( const Vector< real > & aX, Vector< real > & aXi )       = nullptr ;
            bool ( ElementMapper::* mFunGuess )( const Vector< real > & aX, Vector< real > & aXi )       = nullptr ;
            bool ( ElementMapper::* mFunCheck )( const Vector< real > & aXi ) const                       = nullptr ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * constructor; allocates the interpolation-function factory owned by this mapper
             */
            ElementMapper();

//------------------------------------------------------------------------------

            /**
             * destructor; releases cached interpolation functions and the factory
             */
            ~ElementMapper();

//------------------------------------------------------------------------------

            /**
             * non-copyable: the class owns raw pointers (mFactory, cached
             * InterpolationFunction*); copying would double-delete on destruction
             */
            ElementMapper( const ElementMapper & ) = delete ;
            ElementMapper & operator=( const ElementMapper & ) = delete ;

//------------------------------------------------------------------------------

            /**
             * binds the mapper to an element. Caches its node coordinates,
             * pre-computes the bilinear coefficient matrix mA for quads, and
             * selects the evaluator / guess / inside dispatch for the element type.
             *
             * @param[ in ] aElement element to map (non-owning pointer)
             */
            void
            link( mesh::Element * aElement );

//------------------------------------------------------------------------------

            /**
             * solves x(xi) = aX for xi.
             *
             * @param[ in  ] aX  physical coordinate
             *                   ( \<spatial dim\> x 1 )
             * @param[ out ] aXi natural coordinate
             *                   ( \<natural dim\> x 1 )
             * @return true if aXi lies inside the element's reference domain
             *         (within BELFEM_MESH_EPSILON), false otherwise. The Newton
             *         path aborts with BELFEM_ERROR after 100 iterations without
             *         convergence; it does not return false.
             */
            bool
            evaluate( const Vector< real > & aX, Vector< real > & aXi );

//------------------------------------------------------------------------------

            /**
             * sets the spatial dimension of the physical coordinates. Call
             * this before link() to override; if left at 0 (default), the
             * dimension is auto-detected from the first linked element.
             */
            void
            set_dimension( const uint aDim );

//------------------------------------------------------------------------------

            /**
            * computes the weights by evaluating the shape function
            * @param aXi
            * @param aWeights
            */
            void
            weights( const Vector< real > & aXi, Vector< real > & aWeights );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            /**
             * closed-form inverse for the linear triangle (TRI3).
             * Also serves as the Newton seed for higher-order TRI elements:
             * when wired as mFunGuess the boolean return is discarded and
             * only the affine xi is consumed.
             */
            bool
            evaluate_tri3( const Vector< real > & aX, Vector< real > & aXi );

//------------------------------------------------------------------------------

            /**
             * closed-form inverse for the bilinear quadrilateral (QUAD4):
             * uses a 2x2 affine system when the xi*eta coefficients vanish,
             * otherwise solves a quadratic in eta and recovers xi.
             * Also serves as the Newton seed for higher-order QUAD elements.
             */
            bool
            evaluate_quad4( const Vector< real > & aX, Vector< real > & aXi );

//------------------------------------------------------------------------------

            /**
             * closed-form inverse for the linear tetrahedron (TET4)
             */
            bool
            evaluate_tet4( const Vector< real > & aX, Vector< real > & aXi );

//------------------------------------------------------------------------------

            /**
             * Newton iteration  xi := xi - J^{-1} ( N(xi) X - aX ),
             * starting from mFunGuess and validated by mFunCheck.
             * Capped at 100 iterations.
             */
            bool
            evaluate_general( const Vector< real > & aX, Vector< real > & aXi );

//------------------------------------------------------------------------------

            /**
             * seeds aXi with the reference-domain centroid mXi0;
             * used for TET / PENTA / PYRA / HEX where no cheap closed-form
             * seed is provided
             */
            bool
            guess_general( const Vector< real > & aX, Vector< real > & aXi );

//------------------------------------------------------------------------------

            /**
             * inside the reference triangle:
             * { (xi, eta) : xi >= 0, eta >= 0, xi + eta <= 1 }
             */
            bool
            inside_tri( const Vector< real > & aXi ) const ;

//------------------------------------------------------------------------------

            /**
             * inside the reference quadrilateral:
             * { (xi, eta) : |xi| <= 1, |eta| <= 1 }
             */
            bool
            inside_quad( const Vector< real > & aXi ) const ;

//------------------------------------------------------------------------------

            /**
             * inside the reference tetrahedron:
             * { xi, eta, zeta >= 0 ; xi + eta + zeta <= 1 }
             */
            bool
            inside_tet( const Vector< real > & aXi ) const ;

//------------------------------------------------------------------------------

            /**
             * inside the reference pentahedron / triangular prism:
             * triangle in (xi, eta) extruded along zeta in [-1, 1]
             */
            bool
            inside_penta( const Vector< real > & aXi ) const ;

//------------------------------------------------------------------------------

            /**
             * inside the reference pyramid:
             * |xi| <= 1, |eta| <= 1, zeta in [0, 1]
             */
            bool
            inside_pyra( const Vector< real > & aXi ) const ;

//------------------------------------------------------------------------------

            /**
             * inside the reference hexahedron:
             * |xi| <= 1, |eta| <= 1, |zeta| <= 1
             */
            bool
            inside_hex( const Vector< real > & aXi ) const ;

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* end namespace fem */
} /* end namespace belfem */
#endif //BELFEM_CL_IF_ELEMENTMAPPER_HPP
