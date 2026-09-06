//
// Created by christian on 12/1/21.
//

#ifndef BELFEM_CL_EF_EDGEFUNCTION_HPP
#define BELFEM_CL_EF_EDGEFUNCTION_HPP

#include "typedefs.hpp"
#include "Mesh_Enums.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace mesh
    {
        class Element ;
    }

    namespace fem
    {
        class Element ;

        /**
         * the edge function base class
         *
         * @ingroup grp_fem_interpolation
         * @see @ref fem_interpolation_nedelec
         */
        class EdgeFunction
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            //! the Geometry Jacobian (transposed)
            Matrix< real > mJ;

            //! the inverse of the Geometry Jacobian (transposed)
            Matrix< real > mInvJ;

            //! the determinant of the Geometry Jacobian
            real mDetJ = BELFEM_QUIET_NAN ;

            //! the absolute value determinant of the Geometry Jacobian
            real mAbsDetJ = BELFEM_QUIET_NAN ;

            // values for the shape function
            Matrix< real > mE ;

            //! matrix containing the curl for the H-Function
            Matrix< real > mC ;

            //! matrix containing the gradient operator for the H-Function.
            //! Named mGrad rather than mG because seven subclasses already
            //! carry a private shape-coefficient matrix mG ( the mF/mG/mH
            //! triple ) that would shadow a base-class mG.
            //!
            //! Layout ( one convention for every subclass, do not deviate ):
            //! G( aIndex ) returns the ( d*d ) x nDofs matrix of the basis
            //! function gradients at integration point aIndex, d being the
            //! spatial dimension. Column e is the column-major vectorization
            //! of the d x d tensor grad( w_e ) with the convention
            //! ( grad h )_ij = d h_j / d x_i :
            //!
            //!     G( i + d*j , e ) = d ( w_e )_j / d x_i
            //!
            //! shapes and row order:
            //!   3D -> 9 x nDofs, rows: dHx/dx dHx/dy dHx/dz
            //!                          dHy/dx dHy/dy dHy/dz
            //!                          dHz/dx dHz/dy dHz/dz
            //!   2D -> 4 x nDofs, rows: dHx/dx dHx/dy dHy/dx dHy/dy
            //!
            //! derived quantities:
            //!   divergence : div h = sum_i G( i + d*i, : ) * q
            //!   curl tie ( all components explicit — do NOT cycle the row
            //!   integers, cycle (x,y,z) in the formula G(i+3j)-G(j+3i) ):
            //!     3D: ( curl h )_x = [ row(7) - row(5) ] * q
            //!         ( curl h )_y = [ row(2) - row(6) ] * q
            //!         ( curl h )_z = [ row(3) - row(1) ] * q
            //!     2D:   curl h     = [ row(2) - row(1) ] * q
            //!   Must reproduce C() to round-off — this is the anchor test
            //!   for every implementation.
            //! The matrix is sized in the subclass constructors alongside
            //! mE and mC — except EF_LINE3, which deliberately leaves it
            //! unsized until its 1D-manifold gradient contract is decided
            //! ( same pattern as its unsized mC ). The base class does not
            //! allocate it.
            Matrix< real > mGrad ;

            //! sum of all weights
            real mSumW ;
            uint mNumDofs ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            EdgeFunction() = default;

//------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            virtual ~EdgeFunction() = default;

//------------------------------------------------------------------------------

            /**
             * links the shape function with the element and precomputes data
             * @param aElement
             */
             virtual void
             link( Element * aElement ) = 0 ;

//------------------------------------------------------------------------------

             /**
              * only needed for higher order elements
              * @param aXi
              */
             virtual void
             precompute( const Matrix< real > & aXi ) = 0 ;

//------------------------------------------------------------------------------

            // compute the edge function
            virtual const Matrix< real > &
            E( const uint aIndex ) ;

//------------------------------------------------------------------------------

            // compute the curl function
            virtual const Matrix< real > &
            C( const uint aIndex = 0 ) = 0 ;

//------------------------------------------------------------------------------

            // compute the gradient operator of the edge function
            // ( layout contract: see the mGrad member documentation )
            virtual const Matrix< real > &
            G( const uint aIndex = 0 ) = 0 ;

//------------------------------------------------------------------------------

            /**
             * returns the current value of the determinant
             */
            real
            det_J() const ;

//------------------------------------------------------------------------------

            /**
             * returns the current value of the determinant
             */
            real
            abs_det_J() const ;

//------------------------------------------------------------------------------

            /**
             * returns the sum of all integration weights
             */
            real
            sum_w() const ;

//------------------------------------------------------------------------------

            /**
             *
             * @return number of dofs for this element
             */
            uint
            ndofs() const ;

//------------------------------------------------------------------------------

            virtual void
            update_nabla( const uint aIndex );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        inline uint
        EdgeFunction::ndofs() const
        {
            return mNumDofs ;
        }

//------------------------------------------------------------------------------

        inline real
        EdgeFunction::sum_w() const
        {
            return mSumW ;
        }

//------------------------------------------------------------------------------

        inline real
        EdgeFunction::det_J() const
        {
            return mDetJ;
        }

//------------------------------------------------------------------------------

        inline real
        EdgeFunction::abs_det_J() const
        {
            return mAbsDetJ ;
        }

//------------------------------------------------------------------------------

        inline const Matrix <real> &
        EdgeFunction::E( const uint aIndex )
        {
            BELFEM_ERROR( false, "Invalid call to edge function");
            return mE;
        }

//------------------------------------------------------------------------------

        inline const Matrix <real> &
        EdgeFunction::C( const uint aIndex )
        {
            BELFEM_ERROR( false, "Invalid call to curl function");
            return mC;
        }

//------------------------------------------------------------------------------

        inline void
        EdgeFunction::update_nabla( const uint aIndex )
        {
            BELFEM_ERROR( false, "Invalid call to update_nabla()");
        }

//------------------------------------------------------------------------------

    } /* end namespace fem */
}  /* end namespace belfem */


#endif //BELFEM_CL_EF_EDGEFUNCTION_HPP
