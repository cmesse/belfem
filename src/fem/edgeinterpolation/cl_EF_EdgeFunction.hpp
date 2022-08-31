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
#include "cl_FEM_Group.hpp"

namespace belfem
{
    namespace mesh
    {
        class Element ;
    }

    namespace fem
    {
        class Element ;
        class Group ;

        /**
         * the edge function base class
         */
        class EdgeFunction
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            //! the group for this element
            Group * mGroup = nullptr ;

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

            //! matrix containing the grad operator for the phi function
            Matrix< real > mB ;

            //! matrix containing the curl for the matrix
            Matrix< real > mCurlA ;

            //! sum of all weights
            real mSumW ;

            //! vector containing normal
            Vector< real > mNormal ;

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
             * @param aOperatorCurlH    if the B-Operator is used (phi and az)
             * @param aOperatorGradPhi  if the B-A-Operator is used (ax, ay, az)
             * @param aOperatorCurlA    if the curl operator is used
             */
             virtual void
             link( Element * aElement,
                   const bool aOperatorCurlH=true,
                   const bool aOperatorGradPhi=false,
                   const bool aOperatorCurlA=false ) = 0 ;

//------------------------------------------------------------------------------

             /**
              * only needed for higher order elements
              * @param aXi
              */
             virtual void
             precompute( const Matrix< real > & aXi ) = 0 ;

//------------------------------------------------------------------------------

            /**
             * links the group of the element
             */
             void
             link( Group * aGroup );

//------------------------------------------------------------------------------

            // compute the edge function
            virtual const Matrix< real > &
            E( const uint aIndex ) ;

//------------------------------------------------------------------------------

            // compute the N-Operator
            const Matrix< real > &
            N( const uint aIndex ) ;

//------------------------------------------------------------------------------

            // compute the B-Operator
            virtual const Matrix< real > &
            B( const uint aIndex = 0 ) = 0 ;

//------------------------------------------------------------------------------

            // compute the Ca-Operator that is needed for the matrix
            virtual const Matrix< real > &
            CA( const uint aIndex = 0 ) = 0 ;

//------------------------------------------------------------------------------

            // compute the curl function
            virtual const Matrix< real > &
            C( const uint aIndex = 0 ) = 0 ;

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
             * special function for thin shells. Returns the
             * node derivative of the volume element
             */
            virtual const Matrix< real > &
            vol_N_xi( const uint aIndex = 0 ) ;

//------------------------------------------------------------------------------

            /**
             * special function for thin shells. Returns the
             * node derivative of the volume element
             */
            virtual const Vector< real > &
            normal( const uint aIndex , const Matrix< real > & aX ) ;

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        inline void
        EdgeFunction::link( Group * aGroup )
        {
            mGroup = aGroup ;

            switch( mesh::geometry_type( aGroup->element_type() ) )
            {
                case( GeometryType::LINE ) :
                {
                    mSumW = 2.0 ;
                    break ;
                }
                case( GeometryType::TRI ) :
                {
                    mSumW = 0.5 ;
                    break ;
                }
                case( GeometryType::QUAD ) :
                {
                    mSumW = 4.0 ;
                    break ;
                }
                case( GeometryType::TET ) :
                {
                    mSumW = 1.0/6.0 ;
                    break ;
                }
                case( GeometryType::PENTA ) :
                {
                    mSumW = 1.0 ;
                    break ;
                }
                case( GeometryType::HEX ) :
                {
                    mSumW = 8.0 ;
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid geometry type" );
                }
            }
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

        inline const Matrix< real > &
        EdgeFunction::N( const uint aIndex )
        {
            return mGroup->N( aIndex );
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

        inline const Matrix< real > &
        EdgeFunction::vol_N_xi( const uint aIndex  )
        {
            BELFEM_ERROR( false, "Invalid call to EdgeFunction::vol_N_xi");
            return mC;
        }

//-----------------------------------------------------------------------------

        inline const Vector< real > &
        EdgeFunction::normal( const uint aIndex , const Matrix< real > & aX )
        {
            BELFEM_ERROR( false, "Invalid call to EdgeFunction::normal");
            return mNormal ;
        }

//------------------------------------------------------------------------------

    } /* end namespace fem */
}  /* end namespace belfem */


#endif //BELFEM_CL_EF_EDGEFUNCTION_HPP
