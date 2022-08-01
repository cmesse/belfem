//
// Created by Christian Messe on 28.04.20.
//

#ifndef BELFEM_CL_BS_TMATRIX_HPP
#define BELFEM_CL_BS_TMATRIX_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "Mesh_Enums.hpp"

namespace belfem
{
    namespace bspline
    {
        class TMatrix
        {
            const uint mNumberOfDimensions;
            const uint mBsplineOrder ;
            const uint mLagrangeOrder ;
            const uint mNumberOfBasis ;
            const uint mNumberOfNodes ;
            const ElementType mBsplineType ;
            const ElementType mLagrangeType ;
            Matrix< index_t > mBasisIndex;
            Matrix< index_t > mNodeIndex;
            Matrix< real >    mNodeParamCoords;

            Matrix< real > mLagrangeMatrix ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            TMatrix(
                    const uint aNumberOfDimensions,
                    const uint aBsplineOrder,
                    const uint aLangrangeOrder=0 );

//------------------------------------------------------------------------------

            ~TMatrix() = default ;

//------------------------------------------------------------------------------

            inline ElementType
            lagrange_type() const
            {
                return mLagrangeType ;
            }

//------------------------------------------------------------------------------

            inline ElementType
            bspline_type() const
            {
                return mBsplineType ;
            }

//---------------------------------------------------------------------------

            inline const Matrix< index_t > &
            node_index() const
            {
                return mNodeIndex;
            };

//---------------------------------------------------------------------------

            inline const Matrix< index_t > &
            basis_index() const
            {
                return mBasisIndex;
            };


//------------------------------------------------------------------------------

            // the matrix that converts node dofs to basis dofs
            inline const Matrix< real > &
            Lagrange() const
            {
                return mLagrangeMatrix ;
            }

// ------------------------------------------------------------------------------
        private :
//------------------------------------------------------------------------------

            ElementType
            select_type( const uint aNumDimensions, const uint aOrder );

//------------------------------------------------------------------------------

            void
            set_index( const ElementType aType, Matrix< index_t > & aIndex );

//------------------------------------------------------------------------------

            void
            populate_node_param_coords();

//------------------------------------------------------------------------------

            real
            shape_function_1d(
                    const uint & aK,
                    const real aXi ) ;

//------------------------------------------------------------------------------

            void
            shape_function(
                    const real        & aXi,
                    Vector< real >    & aN );

//------------------------------------------------------------------------------

            void
            shape_function(
                    const real        & aXi,
                    const real        & aEta,
                    Vector< real >    & aN );

//------------------------------------------------------------------------------

            void
            shape_function(
                    const real             & aXi,
                    const real             & aEta,
                    const real             & aZeta,
                    Vector< real > & aN );

//------------------------------------------------------------------------------

            void
            compute_lagrange_matrix();

//------------------------------------------------------------------------------

        };
    }
}
#endif //BELFEM_CL_BS_TMATRIX_HPP
