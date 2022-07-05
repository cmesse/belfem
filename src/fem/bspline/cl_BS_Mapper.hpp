//
// Created by Christian Messe on 29.04.20.
//

#ifndef BELFEM_CL_BS_MAPPER_HPP
#define BELFEM_CL_BS_MAPPER_HPP
#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "cl_BS_TMatrix.hpp"
#include "cl_BS_Basis.hpp"
#include "cl_BS_Element.hpp"
#include "cl_SpMatrix.hpp"
#include "cl_IF_InterpolationFunction.hpp"

namespace belfem
{
    namespace bspline
    {
        /**
         * this is a mapper for tensor grids
         * the purpose is to be able to interpolate a dataset using
         * C^(p-1) continuous derivatives.
         */
        class Mapper
        {
            const uint mNumberOfDimensions ;
            const uint mOrder ;

            // number of elements per dimension
            const Vector< index_t > mNumberOfElementsPerDimension ;

            // number of nodes per dimension
            Vector< index_t> mNumberOfNodesPerDimension ;

            index_t mNumberOfNodes ;
            index_t mNumberOfElements ;

            uint mNumberOfNodesPerElement ;
            uint mNumberOfBasisPerElement ;

            // basic gridlines including the padding elements
            Vector< real > mAxisX;
            Vector< real > mAxisY;
            Vector< real > mAxisZ;

            // Length, with and height of an element
            Vector< real > mElementLength ;

            TMatrix * mTMatrix = nullptr ;

            Mesh * mMesh = nullptr ;

            Matrix< real > mIntegrationPoints ;
            Vector< real > mIntegrationWeights ;

            Matrix< real > mIntegrationGrid ;
            index_t mGridSize ;

            Vector< index_t > mNumberOfBasisPerDimension ;
            index_t mNumberOfBasis ;

            Cell< Basis * > mBasis;

            Cell< Element * > mElements ;

            fem::InterpolationFunction * mInterpolationFunction = nullptr ;

            SpMatrix * mJacobian = nullptr ;

            // mass matrix for one lagrange element
            Matrix< real > mMassMatrix ;

            // RHS matrix for integration points
            Matrix< real > mRhsMatrix ;

            Cell< Vector< real > * > mFields ;
            Cell< string > mFieldLables ;

            Vector< real > mRHS ;
            Vector< real > mDOFs ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             *
             * @param aNumberOfDimensions
             * @param aOrder
             * @param aNumberOfElementsPerDiection
             * @param aMinPoint  point with small coordinates of bounding box
             * @param aMaxPoint  point with high coordinates of bounding box
             */
            Mapper( const uint aNumberOfDimensions,
                    const uint aOrder,
                    const Vector< index_t > & aNumberOfElementsPerDimension,
                    const Vector< real >   & aMinPoint,
                    const Vector< real >   & aMaxPoint
                    );

//------------------------------------------------------------------------------

            ~Mapper();

//------------------------------------------------------------------------------

            // expose the mesh
            Mesh *
            mesh()
            {
                return mMesh ;
            }

//------------------------------------------------------------------------------

            const Matrix< real > &
            integration_grid() const;

//------------------------------------------------------------------------------

            Vector< real > &
            create_field( const string & aLabel );

//------------------------------------------------------------------------------

            Vector< real > &
            field( const index_t aIndex );

//------------------------------------------------------------------------------

            index_t
            number_of_fields() const ;

//------------------------------------------------------------------------------

            void
            compute_node_values();

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            create_t_matrix();

//------------------------------------------------------------------------------

            void
            create_axis( const Vector< real >   & aMinPoint,
                         const Vector< real >   & aMaxPoint );

//------------------------------------------------------------------------------

            void
            create_mesh();

//------------------------------------------------------------------------------

            void
            create_nodes();

//------------------------------------------------------------------------------

            void
            create_lagrange_elements();

//------------------------------------------------------------------------------

            void
            create_global_variables();

//------------------------------------------------------------------------------

            void
            create_integration_points();

//------------------------------------------------------------------------------

            void
            create_integration_grid();

//------------------------------------------------------------------------------

            // create the B-Spline basis
            void
            create_basis();

//------------------------------------------------------------------------------

            void
            create_bspline_elements();

//------------------------------------------------------------------------------

            void
            create_basis_adjency();

//------------------------------------------------------------------------------

            void
            create_interpolation_function();

//------------------------------------------------------------------------------

            void
            create_jacobian();



//------------------------------------------------------------------------------

            void
            compute_element_matrices();

//------------------------------------------------------------------------------

            void
            compute_dofs( const index_t aFieldIndex );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        inline index_t
        Mapper::number_of_fields() const
        {
            return mFields.size() ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_BS_MAPPER_HPP
