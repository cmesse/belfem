/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_CL_TENSORMESHCONFIG_HPP
#define BELFEM_CL_TENSORMESHCONFIG_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "Mesh_Enums.hpp"

namespace belfem
{
    class TensorMeshConfig
    {

        const uint mOrder ;
        const uint mDimension ;
        const Vector< index_t > mNumNodesPerDim ;
        const Vector< real > mNodeSteps ;
        const Vector< real > mOrigin ;
        const Vector< index_t > mNumElementsPerDim ;
        const Vector< index_t > mNumControlPointsPerDim ;

        const Vector< real > mLimit ;
        const Vector< real > mElementSteps ;
        const Vector< real > mInvElementSteps ;
        const uint mNumNodes ;
        const uint mNumElements ;
        const uint mNumControlPoints ;

        const ElementType    mElementType ;
    public:

        TensorMeshConfig(
            const uint                aOrder,
            const Vector< index_t > & aNumNodes,
            const Vector< real >    & aStep,
            const Vector< real >      aOrigin = {} );

        ~TensorMeshConfig() = default;

        uint order() const
        {
           return mOrder ;
        }

        uint dimension() const
        {
           return mDimension ;
        }

        const Vector< index_t > & num_nodes_vector() const
        {
           return mNumNodesPerDim ;
        }

        const Vector< index_t > & num_elements_vector() const
        {
            return mNumElementsPerDim ;
        }

        const Vector< real > & step() const
        {
            return mNodeSteps ;
        }

        const Vector< real > & origin() const
        {
            return mOrigin ;
        }

        const Vector< real > & limit() const
        {
            return mLimit ;
        }

        uint
        num_elements( const uint aDimension = BELFEM_UINT_MAX ) const
        {
            switch ( aDimension )
            {
                case 0 : return mNumElementsPerDim( 0 ) ;
                case 1 : return mNumElementsPerDim( 1 ) ;
                case 2 :
                {
                    BELFEM_ERROR( mDimension == 3, "Invalid dimension : %u", ( unsigned int ) aDimension );
                    return mNumElementsPerDim( 2 ) ;
                }
                default : return mNumElements ;
            }
        }

        uint
        num_nodes( const uint aDimension = BELFEM_UINT_MAX ) const
        {
            switch ( aDimension )
            {
                case 0 : return mNumNodesPerDim( 0 ) ;
                case 1 : return mNumNodesPerDim( 1 ) ;
                case 2 :
                {
                    BELFEM_ERROR( mDimension == 3,
                        "Invalid dimension : %u", ( unsigned int ) aDimension );
                    return mNumNodesPerDim( 2 ) ;
                }
                default : return mNumNodes ;
            }
        }

        ElementType element_type() const
        {
            return mElementType ;
        }

        uint
        num_control_points( const uint aDimension = BELFEM_UINT_MAX ) const
        {
            switch ( aDimension )
            {
                case 0 : return mNumControlPointsPerDim( 0 ) ;
                case 1 : return mNumControlPointsPerDim( 1 ) ;
                case 2 :
                {
                    BELFEM_ERROR( mDimension == 3,
                        "Invalid dimension : %u", ( unsigned int ) aDimension );
                    return mNumControlPointsPerDim( 2 ) ;
                }
                default : return mNumControlPoints ;
            }
        }

        real
        min( const uint aDimension ) const
        {
            BELFEM_ASSERT( aDimension < mDimension, "Invalid dimension : %u", ( unsigned int ) aDimension );
            return mOrigin( aDimension );
        }

        real
        max( const uint aDimension ) const
        {
            BELFEM_ASSERT( aDimension < mDimension, "Invalid dimension : %u", ( unsigned int ) aDimension );
            return mLimit( aDimension );
        }

        index_t
        node_index( const index_t i, const index_t j ) const
        {
            BELFEM_ASSERT( mDimension == 2, "Invalid dimension : %u", ( unsigned int ) mDimension );

            BELFEM_ASSERT( i < mNumNodesPerDim( 0 ), "Invalid i : %u", ( unsigned int ) i );
            BELFEM_ASSERT( j < mNumNodesPerDim( 1 ), "Invalid j : %u", ( unsigned int ) j );
            return mNumNodesPerDim( 1 ) * i + j ;
        }

        index_t
        node_index( const index_t i, const index_t j, const index_t k ) const
        {
            BELFEM_ASSERT( mDimension == 3, "Invalid dimension : %u", ( unsigned int ) mDimension );
            BELFEM_ASSERT( i < mNumNodesPerDim( 0 ), "Invalid i : %u", ( unsigned int ) i );
            BELFEM_ASSERT( j < mNumNodesPerDim( 1 ), "Invalid j : %u", ( unsigned int ) j );
            BELFEM_ASSERT( k < mNumNodesPerDim( 2 ), "Invalid k : %u", ( unsigned int ) k );
            return mNumNodesPerDim( 2 ) * ( mNumNodesPerDim( 1 )*i + j ) + k ;
        }

        index_t
        element_index( const index_t i, const index_t j ) const
        {
            BELFEM_ASSERT( mDimension == 2, "Invalid dimension : %u", ( unsigned int ) mDimension );
            BELFEM_ASSERT( i < mNumElementsPerDim( 0 ), "Invalid i : %u", ( unsigned int ) i );
            BELFEM_ASSERT( j < mNumElementsPerDim( 1 ), "Invalid j : %u", ( unsigned int ) j );
            return mNumElementsPerDim( 1 ) * i + j ;
        }

        index_t
        element_index( const index_t i, const index_t j, const index_t k ) const
        {
            BELFEM_ASSERT( mDimension == 3, "Invalid dimension : %u", ( unsigned int ) mDimension );
            BELFEM_ASSERT( i < mNumElementsPerDim( 0 ), "Invalid i : %u", ( unsigned int ) i );
            BELFEM_ASSERT( j < mNumElementsPerDim( 1 ), "Invalid j : %u", ( unsigned int ) j );
            BELFEM_ASSERT( k < mNumElementsPerDim( 2 ), "Invalid k : %u", ( unsigned int ) k );
            return mNumElementsPerDim( 2 ) * ( mNumElementsPerDim( 1 ) * i + j ) + k ;
        }

        real
        element_step( const uint aDimension ) const
        {
            BELFEM_ASSERT( aDimension < mDimension, "Invalid dimension : %u", ( unsigned int ) aDimension );
            return mElementSteps( aDimension );
        }

        real
        inv_element_step( const uint aDimension ) const
        {
            BELFEM_ASSERT( aDimension < mDimension, "Invalid dimension : %u", ( unsigned int ) aDimension );
            return mInvElementSteps( aDimension );
        }

        index_t
        element_ijk( const index_t i, const real x )
        {
            BELFEM_ASSERT( i < mDimension, "Invalid dimension : %u", ( unsigned int ) i );

            return std::clamp( ( index_t ) ( ( x - mOrigin( i ) ) * mInvElementSteps( i ) ),
                ( index_t ) 0, mNumElementsPerDim( i ) -1);
        }
    private:

        Vector< index_t >
        compute_num_elements_per_dim() const ;

        Vector< index_t >
        compute_num_control_points_per_dim() const ;

        Vector< real >
        compute_limit() const ;

        index_t
        compute_num_nodes() const ;

        index_t
        compute_num_elements() const ;

        index_t
        compute_num_control_point() const ;

        ElementType
        determine_element_type() const ;

        Vector< real >
        compute_element_steps() const ;

        Vector< real >
        compute_inv_element_steps() const ;
    };

}
#endif //BELFEM_CL_TENSORMESHCONFIG_HPP