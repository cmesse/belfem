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

#include "cl_TensorMeshConfig.hpp"
#include "assert.hpp"

namespace belfem
{
    TensorMeshConfig::TensorMeshConfig(
        const uint aOrder,
        const Vector< index_t > & aNumNodes,
        const Vector< real >    & aStep,
        const Vector< real >   aOrigin ) :
        mOrder( aOrder ),
        mDimension( aNumNodes.length() ),
        mNumNodesPerDim( aNumNodes ),
        mNodeSteps( aStep ),
        mOrigin( aOrigin.length() == 0 ?
            Vector< real >( aNumNodes.length(), 0.0 ) : aOrigin ),
        mNumElementsPerDim( this->compute_num_elements_per_dim() ),
        mNumControlPointsPerDim( this->compute_num_control_points_per_dim() ),
        mLimit( this->compute_limit() ),
        mElementSteps( this->compute_element_steps() ),
        mInvElementSteps( this->compute_inv_element_steps() ),
        mNumNodes( this->compute_num_nodes() ),
        mNumElements( this->compute_num_elements() ),
        mNumControlPoints( this->compute_num_control_point() ),
        mElementType( this->determine_element_type() )
    {

    }

    Vector< index_t >
    TensorMeshConfig::compute_num_elements_per_dim() const
    {
        Vector< index_t > aNumElements( mDimension, 0 );

        for ( uint i=0; i<mDimension; ++i )
        {
            BELFEM_ERROR( ( mNumNodesPerDim( i ) - 1 ) % mOrder == 0,
                "Invalid number of nodes in direction %u for order %u : %u",
                 ( unsigned int ) i, ( unsigned int ) mOrder , ( unsigned int ) mNumNodesPerDim( i ) );

            aNumElements( i ) = ( mNumNodesPerDim( i ) - 1) / mOrder ;

        }
        return aNumElements;
    }

    Vector< real >
    TensorMeshConfig::compute_limit() const
    {
        Vector< real > aLimit( mDimension, 0 );

        BELFEM_ERROR( mNodeSteps.length() == mDimension, "Invalid step dimension ( is %u, expect %u )",
            ( unsigned int ) mNodeSteps.length(), ( unsigned int ) mDimension );

        for ( uint i=0; i<mDimension; ++i )
        {
            aLimit( i ) = mOrigin( i ) + mOrder * mNumElementsPerDim( i ) * mNodeSteps( i );
        }
        return aLimit;
    }

    ElementType
    TensorMeshConfig::determine_element_type() const
    {
        switch ( mDimension )
        {
            case 2 :
            {
                switch ( mOrder )
                {
                    case 1 : return ElementType::QUAD4 ;
                    case 2 : return ElementType::QUAD9 ;
                    case 3 : return ElementType::QUAD16 ;
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid element order : %u", ( unsigned int ) mOrder );
                        return ElementType::UNDEFINED ;
                    }
                }
            }
            case 3 :
            {
                switch ( mOrder )
                {
                    case 1 : return ElementType::HEX8 ;
                    case 2 : return ElementType::HEX27 ;
                    case 3 : return ElementType::HEX64 ;
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid element order : %u", ( unsigned int ) mOrder );
                        return ElementType::UNDEFINED ;
                    }
                }
            }
            default :
            {
                BELFEM_ERROR( false, "Invalid mesh dimension : %u", ( unsigned int ) mDimension );
                return ElementType::UNDEFINED ;
            }
        }
    }

    index_t
    TensorMeshConfig::compute_num_nodes() const
    {
        index_t aCount = 1 ;
        for ( index_t i=0; i<mDimension; ++i )
        {
            aCount *= mNumNodesPerDim( i );
        }
        return aCount;
    }

    index_t
    TensorMeshConfig::compute_num_elements() const
    {
        index_t aCount = 1 ;
        for ( index_t i=0; i<mDimension; ++i )
        {
            aCount *= mNumElementsPerDim( i );
        }
        return aCount;
    }

    Vector< index_t >
    TensorMeshConfig::compute_num_control_points_per_dim() const
    {
        Vector< index_t > aNumControlPoints( mDimension, 0 );

        for ( uint i=0; i<mDimension; ++i )
        {

            aNumControlPoints( i ) = mNumElementsPerDim( i ) + mOrder ;

        }
        return aNumControlPoints;
    }

    index_t
    TensorMeshConfig::compute_num_control_point() const
    {
        index_t aCount = 1 ;
        for ( index_t i=0; i<mDimension; ++i )
        {
            aCount *= mNumControlPointsPerDim( i );
        }
        return aCount;
    }


    Vector< real >
    TensorMeshConfig::compute_element_steps() const
    {
        Vector< real > aSteps( mDimension, 0 );
        for ( uint32_t i=0; i<mDimension; ++i )
        {
            aSteps( i ) = mNodeSteps( i ) * mOrder ;
        }
        return aSteps;
    }

    Vector< real >
    TensorMeshConfig::compute_inv_element_steps() const
    {
        Vector< real > aInvSteps( mDimension, 0 );
        for ( uint32_t i=0; i<mDimension; ++i )
        {
            aInvSteps( i ) = 1.0 / mElementSteps( i );
        }
        return aInvSteps;
    }
}