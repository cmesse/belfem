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

#include "cl_Gradient.hpp"
#include "../../linalg/lapack/fn_gesv.hpp"
#include "fn_inv.hpp"
namespace belfem
{
    namespace fem
    {
        Gradient::Gradient(
            Kernel * aKernel,
            const Vector< id_t > & aBlocksIDs,
            const bool aFlipSign,
            const uint aDofManagerIndex ) :
                Postprocessor( aKernel ),
                mNumDimensions( aKernel->mesh()->number_of_dimensions() ),
                mFlipSign( aFlipSign ),
                mDofMaganerIndex( aDofManagerIndex )
        {
            mOrder = 0 ;

            mMesh->unflag_all_nodes();

            for ( id_t tID : aBlocksIDs )
            {
                mMesh->block( tID )->flag_nodes();
                mBlocks[ tID ] = aKernel->dofmgr( aDofManagerIndex )->block( tID );

                ElementType tType = mBlocks( tID )->element_type() ;

                uint tOrder = 0 ;

                if ( mesh::geometry_type( tType ) == GeometryType::TRI ||
                     mesh::geometry_type( tType ) == GeometryType::TET )
                {
                    tOrder = mesh::interpolation_order_numeric( tType ) ;
                }
                else
                {
                    tOrder = 2 * mesh::interpolation_order_numeric( tType );
                }
                if ( tOrder > mOrder )
                {
                    mOrder = tOrder;
                }
            }
            mX.set_size( 1, mNumDimensions );

            this->set_order( mOrder );

            mJ.set_size( mNumDimensions, mNumDimensions );

            // collect my nodes
            index_t tCount = 0 ;
            for ( mesh::Node * tNode : mMesh->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    tCount++ ;
                }
            }
            BELFEM_ERROR( tCount > 0, "No nodes selected for gradient computation" );

            mNodes.set_size( tCount, nullptr );
            tCount = 0 ;
            for ( mesh::Node * tNode : mMesh->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    mNodes( tCount++ ) = tNode ;
                }
            }

        }

        void
        Gradient::set_fields( const string & aScalarField, const string & aGradientField )
        {
            mScalarField   = aScalarField;

            BELFEM_ASSERT( mMesh->field_exists( aScalarField ), "Field %s does not exist on mesh", aScalarField.c_str() );

            mGradientFields.push( aGradientField + "x" );
            mGradientFields.push( aGradientField + "y" );
            mGradientFields.push( aGradientField + "z" );

            for ( const string & tGradientField : mGradientFields )
            {
                if ( ! mMesh->field_exists( tGradientField ) )
                {
                    mMesh->create_field( tGradientField );
                }
            }

        }

        void
        Gradient::process_node( mesh::Node * aNode )
        {
            if ( aNode->is_duplicate() ) return ;

            if ( aNode->owner() == mCommRank )
            {
                mV.fill( 0.0 );

                mC.fill( 0.0 );


                uint tElemCount = this->subprocess_node( aNode );
                for ( uint k=0; k<aNode->number_of_duplicates(); ++k )
                {
                    tElemCount += this->subprocess_node( aNode->duplicate( k ) );
                }

                if ( tElemCount == 0 ) return;

                gesv( mV, mC, mPivot );


                for ( int i=0; i<mNumDimensions; ++i )
                {
                    mX( 0, i ) = aNode->x( i );
                }

                this->compute_poly( mX );


                // evaluate field data
                for ( int i=0; i<mNumDimensions; ++i )
                {
                    Vector< real > & tGrad = mMesh->field_data( mGradientFields( i ) );

                    real tVal = dot( mP.col( 0 ), mC.col( i ) );

                    tGrad( aNode->index() ) = tVal ;

                    for ( uint k=0; k<aNode->number_of_duplicates(); ++k )
                    {
                        tGrad( aNode->duplicate( k )->index() ) = tVal ;
                    }
                }
            }
        }

        uint
        Gradient::subprocess_node( mesh::Node * aNode )
        {
            uint aElemCount = 0 ;

            const Vector< real > & tField = mMesh->field_data( mScalarField );

            for ( uint e=0; e<aNode->number_of_elements(); ++e )
            {
                if ( ! aNode->element( e )->is_flagged() ) continue ;

                ++aElemCount ;

                Block * tBlock = mBlocks( aNode->element( e )->block_id() ) ;

                const IntegrationData * tIntegration = tBlock->integration() ;


                mesh::Element * tElement = aNode->element( e ) ;

                mPhi.set_size( tElement->number_of_nodes() );
                for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                {
                    mPhi( k ) = tField( tElement->node( k )->index() ) ;
                }

                // get element coordinates
                mElX.set_size( tElement->number_of_nodes() , mNumDimensions );
                for ( int j=0; j<mNumDimensions; ++j )
                {
                    for ( uint i=0; i<tElement->number_of_nodes(); ++i )
                    {
                        mElX( i, j ) = tElement->node( i )->x( j );
                    }
                }

                for ( uint k=0 ; k<tIntegration->number_of_integration_points(); ++k )
                {
                    mX = tIntegration->N( k ) * mElX ;

                    mJ = tIntegration->dNdXi( k ) * mElX ;

                    mG = inv( mJ.matrix_data() ) * tIntegration->dNdXi( k ).matrix_data() * mPhi.vector_data();

                    this->compute_poly( mX );

                    mV += mP * trans( mP );

                    for ( int j=0; j<mNumDimensions; ++j )
                    {
                        for ( int i=0; i<mN; ++i )
                        {
                            mC( i, j ) += mP( i, 0 ) * mG( j );
                        }
                    }
                }
            }

            return aElemCount ;
        }

        void
        Gradient::set_order( const uint aOrder )
        {
            uint n2D = 0 ;
            uint n3D = 0 ;

            switch ( aOrder )
            {
                case 1 :
                {
                    mFunPoly2D = & Gradient::poly1_2d ;
                    mFunPoly3D = & Gradient::poly1_3d ;
                    n2D = 3 ;
                    n3D = 4 ;
                    break ;
                }
                case 2 :
                {
                    mFunPoly2D = & Gradient::poly2_2d ;
                    mFunPoly3D = & Gradient::poly2_3d ;
                    n2D = 6 ;
                    n3D = 10 ;
                    break ;
                }
                case 3 :
                {
                    mFunPoly2D = & Gradient::poly3_2d ;
                    mFunPoly3D = & Gradient::poly3_3d ;
                    n2D = 10 ;
                    n3D = 20 ;
                    break ;
                }
                case 4 :
                {
                    mFunPoly2D = & Gradient::poly4_2d ;
                    mFunPoly3D = & Gradient::poly4_3d ;
                    n2D = 15 ;
                    n3D = 35 ;
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Invalid order: %u", ( unsigned int ) aOrder );
                }
            }

            mN = mNumDimensions == 2 ? n2D : n3D ;

            if ( mNumDimensions == 2 )
            {
                mFunComputePoly = & Gradient::compute_poly_2d ;
            }
            else
            {
                mFunComputePoly = & Gradient::compute_poly_3d ;
            }
            mP.set_size( mN, 1 );
            mC.set_size( mN, mNumDimensions );
            mPivot.set_size( mN );
            mV.set_size( mN, mN );
        }

        void
        Gradient::run()
        {
            BELFEM_ASSERT( mScalarField.length() > 0, "Fields not defined for gradient processor");

            // make sure that we only use the elements we want to
            mMesh->unflag_all_elements() ;
            for ( auto tPair : mBlocks )
            {
                if ( mMesh->block_exists( tPair.first ) )
                {
                    mMesh->block( tPair.first )->flag_elements() ;
                }
            }

            mKernel->dofmgr( mDofMaganerIndex )->distribute_fields( {mScalarField } );

            for ( mesh::Node * tNode : mNodes )
            {
                this->process_node( tNode );
            }

            if ( mFlipSign )
            {
                for ( const string & tGradientField : mGradientFields )
                {
                    mMesh->field_data( tGradientField ) *= -1.0 ;
                }
            }

            // synch data
            mKernel->dofmgr( mDofMaganerIndex )->synchronize_fields( mGradientFields );

        }

    }
}
