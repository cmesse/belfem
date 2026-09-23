/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 *
 * The history table of the calculator: Calculator::qold( s ) must return,
 * for every local dof, the value of the history field "<label><s>" of that
 * dof's own field. The table is indexed s * stride + field index; with a
 * stride of one less than the largest field index plus one, the slot of
 * ( s, last field ) and ( s + 1, first field ) coincide. Two dof fields at
 * mesh field indices 0 and 1 make that collision observable.
 */

#include <gtest/gtest.h>
#include <memory>
#include <string>

#include "typedefs.hpp"
#include "commtools.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "cl_Mesh_Field.hpp"
#include "cl_Node.hpp"
#include "cl_Block.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_KernelParameters.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Calculator.hpp"
#include "cl_FEM_Dof.hpp"
#include "cl_IWG_Timestep.hpp"
#include "en_IWGs.hpp"
#include "en_SolverEnums.hpp"

using namespace belfem;

namespace
{
    //! a transient IWG with two nodal dof fields; nothing is assembled
    class TwoFieldIwg : public fem::IWG_Timestep
    {
        public:
            TwoFieldIwg() :
                fem::IWG_Timestep( IwgType::Poisson, ModelDimensionality::TwoD, IwgMode::Direct )
            {
                mNumberOfDofsPerNode = 2 ;
                mNumberOfSpatialDimensions = 2 ;
                mNumberOfDerivativeDimensions = 2 ;
                mDofFields = { "phi", "psi" };

                // as IWG_Poisson does: the dof tables must exist before
                // select_blocks() assigns them
                this->initialize();
            }
    };

    //! the value planted in history field "<label><s>" at a node
    real
    sentinel( const uint aStep, const uint aField, const index_t aNode )
    {
        return 1000.0 * ( aStep + 1 ) + 100.0 * aField + static_cast< real >( aNode ) ;
    }

    const char * gLabels[ 2 ] = { "phi", "psi" };

    //! unit square as two TRI3, one block, both dof fields created on the
    //! mesh before anything else so they take field indices 0 and 1
    class HistoryProblem
    {
        public:

            std::unique_ptr< Mesh >                   mMesh ;
            std::unique_ptr< fem::KernelParameters >  mParams ;
            std::unique_ptr< fem::Kernel >            mKernel ;
            TwoFieldIwg            * mIWG ;     // owned by the kernel
            fem::DofManager        * mField ;   // owned by the kernel

            explicit
            HistoryProblem( const EulerMethod aMethod )
            {
                // third argument must stay false: the MeshChecker in the
                // Kernel ctor refuses pre-existing edges
                mMesh.reset( new Mesh( 2, 0, false ) );

                mesh::ElementFactory tFactory ;

                Cell< mesh::Node * > & tNodes = mMesh->nodes();
                tNodes.push( new mesh::Node( 1, 0.0, 0.0 ) );
                tNodes.push( new mesh::Node( 2, 1.0, 0.0 ) );
                tNodes.push( new mesh::Node( 3, 1.0, 1.0 ) );
                tNodes.push( new mesh::Node( 4, 0.0, 1.0 ) );

                mesh::Block * tBlock = new mesh::Block( 1, 2 );
                const uint tTriangles[ 2 ][ 3 ] = { { 0, 1, 2 }, { 0, 2, 3 } };
                for ( uint e = 0; e < 2; ++e )
                {
                    mesh::Element * tTri =
                            tFactory.create_element( ElementType::TRI3, e + 1 );
                    for ( uint k = 0; k < 3; ++k )
                    {
                        tTri->insert_node( tNodes( tTriangles[ e ][ k ] ), k );
                    }
                    tTri->set_block_id( 1 );
                    mMesh->elements().push( tTri );
                    tBlock->insert_element( tTri );
                }
                mMesh->blocks().push( tBlock );
                mMesh->finalize();

                mMesh->create_field( gLabels[ 0 ] );
                mMesh->create_field( gLabels[ 1 ] );

                mParams.reset( new fem::KernelParameters( mMesh.get() ) );
                mKernel.reset( new fem::Kernel( mParams.get() ) );

                mIWG = new TwoFieldIwg();
                mIWG->select_blocks( { 1 } );
                mKernel->add_equation( mIWG );
                mIWG->set_timestepping_method( aMethod );

                mField = mKernel->create_field( mIWG );
                mField->set_solver( gDefaultSolver );
            }

            HistoryProblem( const HistoryProblem & ) = delete ;
            HistoryProblem( HistoryProblem && ) = delete ;
            HistoryProblem & operator=( const HistoryProblem & ) = delete ;
            HistoryProblem & operator=( HistoryProblem && ) = delete ;

            //! plants the sentinels into every history field up to aOrder
            void
            seed( const uint aOrder )
            {
                for ( uint s = 0; s < aOrder; ++s )
                {
                    for ( uint f = 0; f < 2; ++f )
                    {
                        Vector< real > & tData = mMesh->field_data(
                                std::string( gLabels[ f ] ) + std::to_string( s ) );
                        for ( index_t i = 0; i < tData.length(); ++i )
                        {
                            tData( i ) = sentinel( s, f, i );
                        }
                    }
                }
            }

            //! qold( s ) on every element of the block against the sentinels
            void
            expect_history( const uint aOrder )
            {
                fem::Block * tBlock = mField->block( 1 );
                fem::Calculator * tCalc = tBlock->calculator();
                ASSERT_NE( tCalc, nullptr );

                const index_t tPhiIndex = mMesh->field( gLabels[ 0 ] )->index() ;
                const index_t tPsiIndex = mMesh->field( gLabels[ 1 ] )->index() ;

                // an empty walk would pass vacuously
                ASSERT_EQ( tBlock->elements().size(), 2u );

                for ( fem::Element * tElement : tBlock->elements() )
                {
                    tCalc->link( tElement );

                    for ( uint s = 0; s < aOrder; ++s )
                    {
                        const Vector< real > & tQ = tCalc->qold( s );
                        ASSERT_EQ( tElement->number_of_local_dofs(), 6u );
                        ASSERT_EQ( tQ.length(), 6u );

                        for ( uint k = 0; k < tElement->number_of_local_dofs(); ++k )
                        {
                            fem::Dof * tDof = tElement->local_dof( k );
                            const index_t tFieldIndex = tDof->field_index() ;
                            ASSERT_TRUE( tFieldIndex == tPhiIndex || tFieldIndex == tPsiIndex );

                            const uint f = tFieldIndex == tPhiIndex ? 0 : 1 ;

                            EXPECT_EQ( tQ( k ), sentinel( s, f, tDof->dof_index_on_field() ) )
                                << "element " << tElement->id() << ", step " << s
                                << ", local dof " << k << " ( field " << gLabels[ f ] << " )" ;
                        }
                    }
                }
            }
    };
}

//------------------------------------------------------------------------------

TEST( QoldTable, EveryStepAndFieldReadsItsOwnSlot )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    HistoryProblem tProblem( EulerMethod::BackwardDifference2 );

    // the colliding pair exists: the two dof fields are 0 and 1
    ASSERT_EQ( tProblem.mMesh->field( "phi" )->index(), 0u );
    ASSERT_EQ( tProblem.mMesh->field( "psi" )->index(), 1u );

    tProblem.seed( 2 );
    tProblem.mField->initialize();

    // ( step 1, field phi ) must not read ( step 0, field psi )
    tProblem.expect_history( 2 );
}

TEST( QoldTable, OrderRaisedAfterInitializationReadsTheNewLevel )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    HistoryProblem tProblem( EulerMethod::BackwardDifference2 );
    tProblem.mField->initialize();

    // the order change rebuilds the history fields and the calculator's
    // table; the third level must be readable, and read its own values
    tProblem.mIWG->set_timestepping_method( EulerMethod::BackwardDifference3 );

    ASSERT_TRUE( tProblem.mMesh->field_exists( "phi2" ) );
    ASSERT_TRUE( tProblem.mMesh->field_exists( "psi2" ) );

    tProblem.seed( 3 );
    tProblem.expect_history( 3 );
}

TEST( QoldTable, HistoryFieldsAreNotWrittenToFile )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    HistoryProblem tProblem( EulerMethod::BackwardDifference2 );
    tProblem.mField->initialize();

    for ( const char * tLabel : { "phi0", "phi1", "psi0", "psi1" } )
    {
        EXPECT_FALSE( tProblem.mMesh->field( tLabel )->write_field_to_file() ) << tLabel ;
    }
    for ( const char * tLabel : gLabels )
    {
        EXPECT_TRUE( tProblem.mMesh->field( tLabel )->write_field_to_file() ) << tLabel ;
    }
}
