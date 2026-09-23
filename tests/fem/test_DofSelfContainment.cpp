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
 * The flattening invariant of hanging dofs, as a property of the built dof
 * set: a hanging dof never lists itself among its sources, and no source
 * hangs in turn. Dof::set_sources only asserts the second half, and only in
 * an assertions-active build; at that moment the dof being set does not
 * hang yet, so a self-source would pass the assert. The walk below runs in
 * every build.
 */

#include <gtest/gtest.h>
#include <memory>

#include "typedefs.hpp"
#include "commtools.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_Block.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_KernelParameters.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_Dof.hpp"
#include "en_IWGs.hpp"
#include "en_SolverEnums.hpp"

using namespace belfem;

namespace
{
    //! two unit squares that do NOT share their interface nodes; the right
    //! square's interface nodes hang on the left square's, one by a single
    //! unit source and one by two half-weight sources
    //!
    //!   4---3  6---8      node 5 ( 1, 0 ) hangs on node 2
    //!   | L |  | R |      node 6 ( 1, 1 ) hangs on nodes 3 and 2, 0.5 each
    //!   1---2  5---7      TRI3 {1,2,3} {1,3,4} {5,7,8} {5,8,6}, all CCW
    class HangingStrip
    {
        public:

            // owners in construction order: a throw in the constructor
            // destroys what was built, in reverse order
            std::unique_ptr< Mesh >                   mMesh ;
            std::unique_ptr< fem::KernelParameters >  mParams ;
            std::unique_ptr< fem::Kernel >            mKernel ;
            fem::DofManager        * mField ;

            HangingStrip()
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
                tNodes.push( new mesh::Node( 5, 1.0, 0.0 ) );
                tNodes.push( new mesh::Node( 6, 1.0, 1.0 ) );
                tNodes.push( new mesh::Node( 7, 2.0, 0.0 ) );
                tNodes.push( new mesh::Node( 8, 2.0, 1.0 ) );

                // the ties must exist before the Kernel is built: its ctor
                // collects the hanging bases from the mesh
                {
                    Cell< mesh::Node * > tSources = { tNodes( 1 ) };
                    Vector< real >       tWeights( 1, 1.0 );
                    tNodes( 4 )->set_sources( tSources, tWeights );
                }
                {
                    Cell< mesh::Node * > tSources = { tNodes( 2 ), tNodes( 1 ) };
                    Vector< real >       tWeights( 2, 0.5 );
                    tNodes( 5 )->set_sources( tSources, tWeights );
                }

                Cell< mesh::Element * > & tElements = mMesh->elements();
                mesh::Block * tBlock = new mesh::Block( 1, 4 );

                const uint tTriangles[ 4 ][ 3 ] = {
                        { 0, 1, 2 },
                        { 0, 2, 3 },
                        { 4, 6, 7 },
                        { 4, 7, 5 } };

                for ( uint e = 0; e < 4; ++e )
                {
                    mesh::Element * tTri =
                            tFactory.create_element( ElementType::TRI3, e + 1 );
                    for ( uint k = 0; k < 3; ++k )
                    {
                        tTri->insert_node( tNodes( tTriangles[ e ][ k ] ), k );
                    }
                    tTri->set_block_id( 1 );
                    tElements.push( tTri );
                    tBlock->insert_element( tTri );
                }

                mMesh->blocks().push( tBlock );
                mMesh->finalize();

                mParams.reset( new fem::KernelParameters( mMesh.get() ) );
                mKernel.reset( new fem::Kernel( mParams.get() ) );

                fem::IWG * tIWG = mKernel->create_equation(
                        IwgType::Poisson,
                        ModelDimensionality::TwoD,
                        { 1 } );

                mField = mKernel->create_field( tIWG );
                mField->set_solver( gDefaultSolver );
            }

            HangingStrip( const HangingStrip & ) = delete ;
            HangingStrip( HangingStrip && ) = delete ;
            HangingStrip & operator=( const HangingStrip & ) = delete ;
            HangingStrip & operator=( HangingStrip && ) = delete ;

            //! the dof of every mesh node, hanging or not
            Cell< fem::Dof * >
            node_dofs()
            {
                Cell< fem::Dof * > tDofs ;
                for ( mesh::Node * tNode : mMesh->nodes() )
                {
                    tDofs.push( mField->dof( mField->calculate_dof_id( tNode, 0 ) ) );
                }
                return tDofs ;
            }
    };

    //! walks the dof set; aPhase names the moment for the failure message
    void
    expect_flat_and_self_free( HangingStrip & aStrip, const char * aPhase,
                               uint & aNumHanging )
    {
        aNumHanging = 0 ;

        for ( fem::Dof * tDof : aStrip.node_dofs() )
        {
            ASSERT_NE( tDof, nullptr ) << aPhase ;

            if ( ! tDof->is_hanging() ) continue ;
            ++aNumHanging ;

            for ( uint k = 0; k < tDof->number_of_sources(); ++k )
            {
                fem::Dof * tSource = tDof->source( k );

                ASSERT_NE( tSource, nullptr ) << aPhase << ": dof " << tDof->id() ;
                EXPECT_NE( tSource, tDof )
                    << aPhase << ": dof " << tDof->id() << " is its own source " << k ;
                EXPECT_FALSE( tSource->is_hanging() )
                    << aPhase << ": source " << tSource->id() << " of dof "
                    << tDof->id() << " hangs itself" ;
            }
        }
    }
}

//------------------------------------------------------------------------------

TEST( DofSelfContainment, NoHangingDofContainsItselfOrAHangingSource )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    HangingStrip tStrip ;

    // the property must hold where the sources are resolved ( inside
    // create_field ) and after the dof manager is frozen
    uint tNumHanging = 0 ;
    expect_flat_and_self_free( tStrip, "after create_field", tNumHanging );
    EXPECT_EQ( tNumHanging, 2u ) << "the fixture must produce both ties" ;

    tStrip.mField->initialize();

    expect_flat_and_self_free( tStrip, "after initialize", tNumHanging );
    EXPECT_EQ( tNumHanging, 2u );
}

TEST( DofSelfContainment, SolvedDofContainerHoldsNoHangingDof )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    HangingStrip tStrip ;
    tStrip.mField->initialize();

    // the container the solver numbers holds only independent dofs, and
    // every source of a hanging dof is one of them
    Cell< fem::Dof * > & tDofs = tStrip.mField->dofs();
    EXPECT_EQ( tDofs.size(), 6u );

    for ( fem::Dof * tDof : tDofs )
    {
        EXPECT_FALSE( tDof->is_hanging() ) << "dof " << tDof->id() ;
    }

    for ( fem::Dof * tDof : tStrip.node_dofs() )
    {
        for ( uint k = 0; k < tDof->number_of_sources(); ++k )
        {
            bool tFound = false ;
            for ( fem::Dof * tIndependent : tDofs )
            {
                tFound = tFound || tIndependent == tDof->source( k );
            }
            EXPECT_TRUE( tFound ) << "source " << k << " of dof " << tDof->id()
                                  << " is not an independent dof" ;
        }
    }
}
