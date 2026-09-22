/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any
 * required approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include <gtest/gtest.h>
#include <memory>

#include <cmath>
#include <stdexcept>
#include <string>

#include "typedefs.hpp"
#include "commtools.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_Facet.hpp"
#include "cl_Block.hpp"
#include "cl_SideSet.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_KernelParameters.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_SideSet.hpp"
#include "cl_FEM_Dof.hpp"
#include "en_IWGs.hpp"
#include "en_SolverEnums.hpp"

using namespace belfem;

namespace
{
    // the contract under test: a Dirichlet condition imposed on a sideset
    // node whose dof is condensed onto ONE unit-weight source lands on that
    // SOURCE dof ( the alias stays unfixed and reconstructs through the
    // T-matrix ); a non-unit single-source weight and a first free-to-fixed
    // flip after initialize() are refused in every build

    constexpr real tDirichletValue = 7.5 ;

    // runs aCall, requires a std::runtime_error whose what() carries the
    // literal part of the BELFEM_ERROR text
    template< typename F >
    void
    expect_refusal( F && aCall, const char * aSubstring )
    {
        bool tThrown = false ;
        try
        {
            aCall();
        }
        catch( const std::runtime_error & aError )
        {
            tThrown = true ;
            EXPECT_NE( std::string( aError.what() ).find( aSubstring ), std::string::npos )
                << "message was: " << aError.what();
        }
        EXPECT_TRUE( tThrown ) << "expected a refusal carrying \"" << aSubstring << "\"";
    }

    //! two unit squares that do NOT share their interface nodes; the right
    //! square's interface nodes hang on the left square's by single sources,
    //! the shape of a periodic pair
    //!
    //!   4---3  6---8      node 5 ( 1, 0 ) hangs on node 2 ( 1, 0 )
    //!   | L |  | R |      node 6 ( 1, 1 ) hangs on node 3 ( 1, 1 )
    //!   1---2  5---7      TRI3 {1,2,3} {1,3,4} {5,7,8} {5,8,6}, all CCW
    //!
    //! sideset 10 is the top edge of R ( 8 -> 6 ): its node 6 is condensed
    //! onto node 3, which lies OUTSIDE the sideset. Sideset 11 is the left
    //! edge of L ( 4 -> 1 ). Nothing is imposed by the fixture itself.
    class GluedStrip
    {
        public:

            // owners in construction order: a throw in the constructor
            // destroys what was built, in reverse order
            std::unique_ptr< Mesh >                   mMesh ;
            std::unique_ptr< fem::KernelParameters >  mParams ;
            std::unique_ptr< fem::Kernel >            mKernel ;
            fem::DofManager        * mField ;

//------------------------------------------------------------------------------

            //! aWeight: weight of node 6 on node 3 ( node 5 always hangs on
            //! node 2 with unit weight ). With aTieToBothCorners, node 6
            //! hangs on nodes 3 and 2 with weights 0.5 each and aWeight is
            //! ignored
            explicit
            GluedStrip( const real aWeight = 1.0,
                        const bool aTieToBothCorners = false )
            {
                // third argument must stay false: connectivity computation
                // would create edges in finalize() and the MeshChecker in
                // the Kernel ctor refuses pre-existing edges
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

                if ( aTieToBothCorners )
                {
                    Cell< mesh::Node * > tSources = { tNodes( 2 ), tNodes( 1 ) };
                    Vector< real >       tWeights( 2, 0.5 );
                    tNodes( 5 )->set_sources( tSources, tWeights );
                }
                else
                {
                    Cell< mesh::Node * > tSources = { tNodes( 2 ) };
                    Vector< real >       tWeights( 1, aWeight );
                    tNodes( 5 )->set_sources( tSources, tWeights );
                }

                Cell< mesh::Element * > & tElements = mMesh->elements();

                mesh::Block * tBlock = new mesh::Block( 1, 4 );

                // node triples are indices into tNodes
                const uint tTriangles[ 4 ][ 3 ] = {
                        { 0, 1, 2 },     // ids 1, 2, 3
                        { 0, 2, 3 },     // ids 1, 3, 4 — carries the left edge
                        { 4, 6, 7 },     // ids 5, 7, 8
                        { 4, 7, 5 } };   // ids 5, 8, 6 — carries the top edge

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

                // the LINE2 of each facet is owned by the Facet and must
                // not be pushed into elements(); set_master relinks the
                // facet nodes to the canonical face order

                // sideset 10: the top edge of R ( ids 8 -> 6 ) is local
                // edge 1 ( nodes 1 -> 2 ) of triangle 4
                {
                    mesh::Element * tLine =
                            tFactory.create_element( ElementType::LINE2, 101 );
                    tLine->insert_node( tNodes( 7 ), 0 );
                    tLine->insert_node( tNodes( 5 ), 1 );

                    mesh::Facet * tFacet = new mesh::Facet( tLine );
                    tFacet->set_master( tElements( 3 ), 1 );

                    mesh::SideSet * tSideSet = new mesh::SideSet( 10, 1 );
                    tSideSet->insert_facet( tFacet );
                    mMesh->sidesets().push( tSideSet );
                }

                // sideset 11: the left edge of L ( ids 4 -> 1 ) is local
                // edge 2 ( nodes 2 -> 0 ) of triangle 2
                {
                    mesh::Element * tLine =
                            tFactory.create_element( ElementType::LINE2, 102 );
                    tLine->insert_node( tNodes( 3 ), 0 );
                    tLine->insert_node( tNodes( 0 ), 1 );

                    mesh::Facet * tFacet = new mesh::Facet( tLine );
                    tFacet->set_master( tElements( 1 ), 2 );

                    mesh::SideSet * tSideSet = new mesh::SideSet( 11, 1 );
                    tSideSet->insert_facet( tFacet );
                    mMesh->sidesets().push( tSideSet );
                }

                mMesh->finalize();

                mParams.reset( new fem::KernelParameters( mMesh.get() ) );
                mKernel.reset( new fem::Kernel( mParams.get() ) );

                fem::IWG * tIWG = mKernel->create_equation(
                        IwgType::Poisson,
                        ModelDimensionality::TwoD,
                        { 1 },
                        { 10, 11 } );

                // set_equation inside create_field resolves the Dof-level
                // sources, so the alias dof already hangs before initialize()
                mField = mKernel->create_field( tIWG );
                mField->set_solver( gDefaultSolver );
            }

//------------------------------------------------------------------------------

            GluedStrip( const GluedStrip & ) = delete ;
            GluedStrip( GluedStrip && ) = delete ;
            GluedStrip & operator=( const GluedStrip & ) = delete ;
            GluedStrip & operator=( GluedStrip && ) = delete ;

//------------------------------------------------------------------------------

            fem::Dof *
            node_dof( const id_t aNodeID )
            {
                return mField->dof( mField->calculate_dof_id(
                        mMesh->node( aNodeID ), 0 ) );
            }

//------------------------------------------------------------------------------

            void
            impose_top( const real aValue )
            {
                mField->sideset( 10 )->impose_dirichlet( aValue );
            }
    };
}

//------------------------------------------------------------------------------

TEST( DirichletPin, SingleSourcePinFixesTheSource )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    GluedStrip tStrip ;

    fem::Dof * tAlias  = tStrip.node_dof( 6 );
    fem::Dof * tSource = tStrip.node_dof( 3 );

    ASSERT_TRUE( tAlias->is_hanging() );
    ASSERT_EQ( tAlias->number_of_sources(), 1u );
    ASSERT_EQ( tAlias->source( 0 ), tSource );
    EXPECT_EQ( tAlias->weight( 0 ), 1.0 );
    EXPECT_FALSE( tSource->is_fixed() );

    tStrip.impose_top( tDirichletValue );

    // the condition lands on the independent dof; the stored weight is
    // exactly 1.0, so the division lands exactly on the value
    EXPECT_TRUE( tSource->is_fixed() );
    EXPECT_EQ( tSource->value(), tDirichletValue );

    // the alias is left to the T-matrix
    EXPECT_FALSE( tAlias->is_fixed() );

    // the ordinary sideset node is fixed as before
    EXPECT_TRUE( tStrip.node_dof( 8 )->is_fixed() );
    EXPECT_EQ( tStrip.node_dof( 8 )->value(), tDirichletValue );

    // nothing else was touched
    for ( const id_t tNodeID : { 1, 2, 4, 5, 7 } )
    {
        EXPECT_FALSE( tStrip.node_dof( tNodeID )->is_fixed() )
            << "node " << tNodeID ;
    }
}

//------------------------------------------------------------------------------

TEST( DirichletPin, NearUnitWeightIsAcceptedExactly )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    // one ulp above unity passes the unit-weight tolerance, and v / w and
    // v * w differ at that level: this is where the division is observable
    const real tWeight = std::nextafter( 1.0, 2.0 );

    GluedStrip tStrip( tWeight );

    fem::Dof * tAlias  = tStrip.node_dof( 6 );
    fem::Dof * tSource = tStrip.node_dof( 3 );

    ASSERT_EQ( tAlias->number_of_sources(), 1u );
    ASSERT_EQ( tAlias->weight( 0 ), tWeight );

    ASSERT_NO_THROW( tStrip.impose_top( tDirichletValue ) );

    EXPECT_TRUE( tSource->is_fixed() );
    EXPECT_EQ( tSource->value(), tDirichletValue / tWeight );
    EXPECT_NE( tSource->value(), tDirichletValue * tWeight );
    EXPECT_FALSE( tAlias->is_fixed() );
}

//------------------------------------------------------------------------------

TEST( DirichletPin, NonUnitWeightIsRefused )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    GluedStrip tStrip( 0.5 );

    fem::Dof * tAlias  = tStrip.node_dof( 6 );
    fem::Dof * tSource = tStrip.node_dof( 3 );

    ASSERT_EQ( tAlias->number_of_sources(), 1u );
    ASSERT_EQ( tAlias->weight( 0 ), 0.5 );

    expect_refusal( [ & ]() { tStrip.impose_top( tDirichletValue ); },
                    "cannot be pinned" );

    // the refusal precedes any fix
    EXPECT_FALSE( tSource->is_fixed() );
    EXPECT_FALSE( tAlias->is_fixed() );
}

//------------------------------------------------------------------------------

TEST( DirichletPin, FirstFlipAfterInitializationIsRefused )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    GluedStrip tStrip ;

    fem::Dof * tSource = tStrip.node_dof( 3 );

    tStrip.mField->initialize();
    ASSERT_TRUE( tStrip.mField->is_initialized() );

    expect_refusal( [ & ]() { tStrip.impose_top( tDirichletValue ); },
                    "after initialization" );

    EXPECT_FALSE( tSource->is_fixed() );
}

//------------------------------------------------------------------------------

TEST( DirichletPin, ReimpositionAfterInitializationIsAccepted )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    GluedStrip tStrip ;

    fem::Dof * tSource = tStrip.node_dof( 3 );

    // the per-timestep path: the source was fixed before the freeze, so a
    // later imposition only updates its value
    tStrip.impose_top( tDirichletValue );
    tStrip.mField->initialize();
    ASSERT_TRUE( tStrip.mField->is_initialized() );

    // a distinct value, so a reimposition that stops writing is caught
    const real tUpdated = 2.0 * tDirichletValue ;

    ASSERT_NO_THROW( tStrip.impose_top( tUpdated ) );

    EXPECT_TRUE( tSource->is_fixed() );
    EXPECT_EQ( tSource->value(), tUpdated );
    EXPECT_FALSE( tStrip.node_dof( 6 )->is_fixed() );
}

//------------------------------------------------------------------------------

TEST( DirichletPin, TwoSourcesFallBackToTheAliasFix )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    // this locks the helper's documented fallback for the two-source shape:
    // it returns false and the caller keeps fixing the alias ( a single
    // non-unit source is refused instead, see above ). It does NOT show that
    // the multi-source condition is enforced: the alias is eliminated from
    // the system and reconstructed from its ( free ) sources, so with both
    // sources outside the sideset the condition is lost
    GluedStrip tStrip( 1.0, true );

    fem::Dof * tAlias = tStrip.node_dof( 6 );

    ASSERT_TRUE( tAlias->is_hanging() );
    ASSERT_EQ( tAlias->number_of_sources(), 2u );

    ASSERT_NO_THROW( tStrip.impose_top( tDirichletValue ) );

    EXPECT_TRUE( tAlias->is_fixed() );
    EXPECT_EQ( tAlias->value(), tDirichletValue );
    EXPECT_FALSE( tStrip.node_dof( 3 )->is_fixed() );
    EXPECT_FALSE( tStrip.node_dof( 2 )->is_fixed() );
}
