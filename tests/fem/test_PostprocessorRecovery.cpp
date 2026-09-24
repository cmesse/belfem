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

/*
 * The node-local patch recovery of fem::Postprocessor, run for real on
 * code-built meshes through a Kernel and a DofManager. A test postprocessor
 * returns a manufactured field at every integration point, the recovery fits
 * it, and the nodal values are compared with the field at the nodes.
 *
 * The fit is the weighted least-squares patch fit of Zienkiewicz & Zhu 1992.
 * The contracts pinned here are BELFEM's (postprocessor_recovery_theory.md):
 * the basis is centered and scaled at every node (§2.1.1), linear thin-shell
 * layers are fit in the shell coordinates ( s1, s2, eta ) (§2.1.2, BELFEM's
 * extension), and the constant is solved exactly while only the gradient
 * block is truncated (§3.2). A field inside the basis must come back to
 * round-off, whatever the element weights. The intercept contract of §3.2
 * cannot be provoked through valid elements ( a node lies in its own patch ),
 * so the block inverse is also tested directly.
 *
 * Serial only: the fixtures never partition.
 */

#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <utility>

#include "typedefs.hpp"
#include "globals.hpp"
#include "constants.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"
#include "commtools.hpp"
#include "fn_FEM_recovery_inverse.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "fn_trans.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_Edge.hpp"
#include "cl_Facet.hpp"
#include "cl_Block.hpp"
#include "cl_SideSet.hpp"
#include "cl_ThinShell.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_MeshChecker.hpp"
#include "cl_Mesh_ConnectivityCalculator.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_KernelParameters.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_FEM_Calculator.hpp"
#include "cl_FEM_Postprocessor.hpp"
#include "cl_IF_IntegrationData.hpp"
#include "cl_IWG_Maxwell.hpp"
#include "en_Maxwell_Formulations.hpp"
#include "en_SolverEnums.hpp"

using namespace belfem;

namespace
{
    //! field value at a physical point, given the element's reference
    //! thickness coordinate zeta ( last coordinate of the rule ) and block id
    using ManufacturedField = std::function< real( const Vector< real > &, real, id_t ) > ;

//------------------------------------------------------------------------------

    //! recovers one target field "f" from a manufactured field
    class ManufacturedRecovery : public fem::Postprocessor
    {
        ManufacturedField mFunction ;

        // scratch, sized on first use
        Matrix< real > mCoords ;
        Vector< real > mPoint ;
        Vector< real > mValue ;

    public:

        ManufacturedRecovery(
                fem::Kernel     * aKernel,
                fem::DofManager * aField,
                const Vector< id_t > & aBlockIDs,
                ManufacturedField aFunction ) :
            fem::Postprocessor( aKernel ),
            mFunction( std::move( aFunction ) ),
            mValue( 1, 0.0 )
        {
            this->set_field( aField );
            this->set_block_ids( aBlockIDs );
            this->set_target_fields( { "f" } );
        }

    protected:

        void
        update_element_dofs() override
        {
        }

        const Vector< real > &
        compute( const uint aK ) override
        {
            const uint tNumNodes = mElement->element()->number_of_nodes() ;
            const uint tNumDimensions = mMesh->number_of_dimensions() ;
            if ( mCoords.n_rows() != tNumNodes || mCoords.n_cols() != tNumDimensions )
            {
                mCoords.set_size( tNumNodes, tNumDimensions );
                mPoint.set_size( tNumDimensions );
            }
            mElement->get_node_coors( mCoords );
            mPoint = trans( mCoords ) * mCalculator->Nvec( aK );

            const Matrix< real > & tXi = mCalculator->integration()->points() ;
            mValue( 0 ) = mFunction( mPoint, tXi( tXi.n_rows() - 1, aK ), mElement->element()->block_id() );
            return mValue ;
        }
    };

//------------------------------------------------------------------------------

    //! one Edge object per node pair, shared by every element that uses it
    void
    create_unique_edges( Mesh * aMesh, const Cell< mesh::Element * > & aElements )
    {
        Map< luint, mesh::Edge * > tEdges ;
        id_t tEdgeID = 1 ;
        Cell< mesh::Node * > tEdgeNodes ;
        for ( mesh::Element * tElement : aElements )
        {
            tElement->allocate_edge_container();
            for ( uint e = 0; e < tElement->number_of_edges(); ++e )
            {
                tElement->get_nodes_of_edge( e, tEdgeNodes );
                const id_t tA = tEdgeNodes( 0 )->id() ;
                const id_t tB = tEdgeNodes( 1 )->id() ;
                const luint tKey = ( ( luint ) std::min( tA, tB ) << 32 ) | ( luint ) std::max( tA, tB );
                if ( ! tEdges.key_exists( tKey ) )
                {
                    mesh::Edge * tEdge = new mesh::Edge();
                    tEdge->set_id( tEdgeID++ );
                    tEdge->allocate_node_container( 2 );
                    tEdge->insert_node( tEdgeNodes( 0 ), 0 );
                    tEdge->insert_node( tEdgeNodes( 1 ), 1 );
                    tEdge->set_owner( 0 );
                    aMesh->edges().push( tEdge );
                    tEdges[ tKey ] = tEdge ;
                }
                tElement->insert_edge( tEdges( tKey ), e );
            }
        }
        aMesh->finalize_edges();
    }

//------------------------------------------------------------------------------

    //! owns the mesh and the kernel chain of one fixture
    class RecoveryFixture
    {
        //! h_ghost falls back to the nodal "T" field when gTbulk is NaN
        struct TbulkGuard
        {
            const real mSaved ;
            TbulkGuard() : mSaved( gTbulk ) { gTbulk = 77.0 ; }
            ~TbulkGuard() { gTbulk = mSaved ; }
        };

        TbulkGuard mTbulkGuard ;

    public:

        Mesh                  * mMesh   = nullptr ;
        fem::KernelParameters * mParams = nullptr ;
        fem::Kernel           * mKernel = nullptr ;
        fem::DofManager       * mField  = nullptr ;   // owned by the kernel
        fem::IWG_Maxwell      * mIWG    = nullptr ;   // owned by the kernel
        Vector< id_t >          mBlockIDs ;

        RecoveryFixture() = default ;

        ~RecoveryFixture()
        {
            delete mKernel ;
            delete mParams ;
            delete mMesh ;
        }

        RecoveryFixture( const RecoveryFixture & ) = delete ;
        RecoveryFixture & operator=( const RecoveryFixture & ) = delete ;

        //! kernel, equation and dof manager for the given block types
        void
        create_field( const Cell< DomainType > & aBlockTypes )
        {
            mParams = new fem::KernelParameters( mMesh );
            mKernel = new fem::Kernel( mParams );

            mIWG = new fem::IWG_Maxwell(
                    fem::maxwell::Formulation::HPhi,
                    ModelDimensionality::ThreeD,
                    false, false );
            mIWG->set_blocks( mBlockIDs, aBlockTypes );
            mKernel->add_equation( mIWG );

            mField = mKernel->create_field( mIWG );
            mField->create_fields( mField->iwg() );
            if ( mMesh->field_exists( "element_rho" ) )
            {
                mMesh->field_data( "element_rho" ).fill( 1.0e-8 );
            }
            mField->set_solver( gDefaultSolver );
            mField->initialize();
        }

        //! runs the recovery of aFunction and returns the nodal field "f";
        //! aTruncated receives the count of fits that dropped a direction
        const Vector< real > &
        recover( ManufacturedField aFunction, index_t & aTruncated )
        {
            ManufacturedRecovery tRecovery( mKernel, mField, mBlockIDs, std::move( aFunction ) );
            tRecovery.initialize();
            tRecovery.run();
            aTruncated = tRecovery.number_of_truncated_patches() ;
            return mMesh->field_data( "f" );
        }
    };

//------------------------------------------------------------------------------

    constexpr real gRadius     = 1.0 ;
    constexpr real gThickness  = 2.0e-3 ;
    constexpr real gHeight     = 0.1 ;
    constexpr uint gStations   = 5 ;
    constexpr uint gRows       = 3 ;
    constexpr real gSpanDeg    = 40.0 ;

    //! F1: two PENTA6TS layers ( blocks 1, 2 ) on a strip of the cylinder
    //! r = R around the z axis. Node rows at r = R - t, R, R + t, extruded
    //! exactly along the radius. The tape facets sit on the middle row, the
    //! interface of the two layers. The wedges of each layer are inserted in
    //! tape-facet order with consecutive element ids: the kernel measures shell
    //! volumes by the id offset, the dof manager links facets by position
    class CurvedShellFan : public RecoveryFixture
    {
    public:

        //! node ( station, z row, radial level 0 .. 2 )
        Cell< mesh::Node * > mGrid ;

        mesh::Node *
        node( const uint aStation, const uint aRow, const uint aLevel )
        {
            return mGrid( ( aLevel * gRows + aRow ) * gStations + aStation );
        }

        CurvedShellFan()
        {
            mMesh = new Mesh( 3, 0, false );
            mesh::ElementFactory tFactory ;

            Cell< mesh::Node * > & tNodes = mMesh->nodes();
            id_t tNodeID = 1 ;
            const real tRadii[ 3 ] = { gRadius - gThickness, gRadius, gRadius + gThickness };
            for ( uint l = 0; l < 3; ++l )
            {
                for ( uint r = 0; r < gRows; ++r )
                {
                    for ( uint s = 0; s < gStations; ++s )
                    {
                        const real tPhi = ( -0.5 * gSpanDeg + gSpanDeg * s / ( gStations - 1 ) )
                                * constant::pi / 180.0 ;
                        mesh::Node * tNode = new mesh::Node( tNodeID++,
                                tRadii[ l ] * std::cos( tPhi ),
                                tRadii[ l ] * std::sin( tPhi ),
                                gHeight * r );
                        tNodes.push( tNode );
                        mGrid.push( tNode );
                    }
                }
            }

            // mid-surface triangles, counterclockwise seen from outside ( +r ),
            // the diagonal from ( s, r ) to ( s+1, r+1 ): the interior node
            // ( 2, 1 ) touches six of them
            Cell< std::array< uint, 6 > > tTriangles ;
            for ( uint r = 0; r + 1 < gRows; ++r )
            {
                for ( uint s = 0; s + 1 < gStations; ++s )
                {
                    tTriangles.push( { s, r, s + 1, r, s + 1, r + 1 } );
                    tTriangles.push( { s, r, s + 1, r + 1, s, r + 1 } );
                }
            }

            mesh::SideSet * tTapeSet = new mesh::SideSet( 10, tTriangles.size() );
            tTapeSet->set_domain_type( DomainType::ThinShell );

            mesh::Block * tBlocks[ 2 ] = {
                    new mesh::Block( 1, tTriangles.size() ),
                    new mesh::Block( 2, tTriangles.size() ) };

            Cell< mesh::Element * > & tElements = mMesh->elements();
            id_t tFacetID = 1000 ;

            for ( uint l = 0; l < 2; ++l )
            {
                tBlocks[ l ]->set_thickness( gThickness );
                tBlocks[ l ]->set_domain_type( DomainType::ThinShell );
            }

            // layer 0 first, then layer 1: ids are consecutive per block
            id_t tElementID = 1 ;
            for ( uint l = 0; l < 2; ++l )
            {
                for ( const auto & tTri : tTriangles )
                {
                    mesh::Element * tPrism = tFactory.create_element( ElementType::PENTA6TS, tElementID++ );
                    for ( uint k = 0; k < 3; ++k )
                    {
                        tPrism->insert_node( this->node( tTri[ 2 * k ], tTri[ 2 * k + 1 ], l ), k );
                        tPrism->insert_node( this->node( tTri[ 2 * k ], tTri[ 2 * k + 1 ], l + 1 ), k + 3 );
                    }
                    tPrism->set_block_id( l + 1 );
                    tElements.push( tPrism );
                    tBlocks[ l ]->insert_element( tPrism );
                }
            }

            for ( const auto & tTri : tTriangles )
            {
                mesh::Element * tSurf = tFactory.create_element( ElementType::TRI3, tFacetID++ );
                for ( uint k = 0; k < 3; ++k )
                {
                    tSurf->insert_node( this->node( tTri[ 2 * k ], tTri[ 2 * k + 1 ], 1 ), k );
                }
                mesh::Facet * tFacet = new mesh::Facet( tSurf );
                tFacet->set_sideset_id( 10 );
                tTapeSet->insert_facet( tFacet );
            }
            mMesh->sidesets().push( tTapeSet );

            mesh::ThinShell * tShell = new mesh::ThinShell( tTapeSet, nullptr );
            tShell->blocks().push( tBlocks[ 0 ] );
            tShell->blocks().push( tBlocks[ 1 ] );
            {
                Vector< real > tThicknesses( 2, gThickness );
                tShell->set_thicknesses( tThicknesses );
            }
            mMesh->thin_shells().push( tShell );

            mMesh->blocks().push( tBlocks[ 0 ] );
            mMesh->blocks().push( tBlocks[ 1 ] );
            mMesh->finalize();

            {
                MeshChecker tCheck( mMesh );
            }

            // Mesh( 3, 0, false ) skips the connectivity in finalize() so that
            // the mesh checker can run before edges exist; the recovery walks
            // node -> element
            {
                mesh::ConnectivityCalculator tConnect( mMesh );
                tConnect.connect_nodes_to_elements();
            }

            create_unique_edges( mMesh, mMesh->elements() );

            mBlockIDs.set_size( 2 );
            mBlockIDs( 0 ) = 1 ;
            mBlockIDs( 1 ) = 2 ;
            this->create_field( Cell< DomainType >( 2, DomainType::ThinShell ) );
        }
    };

//------------------------------------------------------------------------------

    constexpr real gCell = 1.0e-3 ;
    const real gOrigin[ 3 ] = { 1.0e3, -2.0e3, 5.0e2 };

    //! F2: 2 x 2 x 2 cubes of edge h = 1e-3, six TET4 or TET10 each, one Air
    //! block, 1e3 from the origin. aSplitPlane duplicates the nodes of the
    //! z = h plane for the upper cubes, the pattern of a cut
    class FarTetPatch : public RecoveryFixture
    {
    public:

        Cell< mesh::Node * > mVertices ;   // 3 x 3 x 3 grid, originals

        FarTetPatch( const bool aQuadratic, const bool aSplitPlane = false )
        {
            mMesh = new Mesh( 3, 0, false );
            mesh::ElementFactory tFactory ;
            Cell< mesh::Node * > & tNodes = mMesh->nodes();
            id_t tNodeID = 1 ;

            for ( uint k = 0; k < 3; ++k )
            {
                for ( uint j = 0; j < 3; ++j )
                {
                    for ( uint i = 0; i < 3; ++i )
                    {
                        mesh::Node * tNode = new mesh::Node( tNodeID++,
                                gOrigin[ 0 ] + gCell * i,
                                gOrigin[ 1 ] + gCell * j,
                                gOrigin[ 2 ] + gCell * k );
                        tNodes.push( tNode );
                        mVertices.push( tNode );
                    }
                }
            }

            // the duplicates of the z = h plane, co-located with the originals
            Cell< mesh::Node * > tUpper( 9, nullptr );
            if ( aSplitPlane )
            {
                for ( uint n = 0; n < 9; ++n )
                {
                    mesh::Node * tOriginal = mVertices( 9 + n );
                    mesh::Node * tDuplicate = new mesh::Node( tNodeID++,
                            tOriginal->x(), tOriginal->y(), tOriginal->z() );
                    tOriginal->allocate_duplicate_container( 1 );
                    tOriginal->add_duplicate( tDuplicate );
                    tDuplicate->set_original( tOriginal );
                    tNodes.push( tDuplicate );
                    tUpper( n ) = tDuplicate ;
                }
            }

            // midside nodes, one per vertex pair
            Map< luint, mesh::Node * > tMidNodes ;
            auto tMid = [ & ]( mesh::Node * aA, mesh::Node * aB ) -> mesh::Node *
            {
                const luint tKey = ( ( luint ) std::min( aA->id(), aB->id() ) << 32 ) | ( luint ) std::max( aA->id(), aB->id() );
                if ( tMidNodes.key_exists( tKey ) ) return tMidNodes( tKey );
                mesh::Node * tNode = new mesh::Node( tNodeID++,
                        0.5 * ( aA->x() + aB->x() ), 0.5 * ( aA->y() + aB->y() ), 0.5 * ( aA->z() + aB->z() ) );
                tNodes.push( tNode );
                tMidNodes[ tKey ] = tNode ;
                return tNode ;
            };

            // the six tetrahedra of a cube along its main diagonal ( corner 0 to 7 )
            const uint tKuhn[ 6 ][ 4 ] = {
                    { 0, 1, 3, 7 }, { 0, 1, 5, 7 }, { 0, 2, 3, 7 },
                    { 0, 2, 6, 7 }, { 0, 4, 5, 7 }, { 0, 4, 6, 7 } };

            mesh::Block * tBlock = new mesh::Block( 1, 48 );
            tBlock->set_domain_type( DomainType::Air );
            Cell< mesh::Element * > & tElements = mMesh->elements();
            id_t tElementID = 1 ;

            for ( uint ck = 0; ck < 2; ++ck )
            {
                for ( uint cj = 0; cj < 2; ++cj )
                {
                    for ( uint ci = 0; ci < 2; ++ci )
                    {
                        mesh::Node * tCorner[ 8 ];
                        for ( uint c = 0; c < 8; ++c )
                        {
                            const uint i = ci + ( c & 1 ), j = cj + ( ( c >> 1 ) & 1 ), k = ck + ( ( c >> 2 ) & 1 );
                            mesh::Node * tNode = mVertices( ( k * 3 + j ) * 3 + i );
                            // the upper cubes reference the duplicates of the split plane
                            if ( aSplitPlane && ck == 1 && k == 1 )
                            {
                                tNode = tUpper( j * 3 + i );
                            }
                            tCorner[ c ] = tNode ;
                        }

                        for ( uint t = 0; t < 6; ++t )
                        {
                            mesh::Node * tV[ 4 ] = {
                                    tCorner[ tKuhn[ t ][ 0 ] ], tCorner[ tKuhn[ t ][ 1 ] ],
                                    tCorner[ tKuhn[ t ][ 2 ] ], tCorner[ tKuhn[ t ][ 3 ] ] };

                            // BELFEM's TET reference: node 0 at xi = ( 1, 0, 0 ), node 1 at
                            // ( 0, 0, 1 ), node 2 at ( 0, 1, 0 ), node 3 at the origin; the
                            // mesh checker reorients a negative element before the kernel
                            mesh::Element * tTet = tFactory.create_element(
                                    aQuadratic ? ElementType::TET10 : ElementType::TET4, tElementID++ );
                            for ( uint v = 0; v < 4; ++v )
                            {
                                tTet->insert_node( tV[ v ], v );
                            }
                            if ( aQuadratic )
                            {
                                // midside order from the reference coordinates of TET10
                                const uint tPairs[ 6 ][ 2 ] = {
                                        { 0, 1 }, { 1, 2 }, { 0, 2 }, { 0, 3 }, { 1, 3 }, { 2, 3 } };
                                for ( uint m = 0; m < 6; ++m )
                                {
                                    tTet->insert_node( tMid( tV[ tPairs[ m ][ 0 ] ], tV[ tPairs[ m ][ 1 ] ] ), 4 + m );
                                }
                            }
                            tTet->set_block_id( 1 );
                            tElements.push( tTet );
                            tBlock->insert_element( tTet );
                        }
                    }
                }
            }

            mMesh->blocks().push( tBlock );
            mMesh->finalize();
            {
                MeshChecker tCheck( mMesh );
            }

            // Mesh( 3, 0, false ) skips the connectivity in finalize() so that
            // the mesh checker can run before edges exist; the recovery walks
            // node -> element
            {
                mesh::ConnectivityCalculator tConnect( mMesh );
                tConnect.connect_nodes_to_elements();
            }

            mBlockIDs.set_size( 1, 1 );
            this->create_field( Cell< DomainType >( 1, DomainType::Air ) );
        }
    };

//------------------------------------------------------------------------------

    //! F4: one flat TET4, 1 x 1 x 1e-6: every node's patch is the element,
    //! whose integration points span the thickness direction a millionth of
    //! the in-plane one, so the gradient block's condition number exceeds the
    //! truncation cut
    class FlatTet : public RecoveryFixture
    {
    public:

        FlatTet()
        {
            mMesh = new Mesh( 3, 0, false );
            mesh::ElementFactory tFactory ;
            Cell< mesh::Node * > & tNodes = mMesh->nodes();
            tNodes.push( new mesh::Node( 1, 1.0, 0.0, 0.0 ) );
            tNodes.push( new mesh::Node( 2, 0.0, 0.0, 1.0e-6 ) );
            tNodes.push( new mesh::Node( 3, 0.0, 1.0, 0.0 ) );
            tNodes.push( new mesh::Node( 4, 0.0, 0.0, 0.0 ) );

            mesh::Block * tBlock = new mesh::Block( 1, 1 );
            tBlock->set_domain_type( DomainType::Air );
            mesh::Element * tTet = tFactory.create_element( ElementType::TET4, 1 );
            for ( uint v = 0; v < 4; ++v )
            {
                tTet->insert_node( tNodes( v ), v );
            }
            tTet->set_block_id( 1 );
            mMesh->elements().push( tTet );
            tBlock->insert_element( tTet );
            mMesh->blocks().push( tBlock );
            mMesh->finalize();
            {
                MeshChecker tCheck( mMesh );
            }

            // Mesh( 3, 0, false ) skips the connectivity in finalize() so that
            // the mesh checker can run before edges exist; the recovery walks
            // node -> element
            {
                mesh::ConnectivityCalculator tConnect( mMesh );
                tConnect.connect_nodes_to_elements();
            }

            mBlockIDs.set_size( 1, 1 );
            this->create_field( Cell< DomainType >( 1, DomainType::Air ) );
        }
    };

//------------------------------------------------------------------------------

    real
    radius_xy( const mesh::Node * aNode )
    {
        return std::sqrt( aNode->x() * aNode->x() + aNode->y() * aNode->y() );
    }
}

//------------------------------------------------------------------------------

// the transport current of a superconducting layer makes the tangential
// field change through the layer; on a curved patch a Cartesian linear fit
// spreads that change over both faces. The manufactured field is the ramp
// eta = ( 1 + zeta ) t / 2 + layer * t through the two-layer stack: linear
// in the shell basis of every node, so each node row must come back exactly
TEST( PostprocessorRecovery, ThicknessJumpSurvivesCurvature )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    CurvedShellFan tFan ;

    // why the test discriminates: over the interior node's patch the shell
    // bends away from its chord by several layer thicknesses
    const real tHalfAngle = 0.5 * gSpanDeg / ( gStations - 1 ) * 2.0 * constant::pi / 180.0 ;
    const real tSagitta = gRadius * ( 1.0 - std::cos( tHalfAngle ) );
    ASSERT_GT( tSagitta, 5.0 * gThickness );

    index_t tTruncated = 0 ;
    const Vector< real > & tF = tFan.recover(
            []( const Vector< real > &, const real aZeta, const id_t aBlock ) -> real
            {
                return 0.5 * ( 1.0 + aZeta ) * gThickness + ( aBlock - 1 ) * gThickness ;
            }, tTruncated );

    EXPECT_EQ( tTruncated, 0u );
    for ( uint l = 0; l < 3; ++l )
    {
        for ( uint r = 0; r < gRows; ++r )
        {
            for ( uint s = 0; s < gStations; ++s )
            {
                const mesh::Node * tNode = tFan.node( s, r, l );
                EXPECT_NEAR( tF( tNode->index() ), l * gThickness, 1.0e-12 * gThickness )
                    << "node row " << l << ", station " << s << ", z row " << r ;
            }
        }
    }
}

//------------------------------------------------------------------------------

// the constant is solved through the Schur complement of the patch weight,
// never through the truncation
TEST( PostprocessorRecovery, ConstantReproducedOnCurvedShell )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    CurvedShellFan tFan ;
    index_t tTruncated = 0 ;
    const Vector< real > & tF = tFan.recover(
            []( const Vector< real > &, const real, const id_t ) -> real { return 3.7 ; },
            tTruncated );

    EXPECT_EQ( tTruncated, 0u );
    for ( mesh::Node * tNode : tFan.mGrid )
    {
        EXPECT_NEAR( tF( tNode->index() ), 3.7, 1.0e-13 * 3.7 );
    }
}

//------------------------------------------------------------------------------

// the cylinder axis lies in every node's tangent plane, and eta does not
// depend on z: a field linear in z is inside the shell basis
TEST( PostprocessorRecovery, AxialLinearFieldReproducedOnCurvedShell )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    CurvedShellFan tFan ;
    index_t tTruncated = 0 ;
    const Vector< real > & tF = tFan.recover(
            []( const Vector< real > & aX, const real, const id_t ) -> real { return 2.0 + 5.0 * aX( 2 ); },
            tTruncated );

    EXPECT_EQ( tTruncated, 0u );
    for ( mesh::Node * tNode : tFan.mGrid )
    {
        EXPECT_NEAR( tF( tNode->index() ), 2.0 + 5.0 * tNode->z(), 1.0e-12 * 3.0 )
            << "node at r = " << radius_xy( tNode ) << ", z = " << tNode->z() ;
    }
}

//------------------------------------------------------------------------------

// in absolute coordinates, the Gram matrix of a 1e-3 patch 1e3 from the
// origin has a condition number of about 4e12; centered and scaled at the
// node it is of order one, and a linear field is exact at every node,
// interior and boundary
TEST( PostprocessorRecovery, LinearFieldFarFromOrigin )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    FarTetPatch tPatch( false );
    index_t tTruncated = 0 ;
    auto tLinear = []( const real aX, const real aY, const real aZ ) -> real
    {
        return 1.0 + 2.0e3 * ( aX - gOrigin[ 0 ] ) - 3.0e3 * ( aY - gOrigin[ 1 ] ) + 4.0e3 * ( aZ - gOrigin[ 2 ] );
    };
    const Vector< real > & tF = tPatch.recover(
            [ & ]( const Vector< real > & aX, const real, const id_t ) -> real
            {
                return tLinear( aX( 0 ), aX( 1 ), aX( 2 ) );
            }, tTruncated );

    EXPECT_EQ( tTruncated, 0u );
    for ( mesh::Node * tNode : tPatch.mVertices )
    {
        EXPECT_NEAR( tF( tNode->index() ), tLinear( tNode->x(), tNode->y(), tNode->z() ), 1.0e-10 * 20.0 );
    }
}

//------------------------------------------------------------------------------

// the mesh order sets the Cartesian basis: TET10 gives the complete
// quadratic, which reproduces a quadratic field, far from the origin
TEST( PostprocessorRecovery, QuadraticFieldFarFromOrigin )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    FarTetPatch tPatch( true );
    index_t tTruncated = 0 ;
    auto tQuadratic = []( const real aX, const real aY, const real aZ ) -> real
    {
        const real u = ( aX - gOrigin[ 0 ] ) / gCell ;
        const real v = ( aY - gOrigin[ 1 ] ) / gCell ;
        const real w = ( aZ - gOrigin[ 2 ] ) / gCell ;
        return 1.0 + u - 2.0 * v + 0.5 * w + u * u - 0.7 * u * v + 0.3 * v * w + 1.1 * w * w - 0.4 * u * w + 0.2 * v * v ;
    };
    const Vector< real > & tF = tPatch.recover(
            [ & ]( const Vector< real > & aX, const real, const id_t ) -> real
            {
                return tQuadratic( aX( 0 ), aX( 1 ), aX( 2 ) );
            }, tTruncated );

    EXPECT_EQ( tTruncated, 0u );
    for ( mesh::Node * tNode : tPatch.mMesh->nodes() )
    {
        EXPECT_NEAR( tF( tNode->index() ), tQuadratic( tNode->x(), tNode->y(), tNode->z() ), 1.0e-9 * 10.0 );
    }
}

//------------------------------------------------------------------------------

// co-located siblings ( the cut pattern ) are fit from the union of their
// patches in one frame. The union has the same elements and integration
// points as the unsplit patch, so both siblings must reproduce the unsplit
// recovery. The field is cubic through the split plane, outside the linear
// basis, so a fit of either half alone gives a different value
TEST( PostprocessorRecovery, SiblingsShareTheUnionFit )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    auto tField = []( const real aX, const real aY, const real aZ ) -> real
    {
        const real u = ( aX - gOrigin[ 0 ] ) / gCell ;
        const real v = ( aY - gOrigin[ 1 ] ) / gCell ;
        const real w = ( aZ - gOrigin[ 2 ] ) / gCell ;
        return -2.0 + u + 2.0 * v + w * w * w ;
    };
    auto tManufactured = [ & ]( const Vector< real > & aX, const real, const id_t ) -> real
    {
        return tField( aX( 0 ), aX( 1 ), aX( 2 ) );
    };

    // reference: the same grid without the split
    Vector< real > tReference( 9 );
    {
        FarTetPatch tWhole( false );
        index_t tTruncated = 0 ;
        const Vector< real > & tF = tWhole.recover( tManufactured, tTruncated );
        for ( uint n = 0; n < 9; ++n )
        {
            tReference( n ) = tF( tWhole.mVertices( 9 + n )->index() );
        }
    }

    FarTetPatch tPatch( false, true );
    index_t tTruncated = 0 ;
    const Vector< real > & tF = tPatch.recover( tManufactured, tTruncated );

    for ( uint n = 0; n < 9; ++n )
    {
        mesh::Node * tOriginal = tPatch.mVertices( 9 + n );
        ASSERT_EQ( tOriginal->number_of_duplicates(), 1u );
        const mesh::Node * tSibling = tOriginal->duplicate( 0 );
        EXPECT_NEAR( tF( tOriginal->index() ), tReference( n ), 1.0e-11 * 10.0 );
        EXPECT_NEAR( tF( tSibling->index() ), tReference( n ), 1.0e-11 * 10.0 );
    }
}

//------------------------------------------------------------------------------

// a flat element drops the thickness direction of every fit instead of
// aborting, and keeps the constant and the in-plane gradient. It does not
// separate the Schur inverse from a pseudo-inverse of the whole matrix: the
// node lies in the element's plane, so the constant hardly couples to the
// dropped direction. The coupled case is tested on the inverse directly
TEST( PostprocessorRecovery, DegeneratePatchKeepsConstantAndInPlaneGradient )
{
    if ( comm_size() != 1 )
    {
        GTEST_SKIP() << "serial-only fixture" ;
    }

    FlatTet tTet ;
    index_t tTruncated = 0 ;
    const Vector< real > & tF = tTet.recover(
            []( const Vector< real > & aX, const real, const id_t ) -> real
            {
                return 3.7 + 2.0 * aX( 0 ) - aX( 1 );
            }, tTruncated );

    EXPECT_GE( tTruncated, 1u );
    for ( mesh::Node * tNode : tTet.mMesh->nodes() )
    {
        EXPECT_NEAR( tF( tNode->index() ), 3.7 + 2.0 * tNode->x() - tNode->y(), 1.0e-8 );
    }
}

//------------------------------------------------------------------------------

namespace
{
    //! scratch and call of the block inverse for an m x m Gram matrix
    bool
    block_inverse( const Matrix< real > & aGram, Matrix< real > & aInverse, real & aCondition )
    {
        const uint m = aGram.n_rows();
        aInverse.set_size( m, m, 0.0 );
        Matrix< real > tSchur( m - 1, m - 1 ), tU, tVT ;
        Vector< real > tCoupling( m - 1 ), tSchurCoupling( m - 1 ), tS, tWork ;
        return fem::recovery_block_inverse( aGram, 1.0e10, aInverse, aCondition,
                tSchur, tCoupling, tSchurCoupling, tU, tVT, tS, tWork );
    }
}

//------------------------------------------------------------------------------

// §3.2: on a healthy Gram matrix nothing is dropped and the block inverse is
// the inverse
TEST( PostprocessorRecovery, BlockInverseOfHealthyGramIsTheInverse )
{
    Matrix< real > tGram = {
            { 4.0, 0.3, -0.2, 0.1 },
            { 0.3, 2.0,  0.4, 0.0 },
            {-0.2, 0.4,  3.0, 0.5 },
            { 0.1, 0.0,  0.5, 1.5 } };
    Matrix< real > tInverse ;
    real tCondition = 0.0 ;
    EXPECT_FALSE( block_inverse( tGram, tInverse, tCondition ) );
    EXPECT_LT( tCondition, 1.0e2 );

    for ( uint c = 0; c < 4; ++c )
    {
        for ( uint r = 0; r < 4; ++r )
        {
            real tProduct = 0.0 ;
            for ( uint k = 0; k < 4; ++k )
            {
                tProduct += tGram( r, k ) * tInverse( k, c );
            }
            EXPECT_NEAR( tProduct, r == c ? 1.0 : 0.0, 1.0e-13 );
        }
    }
}

//------------------------------------------------------------------------------

// §3.2: every sample of the patch at one point, off the node: the gradient
// in the x direction is fully coupled to the constant, G = w v v' with
// v = ( 1, 1, 0, 0 ). The constant must come back exactly; a truncated
// pseudo-inverse of the whole matrix returns half of it
TEST( PostprocessorRecovery, BlockInverseKeepsCoupledConstant )
{
    const real v[ 4 ] = { 1.0, 1.0, 0.0, 0.0 };
    const real w = 2.5 ;
    const real f = 3.7 ;
    Matrix< real > tGram( 4, 4, 0.0 );
    Vector< real > tRhs( 4, 0.0 );
    for ( uint c = 0; c < 4; ++c )
    {
        for ( uint r = 0; r < 4; ++r )
        {
            tGram( r, c ) = w * v[ r ] * v[ c ];
        }
        tRhs( c ) = w * v[ c ] * f ;
    }

    Matrix< real > tInverse ;
    real tCondition = 0.0 ;
    EXPECT_TRUE( block_inverse( tGram, tInverse, tCondition ) );

    real tConstant = 0.0 ;
    for ( uint k = 0; k < 4; ++k )
    {
        tConstant += tInverse( 0, k ) * tRhs( k );
    }
    EXPECT_NEAR( tConstant, f, 1.0e-13 * f );
}

//------------------------------------------------------------------------------

// §3.2: samples on a plane through the node ( no z spread ): the z direction
// is dropped, the constant and the in-plane gradient are kept
TEST( PostprocessorRecovery, BlockInverseKeepsInPlaneFitOfCoplanarPatch )
{
    const real tPoints[ 5 ][ 3 ] = {
            { 0.3, 0.1, 0.0 }, { -0.4, 0.2, 0.0 }, { 0.1, -0.5, 0.0 }, { 0.6, 0.6, 0.0 }, { -0.2, -0.3, 0.0 } };
    Matrix< real > tGram( 4, 4, 0.0 );
    Vector< real > tRhs( 4, 0.0 );
    for ( uint p = 0; p < 5; ++p )
    {
        const real tP[ 4 ] = { 1.0, tPoints[ p ][ 0 ], tPoints[ p ][ 1 ], tPoints[ p ][ 2 ] };
        const real tY = 1.5 - 2.0 * tPoints[ p ][ 0 ] + 0.5 * tPoints[ p ][ 1 ] ;
        for ( uint c = 0; c < 4; ++c )
        {
            for ( uint r = 0; r < 4; ++r )
            {
                tGram( r, c ) += tP[ r ] * tP[ c ];
            }
            tRhs( c ) += tP[ c ] * tY ;
        }
    }

    Matrix< real > tInverse ;
    real tCondition = 0.0 ;
    EXPECT_TRUE( block_inverse( tGram, tInverse, tCondition ) );

    const real tExpected[ 3 ] = { 1.5, -2.0, 0.5 };
    for ( uint i = 0; i < 3; ++i )
    {
        real tCoefficient = 0.0 ;
        for ( uint k = 0; k < 4; ++k )
        {
            tCoefficient += tInverse( i, k ) * tRhs( k );
        }
        EXPECT_NEAR( tCoefficient, tExpected[ i ], 1.0e-12 );
    }
}

//------------------------------------------------------------------------------

// §3.2: the same one-point patch with a coupling that is not exactly
// representable: the Schur block is zero in exact arithmetic but the
// subtraction leaves round-off, which must be dropped, not inverted
TEST( PostprocessorRecovery, BlockInverseDropsSchurRoundOff )
{
    const real v[ 4 ] = { 1.0, 0.1, 0.0, 0.0 };
    const real w = 0.1 ;
    const real f = 3.7 ;
    Matrix< real > tGram( 4, 4, 0.0 );
    Vector< real > tRhs( 4, 0.0 );
    for ( uint c = 0; c < 4; ++c )
    {
        for ( uint r = 0; r < 4; ++r )
        {
            tGram( r, c ) = w * v[ r ] * v[ c ];
        }
        tRhs( c ) = w * v[ c ] * f ;
    }

    Matrix< real > tInverse ;
    real tCondition = 0.0 ;
    EXPECT_TRUE( block_inverse( tGram, tInverse, tCondition ) );

    real tConstant = 0.0 ;
    for ( uint k = 0; k < 4; ++k )
    {
        tConstant += tInverse( 0, k ) * tRhs( k );
    }
    EXPECT_NEAR( tConstant, f, 1.0e-12 * f );
}
