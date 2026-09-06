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

#ifndef BELFEM_CL_EF_TESTVOLUME_HPP
#define BELFEM_CL_EF_TESTVOLUME_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "cl_Mesh.hpp"
#include "cl_Node.hpp"
#include "cl_Edge.hpp"
#include "cl_Face.hpp"
#include "cl_Block.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_KernelParameters.hpp"
#include "cl_FEM_DofManagerBase.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_FEM_Element.hpp"

namespace belfem
{
    namespace fem
    {
        namespace test
        {
//------------------------------------------------------------------------------

            /**
             * Single code-built volume element ( TRI3, TRI6, TET4 or TET10 )
             * for the Nedelec edge-function battery.
             *
             * The caller provides the CORNER coordinates ( dim x nCorners,
             * EXODUS node order ); midside nodes of the quadratic elements are
             * placed on the edge midpoints UNLESS aMidOffsets displaces
             * them ( then the quadratic geometry is genuinely curved ). The
             * curved flag is set explicitly from aCurved ( the mesh default
             * is true ) and ONLY selects the evaluation path — it never
             * curves the mesh itself; pair aCurved = true with offsets for
             * a true curved-geometry test.
             *
             * Edge objects are built from the element's own canonical
             * get_nodes_of_edge() order, so all edge directions are positive
             * by construction; aFlipEdge reverses the corner order of one
             * edge to exercise the direction logic. For the TET10, one
             * fixture-owned mesh::Face per element face is inserted with the
             * element as its own master, so the face-ownership branch of
             * EF_TET10::link() resolves to mT = 0.
             *
             * The fem::Element is built with the aura constructor against a
             * minimal real Kernel/DofManagerBase/Block chain, which computes
             * the edge directions the way production does — no seams.
             */
            class EF_TestVolume
            {
                Mesh                   * mMesh ;
                KernelParameters       * mParams ;
                Kernel                 * mKernel ;
                DofManagerBase         * mDofBase ;
                Block                  * mBlock ;

                fem::Element           * mElement = nullptr ;
                Cell< mesh::Edge * >     mEdges ;
                Cell< mesh::Face * >     mFaces ;

                const uint               mDim ;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                EF_TestVolume(
                        const ElementType      aType,
                        const Matrix< real > & aCorners,
                        const int              aFlipEdge = -1,
                        const bool             aCurved   = false,
                        const Matrix< real > * aMidOffsets = nullptr ) :
                        mDim( mesh::dimension( aType ) )
                {
                    const uint tNumCorners = aCorners.n_cols() ;

                    mMesh = new Mesh( mDim, 0, false );

                    mesh::ElementFactory tFactory ;

                    Cell< mesh::Node * > & tNodes = mMesh->nodes();
                    id_t tID = 1 ;

                    // corner nodes
                    for ( uint k = 0; k < tNumCorners; ++k )
                    {
                        tNodes.push( new mesh::Node( tID++,
                                aCorners( 0, k ),
                                mDim > 1 ? aCorners( 1, k ) : 0.0,
                                mDim > 2 ? aCorners( 2, k ) : 0.0 ) );
                    }

                    // midside nodes on the exact edge midpoints, EXODUS order
                    if ( aType == ElementType::TRI6 || aType == ElementType::TET10 )
                    {
                        const uint tNumEdges = ( aType == ElementType::TRI6 ) ? 3 : 6 ;
                        const uint tEdgeCorners[ 6 ][ 2 ] =
                                { { 0, 1 }, { 1, 2 }, { 2, 0 },
                                  { 0, 3 }, { 1, 3 }, { 2, 3 } };

                        for ( uint e = 0; e < tNumEdges; ++e )
                        {
                            const uint a = tEdgeCorners[ e ][ 0 ];
                            const uint b = tEdgeCorners[ e ][ 1 ];

                            // optional offsets ( dim x numEdges ) displace
                            // the midnodes off the midpoints, making the
                            // quadratic geometry GENUINELY curved — the
                            // aCurved flag alone never curves anything
                            real tOx = 0.0, tOy = 0.0, tOz = 0.0 ;
                            if ( aMidOffsets != nullptr )
                            {
                                tOx = ( *aMidOffsets )( 0, e );
                                if ( mDim > 1 ) tOy = ( *aMidOffsets )( 1, e );
                                if ( mDim > 2 ) tOz = ( *aMidOffsets )( 2, e );
                            }

                            tNodes.push( new mesh::Node( tID++,
                                0.5 * ( tNodes( a )->x() + tNodes( b )->x() ) + tOx,
                                0.5 * ( tNodes( a )->y() + tNodes( b )->y() ) + tOy,
                                0.5 * ( tNodes( a )->z() + tNodes( b )->z() ) + tOz ) );
                        }
                    }

                    mesh::Block * tMeshBlock = new mesh::Block( 1, 1 );

                    mesh::Element * tElement = tFactory.create_element( aType, 1 );
                    for ( uint k = 0; k < tNodes.size(); ++k )
                    {
                        tElement->insert_node( tNodes( k ), k );
                    }
                    tElement->set_block_id( 1 );
                    mMesh->elements().push( tElement );
                    tMeshBlock->insert_element( tElement );
                    mMesh->blocks().push( tMeshBlock );

                    mMesh->finalize();

                    // the flag only selects the evaluation path; the
                    // geometry is straight unless aMidOffsets was given
                    if ( aCurved )
                    {
                        tElement->set_curved_flag();
                    }
                    else
                    {
                        tElement->unset_curved_flag();
                    }

                    // minimal real kernel chain
                    mParams  = new KernelParameters( mMesh );
                    mKernel  = new Kernel( mParams );
                    mDofBase = new DofManagerBase( DofManagerType::UNDEFINED, mKernel );
                    mBlock   = new Block( mDofBase, aType );

                    // canonical edges, AFTER the kernel ( the MeshChecker
                    // refuses pre-existing edges )
                    const uint tNumEdges = tElement->number_of_edges();
                    mEdges.set_size( tNumEdges, nullptr );
                    Cell< mesh::Node * > tEdgeNodes ;
                    tElement->allocate_edge_container();
                    for ( uint e = 0; e < tNumEdges; ++e )
                    {
                        tElement->get_nodes_of_edge( e, tEdgeNodes );
                        mesh::Edge * tEdge = new mesh::Edge();
                        tEdge->set_id( e + 1 );
                        tEdge->allocate_node_container( tEdgeNodes.size() );
                        if ( ( int ) e == aFlipEdge )
                        {
                            tEdge->insert_node( tEdgeNodes( 1 ), 0 );
                            tEdge->insert_node( tEdgeNodes( 0 ), 1 );
                        }
                        else
                        {
                            tEdge->insert_node( tEdgeNodes( 0 ), 0 );
                            tEdge->insert_node( tEdgeNodes( 1 ), 1 );
                        }
                        for ( uint k = 2; k < tEdgeNodes.size(); ++k )
                        {
                            tEdge->insert_node( tEdgeNodes( k ), k );
                        }
                        mEdges( e ) = tEdge ;
                        tElement->insert_edge( tEdge, e );
                    }

                    // faces, element is its own master ( mT = 0 )
                    const uint tNumFaces = tElement->number_of_faces();
                    if ( mDim == 3 && tNumFaces > 0 )
                    {
                        mFaces.set_size( tNumFaces, nullptr );
                        tElement->allocate_face_container();
                        for ( uint f = 0; f < tNumFaces; ++f )
                        {
                            mesh::Face * tFace = new mesh::Face(
                                    tElement, f, tElement, f, 0 );
                            tFace->set_id( f + 1 );
                            mFaces( f ) = tFace ;
                            tElement->insert_face( tFace, f );
                        }
                    }

                    mElement = new fem::Element( mBlock, tElement );
                }

//------------------------------------------------------------------------------

                ~EF_TestVolume()
                {
                    delete mElement ;
                    for ( mesh::Face * tFace : mFaces )
                    {
                        delete tFace ;
                    }
                    for ( mesh::Edge * tEdge : mEdges )
                    {
                        delete tEdge ;
                    }
                    delete mBlock ;
                    delete mDofBase ;
                    delete mKernel ;
                    delete mParams ;
                    delete mMesh ;
                }

//------------------------------------------------------------------------------

                fem::Element *
                element() { return mElement ; }

//------------------------------------------------------------------------------

                uint
                dim() const { return mDim ; }
            };

//------------------------------------------------------------------------------
        } /* end namespace test */
    } /* end namespace fem */
} /* end namespace belfem */

#endif // BELFEM_CL_EF_TESTVOLUME_HPP
