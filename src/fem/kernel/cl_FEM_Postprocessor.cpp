/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include "cl_FEM_Postprocessor.hpp"
#include "cl_Logger.hpp"
#include "fn_FEM_recovery_inverse.hpp"
#include "fn_trans.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"

namespace belfem
{
    namespace fem
    {
        // singular values of a gradient block below its largest one over this
        // number are dropped from the node fit
        constexpr real gMaxRecoveryCondition = 1e10 ;

        Postprocessor::Postprocessor( Kernel * aKernel, DofManager * aField ) :
            mCommRank( comm_rank() ),
            mCommSize( comm_size() ),
            mNumDimensions( aKernel->mesh()->number_of_dimensions() ),
            mOrder( this->get_interpolation_order( aKernel->mesh() ) ),
            mKernel( aKernel ),
            mMesh( aKernel->mesh() ),
            mField( aField )
        {
            // when we pass a nullpointer as field, we assume that we are using a recovery method
            // otherwise we expect a lest squares implementation in the field
            if ( aField == nullptr )
            {
                this->set_field( aKernel->dofmgr() );
            }
        }

        uint
        Postprocessor::get_interpolation_order( const Mesh * aMesh ) const
        {
            uint aOrder = 0 ;
            if ( mCommRank == 0 )
            {
                aOrder = aMesh->max_element_order();
            }
            comm_barrier() ;
            broadcast( aOrder );
            return aOrder ;
        }

        Postprocessor::~Postprocessor()
        {
            this->delete_node_matrices();
            if ( mNodeBitset != nullptr )    delete mNodeBitset;
            if ( mEdgeBitset != nullptr )    delete mEdgeBitset;
            if ( mFaceBitset != nullptr )    delete mFaceBitset;
            if ( mElementBitset != nullptr ) delete mElementBitset;
        }

        void
        Postprocessor::set_field( DofManager * aField )
        {
            if ( mNodeBitset != nullptr ) delete mNodeBitset;
            if ( mElementBitset != nullptr ) delete mElementBitset;
            mNodeBitset    = new DynamicBitset( mMesh->number_of_nodes() );
            mElementBitset = new DynamicBitset( mMesh->number_of_elements() );

            mField = aField ;
            mEquation = aField->iwg() ;
            this->select_polynomial();
        }

        void
        Postprocessor::set_block_ids( const Vector< id_t > & aBlockIDs )
        {
            mBlockIDs.set_size( aBlockIDs.length(), 0 );
            index_t tCount = 0 ;
            for ( id_t tID : aBlockIDs )
            {
                mBlockIDs( tCount++ ) = tID ;
            }
        }

        void
        Postprocessor::set_block_ids( const Cell< id_t > & aBlockIDs )
        {
            mBlockIDs = aBlockIDs ;
        }

        void
        Postprocessor::initialize()
        {
            if ( mIsInitialized ) return ;
            BELFEM_ASSERT( mField != nullptr, "No field set for this postprocessor" );

            this->select_elements_and_owned_nodes();
            this->select_all_relevant_nodes();

            delete mNodeBitset ;
            mNodeBitset = nullptr ;

            // check if we need edges
            this->check_if_we_have_edges_and_faces();
            if ( mHaveEdges )
            {
                BELFEM_ASSERT( mMesh->edges_exist() , "Mesh does not contain edges" );
                mEdgeBitset = new DynamicBitset( mMesh->number_of_edges() );
                this->select_edges();
                delete mEdgeBitset ;
                mEdgeBitset = nullptr ;
            }
            if ( mHaveFaces )
            {
                BELFEM_ASSERT( mMesh->faces_exist() , "Mesh does not contain faces" );
                mFaceBitset = new DynamicBitset( mMesh->number_of_faces() );
                this->select_faces();
                delete mFaceBitset ;
                mFaceBitset = nullptr ;
            }

            // tidy up
            delete mElementBitset ;
            mElementBitset = nullptr ;

            this->select_polynomial();
            this->compute_node_matrices();

            this->create_target_fields();

            mIsInitialized = true ;
        }

        void
        Postprocessor::synch_source_fields()
        {
            for ( const string & tField : mSourceFields )
            {
                this->synch_source_field( tField );
            }
        }

        void
        Postprocessor::synch_target_fields( Matrix< real > & aData )
        {
            if ( mCommRank > 0 )
            {
                // we only set the target fields
                // on the main proc
                comm_barrier();
                send( aData );
            }
            else
            {
                comm_barrier() ;
                Cell< Matrix< real > > tAllData ;
                collect( tAllData );

                index_t tCol = 0 ;
                for ( const string & tField : mTargetFields )
                {
                    Vector< real > & tTarget = mMesh->field( tField )->data() ;
                    index_t tCount = 0 ;
                    for ( index_t k : mMyOwnedNodeIndices )
                    {
                        tTarget( k ) = aData( tCount++, tCol ) ;
                    }
                    ++tCol ;
                }

                for ( proc_t p = 0; p < mCommSize; ++p )
                {
                    Matrix< real > & tData = tAllData( p ) ;
                    tCol = 0 ;
                    for ( const string & tField : mTargetFields )
                    {
                        Vector< real > & tTarget = mMesh->field( tField )->data() ;
                        Vector< index_t > & tOtherIndices = mAllOwnedNodeIndices( p );
                        index_t tCount = 0 ;
                        for ( index_t k : tOtherIndices )
                        {
                            tTarget( k ) = tData( tCount++, tCol ) ;
                        }
                        ++tCol ;
                    }
                }
            }
        }

        void
        Postprocessor::check_if_we_have_edges_and_faces()
        {
            mHaveEdges = false ;
            mHaveFaces = false ;

            for ( const string & tField : mSourceFields )
            {
                if ( mMesh->field_exists( tField ) )
                {
                    if ( mMesh->field( tField )->entity_type() == EntityType::EDGE )
                    {
                        mHaveEdges = true ;
                        break ;
                    }
                }
            }
            for ( const string & tField : mSourceFields )
            {
                if ( mMesh->field_exists( tField ) )
                {
                    if ( mMesh->field( tField )->entity_type() == EntityType::FACE )
                    {
                        mHaveFaces = true ;
                        break ;
                    }
                }
            }
        }

        void
        Postprocessor::select_elements_and_owned_nodes()
        {
            mElementBitset->reset();
            mMesh->update_element_indices();
            mMesh->unflag_all_elements() ;

            BELFEM_ASSERT( mField != nullptr, "No field set for this postprocessor" );

            // first round: flag all nodes that belong to this post processor
            if ( mBlockIDs.size() > 0 ) // if block types are given directly we use them
            {
                for ( id_t b : mBlockIDs )
                {
                    if ( mField->block_exists( b ) )
                    {
                        Block * tBlock = mField->block( b );

                        // key on the element type, not only the domain: shell
                        // layers without rho are reclassified to Buffer by
                        // create_buffers(), but their elements still need the
                        // facet link for the thin shell edge functions
                        const ElementType tBlockType = tBlock->element_type() ;
                        const bool tIsThinShell =
                               mMesh->block( b )->domain_type() == DomainType::ThinShell
                            || tBlockType == ElementType::PENTA6TS
                            || tBlockType == ElementType::QUAD4TS
                            || tBlockType == ElementType::HEX8TS
                            || tBlockType == ElementType::HEX8TB ;

                        // the second pass adds the aura, so that the recovery
                        // patches of nodes on partition boundaries span the
                        // full disc rather than only the owned side. thin shells
                        // stay side-local: their compute path reads phi from the
                        // master and slave volume elements, which is not
                        // synchronized for aura copies
                        for ( uint a=0; a<2; ++a )
                        {
                            if ( a == 1 && tIsThinShell ) continue ;

                            Cell< Element * > & tFemElements = ( a == 0 ) ?
                                tBlock->elements() : tBlock->aura_elements() ;

                            for ( Element * tFemElement : tFemElements )
                            {
                                mesh::Element * tElement = tFemElement->element();
                                tElement->flag();
                                for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                                {
                                    // claim only the node this element actually references,
                                    // never the cross-side sibling: recovery is side-local
                                    // so that fields may jump at interface duplicates
                                    mesh::Node * tNode = tElement->node( k );
                                    if ( tNode->owner() == mCommRank )
                                    {
                                        mNodeBitset->set( tNode->index() );
                                    }
                                }

                            }
                        }
                    }
                }
            }
            else // if blocks are not given, we go over the domain type
            {
                for ( mesh::Block * tMeshBlock : mMesh->blocks() )
                {
                    if ( mField->block_exists( tMeshBlock->id() ) && tMeshBlock->domain_type() == mDomainType )
                    {
                        mBlockIDs.push( tMeshBlock->id() );

                        Block * tFemBlock = mField->block( tMeshBlock->id() );

                        // element-type check, see comment in the branch above
                        const ElementType tBlockType = tFemBlock->element_type() ;
                        const bool tIsThinShell =
                               tMeshBlock->domain_type() == DomainType::ThinShell
                            || tBlockType == ElementType::PENTA6TS
                            || tBlockType == ElementType::QUAD4TS
                            || tBlockType == ElementType::HEX8TS
                            || tBlockType == ElementType::HEX8TB ;

                        // second pass adds the aura; thin shells stay side-local,
                        // see comment in the branch above
                        for ( uint a=0; a<2; ++a )
                        {
                            if ( a == 1 && tIsThinShell ) continue ;

                            Cell< Element * > & tFemElements = ( a == 0 ) ?
                                tFemBlock->elements() : tFemBlock->aura_elements() ;

                            for ( Element * tFemElement : tFemElements )
                            {
                                mesh::Element * tElement = tFemElement->element() ;
                                tElement->flag();
                                for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                                {
                                    // side-local claim, see comment in the branch above
                                    mesh::Node * tNode = tElement->node( k );
                                    if ( tNode->owner() == mCommRank )
                                    {
                                        mNodeBitset->set( tNode->index() );
                                    }
                                }
                            }
                        }
                    }
                }

                mBlockIDs.shrink_to_fit();
            }

            mNodeBitset->where( mMyOwnedNodeIndices );

            Cell< mesh::Node * > & tNodes = mMesh->nodes();

            // now we select the elements
            for ( id_t k : mMyOwnedNodeIndices )
            {
                mesh::Node * tNode = tNodes( k )->original() ;

                for ( uint e=0; e<tNode->number_of_elements(); ++e )
                {
                    mesh::Element * tElement = tNode->element( e );

                    // check if element is part of the FEM problem
                    if ( tElement->is_flagged() )
                    {
                        mElementBitset->set( tElement->index() );
                    }
                }

                for ( uint d=0; d<tNode->number_of_duplicates(); ++d )
                {
                    mesh::Node * tDup = tNode->duplicate( d );
                    for ( uint e=0; e<tDup->number_of_elements(); ++e )
                    {
                        mesh::Element * tElement = tDup->element( e );

                        // check if element is part of the FEM problem
                        if ( tElement->is_flagged() )
                        {
                            mElementBitset->set( tElement->index() );
                        }
                    }
                }
            }

            mElementBitset->where( mMyElementIndices );

            this->synch_node_indices( true );

        }

        void
        Postprocessor::create_target_fields()
        {
            if ( mCommRank != 0 ) return ;

            for ( const string & tField : mTargetFields )
            {
                if ( ! mMesh->field_exists( tField ) )
                {
                    mMesh->create_field( tField );
                }
            }
        }

        void
        Postprocessor::select_all_relevant_nodes()
        {
            Cell< mesh::Element * > & tElements = mMesh->elements();
            for ( index_t e : mMyElementIndices )
            {
                mesh::Element * tElement = tElements( e );

                for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                {
                    mNodeBitset->set( tElement->node( k )->index() );
                }
            }

            mNodeBitset->where( mMyNodeIndices );
            this->synch_node_indices( false );
        }

        void
        Postprocessor::synch_node_indices( const bool aOwnedOnly )
        {

            Cell< index_t > & tNodeIndices = aOwnedOnly ? mMyOwnedNodeIndices : mMyNodeIndices;

            // now we need to create the communication table for the master proc
            if ( mCommRank == 0 )
            {
                comm_barrier();
                Cell< Vector< id_t > > tNodeIDs ;
                collect( tNodeIDs );

                Cell< Vector< index_t > > & tAllNodeIndices = aOwnedOnly ? mAllOwnedNodeIndices : mAllNodeIndices;

                tAllNodeIndices.set_size( mCommSize, {} );
                for ( proc_t tProc = 1; tProc < mCommSize; ++tProc )
                {
                    Vector< id_t >    & tIDs = tNodeIDs( tProc );
                    Vector< index_t > & tIndices = tAllNodeIndices( tProc );
                    tIndices.set_size( tIDs.length() );

                    index_t tCount = 0 ;
                    for ( id_t tID : tIDs )
                    {
                        tIndices( tCount++ ) = mMesh->node( tID )->index();
                    }
                }
            }
            else
            {
                Cell< mesh::Node * > & tNodes = mMesh->nodes();
                Vector< id_t > tNodeIDs( tNodeIndices.size() );
                index_t tCount = 0 ;
                for ( index_t tIndex : tNodeIndices )
                {
                    tNodeIDs( tCount++ ) = tNodes( tIndex )->id();
                }

                comm_barrier();
                send( tNodeIDs );
            }
        }

        void
        Postprocessor::select_edges()
        {
            mEdgeBitset->reset();
            mMesh->update_edge_indices();

            Cell< mesh::Element * > & tElements = mMesh->elements();

            for ( index_t e : mMyElementIndices )
            {
                mesh::Element * tElement = tElements( e );

                for ( uint d=0; d<tElement->number_of_edges(); ++d )
                {
                    mEdgeBitset->set( tElement->edge( d )->index() );
                }
            }

            mEdgeBitset->where( mMyEdgeIndices );

            if ( mCommRank == 0 )
            {
                comm_barrier();
                Cell< Vector< id_t > > tEdgeIDs ;
                collect( tEdgeIDs );
                mAllEdgeIndices.set_size( mCommSize, {} );
                for ( proc_t p = 1; p < mCommSize; ++p )
                {
                    Vector< id_t > & tIDs = tEdgeIDs( p );
                    mAllEdgeIndices( p ).set_size( tIDs.length() );
                    index_t tCount = 0 ;
                    for ( id_t tID : tIDs )
                    {
                        mAllEdgeIndices( p )( tCount++ ) = mMesh->edge( tID )->index();
                    }
                }
            }
            else
            {
                Vector< id_t > tEdgeIDs( mMyEdgeIndices.size() );
                Cell< mesh::Edge * > & tEdges = mMesh->edges();
                index_t tCount = 0 ;
                for ( index_t tIndex : mMyEdgeIndices )
                {
                    tEdgeIDs( tCount++ ) = tEdges( tIndex )->id();
                }

                comm_barrier();
                send( tEdgeIDs );
            }
        }

        void
        Postprocessor::select_polynomial()
        {
            Matrix< uint > tSelect = { { 3, 6, 10, 15}, { 4, 10, 20, 35} };

            BELFEM_ASSERT( mOrder > 0, "Element order not set" );

            mNumCoefficients = tSelect( mNumDimensions-2, mOrder-1 );

            if ( mNumDimensions == 2 )
            {
                mFunComputePoly = & Postprocessor::compute_poly_2d ;
                switch (  mOrder  )
                {
                    case 1:
                    {
                        mFunPoly2D = & Postprocessor::poly1_2d ;
                        break;
                    }
                    case 2 :
                    {
                        mFunPoly2D = & Postprocessor::poly2_2d ;
                        break;
                    }
                    case 3 :
                    {
                        mFunPoly2D = & Postprocessor::poly3_2d ;
                        break;
                    }
                    case 4 :
                    {
                        mFunPoly2D = & Postprocessor::poly4_2d ;
                        break;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Order %d not implemented", mOrder );
                    }
                }
            }
            else
            {
                mFunComputePoly = & Postprocessor::compute_poly_3d ;
                switch (  mOrder  )
                {
                case 1:
                {
                    mFunPoly3D = & Postprocessor::poly1_3d ;
                    break;
                }
                case 2 :
                {
                    mFunPoly3D = & Postprocessor::poly2_3d ;
                    break;
                }
                case 3 :
                {
                    mFunPoly3D = & Postprocessor::poly3_3d ;
                    break;
                }
                case 4 :
                {
                    mFunPoly3D = & Postprocessor::poly4_3d ;
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Order %d not implemented", mOrder );
                }
                }
            }
            mPoly.set_size( mNumCoefficients, 1 );
            mVandermonde.set_size( mNumCoefficients, mNumCoefficients );
        }

        void
        Postprocessor::select_faces()
        {
            mFaceBitset->reset();
            mMesh->update_face_indices();

            Cell< mesh::Element * > & tElements = mMesh->elements();

            for ( index_t e : mMyElementIndices )
            {
                mesh::Element * tElement = tElements( e );

                for ( uint f=0; f<tElement->number_of_faces(); ++f )
                {
                    mFaceBitset->set( tElement->face( f )->index() );
                }
            }

            mFaceBitset->where( mMyFaceIndices );

            if ( mCommRank == 0 )
            {
                comm_barrier();
                Cell< Vector< id_t > > tFaceIDs ;
                collect( tFaceIDs );
                mAllFaceIndices.set_size( mCommSize, {} );
                for ( proc_t p = 1; p < mCommSize; ++p )
                {
                    Vector< id_t > & tIDs = tFaceIDs( p );
                    Vector< index_t > & tIndices = mAllFaceIndices( p );
                    tIndices.set_size( tIDs.length() );

                    index_t tCount = 0 ;
                    for ( id_t tID : tIDs )
                    {
                        tIndices( tCount++ ) = mMesh->face( tID )->index();
                    }
                }
            }
            else
            {
                Vector< id_t > tFaceIDs( mMyFaceIndices.size() );
                Cell< mesh::Face * > & tFaces = mMesh->faces();
                index_t tCount = 0 ;
                for ( index_t tIndex : mMyFaceIndices )
                {
                    tFaceIDs( tCount++ ) = tFaces( tIndex )->id();
                }

                comm_barrier();
                send( tFaceIDs );
            }
        }

        const Vector< real > &
        Postprocessor::compute( const uint aK )
        {
            BELFEM_ERROR( false, "compute() is not implemented for prostprocessor main class" );
            return mNull;
        }

        void
        Postprocessor::update_element_dofs()
        {
            BELFEM_ERROR( false, "update_element_dofs() is not implemented for prostprocessor main class" );
        }

        void
        Postprocessor::compute_node_matrices()
        {
            if ( mMyNodeIndices.size() == 0 ) return ;

            Cell< mesh::Element * > & tElements = mMesh->elements();
            Cell< mesh::Node * >    & tNodes = mMesh->nodes();

            mElementType = tElements( mMyElementIndices( 0 ) )->type() ;

            // flag selected elements for the is_flagged() checks in the patch loops
            mMesh->unflag_all_elements();
            mMesh->unflag_all_nodes();
            for ( index_t e : mMyElementIndices )
            {
                BELFEM_ASSERT( tElements( e )->type() == mElementType, "Wrong element type (is %s, expect %s)",
                    to_string( tElements( e )->type() ).c_str(), to_string( mElementType ).c_str() );
                tElements( e )->flag( 0 );
            }

            this->select_recovery_mode();

            // allocate node matrices (must match recover_fields() indexing)
            index_t tCount = 0 ;

            for ( index_t i : mMyNodeIndices )
            {
                tNodes( i )->flag( 0 );       // node is selected
                tNodes( i )->set_index( tCount++ );
            }

            mLastElementID = gNoID ;
            mLastBlockID = gNoID ;
            mNumTruncated = 0 ;

            // a second initialize() after a thrown error must not leak the first attempt
            this->delete_node_matrices();

            for ( index_t i : mMyNodeIndices )
            {
                mesh::Node * tNode = tNodes( i );

                // the maps own the objects from the moment they exist
                NodeFrame * tFrame = new NodeFrame();
                mNodeFrames[ tNode ] = tFrame ;
                this->compute_node_frame( tNode, *tFrame );

                mVandermonde.fill( 0.0 );
                this->add_node_gram( tNode, *tFrame );

                Matrix< real > * tInv = new Matrix< real >( mNumCoefficients, mNumCoefficients, 0.0 );
                mNodeMatrices[ tNode ] = tInv ;
                this->invert_node_matrix( *tInv );
            }

            if ( mNumTruncated > 0 )
            {
                message( InfoLevel::Default,
                    "    ... warning: %lu recovery patches on proc %u are near-singular, a fit direction was dropped",
                    ( long unsigned int ) mNumTruncated,
                    ( unsigned int ) mCommRank );
            }

            mMesh->update_node_indices();
        }

//------------------------------------------------------------------------------

        void
        Postprocessor::delete_node_matrices()
        {
            for ( auto tPair : mNodeMatrices )
            {
                delete tPair.second ;
            }
            mNodeMatrices.clear();

            for ( auto tPair : mNodeFrames )
            {
                delete tPair.second ;
            }
            mNodeFrames.clear();
        }

//------------------------------------------------------------------------------

        void
        Postprocessor::select_recovery_mode()
        {
            mShellLocal = mElementType == ElementType::PENTA6TS
                       || mElementType == ElementType::QUAD4TS ;

            const uint tNumNodes = mesh::number_of_nodes( mElementType );

            if ( mShellLocal )
            {
                // ( s1, s2, eta ) in 3D, ( s, eta ) in 2D, whatever the mesh order
                mNumCoefficients = mNumDimensions + 1 ;

                // the bottom facet is the face at zeta = -1 of the reference element
                const uint tBottomFacet = mElementType == ElementType::PENTA6TS ? 3 : 0 ;

                Cell< mesh::Element * > & tElements = mMesh->elements();
                Cell< mesh::Node * > tFacetNodes ;

                mIsBottomNode.set_size( tNumNodes, 0 );

                bool tIsFirst = true ;
                for ( index_t e : mMyElementIndices )
                {
                    mesh::Element * tElement = tElements( e );
                    tElement->get_nodes_of_facet( tBottomFacet, tFacetNodes );

                    for ( uint k=0; k<tNumNodes; ++k )
                    {
                        uint tIsBottom = 0 ;
                        for ( mesh::Node * tFacetNode : tFacetNodes )
                        {
                            if ( tFacetNode == tElement->node( k ) )
                            {
                                tIsBottom = 1 ;
                            }
                        }

                        if ( tIsFirst )
                        {
                            mIsBottomNode( k ) = tIsBottom ;
                        }
                        else
                        {
                            BELFEM_ERROR( mIsBottomNode( k ) == tIsBottom,
                                "thin shell element %lu does not share the face numbering of its block: local node %u",
                                ( long unsigned int ) tElement->id(), ( unsigned int ) k );
                        }
                    }
                    tIsFirst = false ;
                }

                // the thickness coordinate of the reference element is read as
                // the last row of the integration rule
                const Matrix< real > & tXi = mField->block(
                        tElements( mMyElementIndices( 0 ) )->block_id() )->calculator()->integration()->points() ;
                BELFEM_ERROR( tXi.n_rows() == mNumDimensions,
                    "the integration rule of thin shell block %lu has %lu coordinates, expected %u",
                    ( long unsigned int ) tElements( mMyElementIndices( 0 ) )->block_id(),
                    ( long unsigned int ) tXi.n_rows(),
                    ( unsigned int ) mNumDimensions );

                mZeta.set_size( 0 );
                mShellDirection.set_size( mNumDimensions, 0.0 );
            }

            mPoly.set_size( mNumCoefficients, 1 );
            mVandermonde.set_size( mNumCoefficients, mNumCoefficients );
            mSchur.set_size( mNumCoefficients - 1, mNumCoefficients - 1 );
            mCoupling.set_size( mNumCoefficients - 1 );
            mSchurCoupling.set_size( mNumCoefficients - 1 );
            mNodeCoords.set_size( tNumNodes, mNumDimensions );
            mLocalPoint.set_size( mNumDimensions );
            mOffset.set_size( mNumDimensions );
        }

//------------------------------------------------------------------------------

        void
        Postprocessor::link_patch_element( mesh::Element * aElement )
        {
            if ( mLastBlockID != aElement->block_id() )
            {
                mLastBlockID = aElement->block_id();
                mBlock = mField->block( mLastBlockID );
                mMaterial = mBlock->material();
                mCalculator = mBlock->calculator();
            }

            if ( mLastElementID == aElement->id() ) return ;

            mLastElementID = aElement->id();
            mElement = mBlock->element( mLastElementID );
            mCalculator->link( mElement );
            mElement->get_node_coors( mNodeCoords );

            const uint tNumPoints = mCalculator->num_intpoints() ;
            if ( mIntPoints.n_rows() != tNumPoints )
            {
                mIntPoints.set_size( tNumPoints, mNumDimensions );
            }

            for ( uint k=0; k<tNumPoints; ++k )
            {
                mLocalPoint = trans( mNodeCoords ) * mCalculator->Nvec( k );
                for ( uint d=0; d<mNumDimensions; ++d )
                {
                    mIntPoints( k, d ) = mLocalPoint( d );
                }
            }

            if ( mShellLocal )
            {
                // the thickness coordinate of the reference element is its last one
                const Matrix< real > & tXi = mCalculator->integration()->points() ;
                if ( mZeta.length() != tNumPoints )
                {
                    mZeta.set_size( tNumPoints );
                }
                for ( uint k=0; k<tNumPoints; ++k )
                {
                    mZeta( k ) = tXi( mNumDimensions - 1, k );
                }

                // bottom-face centroid to top-face centroid
                mShellDirection.fill( 0.0 );
                for ( uint i=0; i<mNodeCoords.n_rows(); ++i )
                {
                    const real tSign = mIsBottomNode( i ) == 1 ? -1.0 : 1.0 ;
                    for ( uint d=0; d<mNumDimensions; ++d )
                    {
                        mShellDirection( d ) += tSign * mNodeCoords( i, d );
                    }
                }
                mShellDirection *= 2.0 / ( real ) mNodeCoords.n_rows() ;

                const real tGap = norm( mShellDirection );
                BELFEM_ERROR( tGap > 0.0 && std::isfinite( tGap ),
                    "thin shell element %lu: its bottom and top faces coincide",
                    ( long unsigned int ) aElement->id() );

                // the edge functions blend the faces over the prescribed layer
                // thickness, the geometric one only stands in when none is set
                mElementThickness = mMesh->block( aElement->block_id() )->thickness() ;
                if ( std::isnan( mElementThickness ) )
                {
                    mElementThickness = tGap ;
                }

                mShellDirection /= tGap ;
            }
        }

//------------------------------------------------------------------------------

        uint
        Postprocessor::local_node_index( const mesh::Element * aElement, const mesh::Node * aNode ) const
        {
            for ( uint k=0; k<aElement->number_of_nodes(); ++k )
            {
                if ( aElement->node( k ) == aNode )
                {
                    return k ;
                }
            }
            BELFEM_ERROR( false, "node %lu is not a node of element %lu",
                ( long unsigned int ) aNode->id(), ( long unsigned int ) aElement->id() );
            return 0 ;
        }

//------------------------------------------------------------------------------

        real
        Postprocessor::shell_offset( const uint aLocalNode, const uint aK ) const
        {
            // eta runs from the node through the element along the stack:
            // positive into the element above a bottom-face node, negative
            // into the element below a top-face node
            return mIsBottomNode( aLocalNode ) == 1 ?
                     0.5 * ( 1.0 + mZeta( aK ) ) * mElementThickness :
                    -0.5 * ( 1.0 - mZeta( aK ) ) * mElementThickness ;
        }

//------------------------------------------------------------------------------

        void
        Postprocessor::compute_node_frame( mesh::Node * aNode, NodeFrame & aFrame )
        {
            aFrame.Origin.set_size( mNumDimensions );
            for ( uint d=0; d<mNumDimensions; ++d )
            {
                aFrame.Origin( d ) = aNode->x( d );
            }

            if ( mShellLocal )
            {
                aFrame.Normal.set_size( mNumDimensions, 0.0 );
            }

            real tH = 0.0 ;
            real tTau = 0.0 ;
            bool tFinite = true ;

            // the patch is the same union of selected siblings the fit uses
            mesh::Node * tOrg = aNode->original() ;
            for ( uint d=0; d<=tOrg->number_of_duplicates(); ++d )
            {
                mesh::Node * tSibling = d == 0 ? tOrg : tOrg->duplicate( d-1 );
                if ( ! tSibling->is_flagged( 0 ) ) continue ;

                for ( uint e=0; e<tSibling->number_of_elements(); ++e )
                {
                    mesh::Element * tElement = tSibling->element( e );
                    if ( ! tElement->is_flagged() ) continue ;

                    this->link_patch_element( tElement );
                    const uint tLocalNode = this->local_node_index( tElement, tSibling );

                    for ( uint k=0; k<mIntPoints.n_rows(); ++k )
                    {
                        real tDistance = 0.0 ;
                        for ( uint i=0; i<mNumDimensions; ++i )
                        {
                            const real tDelta = mIntPoints( k, i ) - aFrame.Origin( i );
                            tDistance += tDelta * tDelta ;
                        }
                        // std::max keeps the running value against a NaN, so
                        // finiteness is tracked separately
                        tFinite = tFinite && std::isfinite( tDistance );
                        tH = std::max( tH, std::sqrt( tDistance ) );

                        if ( mShellLocal )
                        {
                            const real tOffset = this->shell_offset( tLocalNode, k );
                            tFinite = tFinite && std::isfinite( tOffset );
                            tTau = std::max( tTau, std::abs( tOffset ) );
                        }
                    }

                    if ( mShellLocal )
                    {
                        aFrame.Normal += mShellDirection ;
                    }
                }
            }

            BELFEM_ERROR( tFinite,
                "recovery patch of node %lu has non-finite coordinates",
                ( long unsigned int ) aNode->id() );

            // a degenerate scale is left at one, invert_node_matrix() catches the patch
            aFrame.H   = tH   > 0.0 ? tH   : 1.0 ;
            aFrame.Tau = tTau > 0.0 ? tTau : 1.0 ;

            if ( mShellLocal )
            {
                const real tLength = norm( aFrame.Normal );
                BELFEM_ERROR( tLength > 0.0,
                    "thin shell elements around node %lu point in opposing stack directions",
                    ( long unsigned int ) aNode->id() );
                aFrame.Normal /= tLength ;

                aFrame.T1.set_size( mNumDimensions, 0.0 );
                if ( mNumDimensions == 2 )
                {
                    aFrame.T1( 0 ) = -aFrame.Normal( 1 );
                    aFrame.T1( 1 ) =  aFrame.Normal( 0 );
                }
                else
                {
                    // the axis least aligned with the normal keeps the cross product regular
                    uint tAxis = 0 ;
                    for ( uint i=1; i<3; ++i )
                    {
                        if ( std::abs( aFrame.Normal( i ) ) < std::abs( aFrame.Normal( tAxis ) ) )
                        {
                            tAxis = i ;
                        }
                    }
                    const uint tA = ( tAxis + 1 ) % 3 ;
                    const uint tB = ( tAxis + 2 ) % 3 ;

                    // T1 = e_axis x n
                    aFrame.T1( tA ) = -aFrame.Normal( tB );
                    aFrame.T1( tB ) =  aFrame.Normal( tA );
                    aFrame.T1 /= norm( aFrame.T1 );

                    // T2 = n x T1
                    aFrame.T2.set_size( 3 );
                    aFrame.T2( 0 ) = aFrame.Normal( 1 ) * aFrame.T1( 2 ) - aFrame.Normal( 2 ) * aFrame.T1( 1 );
                    aFrame.T2( 1 ) = aFrame.Normal( 2 ) * aFrame.T1( 0 ) - aFrame.Normal( 0 ) * aFrame.T1( 2 );
                    aFrame.T2( 2 ) = aFrame.Normal( 0 ) * aFrame.T1( 1 ) - aFrame.Normal( 1 ) * aFrame.T1( 0 );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Postprocessor::add_node_gram( mesh::Node * aNode, const NodeFrame & aFrame )
        {
            // weighted least-squares patch fit ( Zienkiewicz & Zhu 1992 ),
            // written in the node's own frame; the shell coordinates ( s1, s2,
            // eta ) are BELFEM's extension, see postprocessor_recovery_theory.md
            const Vector< real > & tVolumes = mMesh->field_data( "_Volumes" );

            mesh::Node * tOrg = aNode->original() ;
            for ( uint d=0; d<=tOrg->number_of_duplicates(); ++d )
            {
                mesh::Node * tSibling = d == 0 ? tOrg : tOrg->duplicate( d-1 );
                if ( ! tSibling->is_flagged( 0 ) ) continue ;

                for ( uint e=0; e<tSibling->number_of_elements(); ++e )
                {
                    mesh::Element * tElement = tSibling->element( e );
                    if ( ! tElement->is_flagged() ) continue ;

                    this->link_patch_element( tElement );
                    const uint tLocalNode = this->local_node_index( tElement, tSibling );

                    // use element volume as weight
                    const real tW = tVolumes( tElement->index() );

                    for ( uint k=0; k<mIntPoints.n_rows(); ++k )
                    {
                        this->compute_local_poly( aFrame, tLocalNode, k );
                        mVandermonde += tW * mPoly * trans( mPoly );
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Postprocessor::compute_local_poly( const NodeFrame & aFrame, const uint aLocalNode, const uint aK )
        {
            for ( uint d=0; d<mNumDimensions; ++d )
            {
                mOffset( d ) = mIntPoints( aK, d ) - aFrame.Origin( d );
            }

            if ( mShellLocal )
            {
                mPoly( 0, 0 ) = 1.0 ;
                mPoly( 1, 0 ) = dot( mOffset, aFrame.T1 ) / aFrame.H ;
                if ( mNumDimensions == 3 )
                {
                    mPoly( 2, 0 ) = dot( mOffset, aFrame.T2 ) / aFrame.H ;
                }
                mPoly( mNumDimensions, 0 ) = this->shell_offset( aLocalNode, aK ) / aFrame.Tau ;
            }
            else
            {
                mOffset /= aFrame.H ;
                this->compute_poly( mOffset );
            }
        }

//------------------------------------------------------------------------------

        void
        Postprocessor::invert_node_matrix( Matrix< real > & aInverse )
        {
            // with a basis centered at the node, G( 0, 0 ) is the patch weight
            const real tWeightSum = mVandermonde( 0, 0 );
            BELFEM_ERROR( tWeightSum > 0.0 && std::isfinite( tWeightSum ),
                "recovery patch with weight sum %g", tWeightSum );

            real tCondition = 0.0 ;
            if ( recovery_block_inverse( mVandermonde, gMaxRecoveryCondition, aInverse, tCondition,
                    mSchur, mCoupling, mSchurCoupling, mSvdU, mSvdVT, mSvdS, mSvdWork ) )
            {
                ++mNumTruncated ;
            }
        }

//------------------------------------------------------------------------------

        void
        Postprocessor::accumulate_node_rhs(
                mesh::Node * aTarget,
                const uint   aLocalNode,
                const real   aWeight,
                Map< mesh::Node *, Matrix< real > * > & aNodeCoeffs )
        {
            BELFEM_ASSERT( mNodeFrames.key_exists( aTarget ),
                "node %lu has no recovery frame", ( long unsigned int ) aTarget->id() );

            const NodeFrame & tFrame = *mNodeFrames( aTarget );

            mNodeRhs.fill( 0.0 );
            for ( uint k=0; k<mIntPoints.n_rows(); ++k )
            {
                this->compute_local_poly( tFrame, aLocalNode, k );

                for ( uint i=0; i<mNumCoefficients; ++i )
                {
                    for ( uint j=0; j<mNumTargetFields; ++j )
                    {
                        mNodeRhs( i, j ) += aWeight * mPoly( i, 0 ) * mPointValues( k, j );
                    }
                }
            }

            if ( aNodeCoeffs.key_exists( aTarget ) )
            {
                *aNodeCoeffs( aTarget ) += mNodeRhs ;
            }
            else
            {
                aNodeCoeffs[ aTarget ] = new Matrix< real >( mNodeRhs );
            }
        }

//------------------------------------------------------------------------------

        void
        Postprocessor::synch_source_field( const string & aField )
        {
            EntityType tType         = mMesh->field( aField )->entity_type();
            Vector< real > & tSource = mMesh->field( aField )->data();

            BELFEM_ASSERT( tType == EntityType::NODE || tType == EntityType::EDGE || tType == EntityType::FACE,
                "Only node, edge and face fields are supported" );

            // the index lists below are ENTITY indices, but an edge or face field
            // does not have to hold one value per entity: DofManager::create_fields
            // sizes it as multiplicity * number_of_entities, and the quadratic
            // Nedelec fields ( "edge_h", "face_h" ) use multiplicity 2, stored as
            // 2*index and 2*index+1 ( see Calculator::nedelec_data_quadratic_* ).
            // Copying one value per index would silently transfer the wrong slot
            // and never transfer the second dof at all, so recover the stride from
            // the field's own length rather than assuming it is one
            const index_t tNumEntities =
                    tType == EntityType::NODE ?  mMesh->number_of_nodes() :
                    tType == EntityType::EDGE ?  mMesh->number_of_edges() :
                                                 mMesh->number_of_faces() ;

            const index_t tMultiplicity = tNumEntities == 0 ?
                    1 : tSource.length() / tNumEntities ;

            BELFEM_ASSERT( tMultiplicity * tNumEntities == tSource.length(),
                "Field '%s' has length %lu, which is not a multiple of the %lu entities it is indexed by",
                aField.c_str(),
                ( long unsigned int ) tSource.length(),
                ( long unsigned int ) tNumEntities );

            if ( mCommRank == 0 )
            {
                Cell< Vector< real > > tAllData( mCommSize, {} );

                for ( proc_t p=1; p<mCommSize; ++p )
                {
                    Vector< real > & tData = tAllData( p );
                    const Vector< index_t > & tSourceIndices =
                        tType == EntityType::NODE ?  mAllNodeIndices( p ) :
                        tType == EntityType::EDGE ?  mAllEdgeIndices( p ) :
                                                     mAllFaceIndices( p ) ;

                    tData.set_size( tMultiplicity * tSourceIndices.length() );

                    index_t tCount = 0 ;
                    for ( index_t tIndex : tSourceIndices )
                    {
                        for ( index_t k=0; k<tMultiplicity; ++k )
                        {
                            tData( tCount++ ) = tSource( tMultiplicity * tIndex + k );
                        }
                    }
                }

                comm_barrier();
                distribute( tAllData );
            }
            else
            {
                comm_barrier();
                Vector< real > tData ;
                receive( tData );

                const Cell< index_t > & tMyIndices =
                     tType == EntityType::NODE ?  mMyNodeIndices :
                     tType == EntityType::EDGE ?  mMyEdgeIndices :
                                                  mMyFaceIndices ;

                BELFEM_ASSERT( tData.length() == tMultiplicity * tMyIndices.size(),
                    "Received %lu values for field '%s', but expected %lu",
                    ( long unsigned int ) tData.length(),
                    aField.c_str(),
                    ( long unsigned int ) ( tMultiplicity * tMyIndices.size() ) );

                index_t tCount = 0 ;
                for ( index_t tIndex : tMyIndices )
                {
                    for ( index_t k=0; k<tMultiplicity; ++k )
                    {
                        tSource( tMultiplicity * tIndex + k ) = tData( tCount++ );
                    }
                }
            }
        }


        void
        Postprocessor::recover_fields()
        {
            this->synch_source_fields();
            if ( mMyOwnedNodeIndices.size() == 0 && mCommRank > 0 )
            {
                // this early return must mirror the communication pattern of
                // the normal path EXACTLY: one barrier and one matrix send,
                // matching root's single synch_target_fields() collect. An extra
                // send is a stray base-tag message. Its unmatched size header
                // truncates the next one-element collect when a rank owns no
                // postprocessor nodes.
                comm_barrier();
                Matrix< real > tNull ;
                send( tNull );
                return ;
            }
            Cell< mesh::Node * >    & tNodes = mMesh->nodes();
            Cell< mesh::Element * > & tElements = mMesh->elements();
            const Vector< real >    & tVolumes = mMesh->field_data( "_Volumes" );

            index_t tCount = 0 ;

            Matrix< real > tMyData( mMyOwnedNodeIndices.size(), mNumTargetFields );

            Map< mesh::Node *, Matrix< real > * > tNodeCoeffs ;

            if ( mNodeRhs.n_rows() != mNumCoefficients || mNodeRhs.n_cols() != mNumTargetFields )
            {
                mNodeRhs.set_size( mNumCoefficients, mNumTargetFields );
            }

            mLastBlockID = gNoID ;
            mLastElementID = gNoID ;

            // phase 1: accumulate element contributions into node coefficient matrices
            for ( index_t e : mMyElementIndices )
            {
                mesh::Element * tElement = tElements( e );

                this->link_patch_element( tElement );

                // update dof vector
                this->update_element_dofs();

                // get the element volume which we will use as weight
                real tWeight = tVolumes( tElement->index() );

                // [DIAG] a recovery-patch weight must be a real positive volume.
                // A NaN/zero here means an aura element reached the patch without a
                // valid _Volumes entry (missed redistribution / unassigned owner) and
                // silently corrupts the boundary fit -> a parallel-only seam. Remove
                // once the slave-surface air seam is localized.
                BELFEM_ERROR( tWeight > 0.0,
                    "[DIAG] recover_fields: element %lu (block %lu, owner %u) contributes to a patch on rank %u with bad _Volumes weight %g",
                    ( long unsigned int ) tElement->id(),
                    ( long unsigned int ) tElement->block_id(),
                    ( unsigned int ) tElement->owner(),
                    ( unsigned int ) mCommRank,
                    tWeight );

                // the target values of every integration point, once per element
                const uint tNumPoints = mIntPoints.n_rows() ;
                if ( mPointValues.n_rows() != tNumPoints || mPointValues.n_cols() != mNumTargetFields )
                {
                    mPointValues.set_size( tNumPoints, mNumTargetFields );
                }
                for ( uint k=0; k<tNumPoints; ++k )
                {
                    const Vector< real > & tY = this->compute( k );
                    for ( uint j=0; j<mNumTargetFields; ++j )
                    {
                        mPointValues( k, j ) = tY( j );
                    }
                }

                // add element contribution to owned nodes. we credit the node this
                // element references, plus registered siblings that are selected in
                // this postprocessor: same-domain pairs ( cohomology cuts ) recover
                // continuous fields from the full disc, while cross-domain siblings
                // ( material interfaces ) are not selected here and may jump. Each
                // target is fit in its own frame; the face the element sees is the
                // one of the node it references
                for ( uint i=0; i<tElement->number_of_nodes(); ++i )
                {
                    mesh::Node * tNode = tElement->node( i );
                    if ( tNode->owner() == mCommRank )
                    {
                        this->accumulate_node_rhs( tNode, i, tWeight, tNodeCoeffs );
                    }

                    // fast path: node is not part of a registered pair
                    mesh::Node * tOrg = tNode->original();
                    if ( tOrg == tNode && tNode->number_of_duplicates() == 0 )
                    {
                        continue ;
                    }

                    if ( tOrg != tNode
                         && tOrg->owner() == mCommRank
                         && mNodeMatrices.key_exists( tOrg ) )
                    {
                        this->accumulate_node_rhs( tOrg, i, tWeight, tNodeCoeffs );
                    }
                    for ( uint d=0; d<tOrg->number_of_duplicates(); ++d )
                    {
                        mesh::Node * tSib = tOrg->duplicate( d );
                        if ( tSib == tNode ) continue ;
                        if (    tSib->owner() == mCommRank
                             && mNodeMatrices.key_exists( tSib ) )
                        {
                            this->accumulate_node_rhs( tSib, i, tWeight, tNodeCoeffs );
                        }
                    }
                }
            }

            // phase 2: reindex owned nodes and compute recovery
            tCount = 0 ;
            for ( index_t tIndex : mMyOwnedNodeIndices )
            {
                tNodes( tIndex )->set_index( tCount++ );
            }

            for ( auto & tPair : tNodeCoeffs )
            {
                mesh::Node * tNode = tPair.first ;
                const Matrix< real > & tRhs = *tPair.second ;
                const Matrix< real > & tInv = *mNodeMatrices( tNode );

                // the local basis is one at the node in slot 0 and zero in every
                // other slot, so the recovered value is row 0 of G^-1 b
                for ( uint j=0; j<mNumTargetFields; ++j )
                {
                    real tValue = 0.0 ;
                    for ( uint i=0; i<mNumCoefficients; ++i )
                    {
                        tValue += tInv( 0, i ) * tRhs( i, j );
                    }
                    tMyData( tNode->index(), j ) = tValue ;
                }

                delete tPair.second ;
            }
            this->synch_target_fields( tMyData );

            // restore node indices
            for ( index_t tIndex : mMyOwnedNodeIndices )
            {
                tNodes( tIndex )->set_index( tIndex );
            }

        }

        void
        Postprocessor::run()
        {
            if ( ! mIsInitialized )
            {

                this->initialize();
            }
            this->recover_fields();
        }

    }
}