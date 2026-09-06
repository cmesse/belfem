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
#include "fn_posv.hpp"
#include "fn_trans.hpp"

namespace belfem
{
    namespace fem
    {
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
            for ( auto tPair : mNodeMatrices )
            {
                delete tPair.second ;
            }
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

            // free some memory
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

            // flag selected elements for the is_flagged() checks in the node loop
            mMesh->unflag_all_elements();
            mMesh->unflag_all_nodes();
            for ( index_t e : mMyElementIndices )
            {
                BELFEM_ASSERT( tElements( e )->type() == mElementType, "Wrong element type (is %s, expect %s)",
                    to_string( tElements( e )->type() ).c_str(), to_string( mElementType ).c_str() );
                tElements( e )->flag( 0 );
            }

            // allocate node matrices (must match recover_fields() indexing)
            index_t tCount = 0 ;

            for ( index_t i : mMyNodeIndices )
            {
                tNodes( i )->flag( 0 );       // node is selected
                tNodes( i )->unflag( 1 );     // not yet processed
                tNodes( i )->set_index( tCount++ );
            }

            Matrix< real > tNodeCoords(
                mesh::number_of_nodes( mElementType ), mNumDimensions );

            mLastElementID = gNoID ;
            mLastBlockID = gNoID ;

            mVandermonde.set_size( mNumCoefficients, mNumCoefficients );
            Map< mesh::Element *, Matrix< real > * > tElementData ;
            Matrix< real > tEye( mNumCoefficients, mNumCoefficients, 0.0 );
            for ( uint i = 0; i < mNumCoefficients; ++i )
            {
                tEye( i, i ) = 1.0 ;
            }

            tCount = 0 ;
            for ( index_t i : mMyNodeIndices )
            {
                Matrix< real > * tInv = new Matrix< real >( mNumCoefficients, mNumCoefficients, 0.0 );
                mesh::Node * tNode = tNodes( i );

                mVandermonde.fill( 0.0 );

                // the patch spans this node and all registered siblings that are
                // selected in this postprocessor ( the is_flagged(0) guard inside
                // compute_element_coeffs ). same-domain pairs ( cohomology cuts )
                // thus recover from the full disc, while cross-domain siblings
                // ( material interfaces ) are not selected here and stay side-local
                this->compute_element_coeffs( tNode->original(), tNodeCoords, tElementData );
                for ( uint d=0; d<tNode->original()->number_of_duplicates(); ++d )
                {
                    this->compute_element_coeffs( tNode->original()->duplicate( d ), tNodeCoords, tElementData );
                }

                // store the inverse matrix via Cholesky factor/solve

                *tInv = tEye ;
                posv( mVandermonde, *tInv );

                mNodeMatrices[ tNode ] = tInv ;

                // mark node as processed
                tNode->flag( 1 );
            }

            // tidy up any remaining element data
            for ( auto tPair : tElementData )
            {
                delete tPair.second ;
            }

            mMesh->update_node_indices();
        }

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
        Postprocessor::compute_element_coeffs(
            mesh::Node * aNode,
            Matrix< real > & aNodeCoords,
            Map< mesh::Element*, Matrix< real > * > & aElementData )
        {
            if ( ! aNode->is_flagged( 0 ) ) return ;

            const Vector< real > & tVolumes = mMesh->field_data( "_Volumes" );

            // loop over all elements of this node
            for ( uint e = 0; e < aNode->number_of_elements(); ++e )
            {
                mesh::Element * tElement = aNode->element( e );

                // skip if element is not selected
                if ( ! tElement->is_flagged() ) continue ;

                // compute element coefficients if they don't exist yet
                this->compute_element_coeffs( tElement, aNodeCoords, aElementData );

                // use element volume as weight
                real tW = tVolumes( tElement->index() );

                for ( uint k = 0; k < mCalculator->num_intpoints(); ++k )
                {
                    mPoly.set_col( 0, aElementData( tElement )->col( k ) );
                    mVandermonde += tW * mPoly * trans( mPoly );
                }

                // delete element coefficients if they are not needed anymore
                this->check_element_coeffs_done( tElement, aElementData );
            }
        }


        void
        Postprocessor::compute_element_coeffs( mesh::Element * aElement,
                Matrix< real > & aNodeCoords,
                Map< mesh::Element* , Matrix< real > * > & aElementData   )
        {

            if ( aElementData.key_exists( aElement ) )
            {
                if ( mLastElementID == aElement->id() ) return ;

                mLastElementID = aElement->id();

                if ( mLastBlockID != aElement->block_id() )
                {
                    mLastBlockID = aElement->block_id();

                    // set mBlock so that compute_element_coeffs can access it
                    mBlock = mField->block( mLastBlockID );
                }

                mElement = mBlock->element( mLastElementID );

                // get the calculator
                mCalculator = mBlock->calculator();

                // link calculator to element
                mCalculator->link( mElement );

                // node coordinates at this element
                mElement->get_node_coors( aNodeCoords );

                return;
            }
            
            mLastElementID = aElement->id();

            if ( mLastBlockID != aElement->block_id() )
            {
                mLastBlockID = aElement->block_id();

                // set mBlock so that compute_element_coeffs can access it
                mBlock = mField->block( mLastBlockID );
            }

            mElement = mBlock->element( mLastElementID );

            // get the calculator
            mCalculator = mBlock->calculator();

            // link calculator to element
            mCalculator->link( mElement );

            // node coordinates at this element
            mElement->get_node_coors( aNodeCoords );

            Vector< real > tPoint ;

            Matrix< real > * tPointer = new Matrix< real >(
                mNumCoefficients, mCalculator->num_intpoints(), 0. );

            aElementData[ aElement ] = tPointer ;

            Matrix< real > & tMatrix = *tPointer ;

            for ( uint k=0; k<mCalculator->num_intpoints(); ++k )
            {
                // compute point coordinates
                tPoint = trans( aNodeCoords ) * mCalculator->Nvec( k );

                // evaluate polynomial
                this->compute_poly( tPoint );

                // store polynomial
                tMatrix.set_col( k, mPoly.col( 0 ) );
            }
        }

        void
        Postprocessor::check_element_coeffs_done( mesh::Element * aElement, Map< mesh::Element* , Matrix< real > *  > & aElementData  )
        {
            if ( ! aElement->is_flagged( 0 ) ) return ;

            // note: we don't check duplicates, periodic nodes, etc.
            //       because it is faster to just compute these few matrices again
            //       rather than checking all nodes all the time
            for ( uint k=0; k<aElement->number_of_nodes(); ++k )
            {
                if ( aElement->node( k )->is_flagged( 0 ) && ! aElement->node( k )->is_flagged( 1 ) ) return ;
            }

            if ( aElementData.key_exists( aElement ) )
            {
                delete aElementData( aElement );
                aElementData.erase_key( aElement );
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
                // matching root's single synch_target_fields() collect. A
                // second [DIAG] dummy lived here from 2026-07-19 to 2026-08-27,
                // claiming to pair with a "patch_count gather" that never
                // existed in the tree; its unmatched 3-index_t size header
                // ( 12 bytes in default builds ) was a stray base-tag
                // message that truncated the next
                // 1-element collect whenever the partition left a rank with no
                // owned postprocessor nodes ( observed on RLC np=4; gantry
                // np=10 showed the matching signature )
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

            for ( mesh::Element * tElement : tElements )
            {
                tElement->flag( 0 );
            }
            for ( index_t e : mMyElementIndices )
            {
                tElements( e )->unflag(  );
            }

            // point coordinates
            Vector< real > tPoint( mNumDimensions );

            // node coordinates
            Matrix< real > tCoords( mesh::number_of_nodes( mElementType ) , mNumDimensions );

            // work vector
            Matrix< real > tB( mNumCoefficients, mNumTargetFields );

            mLastBlockID = gNoID ;
            mLastElementID = gNoID ;

            // phase 1: accumulate element contributions into node coefficient matrices
            for ( index_t e : mMyElementIndices )
            {
                mesh::Element * tElement = tElements( e );

                if ( mLastBlockID != tElement->block_id() )
                {
                    mLastBlockID = tElement->block_id();
                    mBlock = mField->block( mLastBlockID );

                    mMaterial = mBlock->material();

                    // get the calculator
                    mCalculator = mBlock->calculator();
                }

                if ( mLastElementID != tElement->id() )
                {
                    mLastElementID = tElement->id();

                    mElement = mBlock->element( tElement->id() );

                    // link calculator to element
                    mCalculator->link( mElement );

                    // node coordinates at this element
                    mElement->get_node_coors( tCoords );
                }

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

                // loop over all integration points
                tB.fill( 0.0 );
                for ( uint k=0; k<mCalculator->num_intpoints(); ++k )
                {
                    // compute point coordinates
                    tPoint = trans( tCoords ) * mCalculator->Nvec( k );

                    // evaluate polynomial
                    this->compute_poly( tPoint );

                    // compute work matrix
                    const Vector< real > & tY = this->compute( k );

                    for ( uint i=0; i<mNumCoefficients; ++i )
                    {
                        for ( uint j=0; j<mNumTargetFields; ++j )
                        {
                            tB( i, j ) += tWeight * mPoly( i, 0 ) * tY( j );
                        }
                    }
                }

                // add element contribution to owned nodes. we credit the node this
                // element references, plus registered siblings that are selected in
                // this postprocessor: same-domain pairs ( cohomology cuts ) recover
                // continuous fields from the full disc, while cross-domain siblings
                // ( material interfaces ) are not selected here and may jump
                for ( uint i=0; i<tElement->number_of_nodes(); ++i )
                {
                    mesh::Node * tNode = tElement->node( i );
                    if ( tNode->owner() == mCommRank )
                    {
                        if ( tNodeCoeffs.key_exists( tNode ) )
                        {
                            *tNodeCoeffs( tNode ) += tB ;
                        }
                        else
                        {
                            tNodeCoeffs[ tNode ] = new Matrix< real >( tB );
                        }
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
                        if ( tNodeCoeffs.key_exists( tOrg ) )
                        {
                            *tNodeCoeffs( tOrg ) += tB ;
                        }
                        else
                        {
                            tNodeCoeffs[ tOrg ] = new Matrix< real >( tB );
                        }
                    }
                    for ( uint d=0; d<tOrg->number_of_duplicates(); ++d )
                    {
                        mesh::Node * tSib = tOrg->duplicate( d );
                        if ( tSib == tNode ) continue ;
                        if (    tSib->owner() == mCommRank
                             && mNodeMatrices.key_exists( tSib ) )
                        {
                            if ( tNodeCoeffs.key_exists( tSib ) )
                            {
                                *tNodeCoeffs( tSib ) += tB ;
                            }
                            else
                            {
                                tNodeCoeffs[ tSib ] = new Matrix< real >( tB );
                            }
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

            Matrix< real > tPointData( 1, mNumTargetFields );


            for ( auto & tPair : tNodeCoeffs )
            {
                mesh::Node * tNode = tPair.first ;
                Matrix< real > * tMatrix = tPair.second ;

                for ( uint d=0; d<mNumDimensions; ++d )
                {
                    tPoint( d ) = tNode->x( d );
                }
                this->compute_poly( tPoint );

                tPointData = trans( mPoly ) * (*mNodeMatrices( tNode )) * (*tMatrix) ;

                tMyData.set_row( tNode->index(), tPointData.row( 0 ) );

                delete tMatrix ;
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