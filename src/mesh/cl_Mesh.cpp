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

#include "cl_Block.hpp"
#include "cl_Mesh.hpp"
#include "cl_Mesh_ExodusWriter.hpp"
#include "fn_unique.hpp"
#include "stringtools.hpp"
#include "cl_Mesh_GmshReader.hpp"
#include "cl_Mesh_OrderConverter.hpp"
#include "cl_Mesh_CurvedElementChecker.hpp"
#include "cl_Mesh_Partitioner.hpp"
#include "cl_Mesh_BfmFile.hpp"
#include "cl_Mesh_VtkWriter.hpp"
#include "cl_Timer.hpp"
#include "cl_Logger.hpp"
#include "commtools.hpp"
#include "assert.hpp"
#include "cl_DynamicBitset.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_EdgeFactory.hpp"
#include "cl_FaceFactory.hpp"
#include "fn_max.hpp"
#include "fn_cross.hpp"
#include "fn_norm.hpp"
#include "cl_Mesh_ConnectivityCalculator.hpp"
#include "op_Graph_Vertex_Index.hpp"
#include "cl_Mesh_SourceExpander.hpp"
#include "cl_TensorMeshFactory.hpp"
#include "cl_Mesh_BfmFile.hpp"
#include "hdf5_tools.hpp"

namespace belfem
{
    namespace mesh
    {
        uint
        compute_facet_index( Facet * aFacet, Element * aElement, Cell< Node * > & aNodes )
        {
            uint tNumNodes = aFacet->element()->number_of_corner_nodes() ;
            // loop over all facets of element
            for( uint k=0; k<aElement->number_of_facets(); ++k )
            {
                aFacet->unflag_nodes();

                // get nodes from face
                aElement->get_corner_nodes_of_facet( k, aNodes );

                // flag all nodes from this face
                for ( Node * tNode: aNodes )
                {
                    tNode->flag();
                }

                // count how many nodes are flagged
                uint tCount = 0;

                for ( uint i = 0; i < tNumNodes; ++i )
                {
                    if ( aFacet->node( i )->is_flagged() )
                    {
                        // increment node counter
                        ++tCount;
                    }
                }

                if ( tCount == tNumNodes )
                {
                    return k ;
                }
            }

            return BELFEM_UINT_MAX ;
        }
    }

//------------------------------------------------------------------------------

    Mesh::Mesh(
            const uint aNumberOfDimensions ,
            const proc_t aMasterProc,
            const bool aComputeConnectivities  ) :
        mMasterProc(  aMasterProc < comm_size() ? aMasterProc : 0 ),
        mCommRank( comm_rank() ),
        mCommSize( comm_size() ),
        mNumberOfDimensions( aNumberOfDimensions )
    {
        if ( aComputeConnectivities )
            this->set_connectivity( Connectivity::Compute );
    }

//------------------------------------------------------------------------------

    Mesh::Mesh( const uint aOrder,
          const Vector< index_t > aNumNodes,
          const Vector< real > aStep,
          const Vector< real > aOrigin,
          const proc_t aMasterProc ):
        mMasterProc(  aMasterProc < comm_size() ? aMasterProc : 0 ),
        mCommRank( comm_rank() ),
        mCommSize( comm_size() ),
        mNumberOfDimensions( aNumNodes.length() ),
        mTensorConfig( new TensorMeshConfig( aOrder, aNumNodes, aStep, aOrigin ) )
    {
        TensorMeshFactory tFactory ;
        tFactory.populate_tensor_mesh( mTensorConfig, this, mMasterProc );

        // we have not implemented first order B-Splines because it's quite
        // pointless at this time to do so
        if ( aOrder > 1 )
        {
            tFactory.create_bsplines( this );
        }
    }

//------------------------------------------------------------------------------

    Mesh::Mesh(
            const string & aPath,
            const proc_t aMasterProc,
            const bool aComputeConnectivities,
            const bool aParallelMode ) :
            mMasterProc(  aMasterProc < comm_size() ? aMasterProc : 0 ),
            mCommRank( comm_rank() ),
            mCommSize( comm_size() ),
            mPath( aPath )
    {
        if ( aComputeConnectivities )
            this->set_connectivity( Connectivity::Compute );

        if( mCommRank == aMasterProc )
        {
            // start a timer
            Timer tTimer;

            string tType = string_to_lower( filetype( aPath ));


            message( InfoLevel::Default, "    Reading mesh from %s...",
                     filename( aPath ).c_str() );

            if ( tType == "msh" )
            {
                mesh::GmshReader tReader( aPath, this, aComputeConnectivities );
            }
            else if ( tType == "hdf5" || tType == "bfm" )
            {
                mesh::BfmFile tFile( aPath, this );
                tFile.load();
            }
            else
            {
                BELFEM_ERROR( false, "don't know how to read a mesh of type <%s>.",
                             tType.c_str());
            }

            uint tTime = tTimer.stop();

            message( InfoLevel::Default, "    Nodes   : %lu",
                     ( long unsigned int ) this->number_of_nodes());

            message( InfoLevel::Default, "    Elements: %lu",
                     ( long unsigned int ) this->number_of_elements());

            message( InfoLevel::Default, "    Blocks  : %u",
                     ( unsigned int ) this->number_of_blocks());

            message( InfoLevel::Default, "    Sidesets: %u",
                     ( unsigned int ) this->number_of_sidesets());

            message( InfoLevel::Default, "    Fields  : %u",
                     ( unsigned int ) this->number_of_fields() );

            message( InfoLevel::Default, "    Globals : %u",
                     ( unsigned int ) this->number_of_global_variables() );

            message( InfoLevel::Default, "    Time %u ms.\n",
                     ( unsigned int ) tTime );



            if( gComm.size() > 1 && aParallelMode )
            {
                comm_barrier() ;
                broadcast( mNumberOfDimensions );
            }
        }
        else if ( aParallelMode )
        {
            comm_barrier() ;
            broadcast( mNumberOfDimensions );
        }
    }

//------------------------------------------------------------------------------^

    Mesh::~Mesh()
    {
        if ( mPeriodicity != nullptr )
        {
            delete mPeriodicity;
        }

        // delete global variables
        for( auto tVariable: mGlobalVariables )
        {
            delete tVariable;
        }

        for( auto tBlock : mBlocks )
        {
            delete tBlock;
        }

        for( auto tSideSet : mSideSets )
        {
            delete tSideSet;
        }

        for ( auto tCurve : mCurves )
        {
            delete tCurve;
        }

        for ( auto tThinShell : mThinShells )
        {
            delete tThinShell;
        }

        // delete fields
        for( auto tField : mFields )
        {
            delete tField;
        }

        // delete edges
        for( auto tEdge: mEdges )
        {
            delete tEdge ;
        }

        // delete faces
        for( auto tFace: mFaces )
        {
            delete tFace ;
        }

        // facets are deleted by sideset
        //for ( auto tFacet: mFacets )
        //{
        //    delete tFacet;
        //}

        // elements are deleted by block
        //for ( auto tElement: mElements )
        //{
        //    delete tElement;
        //}

        // delete nodes
        for ( auto tNode: mNodes )
        {
            delete tNode;
        }

        // delete bearings
        for( auto tVertex : mVertices )
        {
            delete tVertex;
        }

        // delete control points
        for ( auto tControlPoint : mControlPoints )
        {
            delete tControlPoint;
        }

        if ( mTensorConfig != nullptr )
        {
            delete mTensorConfig;
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::scale_mesh( const real aFactor )
    {
        uint tNumDim = this->number_of_dimensions() ;

        // temporary container for coords, must always be of dimension 3
        Vector< real > tCoords( 3, 0.0 ) ;

        // loop over all nodes
        for( mesh::Node * tNode : mNodes )
        {
            // poulate coordinate vector
            for( uint k=0; k<tNumDim; ++k )
            {
                tCoords( k ) = tNode->x( k );
            }

            // scale vector
            tCoords *= aFactor ;

            // write coords back into node
            tNode->set_coords( tCoords );
        }

        for ( mesh::ControlPoint * tControlPoint : mControlPoints )
        {
            // poulate coordinate vector
            for( uint k=0; k<tNumDim; ++k )
            {
                tCoords( k ) = tControlPoint->x( k );
            }

            // scale vector
            tCoords *= aFactor ;

            // write coords back into control point
            tControlPoint->set_coords(
                tCoords( 0 ),
                tCoords( 1 ),
                tCoords( 2 ) );
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::save( const string & aFilePath )
    {
        if( mCommRank == mMasterProc )
        {
            BELFEM_ASSERT( mIsFinalized, "Can't save an unfinalized mesh");

            // start the timer
            Timer tTimer;

            message( InfoLevel::Default, "\n    Saving Mesh to %s ...", filename( aFilePath ).c_str() );

            string tType = string_to_lower( filetype( aFilePath ));

            if( tType == "hdf5" || tType == "bfm" )
            {
                mesh::BfmFile tFile( aFilePath , this );
                tFile.save();
            }
            else if( tType == "vtk" && mCommRank == mMasterProc )
            {
                // create writer object and save file
                mesh::VtkWriter tWriter( aFilePath, this );
            }
            else if(  tType == "exo" )
            {
                // create writer object
                 mesh::ExodusWriter tWriter( this );

                // write mesh to file
                tWriter.save( aFilePath );
            }
            else if(  tType == "e-s" )
            {
                // create writer object
                mesh::ExodusWriter tWriter( this );

                string tFilePath;

                if( mTimeStep < 10 )
                {
                    tFilePath = sprint("%s.0000%1u", aFilePath.c_str(), ( unsigned int ) mTimeStep );
                }
                else if ( mTimeStep < 100 )
                {
                    tFilePath = sprint("%s.000%2u", aFilePath.c_str(), ( unsigned int )  mTimeStep );
                }
                else if ( mTimeStep < 1000 )
                {
                    tFilePath = sprint("%s.00%3u", aFilePath.c_str(), ( unsigned int )  mTimeStep );
                }
                else if ( mTimeStep < 10000 )
                {
                    tFilePath = sprint("%s.0%4u", aFilePath.c_str(), ( unsigned int ) mTimeStep );
                }
		        else
		        {

                    tFilePath = sprint("%s.%5u", aFilePath.c_str(), ( unsigned int ) mTimeStep );
		        }

                // write mesh to file
                tWriter.save( tFilePath );
            }
            else
            {
                BELFEM_ERROR( false, "Don't know how to write mesh of type: %s", tType.c_str() );
            }

            uint tTime = tTimer.stop() ;

            if( tTime < 1000 )
            {
                message( InfoLevel::Default, "    Time %u ms.\n",
                         ( unsigned int ) tTime );
            }
            else
            {
                message( InfoLevel::Default, "    Time %4.1f s.\n",
                         ( float ) tTime * 0.001 );
            }

        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::update_node_indices()
    {
        index_t tCount = 0;
        for( mesh::Node * tNode : mNodes )
        {
            tNode->set_index( tCount++ );
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::update_edge_indices()
    {
        index_t tCount = 0;
        for( mesh::Edge * tEdge : mEdges)
        {
            tEdge->set_index( tCount++ );
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::update_face_indices()
    {
        index_t tCount = 0;
        for( mesh::Face * tFace : mFaces )
        {
            tFace->set_index( tCount++ );
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::update_facet_indices()
    {
        index_t tCount = 0;
        for( mesh::Facet * tFacet : mFacets )
        {
            tFacet->element()->set_index( tCount );
            tFacet->set_index( tCount++ );
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::update_element_indices()
    {
        index_t tCount = 0;
        for( mesh::Element * Element : mElements )
        {
            Element->set_index( tCount++ );
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::update_vertex_indices()
    {
        index_t tCount = 0;
        for( mesh::Element * Vertex : mVertices )
        {
            Vertex->set_index( tCount++ );
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::update_control_point_indices()
    {
        index_t tCount = 0;
        for( mesh::ControlPoint * tControlPoint : mControlPoints )
        {
            tControlPoint->set_index( tCount++ );
        }
    }
//------------------------------------------------------------------------------

    Vector< real > &
    Mesh::create_field(
            const string & aLabel,
            const EntityType aEntity,
            const id_t aID )
    {
        BELFEM_ASSERT( ! mFieldMap.key_exists( aLabel ),
            "Field %s already exists on mesh.", aLabel.c_str() );


        // create field object
        mesh::Field * tField;

        if ( aID == 0 )
        {
            tField = new mesh::Field(
                    *this,
                    aLabel,
                    mFields.size(),
                    mFields.size()+1,
                    aEntity,
                    FieldType::SCALAR );
        }
        else
        {
            tField = new mesh::Field(
                    *this,
                    aLabel,
                    mFields.size(),
                    aID,
                    aEntity,
                    FieldType::SCALAR );
        }

        // add entry into map
        mFieldMap[ aLabel ] = tField;

        mFields.push( tField );

        // increment field counter
        ++mNumberOfFields;

        // return ref to data object
        return tField->data();
    }

//------------------------------------------------------------------------------

    real &
    Mesh::create_global_variable(
            const string & aLabel,
            const real aValue,
            const id_t aID )
    {
        // increment variable counter
        ++mNumberOfGlobalVariables;

        mesh::GlobalVariable * tVariable;

        if ( aID == 0 )
        {
            // auto set the id
            tVariable = new mesh::GlobalVariable(
                    aLabel,
                    mNumberOfGlobalVariables,
                    aValue );
        }
        else
        {
            tVariable = new mesh::GlobalVariable(
                    aLabel,
                    aID,
                    aValue );
        }

        // add entry into maps. The ID map must key on the ID the variable
        // actually carries: with auto-ID ( aID == 0 ) the argument is NOT
        // the assigned ID, and keying on it filed every auto-created
        // variable under 0, so any later lookup of a real ID aborted
        // ( found 2026-08-15 )
        mGlobalVariableMap[ aLabel ] = tVariable ;
        mGlobalVariableIDMap[ tVariable->id() ] = tVariable;

        // add entry to container
        mGlobalVariables.push( tVariable );

        return tVariable->value();
    }

//------------------------------------------------------------------------------

    void
    Mesh::collect_elements_from_blocks()
    {

        BELFEM_ERROR( mElements.size() == 0,
            "collect_elements_from_blocks() must not be called if elements container is already filled" );

        // initialize counters
        index_t tCount = 0;

        // count number of elements
        for(  mesh::Block * tBlock : mBlocks )
        {
            tCount += tBlock->number_of_elements();
        }

        if( tCount > 0 )
        {
            // allocate element container
            mElements.set_size( tCount, nullptr );

            // reset element counter
            tCount = 0;

            for ( mesh::Block * tBlock: mBlocks )
            {
                index_t tNumElements = tBlock->number_of_elements();

                for ( index_t e = 0; e < tNumElements; ++e )
                {
                    mElements( tCount++ ) = tBlock->element( e );
                }
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::collect_facets_from_sidesets()
    {
        BELFEM_ERROR( mFacets.size() == 0,
                     "collect_facets_from_sidesets() must not be called if facet container is already filled" );

        index_t tCount = 0;

        // count number of sidesets
        for(  mesh::SideSet * tSideSet : mSideSets )
        {
            tCount += tSideSet->number_of_facets();
        }

        // allocate facet container
        mFacets.set_size( tCount, nullptr );

        // reset facet counter
        tCount = 0;

        for(  mesh::SideSet * tSideSet : mSideSets )
        {
            index_t tNumElements = tSideSet->number_of_facets();

            for( index_t f=0; f<tNumElements; ++f )
            {
                mesh::Facet * tFacet = tSideSet->facet_by_index( f );
                mFacets( tCount++ ) = tFacet;
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::update_facet_nodes()
    {
        for ( mesh::SideSet * tSideSet : mSideSets )
        {
            Cell< mesh::Node * > tNodes ;

            if ( tSideSet->number_of_facets() == 0 ) continue;

            // thin-shell facets ( the tape sidesets and the GeometryOnly aggregate
            // built by the ThinShellFactory ) are extrusion geometry: their node
            // lists are restored by the cut factory and must not be re-derived from
            // the master elements, which may carry cut duplicates or abstract nodes
            // where a cut terminates on a tape
            if ( tSideSet->domain_type() == DomainType::ThinShell ||
                 tSideSet->domain_type() == DomainType::GeometryOnly ) continue;

            if ( tSideSet->facet_by_index( 0 )->master() == nullptr ) continue;

            for ( mesh::Facet * tFacet : tSideSet->facets() )
            {
                tFacet->master()->get_nodes_of_facet( tFacet->index_on_master(), tNodes );
                uint tCount=0;
                for ( mesh::Node * tNode : tNodes )
                {
                    tFacet->element()->insert_node( tNode, tCount++ );
                }
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::collect_elements_from_group(
            Cell< mesh::Element * > & aElements,
            const id_t aGroupId,
            const ElementType aType )
    {

        this->unflag_all_elements();

        // count elements on group
        index_t tCount = 0;

        if ( aType == ElementType::EMPTY  )
        {
            // collect all elements from this groip
            for( mesh::Element * tElement : mElements )
            {
                if( tElement->geometry_tag() == aGroupId )
                {
                    tElement->flag();
                    ++tCount ;
                }
            }
        }
        else
        {
            // only collect elements from this type
            for( mesh::Element * tElement : mElements )
            {
                if( tElement->geometry_tag() == aGroupId )
                {
                    if( tElement->type() == aType )
                    {
                        tElement->flag();
                        ++tCount;
                    }
                }
            }
        }

        aElements.set_size( tCount, nullptr );
        tCount = 0;

        for( mesh::Element * tElement : mElements )
        {
            if( tElement->is_flagged() )
            {
                aElements( tCount++ ) = tElement ;
                tElement->unflag();
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::compute_facet_orientations()
    {
        if ( ! mComputeFacetOrientationsWhenFinalizing ) return;

        for( mesh::Facet * tFacet : mFacets )
        {
            tFacet->compute_orientation();
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::finalize()
    {
        // elements are always connected to nodes, even when unfinalized
        this->set_connectivity( Connectivity::ElementToNode );

        if( mCommRank == mMasterProc && ! mIsFinalized )
        {
            if ( mElements.size() == 0 )
            {
                this->collect_elements_from_blocks();
            }

            this->update_node_indices();
            this->update_element_indices();
            this->update_vertex_indices();
            this->update_control_point_indices();

            if ( mFacets.size() == 0 && mSideSets.size() > 0 )
            {
                this->collect_facets_from_sidesets();
            }

            // update indices of facets
            this->update_facet_indices();
            this->update_facet_nodes();

            this->create_maps();

            if( this->test_connectivity( Connectivity::Compute ) )
            {
                mesh::ConnectivityCalculator tCalc( this );

                tCalc.connect_nodes_to_elements() ;

                this->finalize_edges();
                this->finalize_faces();

                if( ! test_connectivity( Connectivity::FacetToElement ) )
                {
                    tCalc.connect_facets_to_elements();
                }

                // make sure that all facets are flagged
                // this is needed for the node to facets stem to work
                for( mesh::Facet * tFacet : mFacets )
                {
                    tFacet->flag() ;
                }

                tCalc.connect_nodes_to_facets() ;

                //tCalc.connect_facets_to_facets();

                tCalc.connect_nodes_to_nodes() ;

                //tCalc.connect_elements_to_elements() ;

                tCalc.connect_thin_shells_to_thin_shells() ;

                for( mesh::SideSet * tSideset : mSideSets )
                {
                    tSideset->collect_nodes() ;
                }

                if ( this->is_tensormesh())
                {
                    tCalc.connect_control_points_to_elements();
                    tCalc.connect_control_points_to_control_points();
                }

            }

            this->unflag_all_elements() ;

            this->set_vertex_owners();
        }

        mIsFinalized = true ;

        this->set_block_ids() ;
        this->set_block_indices() ;
        this->set_sideset_ids() ;
        this->set_sideset_indices() ;

        this->compute_max_element_order();

        this->compute_facet_orientations();

        // Resize node/element fields whose size has drifted from the
        // current entity count (e.g. a caller added duplicate nodes and
        // then re-finalized without unfinalizing first; see
        // resize_and_reset_fields for the guard semantics).
        this->resize_and_reset_fields( EntityType::NODE,    this->number_of_nodes() );
        this->resize_and_reset_fields( EntityType::ELEMENT, this->number_of_elements() );
    }

//------------------------------------------------------------------------------

    void
    Mesh::unfinalize()
    {
        if( mCommRank == mMasterProc && mIsFinalized )
        {
            if ( this->has_periodicity() )
            {
                mPeriodicity->reset_nodes() ;
                mPeriodicity->reset_facets();
            }

            for( mesh::SideSet * tSideset : mSideSets )
            {
                tSideset->reset_node_container();
            }

            // reset the element-to-element links
            // we don't do this because it is expensive to compute
            // we assume that we don't add elements before next finalize
            // ( except for thin shell elements, which are handled differntly )

            /*for ( mesh::Element * tElement : mElements )
            {
                tElement->reset_element_container();
            }
            this->reset_connectivity( Connectivity::ElementToElement );
            this->reset_connectivity( Connectivity::TsElementToTsElement ); */

            // in the same fashion, we assume that we don't add new facets

            /*for ( mesh::Facet * tFacet : mFacets )
            {
                tFacet->reset_facet_container();
            }
            this->reset_connectivity( Connectivity::FacetToFacet );*/

            mElements.clear() ;
            mFacets.clear() ;

            for ( mesh::Node * tVertex : mNodes )
            {
                tVertex->reset_vertex_container();
                tVertex->reset_node_container();
                tVertex->reset_edge_container();
                tVertex->reset_face_container();
                tVertex->reset_element_container();
                tVertex->reset_facet_container();
            }

            this->reset_connectivity( Connectivity::NodeToVertex );
            this->reset_connectivity( Connectivity::NodeToNode );
            this->reset_connectivity( Connectivity::NodeToEdge );
            this->reset_connectivity( Connectivity::NodeToFace );
            this->reset_connectivity( Connectivity::NodeToElement );
            this->reset_connectivity( Connectivity::NodeToFacet );

            for ( mesh::Edge * tVertex : mEdges )
            {
                tVertex->reset_vertex_container();
                //tVertex->reset_node_container();
                tVertex->reset_edge_container();
                tVertex->reset_face_container();
                tVertex->reset_element_container();
                tVertex->reset_facet_container();
            }

            this->reset_connectivity( Connectivity::NodeToVertex );
            this->reset_connectivity( Connectivity::NodeToNode );
            this->reset_connectivity( Connectivity::NodeToEdge );
            this->reset_connectivity( Connectivity::NodeToFace );
            this->reset_connectivity( Connectivity::NodeToElement );
            this->reset_connectivity( Connectivity::NodeToFacet );

            for ( mesh::Face * tVertex : mFaces )
            {
                tVertex->reset_vertex_container();
                //tVertex->reset_node_container();
                tVertex->reset_edge_container();
                tVertex->reset_face_container();
                tVertex->reset_facet_container();
                tVertex->reset_element_container();
            }

            this->reset_connectivity( Connectivity::EdgeToVertex );
            //this->reset_connectivity( Connectivity::EdgeToNode );
            this->reset_connectivity( Connectivity::EdgeToEdge );
            this->reset_connectivity( Connectivity::EdgeToFace );
            this->reset_connectivity( Connectivity::EdgeToElement );
            this->reset_connectivity( Connectivity::EdgeToFacet );
            this->reset_connectivity( Connectivity::EdgeToFacet );
            this->reset_connectivity( Connectivity::FaceToFacet );
            this->reset_maps();
        }

        mIsFinalized = false ;
    }

//------------------------------------------------------------------------------

    void
    Mesh::finalize_edges()
    {
        this->finalize_edges( mElements );
    }

//------------------------------------------------------------------------------

    void
    Mesh::finalize_edges( Cell< mesh::Element * > & aElements )
    {
        if( mEdges.size() > 0 )
        {
            // flag edge elements
            this->unflag_all_elements() ;

            if( this->test_connectivity( Connectivity::Compute ) )
            {
                mesh::ConnectivityCalculator tCalc( this );
                tCalc.connect_nodes_to_edges();
                tCalc.connect_edges_to_elements( aElements );
                //tCalc.connect_edges_to_ghost_facets();
                tCalc.connect_edges_to_edges();
            }
            this->compute_edge_directions();
            this->create_edge_map();

            this->update_edge_indices() ;
        }

        for ( mesh::Block * tBlock : mBlocks )
        {
            tBlock->set_edges_flag( false );
            if ( tBlock->number_of_elements() > 0 )
            {
                if ( tBlock->element( 0 )->has_edges() )
                {
                    tBlock->set_edges_flag();
                }
            }
        }

        this->resize_and_reset_fields( EntityType::EDGE, this->number_of_edges() );

        mEdgesAreFinalized = true ;
    }

//------------------------------------------------------------------------------

    void
    Mesh::finalize_faces()
    {
        if( mFaces.size() > 0 )
        {
            this->create_face_map();
            if ( this->test_connectivity( Connectivity::Compute ) and mEdges.size() > 0 )
            {
                mesh::ConnectivityCalculator tCalc( this );
                tCalc.connect_faces_to_edges_and_edges_to_faces() ;
                //tCalc.connect_faces_to_ghost_facets();
                tCalc.connect_faces_to_faces();
            }
        }

        this->resize_and_reset_fields( EntityType::FACE, this->number_of_faces() );

        mFacesAreFinalized = true ;
    }

//------------------------------------------------------------------------------

    void
    Mesh::resize_and_reset_fields( const EntityType aEntityType, const index_t aSize )
    {
        // Guard: only touch fields whose size does not already match the
        // current entity count. Keeps finalize() idempotent — a no-op
        // second call (no topology change) leaves field values intact —
        // while still catching the case where a caller adds entities and
        // then re-finalizes without an intervening unfinalize() (e.g.
        // CutProcessor::duplicate_nodes() followed by mMesh->finalize() in
        // cl_CutProcessor.cpp ).
        for ( mesh::Field * tField : mFields )
        {
            if ( tField->entity_type() == aEntityType
                 && tField->data().length() != aSize )
            {
                tField->data().set_size( aSize, 0.0 );
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::unflag_everything( const uint aFlagIndex )
    {
        this->unflag_all_nodes( aFlagIndex ) ;
        this->unflag_all_edges( aFlagIndex ) ;
        this->unflag_all_faces( aFlagIndex ) ;
        this->unflag_all_facets( aFlagIndex ) ;
        this->unflag_all_elements( aFlagIndex ) ;
        this->unflag_all_vertices( aFlagIndex ) ;
        this->unflag_all_control_points( aFlagIndex ) ;
    }


//------------------------------------------------------------------------------

    void
    Mesh::unflag_all_nodes( const uint aFlagIndex )
    {
        if ( mIsFinalized )
        {
            for( mesh::Node * tNode : mNodes )
            {
                tNode->unflag( aFlagIndex );
            }
        }
        else
        {
            for ( mesh::Block * tBlock: mBlocks )
            {
                for ( mesh::Element * tElement : tBlock->elements() )
                {
                    for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                    {
                        tElement->node( k )->unflag( aFlagIndex );
                    }
                }
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::unflag_all_edges( const uint aFlagIndex )
    {
        if ( mIsFinalized )
        {
            for( mesh::Edge * tEdge : mEdges )
            {
                tEdge->unflag( aFlagIndex );
            }
        }
        else
        {
            for ( mesh::Block * tBlock: mBlocks )
            {
                for ( mesh::Element * tElement : tBlock->elements() )
                {
                    // edges may exist on a subset of the blocks only
                    // ( e.g. conductors ); air elements carry no containers
                    if ( ! tElement->has_edges() ) continue;

                    for ( uint e=0; e<tElement->number_of_edges(); ++e )
                    {
                        tElement->edge( e )->unflag( aFlagIndex );
                    }
                }
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::unflag_all_faces( const uint aFlagIndex )
    {
        if ( mIsFinalized )
        {
            for( mesh::Face * tFace : mFaces )
            {
                tFace->unflag( aFlagIndex ) ;
            }
        }
        else
        {
            for ( mesh::Block * tBlock: mBlocks )
            {
                for ( mesh::Element * tElement : tBlock->elements() )
                {
                    // faces may exist on a subset of the blocks only
                    if ( ! tElement->has_faces() ) continue;

                    for ( uint f=0; f<tElement->number_of_faces(); ++f )
                    {
                        tElement->face( f )->unflag( aFlagIndex ) ;
                    }
                }
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::unflag_all_facets( const uint aFlagIndex )
    {
        if ( mIsFinalized )
        {
            for( mesh::Facet * tFacet : mFacets )
            {
                tFacet->unflag( aFlagIndex );
            }
        }
        else
        {
            for ( mesh::SideSet * tSideSet : mSideSets )
            {
                Cell< mesh::Facet * > & tFacets = tSideSet->facets();

                for ( mesh::Facet * tFacet : tFacets )
                {
                    tFacet->unflag( aFlagIndex );
                }
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::unflag_all_elements( const uint aFlagIndex )
    {
        if ( mIsFinalized )
        {
            for( mesh::Element * tElement : mElements )
            {
                tElement->unflag( aFlagIndex );
            }
        }
        else
        {
            for ( mesh::Block * tBlock : mBlocks )
            {
                Cell< mesh::Element * > & tElements = tBlock->elements();

                for ( mesh::Element * tElement : tElements )
                {
                    tElement->unflag( aFlagIndex );
                }
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::unflag_all_vertices( const uint aFlagIndex )
    {
        for( mesh::Element * tVertex : mVertices )
        {
            tVertex->unflag( aFlagIndex );
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::unflag_all_control_points( const uint aFlagIndex )
    {
        for( mesh::ControlPoint * tControlPoint : mControlPoints )
        {
            tControlPoint->unflag( aFlagIndex );
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::partition( const uint & aNumberOfPartitions,
                     const bool aSetProcOwners,
                     const bool aForceContinuousPartitions,
                     const bool aResetVertexContainers )
    {
        if( mCommRank == mMasterProc )
        {
            BELFEM_ERROR( aNumberOfPartitions > 1, "Must have more than one partition" );

            // assume that all elements are part of the mesh
            this->unflag_all_elements() ;

            // flag all blocks
            for( mesh::Block * tBlock : mBlocks )
            {
                tBlock->flag_elements() ;
            }

            // create a partitioner
            mesh::Partitioner(
                this,
                aNumberOfPartitions,
                aSetProcOwners,
                aForceContinuousPartitions,
                aResetVertexContainers );

            // assume that all elements are part of the mesh
            this->unflag_all_elements() ;
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::partition( const uint & aNumberOfPartitions,
                     const Vector< id_t > & aSelectedBlocks,
                     const bool aSetProcOwners,
                     const bool aForceContinuousPartitions,
                     const bool aResetVertexContainers )
    {
        if( mCommRank == mMasterProc )
        {
            BELFEM_ERROR( aNumberOfPartitions > 1, "Must have more than one partition" );

            // assume that all elements are part of the mesh
            this->unflag_all_elements() ;

            if( aSelectedBlocks.length() > 0 )
            {
                // flag elements on selected blocks
                for( id_t tID : aSelectedBlocks )
                {
                    this->block( tID )->flag_elements() ;
                }
            }
            else
            {
                // flag all blocks
                for( mesh::Block * tBlock : mBlocks )
                {
                    tBlock->flag_elements() ;
                }
            }

            // create a partitioner
            mesh::Partitioner( this,
                aNumberOfPartitions,
                aSetProcOwners,
                aForceContinuousPartitions,
                aResetVertexContainers );

            // assume that all elements are part of the mesh
            this->unflag_all_elements() ;
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::partition( const uint & aNumberOfPartitions,
                     const Vector< id_t > & aSelectedBlocks,
                     const Vector< id_t > & aSelectedSideSets,
                     const bool aSetProcOwners,
                     const bool aForceContinuousPartitions,
                     const bool aResetVertexContainers )
    {
        if( mCommRank == mMasterProc )
        {
            BELFEM_ERROR( aNumberOfPartitions > 1, "Must have more than one partition" );

            // assume that all elements are part of the mesh
            this->unflag_all_elements() ;

            // flag elements on selected blocks
            for( id_t tID : aSelectedBlocks )
            {
                this->block( tID )->flag_elements() ;
            }

            for( id_t tID : aSelectedSideSets )
            {
                // loop over all facets
                Cell< mesh::Facet * > & tFacets = this->sideset( tID )->facets() ;
                for( mesh::Facet * tFacet : tFacets )
                {
                    if( tFacet->has_master() )
                    {
                        tFacet->master()->flag() ;
                    }
                    if( tFacet->has_slave() )
                    {
                        tFacet->slave()->flag() ;
                    }
                }
            }

            // create a partitioner
            mesh::Partitioner( this,
                               aNumberOfPartitions,
                               aSetProcOwners,
                               aForceContinuousPartitions,
                               aResetVertexContainers );

            // assume that all elements are part of the mesh
            this->unflag_all_elements() ;
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::create_edges(  const bool aPrint,
                         const Vector< id_t > aNedelecBlocks,
                         const Vector< id_t > aNedelecSideSets,
                         const bool aCreateEdgesOnAllSideSet )
    {
        if( mCommRank == mMasterProc )
        {
            // create the factory
            mesh::EdgeFactory tFactory( this );
            tFactory.create_edges( aNedelecBlocks, aNedelecSideSets, aCreateEdgesOnAllSideSet );

            // print for debugging
            if( aPrint )
            {
                tFactory.print();
            }
            this->update_edge_indices() ;

            this->compute_edge_directions() ;
        }
    }

 //------------------------------------------------------------------------------

    void
    Mesh::reset_edges()
    {
        if( mCommRank == mMasterProc && this->edges_exist() )
        {
            if ( this->has_periodicity() )
            {
                mPeriodicity->reset_edges();
            }

            for( mesh::Node * tNode : mNodes )
            {
                tNode->reset_edge_container();
            }
            this->reset_connectivity( Connectivity::NodeToEdge );

            for ( mesh::Block * tBlock : mBlocks )
            {
                for( mesh::Element * tElement : tBlock->elements() )
                {
                    tElement->reset_edge_container() ;
                }
                tBlock->set_edges_flag( false );
            }
            this->reset_connectivity( Connectivity::ElementToEdge );

            for( mesh::SideSet * tSideSet : mSideSets )
            {
                for( mesh::Facet * tFacet : tSideSet->facets() )
                {
                    tFacet->element()->reset_edge_container() ;
                }
            }
            this->reset_connectivity( Connectivity::FacetToEdge );

            for( mesh::Face * tFace : mFaces )
            {
                tFace->reset_edge_container();
            }
            this->reset_connectivity( Connectivity::FaceToEdge );

            mEdgeMap.clear();
            for ( mesh::Edge * tEdge : mEdges )
            {
                delete tEdge ;
            }
            mEdges.clear();
            mHangingEdges.clear();
        }

        mEdgesAreFinalized = false ;
    }

//------------------------------------------------------------------------------

    void
    Mesh::reset_faces()
    {
        if( mCommRank == mMasterProc && this->faces_exist() )
        {
            if ( this->has_periodicity() )
            {
                mPeriodicity->reset_faces();
            }

            for( mesh::Node * tNode : mNodes )
            {
                tNode->reset_face_container();
            }
            this->reset_connectivity( Connectivity::NodeToFace );

            for ( mesh::Block * tBlock : mBlocks )
            {
                for( mesh::Element * tElement : tBlock->elements() )
                {
                    tElement->reset_face_container() ;
                }
                tBlock->set_faces_flag( false );
            }
            this->reset_connectivity( Connectivity::ElementToFace );

            for( mesh::SideSet * tSideSet : mSideSets )
            {
                for( mesh::Facet * tFacet : tSideSet->facets() )
                {
                    tFacet->element()->reset_face_container() ;
                }
            }
            this->reset_connectivity( Connectivity::FacetToFace );

            for( mesh::Edge * tEdge : mEdges )
            {
                tEdge->reset_face_container();
            }
            this->reset_connectivity( Connectivity::EdgeToFace );


            mFaceMap.clear();
            for ( mesh::Face * tFace : mFaces )
            {
                delete tFace ;
            }
            mFaces.clear();
            mHangingFaces.clear();
        }

        mFacesAreFinalized = false ;
    }

//------------------------------------------------------------------------------

    void
    Mesh::create_faces(  const bool aPrint,
                         const Vector< id_t > aNedelecBlocks,
                         const Vector< id_t > aNedelecSideSets )
    {
        if( mCommRank == mMasterProc )
        {


            // create the factory
            mesh::FaceFactory tFactory( this );
            tFactory.create_faces( aNedelecBlocks, aNedelecSideSets ) ;

            // print for debugging
            if( aPrint )
            {
                tFactory.print();
            }
            this->update_face_indices() ;
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::create_maps()
    {
        // reset node map
        mNodeMap.clear();

        // loop over all nodes
        for( mesh::Node * tNode : mNodes )
        {
            // add node to map
            mNodeMap[ tNode->id() ] = tNode;
        }

        // reset element map
        mElementMap.clear();

        // loop over all elements
        for( mesh::Element * tElement : mElements )
        {
            // add element to map
            mElementMap[ tElement->id() ] = tElement;
        }

        // reset the facet map
        mFacetMap.clear() ;

        // loop over all facets
        for( mesh::Facet * tFacet : mFacets )
        {
            // add element to map
            mFacetMap[ tFacet->id() ] = tFacet;
        }

        // reset group maps
        this->update_block_map();
        this->update_sideset_map();

        // reset vertex map
        mVertexMap.clear();

        for( mesh::Element * tVertex : mVertices )
        {
            mVertexMap[ tVertex->id() ] = tVertex;
        }

        // reset control point map
        mControlPointMap.clear();

        for( mesh::ControlPoint * tControlPoint : mControlPoints )
        {
            mControlPointMap[ tControlPoint->id() ] = tControlPoint;
        }

        // reset curve map
        mCurveMap.clear();
        for ( mesh::Curve * tCurve : mCurves )
        {
            mCurveMap[ tCurve->id() ] = tCurve;
        }
    }

//-----------------------------------------------------------------------------

    void
    Mesh::reset_maps()
    {
        mNodeMap.clear();
        mElementMap.clear();
        mFacetMap.clear();
        mSideSetMap.clear();
        mVertexMap.clear();
        mControlPointMap.clear();
    }

//-----------------------------------------------------------------------------

    void
    Mesh::create_edge_map()
    {
        // reset edge map
        mEdgeMap.clear();

        // loop over all nodes
        for( mesh::Edge * tEdge : mEdges )
        {
            // add edge to map
            mEdgeMap[ tEdge->id() ] = tEdge;
        }
    }

//-----------------------------------------------------------------------------

    void
    Mesh::create_face_map()
    {
        // reset edge map
        mFaceMap.clear();

        // loop over all nodes
        for( mesh::Face * tFace : mFaces )
        {
            // add node to map
            mFaceMap[ tFace->id() ] = tFace;
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::compute_max_element_order()
    {
        mMaxElementOrder = 0;
        for ( mesh::Block * tBlock : mBlocks )
        {
            // get order of elements on block
            uint tOrder = mesh::interpolation_order_numeric ( tBlock->element_type() );

            mMaxElementOrder = tOrder > mMaxElementOrder ? tOrder : mMaxElementOrder ;
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::update_ownerships()
    {
        // now we set the ownerships. each node - original or duplicate - is owned
        // by the smallest rank among its own elements, so that the owner always
        // holds elements referencing the node ( side-local field recovery needs this )
        for ( mesh::Node * tNode : mNodes )
        {
            if ( tNode->is_duplicate() ) continue ;

            proc_t tOwner = mCommSize ;
            for ( uint e=0; e<tNode->number_of_elements(); ++e )
            {
                tOwner = tNode->element(e)->owner() < tOwner ? tNode->element(e)->owner() : tOwner ;
            }
            tNode->set_owner( tOwner ) ;

            for ( uint d=0; d<tNode->number_of_duplicates(); ++d )
            {
                mesh::Node * tDup = tNode->duplicate( d ) ;

                proc_t tDupOwner = mCommSize ;
                for ( uint e=0; e<tDup->number_of_elements(); ++e )
                {
                    tDupOwner = tDup->element(e)->owner() < tDupOwner ? tDup->element(e)->owner() : tDupOwner ;
                }

                // fall back to the original's owner if the duplicate has no elements
                tDup->set_owner( tDupOwner < mCommSize ? tDupOwner : tOwner ) ;
            }
        }

        // set edge ownerships
        for ( mesh::Edge * tEdge : mEdges )
        {
            proc_t tOwner =
                    tEdge->node(0)->owner() < tEdge->node(1)->owner() ?
                tEdge->node(0)->owner() : tEdge->node(1)->owner() ;
            tEdge->set_owner( tOwner ) ;
        }

        // set face ownerships
        for ( mesh::Face * tFace : mFaces )
        {
            proc_t tOwner =  mCommSize ;
            for ( uint k=0; k<tFace->number_of_corner_nodes(); ++k )
            {
                tOwner = tFace->node( k )->owner() < tOwner ? tFace->node( k )->owner() : tOwner ;
            }
            tFace->set_owner( tOwner ) ;
        }

        for ( mesh::ControlPoint * tPoint : mControlPoints )
        {
            proc_t tOwner = mCommSize ;
            for ( uint e=0; e<tPoint->number_of_elements(); ++e )
            {
                tOwner = tPoint->element( e )->owner() < tOwner ? tPoint->element( e )->owner() : tOwner ;
            }
            tPoint->set_owner( tOwner ) ;
        }

        for ( mesh::Element * tElement : mVertices )
        {
            tElement->set_owner( tElement->node( 0 )->owner() );
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::set_vertex_owners()
    {

        // loop over all vertices
        for ( mesh::Element * tVertex : mVertices )
        {
            // set owner to owner of node
            tVertex->set_owner( tVertex->node( 0 )->owner() );
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::set_block_ids()
    {
        for ( mesh::Block * tBlock : mBlocks )
        {
            // get id of block
            id_t tID = tBlock->id() ;

            // grab elements of block
            Cell< mesh::Element * > & tElements = tBlock->elements() ;

            // loop over all elements of block
            for( mesh::Element * tElement : tElements )
            {
                tElement->set_block_id( tID );
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::set_sideset_ids()
    {
        for ( mesh::SideSet * tSideSet : mSideSets )
        {
            // get id of block
            id_t tID = tSideSet->id() ;

            // grab elements of block
            Cell< mesh::Facet * > & tFacets = tSideSet->facets() ;

            // loop over all elements of block
            for( mesh::Facet * tFacet : tFacets )
            {
                tFacet->set_sideset_id( tID );
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::set_sideset_indices()
    {
        index_t tCount = 0 ;
        for ( mesh::SideSet * tSideSet : mSideSets )
        {
            tSideSet->set_index( tCount++ );
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::set_block_indices()
    {
        index_t tCount = 0 ;
        for ( mesh::Block * tBlock : mBlocks )
        {
            tBlock->set_index( tCount++ );
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::flag_curved_elements()
    {
        //Timer tTimer ;

        //proc_t tRank = mCommRank ;

        /*if( tRank == 0 )
        {
            message( InfoLevel::Verbose, "\n    Flagging curved elements ...\n"  );
        }*/

        // check curved elements
        mesh::CurvedElementChecker tChecker( mNumberOfDimensions, mBlocks, mSideSets );
        tChecker.flag_curved_elements();

        /*index_t tCount = tChecker.flag_curved_elements();

        if( tRank == 0 )
        {
            message( InfoLevel::Verbose,
                     "    ... time for searching curved elements : %u ms. Elements found: %lu\n",
                     ( unsigned int ) tTimer.stop(),  ( long unsigned int ) tCount  );
        }*/
    }

//------------------------------------------------------------------------------

    void
    Mesh::compute_edge_directions()
    {
        Cell< mesh::Node * > tNodes ;
        for( mesh::Element * tElement : mElements )
        {
            if( ! tElement->has_edges() )
            {
                continue;
            }
            for ( uint e = 0; e < tElement->number_of_edges(); ++e )
            {
                tElement->get_nodes_of_edge( e, tNodes );

                id_t tA = tNodes( 0 )->original()->id();
                id_t tB = tNodes( 1 )->original()->id();

                id_t tC = tElement->edge( e )->node( 0 )->original()->id();
                id_t tD = tElement->edge( e )->node( 1 )->original()->id();

                if ( tA == tC && tB == tD )
                {
                    tElement->set_edge_direction( e, true );
                }
                else if ( tA == tD && tB == tC )
                {
                    tElement->set_edge_direction( e, false );
                }
                else
                {
                    BELFEM_ERROR( false, "Internal error at element %lu: invalid edge %u",
                                  ( long unsigned int ) tElement->id(),
                                  ( unsigned int ) e );
                }
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::distribute_edge_directions()
    {
        proc_t tRank = mCommRank ;
        proc_t tCommSize = comm_size() ;

        if( tRank == 0 )
        {
            Cell< Vector< id_t > > tAllElementIDs( tCommSize, {} );
            Cell< Vector< short unsigned int > > tAllSigns( tCommSize, {} );

            for( proc_t p = 1; p<tCommSize; ++p )
            {

                index_t tElementCount = 0 ;
                index_t tEdgeCount = 0 ;

                // count elements
                for( mesh::Element * tElement : mElements )
                {
                    if( tElement->owner() == p && tElement->has_edges() )
                    {
                        tElementCount += 1 ;
                        tEdgeCount += tElement->number_of_edges() ;
                    }
                }

                Vector< id_t > & tElementIDs = tAllElementIDs( p );
                tElementIDs.set_size( tElementCount );

                tElementCount = 0 ;

                Vector< short unsigned int > & tSigns = tAllSigns( p );
                tSigns.set_size( tEdgeCount, 0 );

                tEdgeCount = 0 ;

                for( mesh::Element * tElement : mElements )
                {
                    if( tElement->owner() == p && tElement->has_edges() )
                    {
                        tElementIDs( tElementCount++ ) = tElement->id() ;
                        for( uint e=0; e<tElement->number_of_edges(); ++e )
                        {
                            if( tElement->edge_direction( e ) )
                            {
                                tSigns( tEdgeCount ) = 1 ;
                            }
                            tEdgeCount++ ;
                        }
                    }
                }
            }

            comm_barrier() ;

            distribute( tAllElementIDs );
            distribute( tAllSigns  );

            comm_barrier() ;
        }
        else
        {
            comm_barrier() ;

            Vector< id_t > tElementIDs ;
            Vector< short unsigned int > tSigns ;

            receive( tElementIDs );
            receive( tSigns );

            index_t tCount = 0 ;

            for( id_t tID : tElementIDs )
            {
                mesh::Element * tElement = this->element( tID );

                for( uint e=0; e<tElement->number_of_edges(); ++e )
                {
                    tElement->set_edge_direction( e, tSigns( tCount++) == 1 );
                }
            }
            comm_barrier() ;
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::save_faces( const string & aPath )
    {
        if( mCommRank == 0 )
        {
            // flag the nodes that  we need
            this->unflag_all_nodes() ;

            Cell< ElementType > tTypes( mFaces.size(), ElementType::UNDEFINED );

            index_t tCount = 0 ;

            for( mesh::Face * tFace : mFaces )
            {
                for( uint k=0; k<tFace->number_of_nodes(); ++k )
                {
                    tFace->node( k )->flag() ;
                }
                tTypes( tCount++ ) = mesh::element_type_from_numnodes( 2, tFace->number_of_nodes() );
            }
            unique( tTypes );


            // count flagged nodes
            tCount = 0 ;
            for( mesh::Node * tNode : mNodes )
            {
                if( tNode->is_flagged() )
                {
                    tNode->set_index( tCount++ );
                }
            }


            Mesh * tNewMesh = new Mesh( 3, 2, false );

            Cell< mesh::Node * > & tNodes = tNewMesh->nodes() ;
            tNodes.set_size( tCount, nullptr );

            Map< ElementType, id_t > tTypeMap ;
            tCount = 0 ;
            for( ElementType tType : tTypes )
            {
                tTypeMap[ tType ] = tCount++ ;
            }

            Cell< mesh::Block * > & tBlocks = tNewMesh->blocks() ;
            Vector< index_t > tNumElemsPerBlock( tCount, 0 );
            for( mesh::Face * tFace : mFaces )
            {
                ++tNumElemsPerBlock( tTypeMap( mesh::element_type_from_numnodes( 2, tFace->number_of_nodes() ) ) );
            }


            tBlocks.set_size( tCount, nullptr );

            for( uint b=0; b<tCount; ++b )
            {
                tBlocks( b ) = new mesh::Block( b+1, tNumElemsPerBlock(  b ) );
            }


            tCount = 0 ;

            // create copies of the nodes
            for( mesh::Node * tNode : mNodes )
            {
                if( tNode->is_flagged() )
                {
                    tNode->set_index( tCount );

                    mesh::Node * tNewNode = new mesh::Node(
                            tNode->id(),
                            tNode->x(),
                            tNode->y(),
                            tNode->z() );

                    tNodes( tCount ) = tNewNode ;

                    tNewNode->set_index( tCount++ );
                }
            }

            // create elements
            //Cell< mesh::Element * > & tElements = tNewMesh->elements() ;
            //tElements.set_size( mFaces.size(), nullptr );



            tCount = 0 ;

            mesh::ElementFactory tFactory ;

            Vector< real > tA( 3 );
            Vector< real > tB( 3 );
            Vector< real > tN( 3 );

            tNumElemsPerBlock.fill( 0 );

            Vector< real > & tX = tNewMesh->create_field( "nx", EntityType::ELEMENT );
            Vector< real > & tY = tNewMesh->create_field( "ny", EntityType::ELEMENT );
            Vector< real > & tZ = tNewMesh->create_field( "nz", EntityType::ELEMENT );

            tX.set_size( mFaces.size(), 0 );
            tY.set_size( mFaces.size(), 0 );
            tZ.set_size( mFaces.size(), 0 );

            for( mesh::Face * tFace : mFaces )
            {
               // determine type of face
               ElementType tType = mesh::element_type_from_numnodes( 2, tFace->number_of_nodes() );

               // create a new element
               mesh::Element * tElement = tFactory.create_element( tType, tFace->id() );

               // link nodes
               for( uint k=0; k<tFace->number_of_nodes(); ++k )
               {
                   tElement->insert_node( tNodes( tFace->node( k )->index() ), k );
               }

               tElement->set_block_id( tTypeMap( tType ) );



               if( tType == ElementType::TRI3 )
               {
                   tA( 0 ) = tElement->node( 1 )->x() - tElement->node( 0 )->x() ;
                   tA( 1 ) = tElement->node( 1 )->y() - tElement->node( 0 )->y() ;
                   tA( 2 ) = tElement->node( 1 )->z() - tElement->node( 0 )->z() ;

                   tB( 0 ) = tElement->node( 2 )->x() - tElement->node( 0 )->x() ;
                   tB( 1 ) = tElement->node( 2 )->y() - tElement->node( 0 )->y() ;
                   tB( 2 ) = tElement->node( 2 )->z() - tElement->node( 0 )->z() ;

                   tN = cross( tA, tB );
                   tN /= norm( tN );

                   tX( tCount ) = tN( 0 );
                   tY( tCount ) = tN( 1 );
                   tZ( tCount ) = tN( 2 );
               }

                index_t b = tTypeMap( tType );

                tBlocks( b )->elements()( tNumElemsPerBlock( b )++) = tElement ;
                ++tCount ;
            }

            // restore node indices
            tCount = 0 ;
            for( mesh::Node * tNode : mNodes )
            {
                tNode->set_index( tCount++ );
            }

            tNewMesh->finalize() ;
            tNewMesh->save( aPath );

            delete tNewMesh ;
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::collect_hanging_basis()
    {
        mesh::collect_hanging_basis( mNodes, mHangingNodes );
        mesh::collect_hanging_basis( mEdges, mHangingEdges );
        mesh::collect_hanging_basis( mFaces, mHangingFaces );

        mesh::collect_hanging_basis( mFacets, mHangingFacets );
        mesh::collect_hanging_basis( mElements, mHangingElements );
        mesh::collect_hanging_basis( mControlPoints, mHangingControlPoints );
    }

    void
    Mesh::expand_hanging_basis_sources()
    {
        // some sources might be dependent.
        // the expander checks that and flattens the sources out
        mesh::SourceExpander tExpander( this );
        for ( mesh::Node * tNode : mHangingNodes )
        {
            tExpander.expand_sources( tNode );
        }
        // note, we are deliberately not doing this for edges and faces
        // because they need special treatment for the maxwell case
        /*for ( mesh::Edge * tEdge : mHangingEdges )
        {
            tExpander.expand_sources( tEdge );
        }
        for ( mesh::Face * tFace : mHangingFaces )
        {
            tExpander.expand_sources( tFace );
        }*/
        for ( mesh::Element * tElement : mHangingElements )
        {
            tExpander.expand_sources( tElement );
        }
        for ( mesh::Facet * tFacet : mHangingFacets )
        {
            tExpander.expand_sources( tFacet );
        }
        for ( mesh::ControlPoint * tControlPoint : mHangingControlPoints )
        {
            tExpander.expand_sources( tControlPoint );
        }
    }

//------------------------------------------------------------------------------

    id_t
    Mesh::max_node_id()
    {
        id_t aID = 0 ;
        for( mesh::Node * tNode : mNodes )
        {
            if( tNode->id() > aID )
            {
                aID = tNode->id() ;
            }
        }
        return aID ;
    }

//------------------------------------------------------------------------------

    id_t
    Mesh::max_element_id()
    {
        id_t aID = 0 ;

        for( mesh::Edge * tEdge : mEdges )
        {
            if( tEdge->id() > aID )
            {
                aID = tEdge->id() ;
            }
        }
        for ( mesh::Face *tFace : mFaces )
        {
            if ( tFace->id() > aID )
            {
                aID = tFace->id();
            }
        }

        for( mesh::Block * tBlock : mBlocks )
        {
            for( mesh::Element * tElement : tBlock->elements() )
            {
                if( tElement->id() > aID )
                {
                    aID = tElement->id() ;
                }
            }
        }
        for( mesh::SideSet * tSideSet : mSideSets )
        {
            for( mesh::Facet * tFacet : tSideSet->facets() )
            {
                if( tFacet->id() > aID )
                {
                    aID = tFacet->id() ;
                }
            }
        }
        for ( mesh::Curve * tCurve : mCurves )
        {
            for ( mesh::Segment * tSegment : tCurve->segments() )
            {
                if ( tSegment->id() > aID )
                {
                    aID = tSegment->id() ;
                }
            }
        }
        return aID ;
    }

//------------------------------------------------------------------------------

    id_t
    Mesh::max_block_and_sideset_id()
    {
        id_t aID = 0 ;
        for( mesh::Block * tBlock : mBlocks )
        {
            if( tBlock->id() > aID )
            {
                aID = tBlock->id() ;
            }
        }
        for( mesh::SideSet * tSideSet : mSideSets )
        {
            if( tSideSet->id() > aID )
            {
                aID = tSideSet->id() ;
            }
        }
        for ( mesh::Curve * tCurve : mCurves )
        {
            if ( tCurve->id() > aID )
            {
                aID = tCurve->id() ;
            }
        }
        for ( mesh::ThinShell * tThinShell : mThinShells )
        {
            if ( tThinShell->id() > aID )
            {
                aID = tThinShell->id();
            }
        }
        return aID ;
    }

//------------------------------------------------------------------------------

    void
    Mesh::add_block( mesh::Block * aBlock )
    {
        mBlocks.push( aBlock );
        mBlockMap[ aBlock->id() ] = aBlock ;
    }

//------------------------------------------------------------------------------

    void
    Mesh::add_sideset( mesh::SideSet * aSideSet )
    {
        mSideSets.push( aSideSet );
        mSideSetMap[ aSideSet->id() ] = aSideSet ;

        // also register the facets, so that facet lookups work before the
        // next finalize. We deliberately do NOT push them into mFacets here:
        // finalize() only collects the sidesets when mFacets is empty, and a
        // partially filled container would suppress that collection. The
        // next finalize rebuilds this map from the recollected mFacets,
        // which then contains these facets as well
        for ( mesh::Facet * tFacet : aSideSet->facets() )
        {
            mFacetMap[ tFacet->id() ] = tFacet ;
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::set_abstract_nodes( Cell< mesh::Node * > & aNodes )
    {
        BELFEM_ASSERT( mAbstractNodes.size() == 0, "Abstract nodes have already been set");

        index_t k = mNodes.size() ;

        index_t j = 0 ;
        mAbstractNodes.set_size( aNodes.size(), nullptr );

        for ( mesh::Node * tNode : aNodes )
        {
            if( ! mNodeMap.key_exists( tNode->id() ) )
            {
                j++ ;
            }
        }

        mNodes.vector_data().reserve( k + j );

        j = 0 ;

        for ( mesh::Node * tNode : aNodes )
        {
            mAbstractNodes( j++ ) = tNode ;
            if( ! mNodeMap.key_exists( tNode->id() ) )
            {
                tNode->set_index( k++ );
                mNodes.push( tNode );
                mNodeMap[ tNode->id() ] = tNode ;
            }
        }

        // after adjusting the nodes, we need to make sure that the
        // field sizes are correct. Saving the information from the
        // original node fields should not be necessary, but we
        // do it anyways
        for( mesh::Field * tField : mFields )
        {
            if( tField->entity_type() == EntityType::NODE )
            {
                if( tField->data().length() != mNodes.size() )
                {
                    Vector< real > tData( tField->data() );

                    tField->data().set_size( mNodes.size(), 0 );
                    for( index_t i=0; i<tData.length() ; ++i )
                    {
                        tField->data()( i ) = tData( i );
                    }
                }
            }
        }
    }

//------------------------------------------------------------------------------


    // compute and return the checksum
    std::size_t
    Mesh::checksum()
    {
        if ( mHash.value() == 0 )
        {
            this->compute_checksum();
        }
        return mHash.value();
    }

//------------------------------------------------------------------------------

    void
    Mesh::force_checksum( const std::size_t aChecksum )
    {
        mHash.set_value( aChecksum );
    }

//------------------------------------------------------------------------------

    void
    Mesh::set_config_tag( const uint64_t aTag, const string & aText )
    {
        mConfigTag  = aTag ;
        mConfigText = aText ;
    }

//------------------------------------------------------------------------------

    uint64_t
    Mesh::config_tag() const
    {
        return mConfigTag ;
    }

//------------------------------------------------------------------------------

    const string &
    Mesh::config_text() const
    {
        return mConfigText ;
    }

//------------------------------------------------------------------------------

    void
    Mesh::compute_checksum()
    {
        if ( mCommRank != this->master() ) return;

        mHash.reset();

        mHash += this->number_of_dimensions() ;

        // add nodes
        mHash += this->number_of_nodes() ;

        if ( this->number_of_dimensions() == 2 )
        {
            for ( mesh::Node * tNode : mNodes )
            {
                mHash += tNode->x() ;
                mHash += tNode->y() ;
            }
        }
        else if ( this->number_of_dimensions() == 3 )
        {
            for ( mesh::Node * tNode : mNodes )
            {
                mHash += tNode->x() ;
                mHash += tNode->y() ;
                mHash += tNode->z() ;
            }
        }

        // add elements
        mHash += this->number_of_elements() ;

        for ( mesh::Element * tElement : mElements )
        {
            for ( uint k=0; k<tElement->number_of_nodes(); ++k )
            {
                mHash += tElement->node( k )->id() ;
            }
        }
    }

    void
    Mesh::update_element_map()
    {
        mElementMap.clear();
        for ( mesh::Block * tBlock : mBlocks )
        {
            for ( mesh::Element * tElement : tBlock->elements() )
            {
                mElementMap[ tElement->id() ] = tElement;
            }
        }
    }

    void
    Mesh::update_block_map()
    {
        mBlockMap.clear();
        for ( mesh::Block * tBlock : mBlocks )
        {
            mBlockMap[ tBlock->id() ] = tBlock ;
        }
    }

    void
    Mesh::update_sideset_map()
    {
        mSideSetMap.clear();
        for ( mesh::SideSet * tSideSet : mSideSets )
        {
            mSideSetMap[ tSideSet->id() ] = tSideSet ;
        }
    }

    id_t
    Mesh::max_node_id() const
    {
        id_t aID = 0 ;
        for( mesh::Node * tNode : mNodes )
        {
            if( tNode->id() > aID )
            {
                aID = tNode->id() ;
            }
        }
        return aID ;
    }

    id_t
    Mesh::max_element_id() const
    {
        id_t aID = 0 ;
        for ( mesh::Block * tBlock : mBlocks )
        {
            for ( mesh::Element * tElement : tBlock->elements() )
            {
                if ( tElement->id() > aID )
                {
                    aID = tElement->id() ;
                }
            }
        }
        for ( mesh::SideSet * tSideSet : mSideSets )
        {
            for ( mesh::Facet * tFacet : tSideSet->facets() )
            {
                if ( tFacet->id() > aID )
                {
                    aID = tFacet->id() ;
                }
            }
        }
        for ( mesh::Curve * tCurve : mCurves )
        {
            for ( mesh::Segment * tSegment : tCurve->segments() )
            {
                if ( tSegment->id() > aID )
                {
                    aID = tSegment->id() ;
                }
            }
        }
        for ( mesh::Edge * tEdge : mEdges )
        {
            if ( tEdge->id() > aID )
            {
                aID = tEdge->id() ;
            }
        }
        for ( mesh::Face * tFace : mFaces )
        {
            if ( tFace->id() > aID )
            {
                aID = tFace->id() ;
            }
        }
        for ( mesh::Element * tElement : mBoundaryEdges )
        {
            if ( tElement->id() > aID )
            {
                aID = tElement->id() ;
            }
        }
        for ( mesh::Element * tElement : mVertices )
        {
            if ( tElement->id() > aID )
            {
                aID = tElement->id() ;
            }
        }
        return aID ;
    }

    Mesh *
    Mesh::extract_thin_shell_mesh()
    {
        BELFEM_ERROR( this->thin_shells().size() > 0 , "Can't extract thin shells of a mesh that doesn't have any." );

        // make sure that data is healthy
        this->update_node_indices() ;

        DynamicBitset tNodeFlags( this->number_of_nodes() ) ;

        // count elements and flag nodes
        index_t tNumBlocks = 0 ;

        for ( mesh::ThinShell * tShell : this->thin_shells() )
        {
            tNumBlocks += tShell->blocks().size() ;

            for ( mesh::Block * tBlock : tShell->blocks() )
            {
                for ( mesh::Element * tElement : tBlock->elements() )
                {
                    for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                    {
                        tNodeFlags.set( tElement->node( k )->index() );
                    }
                }
            }
        }

        // creating the output mesh
        Mesh * aMesh = new Mesh( this->number_of_dimensions(), this->master() ) ;

        // grab the nodes indices that are relevant
        Cell< index_t > tIndices ;
        tNodeFlags.where( tIndices );
        key_t tNumNodes = tIndices.size() ;

        // temporary node map for output mesh
        Map< id_t, mesh::Node * > tNodeMap ;

        Cell< mesh::Node * > & tNodes = aMesh->nodes() ;
        tNodes.set_size( tNumNodes, nullptr );
        tNumNodes = 0 ;
        for ( index_t tIndex : tIndices )
        {
            mesh::Node * tOrg = mNodes( tIndex ) ;
            mesh::Node * tDup = new mesh::Node( tOrg->id(), tOrg->x(), tOrg->y(), tOrg->z() ) ;
            tNodeMap[ tOrg->id() ] = tDup ;

            tNodes( tNumNodes++ ) = tDup ;
        }

        aMesh->blocks().reserve( tNumBlocks );

        mesh::ElementFactory tFactory ;

        for ( mesh::ThinShell * tShell : this->thin_shells() )
        {
            for ( mesh::Block * tBlock : tShell->blocks() )
            {
                ElementType tType = ElementType::EMPTY ;
                if ( tBlock->element_type() == ElementType::PENTA6TS )
                {
                    tType = ElementType::PENTA6 ;
                }
                else if ( tBlock->element_type() == ElementType::PENTA18TS )
                {
                    tType = ElementType::PENTA18 ;
                }
                else if (tBlock->element_type() == ElementType::QUAD4TS)
                {
                    tType = ElementType::QUAD4 ;
                }
                else if (tBlock->element_type() == ElementType::QUAD9TS)
                {
                    tType = ElementType::QUAD9 ;
                }
                else
                {
                    BELFEM_ERROR( false, "Invalid element type");
                }

                mesh::Block * tNewBlock = new mesh::Block( tBlock->id(), tBlock->number_of_elements() );

                for ( mesh::Element * tOrg : tBlock->elements() )
                {
                    mesh::Element * tDup = tFactory.create_element( tType, tOrg->id() );
                    for ( uint k=0; k<tOrg->number_of_nodes(); ++k )
                    {
                        tDup->insert_node( tNodeMap( tOrg->node( k )->id() ), k );
                    }
                    tNewBlock->insert_element( tDup );
                }
                aMesh->blocks().push( tNewBlock );
            }
        }

        aMesh->finalize() ;
        return aMesh ;
    }

    void
    Mesh::populate_element_neighbors()
    {
        this->update_node_indices() ;
        key128_t tNumNodes = this->number_of_nodes() ;

        // special case if faces have been created
        if ( this->number_of_dimensions() == 3 && this->faces_exist() )
        {
            for ( mesh::Block * tBlock : mBlocks )
            {
                for ( mesh::Element * tElement : tBlock->elements() )
                {
                    tElement->allocate_neighbor_container() ;
                    for ( uint f=0; f<tElement->number_of_faces(); ++f )
                    {
                        mesh::Face * tFace = tElement->face( f );

                        if ( tFace->slave() == nullptr ) continue ;

                        tElement->insert_neighbor(
                            tFace->master()->id() == tElement->id() ? tFace->slave() : tFace->master() , f );
                    }
                }
            }
            return ;
        }

        // count possible keys
        index_t tCount = 0 ;
        Cell< key128_t > tKeys ;

        key128_t tA ;
        key128_t tB ;
        key128_t tC ;
        Cell< mesh::Node * > tWork ;
        for ( mesh::Block * tBlock : mBlocks )
        {
            tCount += tBlock->number_of_elements() * mesh::number_of_facets( tBlock->element_type() ) ;
        }
        tKeys.set_size( tCount, 0 );
        tCount = 0 ;
        if ( this->number_of_dimensions() == 2 )
        {
            for ( mesh::Block * tBlock : mBlocks )
            {
                for ( mesh::Element * tElement : tBlock->elements() )
                {
                    for ( uint f=0; f<tElement->number_of_facets(); ++f )
                    {
                        tElement->get_corner_nodes_of_facet( f, tWork );

                        if ( tWork( 0 )->index() > tWork( 1 )->index() )
                        {
                            tA = tWork( 0 )->index() ;
                            tB = tWork( 1 )->index() ;
                        }
                        else
                        {
                            tA = tWork( 1 )->index() ;
                            tB = tWork( 0 )->index() ;
                        }
                        tKeys( tCount++ ) = tA * tNumNodes + tB ;
                    }
                }
            }
        }
        else if ( this->number_of_dimensions() == 3 )
        {
            for ( mesh::Block * tBlock : mBlocks )
            {
                for ( mesh::Element * tElement : tBlock->elements() )
                {
                    for ( uint f=0; f<tElement->number_of_facets(); ++f )
                    {
                        tElement->get_corner_nodes_of_facet( f, tWork );

                        sort( tWork, opVertexIndex );
                        tA = tWork( 2 )->index() ;
                        tB = tWork( 1 )->index() ;
                        tC = tWork( 0 )->index() ;
                        tKeys( tCount++ ) = ( tA * tNumNodes + tB ) * tNumNodes + tC ;
                    }
                }
            }
        }

        Cell< key128_t > tUniqueKeys( tKeys );
        unique( tUniqueKeys ) ;

        // one graph vertex per unique facet key, in sorted key order,
        // so that find_index_in_unique_cell() replaces a map
        index_t tNumFacets = tUniqueKeys.size() ;
        Cell< graph::Vertex * > tFacets( tNumFacets, nullptr );
        for ( index_t k=0; k<tNumFacets; ++k )
        {
            tFacets( k ) = new graph::Vertex() ;
            tFacets( k )->init_vertex_container( 2 );
        }

        tCount = 0 ;

        for ( mesh::Block * tBlock : mBlocks )
        {
            for ( mesh::Element * tElement : tBlock->elements() )
            {
                for ( uint f=0; f<tElement->number_of_facets(); ++f )
                {
                    // get the facet
                    graph::Vertex * tFacet = tFacets(
                        find_index_in_unique_cell( tUniqueKeys, tKeys( tCount++ ) ) );
                    tFacet->insert_vertex( tElement );
                }
            }
        }

        tCount = 0 ;
        for ( mesh::Block * tBlock : mBlocks )
        {
            for ( mesh::Element * tElement : tBlock->elements() )
            {
                tElement->allocate_neighbor_container()  ;

                for ( uint f=0; f<tElement->number_of_facets(); ++f )
                {
                    // get the facet
                    graph::Vertex * tFacet = tFacets(
                        find_index_in_unique_cell( tUniqueKeys, tKeys( tCount++ ) ) );

                    if ( tFacet->number_of_vertices() < 2 ) continue ;
                    if ( tFacet->vertex( 0 )->id() == tElement->id() )
                    {
                        tElement->insert_neighbor( reinterpret_cast< mesh::Element * >( tFacet->vertex( 1 ) ), f );
                    }
                    else
                    {
                        tElement->insert_neighbor( reinterpret_cast< mesh::Element * >( tFacet->vertex( 0 ) ), f );
                    }
                }
            }
        }

        // tidy up
        for ( graph::Vertex * tFacet : tFacets )
        {
            delete tFacet ;
        }

    }

    void
    Mesh::create_curve_map()
    {
        mCurveMap.clear() ;
        for ( mesh::Curve * tCurve : mCurves )
        {
            mCurveMap[ tCurve->id() ] = tCurve ;
        }
    }

//------------------------------------------------------------------------------

    size_t
    Mesh::memory() const
    {
        size_t aMem = sizeof( Mesh );

        // Path string
        aMem += mPath.capacity() * sizeof( char );
        aMem += mConfigText.capacity() * sizeof( char );

        // Nodes (owned objects)
        aMem += mNodes.size() * sizeof( mesh::Node * );
        for ( const mesh::Node * tNode : mNodes )
        {
            aMem += tNode->memory();
        }

        // Elements (owned via Blocks, just pointer storage here)
        aMem += mElements.size() * sizeof( mesh::Element * );
        for ( const mesh::Element * tElement : mElements )
        {
            aMem += tElement->memory();
        }

        // Facets (owned objects, each facet owns its element)
        aMem += mFacets.size() * sizeof( mesh::Facet * );
        for ( const mesh::Facet * tFacet : mFacets )
        {
            aMem += tFacet->memory();
        }

        // Edges (owned objects)
        aMem += mEdges.size() * sizeof( mesh::Edge * );
        for ( const mesh::Edge * tEdge : mEdges )
        {
            aMem += tEdge->memory();
        }

        // Faces (owned objects)
        aMem += mFaces.size() * sizeof( mesh::Face * );
        for ( const mesh::Face * tFace : mFaces )
        {
            aMem += tFace->memory();
        }

        // ControlPoints (owned objects)
        aMem += mControlPoints.size() * sizeof( mesh::ControlPoint * );
        for ( const mesh::ControlPoint * tControlPoint : mControlPoints )
        {
            aMem += tControlPoint->memory();
        }

        // Vertices (owned VERTEX elements, deleted in ~Mesh)
        aMem += mVertices.size() * sizeof( mesh::Element * );
        for ( const mesh::Element * tVertex : mVertices )
        {
            aMem += tVertex->memory();
        }

        // Boundary edges (references, just pointer storage)
        aMem += mBoundaryEdges.size() * sizeof( mesh::Element * );

        // Hanging containers (references, just pointer storage)
        aMem += mHangingNodes.size() * sizeof( mesh::Node * );
        aMem += mHangingEdges.size() * sizeof( mesh::Edge * );
        aMem += mHangingFaces.size() * sizeof( mesh::Face * );
        aMem += mHangingFacets.size() * sizeof( mesh::Facet * );
        aMem += mHangingElements.size() * sizeof( mesh::Element * );
        aMem += mHangingControlPoints.size() * sizeof( mesh::ControlPoint * );

        // TensorMeshConfig (optional, owned object)
        if ( mTensorConfig != nullptr )
        {
            aMem += sizeof( TensorMeshConfig );
        }

        // Periodicity (optional, owned object with non-owning pointer arrays)
        if ( mPeriodicity != nullptr )
        {
            aMem += sizeof( mesh::Periodicity );
            aMem += mPeriodicity->master_nodes().size()     * sizeof( mesh::Node * );
            aMem += mPeriodicity->slave_nodes().size()      * sizeof( mesh::Node * );
            aMem += mPeriodicity->master_edges().size()     * sizeof( mesh::Edge * );
            aMem += mPeriodicity->slave_edges().size()      * sizeof( mesh::Edge * );
            aMem += mPeriodicity->master_faces().size()     * sizeof( mesh::Face * );
            aMem += mPeriodicity->slave_faces().size()      * sizeof( mesh::Face * );
            aMem += mPeriodicity->master_facets().size()     * sizeof( mesh::Facet * );
            aMem += mPeriodicity->slave_facets().size()      * sizeof( mesh::Facet * );
        }

        // Blocks (owned objects)
        aMem += mBlocks.size() * sizeof( mesh::Block * );
        for ( const mesh::Block * tBlock : mBlocks )
        {
            aMem += tBlock->memory();
        }

        // SideSets (owned objects)
        aMem += mSideSets.size() * sizeof( mesh::SideSet * );
        for ( const mesh::SideSet * tSideSet : mSideSets )
        {
            aMem += tSideSet->memory();
        }

        // Curves (owned objects)
        aMem += mCurves.size() * sizeof( mesh::Curve * );
        for ( const mesh::Curve * tCurve : mCurves )
        {
            aMem += tCurve->memory();
        }

        // ThinShells (owned objects)
        aMem += mThinShells.size() * sizeof( mesh::ThinShell * );
        for ( const mesh::ThinShell * tThinShell : mThinShells )
        {
            aMem += tThinShell->memory();
        }

        // Fields (owned objects)
        aMem += mFields.size() * sizeof( mesh::Field * );
        for ( const mesh::Field * tField : mFields )
        {
            aMem += tField->memory();
        }

        // GlobalVariables (owned objects)
        aMem += mGlobalVariables.size() * sizeof( mesh::GlobalVariable * );
        for ( const mesh::GlobalVariable * tVariable : mGlobalVariables )
        {
            aMem += tVariable->memory();
        }

        // Abstract and orphaned nodes (references, just pointer storage)
        aMem += mAbstractNodes.size() * sizeof( mesh::Node * );
        aMem += mOrphanedNodes.size() * sizeof( mesh::Node * );

        // Maps (unordered_map: pair + hash-node overhead per entry)
        const size_t tHashNodeOverhead = sizeof( void * ) + sizeof( size_t );
        aMem += mNodeMap.size() * ( sizeof( std::pair< id_t, mesh::Node * > ) + tHashNodeOverhead );
        aMem += mElementMap.size() * ( sizeof( std::pair< id_t, mesh::Element * > ) + tHashNodeOverhead );
        aMem += mFacetMap.size() * ( sizeof( std::pair< id_t, mesh::Facet * > ) + tHashNodeOverhead );
        aMem += mBlockMap.size() * ( sizeof( std::pair< id_t, mesh::Block * > ) + tHashNodeOverhead );
        aMem += mSideSetMap.size() * ( sizeof( std::pair< id_t, mesh::SideSet * > ) + tHashNodeOverhead );
        aMem += mVertexMap.size() * ( sizeof( std::pair< id_t, mesh::Element * > ) + tHashNodeOverhead );
        aMem += mEdgeMap.size() * ( sizeof( std::pair< id_t, mesh::Edge * > ) + tHashNodeOverhead );
        aMem += mFaceMap.size() * ( sizeof( std::pair< id_t, mesh::Face * > ) + tHashNodeOverhead );
        aMem += mControlPointMap.size() * ( sizeof( std::pair< id_t, mesh::ControlPoint * > ) + tHashNodeOverhead );
        aMem += mCurveMap.size() * ( sizeof( std::pair< id_t, mesh::Curve * > ) + tHashNodeOverhead );

        // Field and global variable string maps
        aMem += mFieldMap.size() * ( sizeof( std::pair< string, mesh::Field * > ) + tHashNodeOverhead );
        aMem += mGlobalVariableMap.size() * ( sizeof( std::pair< string, mesh::GlobalVariable * > ) + tHashNodeOverhead );

        // Account for string keys in maps (approximate)
        for ( const auto & tPair : mFieldMap )
        {
            aMem += tPair.first.capacity() * sizeof( char );
        }
        for ( const auto & tPair : mGlobalVariableMap )
        {
            aMem += tPair.first.capacity() * sizeof( char );
        }

        return aMem;
    }

    /*void
    Mesh::save_fields( const string & aFilename, const uint aRunningTimestep )
    {
        BELFEM_ERROR( comm_rank() == 0 , "save_fields can only be called by rank 0");

        HDF5 tFile( aFilename, FileMode::NEW );

        tFile.create_group( "meta" );

        tFile.save_data( "timestep", mTimeStep );
        tFile.save_data( "timestamp", mTimeStamp );
        tFile.save_data( "running_timestep", aRunningTimestep );

        tFile.save_data( "checksum", this->checksum() );

        tFile.close_active_group();

        tFile.create_group( "fields" );


        Cell< string > tLabels( mFields.size() );
        Cell< uchar > tTypes( mFields.size() );
        for ( mesh::Field * tField : mFields )
        {
            tLabels.push( tField->label() );
            tTypes.push( static_cast< uchar >( tField->entity_type() ) );

            tFile.save_data( tField->label(), tField->data() );
        }

        tFile.save_data( "labels", tLabels );
        tFile.save_data( "types", tTypes );
        tFile.close_active_group();
        tFile.close();
    }*/

    void
    Mesh::save_meta( hid_t aFile, const uint aRunningTimestep )
    {
        herr_t tStatus = 0 ;
        hdf5::save_scalar_to_file( aFile,"timestep", mTimeStep, tStatus );
        hdf5::save_scalar_to_file(aFile, "timestamp", mTimeStamp, tStatus );
        hdf5::save_scalar_to_file( aFile, "running_timestep", aRunningTimestep, tStatus );
        hdf5::save_scalar_to_file( aFile, "checksum", this->checksum(), tStatus );

        // the checksum hashes nodes and element connectivity only, so it
        // cannot see a different edge / face layout on the same nodes --
        // thin-shell interfaces with or without duplicate dofs hash the
        // same. The entity counts do see it
        hdf5::save_scalar_to_file( aFile, "num_edges", ( size_t ) this->number_of_edges(), tStatus );
        hdf5::save_scalar_to_file( aFile, "num_faces", ( size_t ) this->number_of_faces(), tStatus );
    }

    uint
    Mesh::load_meta( hid_t aFile )
    {
        herr_t tStatus = 0 ;
        hdf5::load_scalar_from_file( aFile,"timestep", mTimeStep, tStatus );
        hdf5::load_scalar_from_file(aFile, "timestamp", mTimeStamp, tStatus );

        uint aRunningTimestep = 0 ;
        hdf5::load_scalar_from_file( aFile, "running_timestep", aRunningTimestep, tStatus );

        size_t tChecksum = 0 ;
        hdf5::load_scalar_from_file( aFile, "checksum",tChecksum , tStatus );

        BELFEM_ERROR( tChecksum == this->checksum(), "Memdump checksum mismatch" );

        // dumps written before 2026-09-01 carry no counts; those are accepted
        // on the checksum alone, as before
        if (    hdf5::dataset_exists( aFile, "num_edges" )
             && hdf5::dataset_exists( aFile, "num_faces" ) )
        {
            size_t tNumEdges = 0 ;
            size_t tNumFaces = 0 ;
            hdf5::load_scalar_from_file( aFile, "num_edges", tNumEdges, tStatus );
            hdf5::load_scalar_from_file( aFile, "num_faces", tNumFaces, tStatus );

            BELFEM_ERROR(    tNumEdges == ( size_t ) this->number_of_edges()
                          && tNumFaces == ( size_t ) this->number_of_faces(),
                "Memdump was written for a mesh with %lu edges and %lu faces; this mesh has %lu and %lu.\n"
                "       Same nodes, different discretization -- most likely the thin-shell ghost switch\n"
                "       ( nitsche ghost penalty { eta } ) differs between the dump and this deck.\n"
                "       Start fresh ( timestep { restart : false ; } ) or match the deck to the dump.",
                ( long unsigned int ) tNumEdges, ( long unsigned int ) tNumFaces,
                ( long unsigned int ) this->number_of_edges(),
                ( long unsigned int ) this->number_of_faces() );
        }

        return aRunningTimestep ;
    }

    void
    Mesh::save_fields( hid_t aFile )
    {
        Cell< string > tLabels( mFields.size() );
        Cell< uchar > tTypes( mFields.size() );

        herr_t tStatus = 0 ;
        for ( mesh::Field * tField : mFields )
        {
            tLabels.push( tField->label() );
            tTypes.push( static_cast< uchar >( tField->entity_type() ) );

            hdf5::save_vector_to_file( aFile, tField->label(), tField->data(), tStatus );
        }
        hdf5::save_strings_to_file( aFile, "labels", tLabels, tStatus );
        hdf5::save_array_to_file( aFile, "types", tTypes.data(), tTypes.size(), tStatus );
    }

    void
    Mesh::save_globals( hid_t aFile )
    {
        uint n = mGlobalVariables.size();

        Vector< id_t > tIds( n, 0 );
        Cell< string > tLabels( n, "" );
        Vector< real > tValues( n, 0 );

        uint tCount = 0 ;

        for ( auto tVar : mGlobalVariables )
        {
            tIds( tCount ) = tVar->id() ;
            tLabels( tCount ) = tVar->label() ;
            tValues( tCount ) = tVar->value() ;
            ++tCount;
        }
        herr_t tStatus = 0 ;
        hdf5::save_vector_to_file( aFile, "ids", tIds, tStatus );
        hdf5::save_strings_to_file( aFile, "labels", tLabels, tStatus );
        hdf5::save_vector_to_file( aFile, "values", tValues, tStatus );
    }

    void
    Mesh::load_globals( hid_t aFile )
    {
        Vector< id_t > tIDs;
        Cell< string > tLabels ;
        Vector< real > tValues ;
        herr_t tStatus = 0 ;
        hdf5::load_vector_from_file( aFile, "ids", tIDs, tStatus );
        hdf5::load_strings_from_file( aFile, "labels", tLabels, tStatus );
        hdf5::load_vector_from_file( aFile, "values", tValues, tStatus );

        uint tCount = 0 ;

        for ( const string & tLabel : tLabels )
        {
            real tVal = tValues( tCount );
            id_t tID = tIDs( tCount );
            mesh::GlobalVariable * tVar = nullptr ;
            if ( mGlobalVariableMap.key_exists( tLabel ) )
            {
                tVar = mGlobalVariableMap( tLabel );
                tVar->value() = tVal ;
                BELFEM_ERROR( tVar->id() == tID, "ID mismatch for global variable %s. ( is %lu, expect %lu )",
                    tLabel.c_str(), ( long unsigned int ) tVar->id(), ( long unsigned int ) tID );
            }
            else
            {
                tVar = new mesh::GlobalVariable( tLabel, tID, tVal );
                mGlobalVariableMap[ tLabel ] = tVar ;
                mGlobalVariableIDMap[ tID ] = tVar ;
                mGlobalVariables.push( tVar );

                // keep the auto-ID counter in sync: a later
                // create_global_variable() with auto-ID must not reuse an
                // ID this dump just occupied ( latent until a caller
                // creates globals after a load — audited 2026-08-15 )
                mNumberOfGlobalVariables =
                    mNumberOfGlobalVariables < tID ? tID : mNumberOfGlobalVariables ;
            }
            ++tCount ;
        }

        mGlobalVariables.shrink_to_fit();
    }

    /*uint
    Mesh::load_fields( const string & aFilename )
    {
        BELFEM_ERROR( comm_rank() == 0 , "load_fields can only be called by rank 0");

        HDF5 tFile( aFilename, FileMode::OPEN_RDONLY );

        tFile.select_group( "meta" );
        size_t tChecksum ;
        tFile.load_data( "checksum", tChecksum );
        BELFEM_ERROR( tChecksum == this->checksum(), "Checksum mismatch" );

        tFile.load_data( "timestep", mTimeStep );
        tFile.load_data( "timestamp", mTimeStamp );

        uint aRunningTimestep = 0 ;
        tFile.load_data( "running_timestep", aRunningTimestep );

        tFile.close_active_group();
        tFile.select_group( "fields" );




        tFile.close_active_group();
        tFile.close();

        return aRunningTimestep ;
    }*/

    void
    Mesh::load_fields( hid_t aFile )
    {
        Cell< string > tLabels ;
        Cell< uchar > tTypes ;

        herr_t tStatus = 0 ;
        hdf5::load_strings_from_file( aFile, "labels", tLabels, tStatus );

        hsize_t tSize = hdf5::get_array_size( aFile,"types" );
        tTypes.set_size( tSize,  0 );
        hdf5::load_array_from_file( aFile, "types", tTypes.data(), tSize, tStatus );


        uint f = 0 ;
        for ( const string & tLabel : tLabels )
        {
            if ( this->field_exists( tLabel ) )
            {
                // a sized target is authoritative ( node fields from the Field
                // ctor, edge / face fields from create_fields with the entity
                // multiplicity ), so a mismatch is a dump of a different
                // discretization. An empty target is an edge / face shell the
                // Field ctor left unsized -- the BDF history levels -- and the
                // loader sizes it before it reads
                const hsize_t tLength = hdf5::get_array_size( aFile, tLabel );
                Vector< real > & tTarget = this->field( tLabel )->data() ;

                BELFEM_ERROR( tTarget.length() == 0
                              || tLength == ( hsize_t ) tTarget.length(),
                    "Memdump field '%s' has %lu entries, the mesh field has %lu: the dump was written\n"
                    "       for a different discretization. Start fresh ( timestep { restart : false ; } ).",
                    tLabel.c_str(),
                    ( long unsigned int ) tLength,
                    ( long unsigned int ) tTarget.length() );

                hdf5::load_vector_from_file(  aFile, tLabel, tTarget , tStatus );
            }
            else
            {
                EntityType tType = static_cast< EntityType >( tTypes( f ) );

                Vector< real > & tData = this->create_field( tLabel, tType );
                hdf5::load_vector_from_file(  aFile, tLabel, tData , tStatus );
            }

            ++f ;
        }
    }

//------------------------------------------------------------------------------

}
