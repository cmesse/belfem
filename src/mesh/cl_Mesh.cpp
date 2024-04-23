//
// Created by Christian Messe on 2019-07-25.
//
#include "cl_Block.hpp"
#include "cl_Mesh.hpp"
#include "cl_Mesh_ExodusWriter.hpp"
#include "fn_unique.hpp"
#include "stringtools.hpp"
#include "cl_Mesh_GmshReader.hpp"
#include "cl_Mesh_OrderConverter.hpp"
#include "cl_Mesh_CurvedElementChecker.hpp"
#include "cl_Mesh_Partitioner.hpp"
#include "cl_Mesh_HDF5Reader.hpp"
#include "cl_Mesh_HDF5Writer.hpp"
#include "cl_Mesh_VtkWriter.hpp"
#include "cl_Timer.hpp"
#include "cl_Logger.hpp"
#include "commtools.hpp"
#include "assert.hpp"
#include "cl_Element_Factory.hpp"
#include "cl_EdgeFactory.hpp"
#include "cl_FaceFactory.hpp"
#include "fn_max.hpp"
#include "cl_Mesh_OrientationChecker.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    Mesh::Mesh( const proc_t aMasterProc, const bool aComputeConnectivities ) :
        mMasterProc(  aMasterProc < comm_size() ? aMasterProc : 0 ),
        mComputeConnectivities( aComputeConnectivities )
    {

    }

//------------------------------------------------------------------------------

    Mesh::Mesh(
            const uint & aNumberOfDimensions ,
            const proc_t aMasterProc,
            const bool aComputeConnectivities  ) :
        mMasterProc(  aMasterProc < comm_size() ? aMasterProc : 0 ),
        mNumberOfDimensions( aNumberOfDimensions ),
        mComputeConnectivities( aComputeConnectivities )
    {

    }

//------------------------------------------------------------------------------

    Mesh::Mesh(
            const string & aPath,
            const proc_t aMasterProc,
            const bool aComputeConnectivities ) :
            mMasterProc(  aMasterProc < comm_size() ? aMasterProc : 0 ),
            mComputeConnectivities( aComputeConnectivities )
    {
        if( comm_rank() == aMasterProc )
        {
            // start a timer
            Timer tTimer;

            string tType = string_to_lower( filetype( aPath ));


            message( 2, "    Reading mesh from %s...",
                     filename( aPath ).c_str() );

            if ( tType == "msh" )
            {
                mesh::GmshReader tReader( aPath, this );
            }
            else if ( tType == "hdf5" )
            {
                // read the mesh from the file
                mesh::HDF5Reader tReader( aPath, this );
            }
            else
            {
                BELFEM_ERROR( false, "don't know how to read a mesh of type <%s>.",
                             tType.c_str());
            }

            uint tTime = tTimer.stop();

            message( 2, "    Nodes   : %lu",
                     ( long unsigned int ) this->number_of_nodes());

            message( 2, "    Elements: %lu",
                     ( long unsigned int ) this->number_of_elements());

            message( 2, "    Blocks  : %u",
                     ( unsigned int ) this->number_of_blocks());

            message( 2, "    Sidesets: %u",
                     ( unsigned int ) this->number_of_sidesets());

            message( 2, "    Fields  : %u",
                     ( unsigned int ) this->number_of_fields() );

            message( 2, "    Globals : %u",
                     ( unsigned int ) this->number_of_global_variables() );

            message( 2, "    Time %u ms.\n",
                     ( unsigned int ) tTime );



            if( gComm.size() > 1 )
            {
                Vector< proc_t > tCommTable( gComm.size() );
                uint tCount = 0 ;
                for( proc_t p=0; p<gComm.size(); ++p )
                {
                    tCommTable( tCount++ ) = p ;
                }

                comm_barrier() ;

                send( tCommTable, mNumberOfDimensions );
            }
        }
        else
        {
            comm_barrier() ;

            receive( aMasterProc, mNumberOfDimensions );
        }
    }

//------------------------------------------------------------------------------^

    Mesh::~Mesh()
    {
        // delete global variables
        for( auto tVariable: mGlobalVariables )
        {
            delete tVariable;
        }

        for( auto tBlock : mBlocks )
        {
            delete tBlock;
        }

        for( auto tCut : mCuts )
        {
            delete tCut;
        }

        for( auto tSideSet : mSideSets )
        {
            delete tSideSet;
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

        // delete facets
        for ( auto tFacet: mFacets )
        {
            delete tFacet;
        }

        // delete elements
        for ( auto tElement: mElements )
        {
            delete tElement;
        }

        // delete connectors
        for ( auto tConnector: mConnectors )
        {
            delete tConnector;
        }

        // delete nodes
        for ( auto tNode: mNodes )
        {
            delete tNode;
        }

        // delete edges along boundary
        for ( auto tEdge: mBoundaryEdges )
        {
            delete tEdge;
        }

        // delete bearings
        for( auto tVertex : mVertices )
        {
            delete tVertex;
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
    }

//------------------------------------------------------------------------------

    void
    Mesh::save( const string & aFilePath )
    {
        if( comm_rank() == mMasterProc )
        {
            // start the timer
            Timer tTimer;

            message( 2, "\n    Saving Mesh to %s ...", filename( aFilePath ).c_str() );

            string tType = string_to_lower( filetype( aFilePath ));


            if( tType == "hdf5" )
            {
                // create writer object and save file
                mesh::HDF5Writer tWriter( aFilePath, this );
            }
            else if( tType == "vtk" && comm_rank() == mMasterProc )
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
                message( 2, "    Time %u ms.\n",
                         ( unsigned int ) tTime );
            }
            else
            {
                message( 2, "    Time %4.1f s.\n",
                         ( float ) tTime * 0.001 );
            }

        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::connect_nodes_to_elements()
    {
        for( mesh::Node * tNode: mNodes )
        {
            tNode->reset_element_container();
        }

        // loop over all elements
        for( mesh::Element * tElement: mElements )
        {
            for( uint k=0; k<tElement->number_of_nodes(); ++k )
            {
                tElement->node( k )->increment_element_counter();
            }
        }
        for( mesh::Facet * tFacet: mConnectors )
        {
            // get element of facet
            mesh::Element * tElement = tFacet->element() ;

            for( uint k=0; k<tElement->number_of_nodes(); ++k )
            {
                tElement->node( k )->increment_element_counter();
            }
        }

        // allocate container for nodes
        for( mesh::Node * tNode: mNodes )
        {
            tNode->allocate_element_container();
        }

        for( mesh::Element * tElement: mElements )
        {
            for( uint k=0; k<tElement->number_of_nodes(); ++k )
            {
                tElement->node( k )->add_element( tElement );
            }
        }
        for( mesh::Facet * tFacet: mConnectors )
        {
            // get element of facet
            mesh::Element * tElement = tFacet->element() ;

            for( uint k=0; k<tElement->number_of_nodes(); ++k )
            {
                tElement->node( k )->add_element( tElement );
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::connect_nodes_to_edges()
    {
        for( mesh::Node * tNode : mNodes )
        {
            tNode->reset_edge_container();
        }

        for( mesh::Edge * tEdge : mEdges )
        {
            for( uint k=0; k<tEdge->number_of_nodes(); ++k )
            {
                tEdge->node( k )->increment_edge_counter();
            }
        }

        for( mesh::Node * tNode : mNodes )
        {
            tNode->allocate_edge_container();
        }

        for( mesh::Edge * tEdge : mEdges )
        {
            for( uint k=0; k<tEdge->number_of_nodes(); ++k )
            {
                tEdge->node( k )->add_edge( tEdge );
            }
        }

    }

//------------------------------------------------------------------------------

    void
    Mesh::connect_edges_to_elements( Cell< mesh::Element * > & aElements )
    {
        for( mesh::Edge * tEdge: mEdges )
        {
            tEdge->reset_element_container();
        }

        // loop over all elements
        for( mesh::Element * tElement: aElements )
        {
            if( tElement->has_edges() )
            {
                for ( uint k = 0; k < tElement->number_of_edges(); ++k )
                {
                    tElement->edge( k )->increment_element_counter();
                }
            }
        }

        // allocate container for nodes
        for( mesh::Edge * tEdge: mEdges )
        {
            tEdge->allocate_element_container();
        }

        for( mesh::Element * tElement: aElements )
        {
            if( tElement->has_edges() )
            {
                for ( uint k = 0; k < tElement->number_of_edges(); ++k )
                {
                    tElement->edge( k )->add_element( tElement );
                }
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::connect_nodes_to_facets()
    {
        for( mesh::Node * tNode: mNodes )
        {
            tNode->reset_facet_container();
        }

        // loop over all elements
        for( mesh::Facet * tFacet: mFacets )
        {
            if( tFacet->is_flagged() )
            {
                for ( uint k = 0; k < tFacet->element()->number_of_nodes(); ++k )
                {
                    tFacet->element()->node( k )->increment_facet_counter();
                }
            }
        }

        // allocate container for nodes
        for( mesh::Node * tNode: mNodes )
        {
            tNode->allocate_facet_container();
        }

        for( mesh::Facet * tFacet: mFacets )
        {
            if( tFacet->is_flagged() )
            {
                for ( uint k = 0; k < tFacet->element()->number_of_nodes(); ++k )
                {
                    tFacet->element()->node( k )->add_facet( tFacet );
                }
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::connect_facets_to_elements()
    {
        // loop over all facets
        for ( mesh::Facet * tFacet : mFacets )
        {

            // get number of nodes from this facet
            uint tNumNodes = tFacet->element()->number_of_corner_nodes();

            // count max size of candidates
            uint tCount = 0;
            for ( uint k = 0; k < tNumNodes; ++k )
            {
                mesh::Node * tNode = tFacet->node( k );

                tCount += tNode->number_of_elements();
            }

            // create list of Element candidates
            Cell<mesh::Element *> tCandidates( tCount, nullptr );

            // reset counter
            tCount = 0;

            // populate candidates
            for ( uint k = 0; k < tNumNodes; ++k )
            {
                mesh::Node * tNode = tFacet->node( k );

                for ( uint e = 0; e < tNode->number_of_elements(); ++e )
                {
                    tCandidates( tCount++ ) = tNode->element( e );
                }
            }

            // make result unique
            unique( tCandidates );


            // loop over all candidates
            Cell<mesh::Node *> tNodes;

            mesh::Element * tMaster = nullptr;
            uint tMasterFaceIndex = BELFEM_UINT_MAX;
            uint tSlaveFaceIndex = BELFEM_UINT_MAX;

            mesh::Element * tSlave = nullptr;

            for ( mesh::Element * tElement : tCandidates )
            {
                // get number of faces from element
                uint tNumFaces = tElement->number_of_facets();

                // loop over all facets
                for ( uint k = 0; k < tNumFaces; ++k )
                {
                    // unflag my nodes
                    tFacet->unflag_nodes();

                    // get nodes from face
                    tElement->get_corner_nodes_of_facet( k, tNodes );

                    // flag all nodes from this face
                    for ( mesh::Node * tNode: tNodes )
                    {
                        tNode->flag();
                    }

                    // count how many nodes are flagged
                    uint tFCount = 0;

                    for ( uint i = 0; i < tNumNodes; ++i )
                    {
                        if ( tFacet->node( i )->is_flagged() )
                        {
                            // increment node counter
                            ++tFCount;
                        }
                    }

                    if ( tFCount == tNumNodes )
                    {
                        if ( tMaster == nullptr )
                        {
                            tMaster = tElement;
                            tMasterFaceIndex = k;
                        }
                        else if ( tSlave == nullptr )
                        {
                            tSlave = tElement;
                            tSlaveFaceIndex = k;
                        }
                        else
                        {
                            break;
                        }
                    }
                }
                //  break if slave has been found
                if ( tSlave != nullptr )
                {
                    break;
                }
            } // end loop over candidates

            BELFEM_ASSERT( tMaster != nullptr,
                "Could not find master for facet %u",
                          ( unsigned int ) tFacet->id() );

            // check if slave is found
            if ( tSlave == nullptr )
            {
                tFacet->set_master( tMaster, tMasterFaceIndex );
                tFacet->element()->allocate_element_container( 1 );
                tFacet->element()->insert_element( tMaster );
            }
            else
            {
                if ( tMaster->id() < tSlave->id() )
                {
                    tFacet->set_master( tMaster, tMasterFaceIndex );
                    tFacet->set_slave( tSlave, tSlaveFaceIndex );
                }
                else
                {
                    tFacet->set_master( tSlave, tSlaveFaceIndex );
                    tFacet->set_slave( tMaster, tMasterFaceIndex );
                }
                tFacet->element()->allocate_element_container( 2 );
                tFacet->element()->insert_element( tMaster );
                tFacet->element()->insert_element( tSlave );
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

        //tCount = 0;
        for( mesh::Facet * tFacet : mConnectors )
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
    Mesh::connect_elements_to_elements()
    {
        // exit the function if there are no elements on this mesh
        if( mElements.size() == 0 )
        {
            return ;
        }

        index_t tCount = 0 ;

        Vector< uint > tNumFacetsPerElement( mElements.size(), 0 );

        // count facets per element, since some elements are only connected via facets
        for( mesh::Facet * tFacet : mFacets )
        {
            if( tFacet->has_master()  )
            {
                tNumFacetsPerElement( tFacet->master()->index())++;
            }
            if( tFacet->has_slave() )
            {
                tNumFacetsPerElement( tFacet->slave()->index())++;
            }
        }

        // allocate memory
        Matrix< index_t > tFacetsPerElement( max(tNumFacetsPerElement), mElements.size(), 0 );

        // reset counter
        tNumFacetsPerElement.fill( 0 );

        // count facets per element, since some elements are only connected via facets
        for( mesh::Facet * tFacet : mFacets )
        {
            if( tFacet->has_master()  )
            {
                tFacetsPerElement( tNumFacetsPerElement( tFacet->master()->index() )++, tFacet->master()->index() ) = tCount ;
            }
            if( tFacet->has_slave() )
            {
                tFacetsPerElement( tNumFacetsPerElement( tFacet->slave()->index() )++, tFacet->slave()->index() ) = tCount ;
            }
            tFacet->set_index( tCount++ );
        }

        Vector< index_t > tIndices;

        // count elements per element
        for( mesh::Element * tElement : mElements )
        {
            uint tN = tElement->number_of_nodes();

            // get number of dimensions of this element
            int tNumDim = mesh::dimension( tElement->type() );

            // get id of this element
            const id_t tID = tElement->id();

            // count space needed for index vector
            tCount = 0 ;

            // loop over all nodes and count connected elements
            for( uint k=0; k<tN; ++k )
            {
                mesh::Node * tNode = tElement->node( k );

                // loop over all elements of this node
                for( uint e=0; e<tNode->number_of_elements(); ++e )
                {
                    // the connected element must have the same dimension
                    // but must not be the element itself
                    if ( mesh::dimension( tNode->element( e )->type() ) == tNumDim )
                    {
                        ++tCount ;
                    }
                }
            }

            for( uint f=0; f<tNumFacetsPerElement( tElement->index()) ; ++f )
            {
                // get facet
                mesh::Facet * tFacet = mFacets( tFacetsPerElement( f, tElement->index() ) );

                if( tFacet->has_master() )
                {
                    ++tCount ;
                }
                if( tFacet->has_slave() )
                {
                    ++tCount ;
                }
            }

            // allocate index vector
            tIndices.set_size( tCount, tID );

            // reset counter
            tCount = 0 ;

            // loop over all nodes and count connected elements
            for( uint k=0; k<tN; ++k )
            {
                // get node pointer
                mesh::Node * tNode = tElement->node( k );

                // loop over all elements of this node
                for( uint e=0; e<tNode->number_of_elements(); ++e )
                {
                    // the connected element must have the same dimension
                    // but must not be the element itself
                    if ( mesh::dimension( tNode->element( e )->type()) == tNumDim )
                    {
                        tIndices( tCount++ ) = tNode->element( e )->index();
                    }
                }
            }

            // check connections that are only over facets
            for( index_t f=0; f<tNumFacetsPerElement( tElement->index() ) ; ++f )
            {
                // get facet
                mesh::Facet * tFacet = mFacets( tFacetsPerElement( f, tElement->index() ) );

                if( tFacet->has_master()  )
                {
                    tIndices( tCount++ ) = tFacet->master()->index() ;
                }
                if( tFacet->has_slave() )
                {
                    tIndices( tCount++ ) = tFacet->slave()->index() ;
                }
            }

            // make indices unique
            unique(tIndices );

            uint tNumElements = tIndices.length();

            // test if this element has neighbors
            if( tNumElements > 1 )
            {
                // allocate element container
                tElement->allocate_element_container( tNumElements - 1 );

                // loop over all connected elements
                for( uint k=0; k<tNumElements; ++k )
                {
                    // make sure that this is not the element itself
                    if( tIndices( k )  != tElement->index() )
                    {
                        // add element to container
                        tElement->insert_element( mElements( tIndices( k ) ) );
                    }
                }
            }
        }
    }


//------------------------------------------------------------------------------

    void
    Mesh::connect_nodes_to_nodes()
    {
        // loop over all nodes
        for( mesh::Node * tNode : mNodes )
        {
            // get number of elements that are connected to this node
            uint tNumElements = tNode->number_of_elements();

            // count maximum number of nodes
            uint tCount = 0;
            for( uint e=0; e<tNumElements; ++e )
            {
                tCount += tNode->element( e )->number_of_nodes();
            }

            // Allocate cell
            Cell< mesh::Node * > tNodes( tCount, nullptr );

            // reset counter
            tCount = 0;

            // populate cell
            for( uint e=0; e<tNumElements; ++e )
            {
                // get element
                mesh::Element * tElement = tNode->element( e );

                // get number of nodes connected to this element
                uint tNumNodes = tElement->number_of_nodes();

                // loop over all nodes that are connected to this element
                for ( uint k = 0; k <tNumNodes; ++k )
                {
                    tNodes( tCount++ ) = tElement->node( k );
                }
            }

            // note: unique might also work here

            // unflag all nodes on this list
            for( mesh::Node * tOtherNode : tNodes )
            {
                tOtherNode->unflag();
            }

            // flag myself
            tNode->flag();

            // count actual nodes
            tCount = 0;
            for( mesh::Node * tOtherNode : tNodes )
            {
                if( ! tOtherNode->is_flagged() )
                {
                    // flag this node
                    tOtherNode->flag();

                    // increment counter
                    ++tCount;
                }
            }

            // allocate node container
            tNode->allocate_node_container( tCount );

            // unflag all nodes on this list
            for( mesh::Node * tOtherNode : tNodes )
            {
                tOtherNode->unflag();
            }

            // flag myself
            tNode->flag();

            // reset counter
            tCount = 0;

            for( mesh::Node * tOtherNode : tNodes )
            {
                if( ! tOtherNode->is_flagged() )
                {
                    // flag this node
                    tOtherNode->flag();

                    // add this node to the list
                    tNode->insert_node( tOtherNode, tCount++ );
                }
            }
        }

        // tidy up
        this->unflag_all_nodes();

    }

//------------------------------------------------------------------------------

    void
    Mesh::connect_edges_to_edges()
    {
        // loop over all edges
        for( mesh::Edge * tEdge : mEdges )
        {
            // get number of elements that are connected to this edge
            uint tNumElements = tEdge->number_of_elements();

            // count maximum number of edges
            uint tCount = 0;
            for( uint e=0; e<tNumElements; ++e )
            {
                tCount += tEdge->element( e )->number_of_edges();
            }

            // Allocate cell
            Cell< mesh::Edge * > tEdges( tCount, nullptr );

            // reset counter
            tCount = 0;

            // populate cell
            for( uint e=0; e<tNumElements; ++e )
            {
                // get element
                mesh::Element * tElement = tEdge->element( e );

                // get number of edges connected to this element
                uint tNumEdges = tElement->number_of_edges();

                // loop over all edges that are connected to this element
                for ( uint k = 0; k <tNumEdges; ++k )
                {
                    tEdges( tCount++ ) = tElement->edge( k );
                }
            }

            // note: unique might also work here

            // unflag all edges on this list
            for( mesh::Edge * tOtherEdge : tEdges )
            {
                tOtherEdge->unflag();
            }

            // flag myself
            tEdge->flag();

            // count actual edges
            tCount = 0;
            for( mesh::Edge * tOtherEdge : tEdges )
            {
                if( ! tOtherEdge->is_flagged() )
                {
                    // flag this edge
                    tOtherEdge->flag();

                    // increment counter
                    ++tCount;
                }
            }

            // make sure that container is empty
            tEdge->reset_edge_container() ;

            // allocate edge container
            tEdge->allocate_edge_container( tCount );

            // unflag all edges on this list
            for( mesh::Edge * tOtherEdge : tEdges )
            {
                tOtherEdge->unflag();
            }

            // flag myself
            tEdge->flag();

            // reset counter
            tCount = 0;

            for( mesh::Edge * tOtherEdge : tEdges )
            {
                if( ! tOtherEdge->is_flagged() )
                {
                    // flag this edge
                    tOtherEdge->flag();

                    // add this edge to the list
                    tEdge->insert_edge( tOtherEdge, tCount++ );
                }
            }
        }

        // tidy up
        this->unflag_all_edges();

    }

//------------------------------------------------------------------------------

    void
    Mesh::connect_faces_to_edges_and_edges_to_faces()
    {
        if( mNumberOfDimensions == 2 )
        {
            return ;
        }

        this->unflag_all_nodes() ;
        this->unflag_all_edges();
        this->unflag_all_faces();
        this->unflag_all_elements() ;

        Cell< mesh::Edge * > tEdges ;

        for( mesh::Face * tFace : mFaces )
        {
            tFace->master()->get_edges_of_facet( tFace->index_on_master(),
                                                 tEdges );

            tFace->allocate_edge_container( tEdges.size() ) ;

            uint tCount = 0 ;
            for( mesh::Edge * tEdge : tEdges )
            {
                tFace->insert_edge( tEdge, tCount++ );
                tEdge->increment_face_counter();
            }
        }

        for( mesh::Edge * tEdge : mEdges )
        {
            tEdge->allocate_face_container();
        }

        for( mesh::Face * tFace : mFaces )
        {
            tFace->master()->get_edges_of_facet( tFace->index_on_master(),
                                                 tEdges );
            for( mesh::Edge * tEdge : tEdges )
            {
                tEdge->add_face( tFace );
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::connect_vertices_to_vertices()
    {


        index_t tNumVertices = this->number_of_nodes() + this->number_of_edges() ;

        // loop over all edges
        for( index_t v=0; v<tNumVertices; ++v )
        {
            mesh::Vertex * tVertex =  v < this->number_of_nodes() ?
                    reinterpret_cast<mesh::Vertex *>(mNodes( v )) :
                    reinterpret_cast<mesh::Vertex *>(mEdges( v - this->number_of_nodes() ) );

            // get number of elements that are connected to this edge
            uint tNumElements = tVertex->number_of_elements();

            // count maximum number of edges
            uint tCount = 0;
            for( uint e=0; e<tNumElements; ++e )
            {
                tCount   += tVertex->element( e )->number_of_nodes()
                          + tVertex->element( e )->number_of_edges() ;
            }

            // Allocate temporary cell
            Cell< mesh::Vertex * > tVertices( tCount, nullptr );

            // reset counter
            tCount = 0;

            // populate cell
            for( uint e=0; e<tNumElements; ++e )
            {
                // get element
                mesh::Element * tElement = tVertex->element( e );

                // get number of nodes connected to this element
                uint tNumNodes = tElement->number_of_nodes();

                // loop over all edges that are connected to this element
                for ( uint k = 0; k <tNumNodes; ++k )
                {
                    tVertices( tCount++ ) = tElement->node( k );
                }

                // get number of edges connected to this element
                uint tNumEdges = tElement->number_of_edges();

                // loop over all edges that are connected to this element
                for ( uint k = 0; k <tNumEdges; ++k )
                {
                    tVertices( tCount++ ) = tElement->edge( k );
                }
            }

            // note: unique might also work here

            // unflag all edges on this list
            for( mesh::Vertex * tOtherVertex : tVertices )
            {
                tOtherVertex->unflag();
            }

            // flag myself
            tVertex->flag();

            // count actual edges
            tCount = 0;
            for( mesh::Vertex * tOtherVertex : tVertices )
            {
                if( ! tOtherVertex->is_flagged() )
                {
                    // flag this vertex
                    tOtherVertex->flag();

                    // increment counter
                    ++tCount;
                }
            }

            // allocate edge container
            tVertex->init_vertex_container( tCount );

            // unflag all edges on this list
            for( mesh::Vertex * tOtherVertex : tVertices )
            {
                tOtherVertex->unflag();
            }

            // flag myself
            tVertex->flag();

            // reset counter
            tCount = 0;

            for( mesh::Vertex * tOtherVertex : tVertices )
            {
                if( ! tOtherVertex->is_flagged() )
                {
                    // flag this edge
                    tOtherVertex->flag();

                    // add this edge to the list
                    tVertex->insert_vertex( tOtherVertex );
                }
            }
        }

        // tidy up
        this->unflag_all_nodes();
        this->unflag_all_edges();

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

        // add entry into map
        mGlobalVariableMap[ aLabel ] = tVariable ;

        // add entro to container
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
                mFacets( tCount++ ) = tSideSet->facet_by_index( f );
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::collect_facets_from_cuts()
    {
        BELFEM_ERROR( mConnectors.size() == 0,
                     "collect_facets_from_sidesets() must not be called if facet container is already filled" );

        index_t tCount = 0;

        // count number of sidesets
        for(  mesh::SideSet * tSideSet : mCuts )
        {
            tCount += tSideSet->number_of_facets();
        }

        // allocate facet container
        mConnectors.set_size( tCount, nullptr );


        // reset facet counter
        tCount = 0;

        for(  mesh::SideSet * tSideSet : mCuts )
        {
            index_t tNumElements = tSideSet->number_of_facets();

            for( index_t f=0; f<tNumElements; ++f )
            {
                mConnectors( tCount++ ) = tSideSet->facet_by_index( f );
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
    Mesh::collect_edges_from_group(
            Cell< mesh::Element * > & aElements,
            const id_t aGroupId )
    {

        // count elemenents on group
        index_t tCount = 0;

        // only collect elements from this type
        for( mesh::Element * tElement : mBoundaryEdges )
        {

            if( tElement->geometry_tag() == aGroupId )
            {
                ++tCount;
            }
        }


        aElements.set_size( tCount, nullptr );
        tCount = 0;

        for( mesh::Element * tElement : mBoundaryEdges )
        {
            if( tElement->geometry_tag() == aGroupId )
            {
                aElements( tCount++ ) = tElement;
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::compute_facet_orientations()
    {
        for( mesh::Facet * tFacet : mFacets )
        {
            tFacet->compute_orientation();
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::finalize()
    {
        if( comm_rank() == mMasterProc && ! mIsFinalized )
        {
            if ( mElements.size() == 0 )
            {
                this->collect_elements_from_blocks();
            }

            this->update_node_indices();
            this->update_element_indices();

            if ( mFacets.size() == 0 && mSideSets.size() > 0 )
            {
                this->collect_facets_from_sidesets();
            }

            if( mConnectors.size() == 0 && mCuts.size() > 0 )
            {
                this->collect_facets_from_cuts();
            }

            // update indices of facets and connectors
            this->update_facet_indices();

            if( mComputeConnectivities )
            {
                this->connect_nodes_to_elements();

                if( ! mFacetsAreLinked )
                {
                    this->connect_facets_to_elements();
                }

                // make sure that all facets are flagged
                // this is needed for the node to facets stem to work
                for( mesh::Facet * tFacet : mFacets )
                {
                    tFacet->flag() ;
                }

                this->connect_nodes_to_facets();

                this->connect_elements_to_elements();

                this->connect_nodes_to_nodes();

                for( mesh::SideSet * tSideset : mSideSets )
                {
                    tSideset->collect_nodes() ;
                }
            }

            this->create_maps();

            this->set_node_owners();

            this->unflag_all_elements() ;

            this->set_connector_owners();

            this->set_vertex_owners();
        }

        mIsFinalized = true ;

        this->set_block_ids() ;

        this->compute_max_element_order();

        this->compute_facet_orientations();

        this->link_ghost_elements() ;

    }

//------------------------------------------------------------------------------

    void
    Mesh::unfinalize()
    {
        if( comm_rank() == mMasterProc && mIsFinalized )
        {

            for( mesh::SideSet * tSideset : mSideSets )
            {
                tSideset->reset_node_container();
            }
            for( mesh::SideSet * tCut : mCuts )
            {
                tCut->reset_node_container();
            }

            // reset the element-to-element links
            for ( mesh::Element * tElement : mElements )
            {
                tElement->reset_element_container();
            }
            for ( mesh::Facet * tFacet : mConnectors )
            {
                tFacet->element()->reset_element_container();
            }

            mElements.clear() ;
            mConnectors.clear() ;
            mFacets.clear() ;

            for ( mesh::Node * tVertex : mNodes )
            {
                tVertex->reset_vertex_container();
                tVertex->reset_node_container();
                tVertex->reset_edge_container();
                tVertex->reset_face_container();
                tVertex->reset_facet_container();
                tVertex->reset_element_container();
            }

            for ( mesh::Edge * tVertex : mEdges )
            {
                tVertex->reset_vertex_container();
                tVertex->reset_node_container();
                tVertex->reset_edge_container();
                tVertex->reset_face_container();
                tVertex->reset_facet_container();
                tVertex->reset_element_container();
            }

            for ( mesh::Face * tVertex : mFaces )
            {
                tVertex->reset_vertex_container();
                tVertex->reset_node_container();
                tVertex->reset_edge_container();
                tVertex->reset_face_container();
                tVertex->reset_facet_container();
                tVertex->reset_element_container();
            }

            // we do NOT reset these links!
            mFacetsAreLinked = true ;

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
        std::cout << "finalize edges" << std::endl ;
        if( mEdges.size() > 0 )
        {
            // flag edge elements
            this->unflag_all_elements() ;

            if( mComputeConnectivities )
            {
                this->connect_nodes_to_edges();
                this->connect_edges_to_elements( aElements );
                this->connect_edges_to_edges();
            }
            this->compute_edge_orientations();
            this->create_edge_map();
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::finalize_faces()
    {
        std::cout << "finalize faces" << std::endl ;
        if( mFaces.size() > 0 )
        {
            this->create_face_map();
            if ( mComputeConnectivities and mEdges.size() > 0 )
            {
                this->connect_faces_to_edges_and_edges_to_faces();
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::unflag_everything()
    {
        this->unflag_all_nodes() ;
        this->unflag_all_edges() ;
        this->unflag_all_faces() ;
        this->unflag_all_facets() ;
        this->unflag_all_connectors() ;
        this->unflag_all_vertices() ;
    }


//------------------------------------------------------------------------------

    void
    Mesh::unflag_all_nodes()
    {
        for( mesh::Node * tNode : mNodes )
        {
            tNode->unflag();
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::unflag_all_edges()
    {
        for( mesh::Edge * tEdge : mEdges )
        {
            tEdge->unflag();
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::unflag_all_faces()
    {
        for( mesh::Face * tFace : mFaces )
        {
            tFace->unflag() ;
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::unflag_all_facets()
    {
        for( mesh::Facet * tFacet : mFacets )
        {
            tFacet->unflag();
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::unflag_all_connectors()
    {
        for( mesh::Facet * tFacet : mConnectors )
        {
            tFacet->unflag();
        }
    }


//------------------------------------------------------------------------------

    void
    Mesh::unflag_all_elements()
    {
        for( mesh::Element * tElement : mElements )
        {
            tElement->unflag();
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::unflag_all_vertices()
    {
        for( mesh::Element * tVertex : mVertices )
        {
            tVertex->unflag();
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::partition( const uint & aNumberOfPartitions,
                     const bool aSetProcOwners )
    {
        if( comm_rank() == mMasterProc )
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
            mesh::Partitioner( this, aNumberOfPartitions, aSetProcOwners );

            // assume that all elements are part of the mesh
            this->unflag_all_elements() ;
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::partition( const uint & aNumberOfPartitions,
                     const Vector< id_t > & aSelectedBlocks,
                     const bool aSetProcOwners,
                     const bool aForceContinuousPartition )
    {
        if( comm_rank() == mMasterProc )
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

            // unflag elements connected to ghosts, since the kernel handles that
            /*if( this->ghost_block_ids().length() > 0 )
            {
                for( id_t tID : mGhostSideSetIDs )
                {
                    Cell< mesh::Facet * > & tFacets = this->sideset( tID )->facets() ;

                    for( mesh::Facet * tFacet : tFacets )
                    {
                        if( tFacet->has_master() )
                        {
                            tFacet->master()->unflag() ;
                        }
                        if( tFacet->has_slave() )
                        {
                            tFacet->slave()->unflag() ;
                        }
                    }
                }
            }*/

            // create a partitioner
            mesh::Partitioner( this, aNumberOfPartitions, aSetProcOwners, aForceContinuousPartition  );

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
                     const bool aForceContinuousPartition  )
    {
        if( comm_rank() == mMasterProc )
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
                               aForceContinuousPartition );

            // assume that all elements are part of the mesh
            this->unflag_all_elements() ;
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::create_edges(  const bool aPrint,
                         const Vector< id_t > aNedelecBlocks,
                         const Vector< id_t > aNedelecSideSets )
    {
        if( comm_rank() == mMasterProc )
        {
            // create the factory
            mesh::EdgeFactory tFactory( this );
            tFactory.create_edges( aNedelecBlocks, aNedelecSideSets );

            // print for debugging
            if( aPrint )
            {
                tFactory.print();
            }
            this->update_edge_indices() ;

            // todo: make this available in parallel
            this->compute_edge_orientations() ;
        }
    }

 //------------------------------------------------------------------------------

    void
    Mesh::reset_edges()
    {
        if( comm_rank() == mMasterProc &&  this->edges_exist() )
        {
            mEdgeMap.clear();
            for ( mesh::Edge * tEdge : mEdges )
            {
                delete tEdge ;
            }
            mEdges.clear();
        }
    }
//------------------------------------------------------------------------------

    void
    Mesh::create_faces(  const bool aPrint,
                         const Vector< id_t > aNedelecBlocks,
                         const Vector< id_t > aNedelecSideSets )
    {
        if( comm_rank() == mMasterProc )
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
    Mesh::reset_faces()
    {
        if( comm_rank() == mMasterProc && this->faces_exist() )
        {
           mFacetMap.clear() ;
           for( mesh::Face * tFace : mFaces )
           {
               delete tFace ;
           }
           mFaces.clear();
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
        for( mesh::Facet * tFacet : mConnectors )
        {
            // add element to map
            mFacetMap[ tFacet->id() ] = tFacet;
        }

        // reset block map
        mBlockMap.clear() ;
        for( mesh::Block * tBlock : mBlocks )
        {
            // add block to map
            mBlockMap[ tBlock->id() ] = tBlock;
        }

        // reset sideset map
        mSideSetMap.clear();

        // loop over all sidesets
        for( mesh::SideSet * tSideSet : mSideSets )
        {

            mSideSetMap[ tSideSet->id() ] = tSideSet;
        }

        // reset cut map
        mCutMap.clear() ;

        // loop over all cuts
        for( mesh::SideSet * tCut : mCuts )
        {
            mCutMap[ tCut->id() ] = tCut;
            mSideSetMap[ tCut->id() ] = tCut ;
        }

        // reset vertex map
        mVertexMap.clear();

        for( mesh::Element * tVertex : mVertices )
        {
            mVertexMap[ tVertex->id() ] = tVertex;
        }


        // map for ghost facets
        mTapeFacetMap.clear() ;
        index_t tNumFacets = mTapeFacetTable.n_cols() ;
        {
            for ( index_t f=0; f<tNumFacets; ++f )
            {
                mTapeFacetMap[ mTapeFacetTable( 0, f ) ]
                    = mFacetMap[ mTapeFacetTable( 1, f ) ];
            }
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
    Mesh::set_node_owners()
    {
        this->unflag_all_nodes() ;
        this->unflag_all_facets();
        this->unflag_all_connectors();

        Cell< mesh::Node * > tNodes ;

        // reset node owners
        for ( mesh::Node * tNode: mNodes )
        {
            tNode->set_owner( mNumberOfPartitions );
        }


        // loop over all nodes
        for ( mesh::Node * tNode: mNodes )
        {
            proc_t tOwner = mNumberOfPartitions;

            // loop over all elements of this node
            for ( uint e = 0; e < tNode->number_of_elements(); ++e )
            {
                tOwner = ( tNode->element( e )->owner() < tOwner ) ?
                         tNode->element( e )->owner() : tOwner;
            }

            if( tNode->number_of_duplicates() == 0 )
            {
                tNode->set_owner( tOwner );
            }
            else if ( ! tNode->is_flagged() )
            {
                for( uint d=0; d<tNode->number_of_duplicates(); ++d )
                {
                    mesh::Node * tDup = tNode->duplicate( d );

                    // loop over all elements of this node
                    for ( uint e = 0; e < tDup->number_of_elements(); ++e )
                    {
                        tOwner = ( tDup->element( e )->owner() < tOwner ) ?
                                 tDup->element( e )->owner() : tOwner;
                    }
                }

                tNode->set_owner( tOwner );
                tNode->flag() ;
                for( uint d=0; d<tNode->number_of_duplicates(); ++d )
                {
                    mesh::Node * tDup = tNode->duplicate( d );
                    tDup->set_owner( tOwner );
                    tDup->flag() ;
                }
            }
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
    Mesh::set_connector_owners()
    {
        for( mesh::Facet * tConnector : mConnectors )
        {
            proc_t tOwner =
                    tConnector->node( 0 )->owner() < tConnector->node( 1 )->owner() ?
                    tConnector->node( 0 )->owner() : tConnector->node( 1 )->owner() ;

            tConnector->set_owner( tOwner );
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
    Mesh::flag_curved_elements()
    {
        Timer tTimer ;

        proc_t tRank = comm_rank() ;

        if( tRank == 0 )
        {
            message( 4, "\n    Flagging curved elements ...\n"  );
        }
        comm_barrier() ;

        // check curved elements
        mesh::CurvedElementChecker tChecker( mNumberOfDimensions, mBlocks, mSideSets );

        index_t tCount = tChecker.flag_curved_elements();

        if( tRank == 0 )
        {
            message( 4,
                     "    ... time for searching curved elements : %u ms. Elements found: %lu\n",
                     ( unsigned int ) tTimer.stop(),  ( long unsigned int ) tCount  );
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::create_ghost_sidesets(
            const Vector< id_t >    & aGhostBlockIDs,
            const Vector< id_t >    & aGhostSideSetIDs,
            const Vector< id_t >    & aElementIDs,
            Cell< mesh::Layer  * >  & aLayers )
    {
        BELFEM_ASSERT( this->master() == comm_rank(),
                       "this function mut be called by master proc only" );

        // initialize counter
        index_t tCount ;

        // - - - - - - - - - - - - - - - - - - - - - -
        // step 1: reorganize the nodes
        // - - - - - - - - - - - - - - - - - - - - - -

        // count entities
        tCount = mNodes.size() ;

        for( mesh::Layer * tLayer : aLayers )
        {
            tCount += tLayer->Nodes.size() ;
        }

        // backup nodes
        Cell< mesh::Node * > tNodes ;
        tNodes.vector_data() = std::move(mNodes.vector_data() );
        mNodes.set_size( tCount, nullptr );

        // repopulate node list
        // repopulate nodes
        tCount = 0 ;
        for( mesh::Node * tNode : tNodes )
        {
            mNodes( tCount ) = tNode ;
            tNode->set_index( tCount++ ) ;
        }

        tNodes.clear() ; // tidy up memory
        for( mesh::Layer * tLayer : aLayers )
        {
            for ( mesh::Node * tNode: tLayer->Nodes )
            {
                mNodes( tCount ) = tNode;
                tNode->set_index( tCount++ );
            }
        }


        // - - - - - - - - - - - - - - - - - - - - - -
        // step 2: reorganize edges
        // - - - - - - - - - - - - - - - - - - - - - -

        // count new edges
        tCount = 0 ;
        for( mesh::Layer * tLayer : aLayers )
        {
            tCount += tLayer->Edges.size() ;
        }

        // backup edges
        if( mEdges.size() > 0 )
        {
            tCount += mEdges.size() ;
            Cell< mesh::Edge * > tEdges ;
            tEdges.vector_data() = std::move(mEdges.vector_data() );
            mEdges.set_size( tCount, nullptr );

            // repopulate edges
            tCount = 0 ;
            for( mesh::Edge * tEdge : tEdges )
            {
                mEdges( tCount ) = tEdge ;
                tEdge->set_index( tCount++ ) ;
            }
        }
        else
        {
            mEdges.set_size( tCount, nullptr );
            tCount = 0 ;
        }

        // add new edges to container
        for( mesh::Layer * tLayer : aLayers )
        {
            for ( mesh::Edge * tEdge: tLayer->Edges )
            {
                mEdges( tCount ) = tEdge;
                tEdge->set_index( tCount++ );
            }
        }

        // - - - - - - - - - - - - - - - - - - - - - -
        // step 3: reorganize faces
        // - - - - - - - - - - - - - - - - - - - - - -
        if( this->number_of_dimensions() == 3 )
        {
            // count new faces
            tCount = 0 ;
            for( mesh::Layer * tLayer : aLayers )
            {
                tCount += tLayer->Faces.size() ;
            }

            if( tCount > 0 )
            {
                if ( mFaces.size() > 0 )
                {
                    tCount += mFaces.size();
                    Cell< mesh::Face * > tFaces;
                    tFaces.vector_data() = std::move( mFaces.vector_data());
                    mFaces.set_size( tCount, nullptr );
                    // repopulate faces
                    tCount = 0;
                    for ( mesh::Face * tFace: tFaces )
                    {
                        mFaces( tCount ) = tFace;
                        tFace->set_index( tCount++ );
                    }
                }
                else
                {
                    mFaces.set_size( tCount, nullptr );
                    tCount = 0;
                }
                for ( mesh::Layer * tLayer: aLayers )
                {
                    for ( mesh::Face * tFace: tLayer->Faces )
                    {
                        mFaces( tCount ) = tFace;
                        tFace->set_index( tCount++ );
                    }
                }
            } // else new faces  exist
        } // end if 3D

        // remember IDs
        mGhostBlockIDs   = aGhostBlockIDs ;
        mGhostSideSetIDs = aGhostSideSetIDs ;
        mGhostFacetIDs   = aElementIDs ;

        // loop over all layers
        for( uint l=0; l<aGhostSideSetIDs.length(); ++l )
        {
            // grab element container
            Cell< mesh::Element * > & tElements = aLayers( l )->Elements ;

            index_t tNumFacets = tElements.size() ;

            // create the new sideset
            mesh::SideSet * tSideSet = new mesh::SideSet( aGhostSideSetIDs( l ),
                                                          tNumFacets );

            // grab the facet container
            Cell< mesh::Facet * > & tFacets = tSideSet->facets();

            // loop over all elements
            for( index_t e = 0; e<tNumFacets; ++e )
            {
                // grab original element
                mesh::Facet * tOriginal = this->facet( aElementIDs( e ) );

                // create the new facet
                mesh::Facet * tFacet = new mesh::Facet( tElements( e ) );

                if( tOriginal->has_master() )
                {
                    tFacet->set_master( tOriginal->master(), tOriginal->master_index() );
                }
                if( tOriginal->has_slave() )
                {
                    tFacet->set_slave( tOriginal->slave(), tOriginal->slave_index() );
                }
                tFacets( e ) = tFacet ;

            } // end loop over all facets

            // hide this sideset from exodus
            tSideSet->hide();

            mSideSets.push( tSideSet );

        } // end loop over all layers

        this->unfinalize() ;
        this->finalize() ;

        if( this->edges_exist() )
        {
            this->finalize_edges( );
        }
        if( this->faces_exist() )
        {
            this->finalize_faces() ;
        }

        // create the map
        tCount = 0 ;
        mGhostFacetMap.clear() ;
        for( id_t tElementID : aElementIDs )
        {
            mGhostFacetMap[ tElementID ] = tCount++ ;
        }

        this->link_ghost_elements();
    }

//------------------------------------------------------------------------------

    void
    Mesh::distribute_ghost_sidesets( const proc_t aTarget, const proc_t aMasterProc )
    {
        proc_t tNumProcs = comm_size() ;

        if( tNumProcs > 1 && aTarget != aMasterProc )
        {
            if ( comm_rank() == aMasterProc  )
            {
                send( aTarget, mGhostBlockIDs );
                send( aTarget, mGhostSideSetIDs ) ;

                if( mGhostSideSetIDs.length() > 0 )
                {

                    // count facets per proc
                    index_t tCount = 0 ;

                    for ( id_t tID: mGhostFacetIDs )
                    {
                        if( this->facet( tID )->owner() == aTarget )
                        {
                            ++tCount ;
                        }
                    }

                    // create data that need to be sent
                    Vector< id_t > tFacetIDs( tCount );


                    // reset the counter
                    tCount = 0 ;

                    for ( id_t tID: mGhostFacetIDs )
                    {
                        if( this->facet( tID )->owner() == aTarget )
                        {
                            tFacetIDs( tCount++ ) = tID ;
                        }
                    }

                    uint tNumLayers = this->number_of_thin_shell_layers();

                    Matrix< id_t > tTapeFacetTable( 2, tNumLayers * tCount );

                    tCount = 0 ;
                    for( uint l=0; l<tNumLayers; ++l )
                    {
                        for ( id_t tID: tFacetIDs )
                        {
                            tTapeFacetTable( 0, tCount )
                                = this->ghost_facet( tID, l )->id() ;
                            tTapeFacetTable( 1, tCount++ ) = tID ;
                        }
                    }

                    // send data
                    send( aTarget, tFacetIDs );
                    send( aTarget, tTapeFacetTable );
                }
            }
            else
            {
                receive( aMasterProc, mGhostBlockIDs );
                receive( aMasterProc, mGhostSideSetIDs );

                // create the map
                mGhostFacetMap.clear() ;

                if( mGhostSideSetIDs.length() > 0 )
                {
                    receive( aMasterProc, mGhostFacetIDs );
                    receive( aMasterProc, mTapeFacetTable );

                    index_t tCount = 0;
                    for ( id_t tID: mGhostFacetIDs )
                    {
                        mGhostFacetMap[ tID ] = tCount++;
                    }
                }
            }
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::set_ghost_information(
            const Vector< id_t > & aGhostBlocks,
            const Vector< id_t > & aGhostSideSets )
    {
        mGhostBlockIDs   = aGhostBlocks ;
        mGhostSideSetIDs = aGhostSideSets ;
    }

//------------------------------------------------------------------------------

    void
    Mesh::link_ghost_elements()
    {
        // check if ghosts exist
        if( mGhostBlockIDs.length() > 0 )
        {
            // check if this proc has ghosts
            if ( this->block_exists( mGhostBlockIDs( 0 ) ) )
            {
                // get order in thickness direction
                uint tOrder = ( mGhostSideSetIDs.length() - 1 ) / mGhostBlockIDs.length() ;

                uint tNumLayersPerElement = tOrder + 1 ;
                uint tOff = 0 ;

                // loop over all blocks
                for( id_t tBlockID : mGhostBlockIDs )
                {
                    // get the elements on this block
                    Cell< mesh::Element * > & tElements = this->block( tBlockID )->elements() ;

                    // loop over all elements and allocate the vertex container
                    for( mesh::Element * tElement : tElements )
                    {
                        tElement->allocate_facet_container( tNumLayersPerElement );
                    }

                    // loop over all sidesets
                    for( uint l=0; l<tNumLayersPerElement; ++l )
                    {
                        // get  the sideset
                        Cell< mesh::Facet * > & tFacets = this->sideset(
                                mGhostSideSetIDs( tOff + l ) )->facets() ;

                        BELFEM_ASSERT( tElements.size() == tFacets.size(),
                                       "Length of ghost sideset and ghost block does not match");

                        // initialize counter
                        index_t tCount = 0 ;

                        // loop over all elements and add the facets
                        for( mesh::Element * tElement : tElements )
                        {
                            tElement->insert_facet( tFacets( tCount++ ) , l );
                        }
                    }

                    // increment offset
                    tOff += tOrder ;
                }
            }

            this->fix_ghost_adjacency();
        }
    }

//------------------------------------------------------------------------------

    void
    Mesh::fix_ghost_adjacency()
    {
        this->unflag_all_elements();
        this->unflag_all_nodes();

        // collect elements
        index_t tNumElems = 0 ;
        for( id_t tID : mGhostBlockIDs )
        {
            if( this->block_exists( tID ) )
            {
                tNumElems += this->block( tID )->number_of_elements() ;
            }
        }

        if( tNumElems == 0 )
        {
            return;
        }

        // collect all elements
        Cell< mesh::Element * > tElements( tNumElems, nullptr );

        index_t tCount = 0 ;
        for( id_t tID : mGhostBlockIDs )
        {
            if( this->block_exists( tID ) )
            {
               for( mesh::Element * tElement : this->block( tID )->elements() )
               {
                   tElement->flag() ;
                   tElement->flag_nodes();
                   tElement->set_index( tCount );
                   tElements( tCount++ ) = tElement ;
               }
            }
        }

        // count nodes
        tCount = 0 ;
        for( mesh::Node * tNode : mNodes )
        {
            if( tNode->is_flagged() )
            {
                tNode->set_index( tCount++ );
            }
        }

        // count nodes per elements
        Vector< uint > tNumElemsPerNode( tCount, 0 );

        for( mesh::Element * tElement : tElements )
        {
            for( uint k=0; k<tElement->number_of_nodes(); ++k )
            {
                ++tNumElemsPerNode( tElement->node( k )->index() );
            }
        }

        uint tMaxNumElems = max( tNumElemsPerNode );

        // allocate node index container
        Matrix< index_t > tElemsPerNode( tMaxNumElems, tCount , 0 );

        tNumElemsPerNode.fill( 0 );
        for( mesh::Element * tElement : tElements )
        {
            for( uint k=0; k<tElement->number_of_nodes(); ++k )
            {
                index_t tIndex = tElement->node( k )->index() ;

                tElemsPerNode( tNumElemsPerNode( tIndex )++, tIndex ) = tElement->index();
            }
        }

        for( mesh::Element * tElement : tElements )
        {
            // loop over all nodes and unflag connected elements
            for( uint k=0; k<tElement->number_of_nodes(); ++k )
            {
                index_t tIndex = tElement->node( k )->index() ;

                uint tN = tNumElemsPerNode( tIndex );
                for( uint e=0; e<tN; ++e )
                {
                    tElements( tElemsPerNode( e, tIndex ) )->unflag() ;
                }
            }

            // flag myself
            tElement->flag() ;

            // count connected elements
            tCount = 0 ;
            for( uint k=0; k<tElement->number_of_nodes(); ++k )
            {
                index_t tIndex = tElement->node( k )->index() ;

                uint tN = tNumElemsPerNode( tIndex );
                for( uint e=0; e<tN; ++e )
                {
                    mesh::Element * tNeighbor = tElements( tElemsPerNode( e, tIndex ) );
                    if( ! tNeighbor->is_flagged() )
                    {
                        ++tCount ;
                        tNeighbor->flag() ;
                    }
                }
            }

            // allocate element container
            tElement->allocate_element_container( tCount );


            // unflag myself
            tElement->unflag() ;

            for( uint k=0; k<tElement->number_of_nodes(); ++k )
            {
                index_t tIndex = tElement->node( k )->index() ;

                uint tN = tNumElemsPerNode( tIndex );
                for( uint e=0; e<tN; ++e )
                {
                    mesh::Element * tNeighbor = tElements( tElemsPerNode( e, tIndex ) );
                    if( tNeighbor->is_flagged() )
                    {
                        tElement->insert_element( tNeighbor );
                        tNeighbor->unflag();
                    }
                }
            }
        }

        // tidy up
        this->unflag_all_elements();
        this->unflag_all_nodes();
        this->update_node_indices();
        this->update_element_indices();
    }

//------------------------------------------------------------------------------

    void
    Mesh::compute_edge_orientations()
    {
        Cell< mesh::Node * > tNodes ;
        for( mesh::Element * tElement : mElements )
        {
            for( uint e=0; e<tElement->number_of_edges(); ++e )
            {
                tElement->get_nodes_of_edge( e, tNodes );
                id_t tA = tNodes( 0 )->id() ;
                id_t tB = tNodes( 1 )->id() ;

                id_t tC = tElement->edge( e )->node( 0 )->id() ;
                id_t tD = tElement->edge( e )->node( 1 )->id() ;

                if( tA == tC && tB == tD )
                {
                    tElement->set_edge_orientation( e, true );
                }
                else if ( tA == tD && tB == tC )
                {
                    tElement->set_edge_orientation( e, false );
                }
                else
                {
                    BELFEM_ERROR( false, "Internal error at element %lu: invalid edge %u",
                                  ( long unsigned int ) tElement->id(),
                                  ( unsigned  int ) e );
                }
            }
        }
    }

//------------------------------------------------------------------------------

    index_t
    Mesh::check()
    {
        index_t aCount = 0;
        for( mesh::Block * tBlock : mBlocks )
        {

            mesh::OrientationChecker * tChecker = new mesh::OrientationChecker();
            tChecker->set_element_type( tBlock->element_type() );

            for ( mesh::Element * tElement: mElements )
            {
                aCount += tChecker->process_element( tElement );
            }

            delete tChecker ;
        }

        return aCount ;
    }

}
