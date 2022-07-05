//
// Created by Christian Messe on 27.10.19.
//

#include "HDF5_Tools.hpp"
#include "meshtools.hpp"
#include "cl_Mesh_HDF5Writer.hpp"
#include "cl_Timer.hpp"
#include "cl_Logger.hpp"
namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        HDF5Writer::HDF5Writer( const string & aFilePath, Mesh * aMesh   ) :
        mFilePath( aFilePath ),
        mFile( aFilePath, FileMode::NEW ),
        mMesh( aMesh )
        {


            mMesh->update_node_indices();
            mMesh->update_element_indices();

            this->save_nodes();
            this->save_edges();
            this->save_elements();
            this->save_sidesets();
            this->save_vertices();
            this->save_globals();
            this->save_fields();

            mFile.close();
        }

//------------------------------------------------------------------------------

        void
        HDF5Writer::save_nodes()
        {
#ifdef BELFEM_HDF5
            index_t tNumberOfNodes = mMesh->number_of_nodes();

            mFile.create_group( "Nodes");

            mFile.save_data( "NumberOfNodes", tNumberOfNodes );

            // populate node coords
            Matrix< real > tCoords( 3, tNumberOfNodes );

            index_t tCount = 0;

            for( mesh::Node * tNode : mMesh->nodes() )
            {
                tCoords( 0, tCount ) = tNode->x();
                tCoords( 1, tCount ) = tNode->y();
                tCoords( 2, tCount ) = tNode->z();
                ++tCount;
            }

            mFile.save_data( "NodeCoords", tCoords );


            // populate node IDs
            Vector< id_t > tNodeIDs( tNumberOfNodes );

            tCount = 0;
            for( Node * tNode : mMesh->nodes() )
            {
                tNodeIDs( tCount++ ) = tNode->id();
            }

            mFile.save_data( "NodeID", tNodeIDs );


            Vector< uint > tNodeOwners( tNumberOfNodes );

            tCount = 0;
            for( Node * tNode : mMesh->nodes() )
            {
                tNodeOwners( tCount++ ) = tNode->owner();
            }

            mFile.save_data( "NodeOwner", tNodeOwners );

            mFile.close_active_group();
#endif
        }

//------------------------------------------------------------------------------

        void
        HDF5Writer::save_edges()
        {
#ifdef BELFEM_HDF5
            if( mMesh->edges_exist() )
            {
                index_t tNumberOfEdges = mMesh->number_of_edges();

                mFile.create_group( "Edges" );
                mFile.save_data( "NumberOfEdges", tNumberOfEdges );

                // populate node IDs
                Vector< id_t > tEdgeIDs( tNumberOfEdges );

                index_t tCount = 0;
                for( Edge * tEdge : mMesh->edges() )
                {
                    tEdgeIDs( tCount++ ) = tEdge->id() ;
                }
                mFile.save_data( "EdgeID", tEdgeIDs );

                Vector< uint > tEdgeOwners( tNumberOfEdges );

                tCount = 0;
                for( Edge * tEdge : mMesh->edges() )
                {
                    tEdgeOwners( tCount++ ) = tEdge->owner();
                }
                mFile.save_data( "EdgeOwner", tEdgeOwners );

                uint tNumberOfNodes = mMesh->max_element_order() + 1 ;
                Matrix< index_t > tNodeIndices( tNumberOfNodes, tNumberOfEdges );

                tCount = 0;
                for( Edge * tEdge : mMesh->edges() )
                {
                    for( uint k=0; k<tNumberOfNodes; ++k )
                    {
                        tNodeIndices( k, tCount ) = tEdge->node( k )->index() ;
                    }
                    ++tCount ;
                }
                mFile.save_data( "EdgeNodeConnectivity", tNodeIndices );

                mFile.close_active_group();
            }
#endif
        }

//------------------------------------------------------------------------------

        void
        HDF5Writer::save_faces()
        {
#ifdef BELFEM_HDF5
            if( mMesh->faces_exist() )
            {
                index_t tNumberOfFaces = mMesh->number_of_faces();

                mFile.create_group( "Faces" );
                mFile.save_data( "NumberOfFaces", tNumberOfFaces );

                // populate node IDs
                Vector< id_t > tFaceIDs( tNumberOfFaces );

                index_t tCount = 0;
                for( Face * tFace : mMesh->faces() )
                {
                    tFaceIDs( tCount++ ) = tFace->id() ;
                }
                mFile.save_data( "FaceID", tFaceIDs );

                Matrix< index_t > tElement( 2, tNumberOfFaces );

                tCount = 0;
                for( Face * tFace : mMesh->faces() )
                {
                    if( tFace->master() == nullptr )
                    {
                        tElement( 0, tCount ) = gNoIndex ;
                    }
                    else
                    {
                        tElement( 0, tCount ) = tFace->master()->index();
                    }
                    if( tFace->slave() == nullptr )
                    {
                        tElement( 1, tCount++ ) = gNoIndex ;
                    }
                    else
                    {
                        tElement( 1, tCount++ ) = tFace->slave()->index();
                    }
                }
                mFile.save_data( "FaceElementConnectivity", tElement );


                Matrix< uint > tIndex( 2, tNumberOfFaces );

                tCount = 0;
                for( Face * tFace : mMesh->faces() )
                {
                    if( tFace->master() == nullptr )
                    {
                        tIndex( 0, tCount ) =  BELFEM_UINT_MAX ;
                    }
                    else
                    {
                        tIndex( 0, tCount ) = tFace->index_on_master() ;
                    }
                    if( tFace->slave() == nullptr )
                    {
                        tIndex( 1, tCount++ ) =  BELFEM_UINT_MAX ;
                    }
                    else
                    {
                        tIndex( 1, tCount++ ) = tFace->index_on_slave() ;
                    }
                }

                mFile.save_data( "FaceIndexOnElement", tIndex );

                mFile.close_active_group();
            }
#endif
        }

//------------------------------------------------------------------------------

        void
        HDF5Writer::save_elements()
        {
#ifdef BELFEM_HDF5
            uint tNumberOfBlocks = mMesh->number_of_blocks();

            mFile.save_data( "NumberOfBlocks", tNumberOfBlocks );

            Cell< mesh::Block * > & tBlocks = mMesh->blocks();

            for( uint b=0; b<tNumberOfBlocks; ++b )
            {
                Block * tBlock = tBlocks( b );

                mFile.create_group( sprint( "Block_%u", b+1 ) );

                mFile.save_data( "Label", tBlock->label() );
                mFile.save_data( "BlockID", tBlock->id() );

                index_t tNumberOfElements = tBlock->number_of_elements();

                mFile.save_data( "NumberOfElements", tNumberOfElements );

                mFile.save_data( "ElementType", ( uint ) tBlock->element_type() );

                // get number of nodes per element
                index_t tNumNodesPerElement =number_of_nodes( tBlock->element_type() );

                Matrix< index_t > tNodeConnectivity( tNumNodesPerElement, tNumberOfElements );

                for( index_t e=0; e<tNumberOfElements; ++e )
                {
                    // get pointer to element
                    Element * tElement = tBlock->element( e );

                    for( index_t k=0; k<tNumNodesPerElement; ++k )
                    {
                        tNodeConnectivity( k, e ) = tElement->node( k )->index();
                    }
                }

                mFile.save_data( "ElementNodeConnectivity", tNodeConnectivity );

                tBlock->element( 0 )->has_edges() ;

                // check if edges exist on this block
                bool tEdgesExist = false ;
                if( tBlock->number_of_elements() > 0 )
                {
                    tEdgesExist = tBlock->element( 0 )->has_edges() ;
                }

                if( tEdgesExist )
                {
                    uint tNumEdgesPerElement = mesh::number_of_edges( tBlock->element_type() );
                    Matrix< index_t > tEdgeConnectivity( tNumEdgesPerElement, tNumberOfElements );

                    for( index_t e=0; e<tNumberOfElements; ++e )
                    {
                        // get pointer to element
                        Element * tElement = tBlock->element( e );

                        for( index_t k=0; k<tNumEdgesPerElement; ++k )
                        {
                            tEdgeConnectivity( k, e ) = tElement->edge( k )->index();
                        }
                    }

                    mFile.save_data( "ElementEdgeConnectivity", tEdgeConnectivity );
                }

                // check if faces exist on this block
                bool tFacesExist = false ;
                if( tBlock->number_of_elements() > 0 )
                {
                    tFacesExist = tBlock->element( 0 )->has_faces() ;
                }

                if( tFacesExist )
                {
                    uint tNumFacesPerElement = mesh::number_of_faces( tBlock->element_type() );
                    Matrix< index_t > tFaceConnectivity( tNumFacesPerElement, tNumberOfElements );

                    for( index_t e=0; e<tNumberOfElements; ++e )
                    {
                        // get pointer to element
                        Element * tElement = tBlock->element( e );

                        for( index_t k=0; k<tNumFacesPerElement; ++k )
                        {
                            tFaceConnectivity( k, e ) = tElement->face( k )->index();
                        }
                    }
                    mFile.save_data( "ElementFaceConnectivity", tFaceConnectivity );
                }

                Vector< uint > tData( tNumberOfElements );
                Vector< id_t > tIDs( tNumberOfElements );

                // save geometry tag
                for( index_t e=0; e<tNumberOfElements; ++e )
                {
                    tData( e ) = tBlock->element( e )->geometry_tag();
                }
                mFile.save_data( "GeometryTag", tData );

                // save physical tag
                for( index_t e=0; e<tNumberOfElements; ++e )
                {
                    tData( e ) = tBlock->element( e )->physical_tag();
                }
                mFile.save_data( "PhysicalTag", tData );

                // save owners
                for( index_t e=0; e<tNumberOfElements; ++e )
                {
                    tData( e ) = tBlock->element( e )->owner();
                }
                mFile.save_data( "ElementOwner", tData );

                // save ids
                for( index_t e=0; e<tNumberOfElements; ++e )
                {
                    tIDs( e ) = tBlock->element( e )->id();
                }
                mFile.save_data( "ElementID", tIDs );

                mFile.close_active_group();
            }
        }
#endif

//------------------------------------------------------------------------------

        void
        HDF5Writer::save_sidesets()
        {
#ifdef BELFEM_HDF5
            uint tNumberOfSidesets = mMesh->number_of_sidesets();
            index_t tNumberOfElements = mMesh->number_of_elements();

            mFile.save_data( "NumberOfSideSets", tNumberOfSidesets );

            Cell< SideSet * > & tSideSets = mMesh->sidesets() ;

            for( uint s=0; s<tNumberOfSidesets; ++s )
            {
                SideSet * tSideSet = tSideSets( s );

                mFile.create_group( sprint( "SideSet_%u", s+1 ) );

                index_t tNumberOfFacets = tSideSet->number_of_facets();

                mFile.save_data( "Label", tSideSet->label() );
                mFile.save_data( "SideSetID", tSideSet->id() );
                mFile.save_data( "NumberOfFacets", tNumberOfFacets );

                Vector< id_t > tIDs( tNumberOfFacets );
                Vector< index_t > tElements( tNumberOfFacets );

                // save IDs
                for( index_t f=0; f<tNumberOfFacets; ++f )
                {
                    tIDs( f ) = tSideSet->facet_by_index( f )->id();
                }

                mFile.save_data( "FacetID", tIDs );

                for( index_t f=0; f<tNumberOfFacets; ++f )
                {
                    tElements( f ) = tSideSet->facet_by_index( f )->master()->index();
                }

                mFile.save_data( "FacetOwner", tElements );

                for( index_t f=0; f<tNumberOfFacets; ++f )
                {
                    if( tSideSet->facet_by_index( f )->has_slave() )
                    {
                        tElements( f ) = tSideSet->facet_by_index( f )->slave()->index();
                    }
                    else
                    {
                        tElements( f ) = tNumberOfElements;
                    }
                }
                mFile.save_data( "FacetNeighbor", tElements );

                Vector< uint > tData( tNumberOfFacets );

                for( index_t f=0; f<tNumberOfFacets; ++f )
                {
                    tData( f ) = tSideSet->facet_by_index( f )->master_index();
                }
                mFile.save_data( "FacetOwnerIndex", tData );

                for( index_t f=0; f<tNumberOfFacets; ++f )
                {
                    tData( f ) = tSideSet->facet_by_index( f )->slave_index();
                }
                mFile.save_data( "FacetNeighborIndex", tData );
            }

            mFile.close_active_group();
#endif
        }
//------------------------------------------------------------------------------

        void
        HDF5Writer::save_vertices()
        {
#ifdef BELFEM_HDF5
            // get link to vertex container
            Cell< mesh::Element * > & tVertices = mMesh->vertices();

            // get number of vertices
            index_t tNumberOfVertices = tVertices.size();

            mFile.save_data( "NumberOfVertices", tNumberOfVertices );

            // test if any vertices are defined
            if( tNumberOfVertices > 0 )
            {
                // create a group

                // create a container for the element ids
                Vector< id_t > tIDs( tNumberOfVertices );

                // create a contaner for element indices
                Vector< index_t > tIndices( tNumberOfVertices );

                // Tags
                Vector< uint > tGeometryTags( tNumberOfVertices );
                Vector< uint > tPhysicalTags( tNumberOfVertices );
                Vector< uint > tProcOwners( tNumberOfVertices );

                // write element IDs ito vertex
                for( index_t k=0; k<tNumberOfVertices; ++k )
                {
                    tIDs( k )          = tVertices( k )->id();
                    tIndices( k )      = tVertices( k )->node( 0 )->index();
                    tGeometryTags( k ) = tVertices( k )->geometry_tag();
                    tPhysicalTags( k ) = tVertices( k )->physical_tag();
                    tProcOwners( k )   = tVertices( k )->owner();
                }

                mFile.create_group( "Vertices" );

                // save ids into file
                mFile.save_data( "VertexID", tIDs );

                // save node index into file
                mFile.save_data( "NodeIndex", tIndices );

                // save physical Tags
                mFile.save_data( "GeometryTag", tGeometryTags );

                // save node index into file
                mFile.save_data( "PhysicalTag", tPhysicalTags );

                // save owners
                mFile.save_data( "ElementOwner", tProcOwners );

                // loop over all vertices
                mFile.close_active_group();
            }
#endif
        }

//------------------------------------------------------------------------------
        void
        HDF5Writer::save_globals()
        {
#ifdef BELFEM_HDF5
            mFile.save_data( "TimeStamp", mMesh->time_stamp() );
            mFile.save_data( "TimeStep", mMesh->time_step() );
            mFile.save_data( "NumberOfDimensions", mMesh->number_of_dimensions() );

            uint mNumberOfGlobals = mMesh->number_of_global_variables();

            mFile.save_data( "NumberOfGlobals", mNumberOfGlobals );

            for( uint g=0; g<mNumberOfGlobals; ++g )
            {
                mFile.create_group( sprint( "Global_%u", g+1 ) );
                mFile.save_data( "Label", mMesh->global_variable( g )->label() );
                mFile.save_data( "Value", mMesh->global_variable( g )->value() );
                mFile.save_data( "GlobalID", mMesh->global_variable( g )->id() );
            }

            mFile.close_active_group();
#endif
        }

//------------------------------------------------------------------------------

        void
        HDF5Writer::save_fields()
        {
#ifdef BELFEM_HDF5
            uint tNumberOfFields = mMesh->number_of_fields();

            mFile.save_data( "NumberOfFields", tNumberOfFields );

            for( uint f=0; f<tNumberOfFields; ++f )
            {
                mFile.create_group( sprint( "Field_%u", f+1 ) );
                mFile.save_data( "FieldID", mMesh->field( f )->id() );
                mFile.save_data( "Label", mMesh->field( f )->label() );
                mFile.save_data( "Data", mMesh->field( f )->data() );
            }
            mFile.close_active_group();
#endif
        }
    }
}