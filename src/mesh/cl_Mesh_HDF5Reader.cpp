//
// Created by Christian Messe on 27.10.19.
//

#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "HDF5_Tools.hpp"
#include "meshtools.hpp"
#include "cl_Mesh_HDF5Reader.hpp"
#include "cl_Cell.hpp"
#include "cl_Element_Factory.hpp"
#include "commtools.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        HDF5Reader::HDF5Reader( const string & aFilePath, Mesh * aMesh ) :
                mFilePath( aFilePath ),
                mFile( aFilePath, FileMode::OPEN_RDONLY )
        {
            this->read_header();

            if( aMesh == nullptr )
            {
                mOwnMesh = true;
                mMesh = new Mesh( mNumberOfDimensions, comm_rank() );
            }
            else
            {
                mMesh = aMesh;
                mMesh->set_number_of_dimensions( mNumberOfDimensions );
            }

            this->read_nodes();
            this->read_blocks();
            this->read_sidesets();
            this->read_vertices();
            this->read_globals();
            this->read_fields();

            // finalize mesh
            mMesh->finalize();
        }

//------------------------------------------------------------------------------

        HDF5Reader::~HDF5Reader()
        {
            if( mOwnMesh )
            {
                delete mMesh;
            }
        }

//------------------------------------------------------------------------------

        Mesh *
        HDF5Reader::mesh()
        {
            mOwnMesh = false;
            return mMesh;
        }

//------------------------------------------------------------------------------

        void
        HDF5Reader::read_header()
        {
#ifdef BELFEM_HDF5
            mFile.load_data( "NumberOfDimensions", mNumberOfDimensions );
            mFile.load_data( "NumberOfBlocks", mNumberOfBlocks );
            mFile.load_data( "NumberOfSideSets", mNumberOfSideSets );
            mFile.load_data( "NumberOfGlobals", mNumberOfGlobals );
            mFile.load_data( "NumberOfFields", mNumberOfFields );
            mFile.load_data( "TimeStep", mTimeStep );
            mFile.load_data( "TimeStamp", mTimeStamp );
#endif
        }

//------------------------------------------------------------------------------

        void
        HDF5Reader::read_nodes()
        {
#ifdef BELFEM_HDF5
            mFile.select_group( "Nodes" );
            mFile.load_data( "NumberOfNodes", mNumberOfNodes );

            // create the coordinate array
            Matrix< real > tCoords;

            // data container for various stuff
            Vector< index_t > tData;

            // read coords
            mFile.load_data( "NodeCoords", tCoords );

            // read node IDs
            mFile.load_data( "NodeID", tData );

            // get the node container
            Cell< Node * > & tNodes = mMesh->nodes();

            // allocate memory
            tNodes.set_size( mNumberOfNodes, nullptr );

            // create nodes
            for( index_t k=0; k<mNumberOfNodes; ++k )
            {
                tNodes( k ) = new Node(
                        tData( k ),
                        tCoords( 0, k ),
                        tCoords( 1, k ),
                        tCoords( 2, k ) );
            }

            // read node owners
            mFile.load_data( "NodeOwner", tData );
            for( index_t k=0; k<mNumberOfNodes; ++k )
            {
                tNodes( k )->set_owner( tData( k ) );
            }

            // set the indices of the nodes in the memory
            mMesh->update_node_indices();
#endif
        }

//------------------------------------------------------------------------------

        void
        HDF5Reader::read_blocks()
        {
#ifdef BELFEM_HDF5
            // nodes of the mesh
            Cell < Node * > & tNodes = mMesh->nodes();

            // the factory that produces elements
            ElementFactory tFactory;

            // get the block container from the mesh
            Cell< Block * > & tBlocks = mMesh->blocks();

            // allocate memory
            tBlocks.set_size( mNumberOfBlocks, nullptr );

            // loop over all blocks
            for( uint b=0; b<mNumberOfBlocks; ++b )
            {
                // select the group
                mFile.select_group( sprint("Block_%u", b+1 ) );

                // read the id of the block
                id_t tID;
                mFile.load_data( "BlockID", tID );

                // read the number of elements of this block
                index_t tNumberOfElements;
                mFile.load_data("NumberOfElements", tNumberOfElements );

                // read the element type of this block
                uint tType;
                mFile.load_data( "ElementType", tType );
                ElementType tElementType = ( ElementType ) tType;

                // create the block
                Block * tBlock = new Block( tID, tNumberOfElements );

                // read label
                mFile.load_data( "Label", tBlock->label() );

                // read element ids
                Vector< id_t > tIDs;
                mFile.load_data( "ElementID", tIDs );

                // read the connectivity
                Matrix< index_t > tConnectivity;
                mFile.load_data( "ElementNodeConnectivity", tConnectivity );

                // how many nodes does an element have
                uint tNumNodesPerElement = number_of_nodes( tElementType );

                for( index_t e=0; e<tNumberOfElements; ++e )
                {
                    // create the new element
                    Element * tElement
                        = tFactory.create_element( tElementType, tIDs( e ) );

                    // link element with nodes
                    for( uint k=0; k<tNumNodesPerElement; ++k )
                    {
                        tElement->insert_node( tNodes( tConnectivity( k, e ) ), k );
                    }

                    tBlock->insert_element( tElement );
                }

                Vector< uint > tData;

                // write geometry tag
                mFile.load_data("GeometryTag", tData );
                for( index_t e=0; e<tNumberOfElements; ++e )
                {
                    tBlock->element( e )->set_geometry_tag( tData( e ) );
                }

                // write physical tag
                mFile.load_data( "PhysicalTag", tData );
                for( index_t e=0; e<tNumberOfElements; ++e )
                {
                    tBlock->element( e )->set_physical_tag( tData( e ) );
                }

                // write owner
                mFile.load_data( "ElementOwner", tData );
                for( index_t e=0; e<tNumberOfElements; ++e )
                {
                    tBlock->element( e )->set_owner( tData( e ) );
                }

                // add block to container
                tBlocks( b ) = tBlock;
            }

            mFile.close_active_group();

            // with the blocks read, we collocate the elements and index them
            // ( this is needed for the sidesets to work )
            mMesh->collect_elements_from_blocks();
            mMesh->update_element_indices();
            mMesh->update_facet_indices();
#endif
        }

//------------------------------------------------------------------------------

        void
        HDF5Reader::read_sidesets()
        {
#ifdef BELFEM_HDF5
            // the factory that produces elements
            ElementFactory tFactory;

            // get the sideset container
            Cell< SideSet * > & tSideSets = mMesh->sidesets();

            // list with all elements on mesh
            Cell< Element * > & tElements = mMesh->elements();

            index_t tNumberOfElements = tElements.size();

            // reserve the memory
            tSideSets.set_size( mNumberOfSideSets, nullptr );

            for( uint s=0; s<mNumberOfSideSets; ++s )
            {
                // select the group
                mFile.select_group( sprint( "SideSet_%u", s+1 ) );

                // read the id
                id_t tID;
                mFile.load_data( "SideSetID", tID );

                // get number of facets
                index_t tNumberOfFacets;
                mFile.load_data( "NumberOfFacets", tNumberOfFacets );

                // create the sideset
                SideSet * tSideSet = new SideSet( tID, tNumberOfFacets );

                // read the label
                mFile.load_data( "Label", tSideSet->label() );

                // read the ids
                Vector< id_t > tIDs;
                mFile.load_data( "FacetID", tIDs );


                // read the indices of the master element
                Vector< index_t > tElementIndex;
                mFile.load_data( "FacetOwner", tElementIndex );

                // read the facet indices
                Vector< uint > tFacetIndex;
                mFile.load_data( "FacetOwnerIndex", tFacetIndex );

                // loop over all facets
                for( index_t f=0; f<tNumberOfFacets; ++f )
                {
                    // get pointer to master
                    Element * tMaster = tElements( tElementIndex( f ) );

                    // create the Background Element
                    Element * tElement = tFactory.create_element(
                            element_type_of_facet(
                                    tMaster->type(),
                                    tFacetIndex( f ) ), tIDs( f ) );

                    // create the facet
                    Facet * tFacet = new Facet( tElement );

                    tFacet->set_master( tMaster, tFacetIndex( f ) );

                    // add facet to list
                    tSideSet->insert_facet( tFacet );
                }

                // read neighbors
                mFile.load_data( "FacetNeighbor", tElementIndex );

                // read the facet indices
                mFile.load_data( "FacetNeighborIndex", tFacetIndex );

                // loop over all facets
                for( index_t f=0; f<tNumberOfFacets; ++f )
                {
                    // test if a slave exists
                    if( tElementIndex( f ) < tNumberOfElements )
                    {
                        // link facet with slave
                        tSideSet->facet_by_index( f )->set_slave(
                                tElements( tElementIndex( f ) ),
                                tFacetIndex( f ) );
                    }
                }

                // add sideset to list
                tSideSets( s ) = tSideSet;
            }

            mFile.close_active_group();

            // tell the mesh that the faces have been linked
            mMesh->set_facets_are_linked_flag();
#endif
        }

//------------------------------------------------------------------------------

        void
        HDF5Reader::read_vertices()
        {
            index_t tNumberOfVertices;
            mFile.load_data( "NumberOfVertices", tNumberOfVertices );

            if( tNumberOfVertices > 0 )
            {
                // select the group
                mFile.select_group( "Vertices" );

                // create a container for the element ids
                Vector< id_t > tIDs;
                mFile.load_data( "VertexId", tIDs );

                // create a contaner for element indices
                Vector< index_t > tIndices;
                mFile.load_data( "NodeIndex", tIndices );

                // Tags
                Vector< uint > tGeometryTags;
                mFile.load_data( "GeometryTag",  tGeometryTags );

                Vector< uint > tPhysicalTags;
                mFile.load_data( "PhysicalTag",  tPhysicalTags );

                Vector< uint > tProcOwners;
                mFile.load_data( "ElementOwner", tProcOwners );

                // create an element factory
                ElementFactory tFactory;

                // get vertex container of mesh
                Cell< Element * > & tVertices = mMesh->vertices();

                // allocate container
                tVertices.set_size( tNumberOfVertices, nullptr );

                // get nodes
                Cell< Node * > & tNodes = mMesh->nodes();

                // loop over all indices
                for( index_t k=0; k<tNumberOfVertices; ++k )
                {
                    // create a new vertex
                    Element * tVertex = tFactory.create_element( ElementType::VERTEX, tIDs( k ) );

                    // set the node
                    tVertex->insert_node( tNodes( tIndices( k ) ), 0 );

                    // set the tags
                    tVertex->set_physical_tag( tPhysicalTags( k ) );
                    tVertex->set_geometry_tag( tGeometryTags( k ) );
                    tVertex->set_owner( tProcOwners( k ) );

                    // add vertex to container
                    tVertices( k ) = tVertex;
                }

                mFile.close_active_group();
            }
        }

//------------------------------------------------------------------------------

        void
        HDF5Reader::read_globals()
        {
#ifdef BELFEM_HDF5
            // also write time stuff
            mMesh->time_stamp() = mTimeStamp;
            mMesh->set_time_step( mTimeStep );

            for( uint g=0; g<mNumberOfGlobals; ++g )
            {
                // select the group
                mFile.select_group( sprint( "Global_%u", g+1 ) );

                // read the id
                id_t tID;
                mFile.load_data( "GlobalID", tID );

                // read the value
                real tValue;
                mFile.load_data( "Value", tValue );

                // read the label
                string tLabel;
                mFile.load_data( "Label", tLabel );

                // create the global
                mMesh->create_global_variable( tLabel, tValue, tID );
            }

            mFile.close_active_group();
#endif
        }

//------------------------------------------------------------------------------

        void
        HDF5Reader::read_fields()
        {
#ifdef BELFEM_HDF5
            index_t tNumElements = mMesh->number_of_elements();

            Vector< real > tData;

            for( uint f=0; f<mNumberOfFields; ++f )
            {
                // select the group
                mFile.select_group( sprint( "Field_%u", f+1 ) );

                // read the id
                id_t tID;
                mFile.load_data( "FieldID", tID );

                // read the value
                Vector< real > tData;
                mFile.load_data( "Data", tData );

                // read the label
                string tLabel;
                mFile.load_data( "Label", tLabel );

                if( tData.length() == tNumElements )
                {
                    Vector< real > & tValues
                        = mMesh->create_field( tLabel, EntityType::ELEMENT, tID );

                    tValues = tData;
                }
                else if ( tData.length() == mNumberOfNodes )
                {
                    Vector< real > & tValues
                            = mMesh->create_field( tLabel, EntityType::NODE, tID );
                    tValues = tData;
                }
                else
                {
                    BELFEM_ERROR( false, "don't know how to interpret field %s", tLabel.c_str() );
                }
            }

            mFile.close_active_group();
#endif
        }

//------------------------------------------------------------------------------
    }
}