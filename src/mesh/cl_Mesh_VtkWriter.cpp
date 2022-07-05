//
// Created by Christian Messe on 22.06.20.
//

#include "cl_Mesh_VtkWriter.hpp"
#include "cl_Vector.hpp"
#include "cl_Node.hpp"
#include "meshtools.hpp"
#include "assert.hpp"
#include "vtktools.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        VtkWriter::VtkWriter( const string & aFilePath, Mesh * aMesh ) :
                mFilePath( aFilePath ),
                mMesh( aMesh ),
                mFile( aFilePath, std::ios::binary )
        {
            this->write_header();
            // causes a bug in new ParaView
            //this->write_time();
            this->write_nodes();
            this->write_elements();
            this->write_element_fields() ;
            this->write_node_fields() ;
            mFile.close();
        }

//------------------------------------------------------------------------------

        void
        VtkWriter::write_header()
        {

            // dump header information
            mFile << "# vtk DataFile Version 3.0" << std::endl;
            mFile << "GO BUFFS!" << std::endl;
            mFile << "BINARY" << std::endl;

        }

//------------------------------------------------------------------------------

        void
        VtkWriter::write_time()
        {

            // dump time information
            mFile << "DATASET POLYDATA" << std::endl;
            mFile << "FIELD FieldData 2" << std::endl;
            mFile << "TIME 1 1 double" << std::endl;

            double tDChar = vtk::swap_byte_endian( ( double ) mMesh->time_stamp() );
            mFile.write(( char * ) &tDChar, sizeof( double ));
            mFile << std::endl;

            mFile << "CYCLE 1 1 int" << std::endl;
            int tIChar = vtk::swap_byte_endian( ( int ) mMesh->time_step() );
            mFile.write(( char * ) &tIChar, sizeof( int ));
            mFile << std::endl;
        }

//------------------------------------------------------------------------------

        void
        VtkWriter::write_nodes()
        {
            // only use nodes that are connected to elements
            mMesh->unflag_all_nodes();

            // get elements from mesh
            Cell< mesh::Element * > & tElements = mMesh->elements();

            // loop over elements and flag used nodes
            for ( mesh::Element * tElement : tElements )
            {
                tElement->flag_nodes();
            }

            // grab node container
            Cell< mesh::Node * > & tNodes = mMesh->nodes();

            // node counter
            mNumberOfNodes = 0;

            // backup original node indices
            for ( mesh::Node * tNode : tNodes )
            {
                // check if node is used
                if ( tNode->is_flagged() )
                {
                    // increment counter
                    ++mNumberOfNodes;
                }
            }

            // specify grid type
            mFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

            // write number of nodes to file
            mFile << "POINTS " << mNumberOfNodes << " float" << std::endl;

            // float container
            float tFChar;

            // reset counter
            mNumberOfNodes = 0;

            // loop over all nodes
            for ( mesh::Node * tNode : mMesh->nodes() )
            {
                // check if node is used
                if ( tNode->is_flagged() )
                {
                    // remember index of this node
                    mNodeMap[ tNode->id() ] = mNumberOfNodes++;

                    // write node coordinates
                    tFChar = vtk::swap_byte_endian(( float ) tNode->x());
                    mFile.write(( char * ) &tFChar, sizeof( float ));
                    tFChar = vtk::swap_byte_endian(( float ) tNode->y());
                    mFile.write(( char * ) &tFChar, sizeof( float ));
                    tFChar = vtk::swap_byte_endian(( float ) tNode->z());
                    mFile.write(( char * ) &tFChar, sizeof( float ));
                }
            }

            // create new line
            mFile << std::endl;

        }

//------------------------------------------------------------------------------

        void
        VtkWriter::write_elements()
        {
            // get block container of mesh
            Cell< mesh::Block * > & tBlocks = mMesh->blocks();

            // count memory index that is required by vtk
            int tCount = 0;

            // number of elements in total
            mNumberOfElements = 0 ;

            // loop over all blocks
            for ( mesh::Block * tBlock : tBlocks )
            {
                // get number of nodes from this element
                index_t tNumNodesPerElement = mesh::number_of_nodes( tBlock->element_type() );

                // add elements of this block to counter
                mNumberOfElements += tBlock->number_of_elements();

                // add to memory counter
                tCount += tNumNodesPerElement * tBlock->number_of_elements() + tBlock->number_of_elements();
            }

            // get cut container of mesh
            Cell< mesh::SideSet * > & tCuts = mMesh->cuts();

            for ( mesh::SideSet * tCut : tCuts )
            {
                // get number of nodes from this element
                index_t tNumNodesPerElement = 2 ;

                // add elements of this block to counter
                mNumberOfElements += tCut->number_of_facets();

                // add to memory counter
                tCount += tNumNodesPerElement * tCut->number_of_facets() + tCut->number_of_facets();
            }

            // write header for cells
            mFile << "CELLS " << mNumberOfElements << " " << tCount << std::endl;

            // integer container
            int tIChar;

            // loop over all blocks
            for ( mesh::Block * tBlock : tBlocks )
            {
                // get element type from block
                ElementType tType = tBlock->element_type();

                // get number of nodes from this element
                uint tNumNodesPerElement = mesh::number_of_nodes( tType );

                // get element container from this block
                Cell< mesh::Element * > & tElements = tBlock->elements();

                // container with node IDS
                Vector< id_t > tNodeIDs( tNumNodesPerElement );
                Vector< uint > tNodeIndices( tNumNodesPerElement );

                // loop over all elements from this block
                for ( mesh::Element * tElement : tElements )
                {
                    // write number of nodes
                    tIChar = vtk::swap_byte_endian(( int ) tNumNodesPerElement );
                    mFile.write(( char * ) &tIChar, sizeof( int ));

                    // populate node IDs
                    vtk::get_node_ids( tElement, tNodeIDs );

                    // write node ID
                    for ( uint k = 0; k < tNumNodesPerElement; ++k )
                    {
                        tIChar = vtk::swap_byte_endian(( int ) mNodeMap( tNodeIDs( k ) ) );
                        mFile.write( ( char * ) &tIChar, sizeof( int ));
                    }
                }
            }

            // loop over all cuts
            for ( mesh::SideSet * tCut : tCuts )
            {
                // get element type from block
                ElementType tType = ElementType::LINE2 ;

                // get number of nodes from this element
                uint tNumNodesPerElement = mesh::number_of_nodes( tType );

                // get element container from this block
                Cell< mesh::Facet * > & tFacets = tCut->facets();

                // container with node IDS
                Vector< id_t > tNodeIDs( tNumNodesPerElement );
                Vector< uint > tNodeIndices( tNumNodesPerElement );

                // loop over all elements from this block
                for ( mesh::Facet * tFacet : tFacets )
                {
                    // get element
                    mesh::Element * tElement = tFacet->element() ;

                    // write number of nodes
                    tIChar = vtk::swap_byte_endian(( int ) tNumNodesPerElement );
                    mFile.write(( char * ) &tIChar, sizeof( int ));

                    // populate node IDs
                    vtk::get_node_ids( tElement, tNodeIDs );

                    // write node ID
                    for ( uint k = 0; k < tNumNodesPerElement; ++k )
                    {
                        tIChar = vtk::swap_byte_endian(( int ) mNodeMap( tNodeIDs( k ) ) );
                        mFile.write( ( char * ) &tIChar, sizeof( int ));
                    }
                }
            }

            // create new line
            mFile << std::endl;

        }

//------------------------------------------------------------------------------

        void
        VtkWriter::write_element_fields()
        {
            // get block container of mesh
            Cell< mesh::Block * > & tBlocks = mMesh->blocks();
            Cell< mesh::SideSet * > & tCuts = mMesh->cuts();

            mFile << "CELL_TYPES " << mNumberOfElements << std::endl;

            // integer container
            int tIChar ;
            float tFChar ;

            // write cell types
            for ( mesh::Block * tBlock : tBlocks )
            {
                // get type
                tIChar = vtk::swap_byte_endian(( int ) vtk::vtk_type( tBlock->element_type()));

                // populate data
                for ( index_t k = 0; k < tBlock->number_of_elements(); ++k )
                {
                    mFile.write(( char * ) &tIChar, sizeof( int ));
                }
            }
            for ( mesh::SideSet * tCut : tCuts )
            {
                // get type
                tIChar = vtk::swap_byte_endian(( int ) vtk::vtk_type( ElementType::LINE2 ) );

                // populate data
                for ( index_t k = 0; k < tCut->number_of_facets(); ++k )
                {
                    mFile.write(( char * ) &tIChar, sizeof( int ));
                }
            }

            // create new line
            mFile << std::endl;

            // write element data
            mFile << "CELL_DATA " << mNumberOfElements << std::endl;

            // write element ID
            mFile << "SCALARS ELEMENT_ID int" << std::endl;
            mFile << "LOOKUP_TABLE default" << std::endl;

            for ( mesh::Block * tBlock : tBlocks )
            {
                // get element container from this block
                Cell< mesh::Element * > & tElements = tBlock->elements();

                // populate data
                for ( mesh::Element * tElement : tElements )
                {
                    tIChar =  vtk::swap_byte_endian( ( int ) tElement->id() );
                    mFile.write(( char * ) &tIChar, sizeof( int ));
                }
            }
            for ( mesh::SideSet * tCut : tCuts )
            {
                // get element container from this block
                Cell< mesh::Facet * > & tFacets = tCut->facets();

                // populate data
                for ( mesh::Facet * tFacet : tFacets )
                {
                    tIChar =  vtk::swap_byte_endian( ( int ) tFacet->id() );
                    mFile.write(( char * ) &tIChar, sizeof( int ));
                }
            }

            // create new line
            mFile << std::endl;

            // write element ID
            mFile << "SCALARS ELEMENT_OWNER int" << std::endl;
            mFile << "LOOKUP_TABLE default" << std::endl;

            for ( mesh::Block * tBlock : tBlocks )
            {
                // get element container from this block
                Cell< mesh::Element * > & tElements = tBlock->elements();

                // populate data
                for ( mesh::Element * tElement : tElements )
                {
                    tIChar =  vtk::swap_byte_endian( ( int ) tElement->owner() );
                    mFile.write(( char * ) &tIChar, sizeof( int ));
                }
            }
            for ( mesh::SideSet * tCut : tCuts )
            {
                // get element container from this block
                Cell< mesh::Facet * > & tFacets = tCut->facets();

                // populate data
                for ( mesh::Facet * tFacet : tFacets )
                {
                    tIChar =  vtk::swap_byte_endian( ( int ) tFacet->element()->owner() );
                    mFile.write(( char * ) &tIChar, sizeof( int ));
                }
            }

            // create new line
            mFile << std::endl;


            // write block ID
            mFile << "SCALARS BLOCK_ID int" << std::endl;
            mFile << "LOOKUP_TABLE default" << std::endl;

            for ( mesh::Block * tBlock : tBlocks )
            {
                // get ID of Block
                tIChar =  vtk::swap_byte_endian( ( int ) tBlock->id() );

                // populate data
                for ( index_t e=0; e<tBlock->number_of_elements(); ++e )
                {
                    mFile.write(( char * ) &tIChar, sizeof( int ));
                }
            }
            for ( mesh::SideSet * tCut : tCuts )
            {
                // get ID of Block
                tIChar =  vtk::swap_byte_endian( ( int ) tCut->id() );

                // populate data
                for ( index_t f=0; f<tCut->number_of_facets(); ++f )
                {
                    mFile.write(( char * ) &tIChar, sizeof( int ));
                }
            }

            // create new line
            mFile << std::endl;

            // get number of fields
            uint tNumFields = mMesh->number_of_fields() ;

            for( uint k=0; k<tNumFields; ++k )
            {
                // get field
                mesh::Field * tField = mMesh->field( k ) ;

                if( tField->entity_type() == EntityType::ELEMENT )
                {
                    string tLabel = search_and_replace( mMesh->field( k )->label()," ", "_" );

                    mFile << "SCALARS " << tLabel << " float" << std::endl;
                    mFile << "LOOKUP_TABLE default" << std::endl;

                    // get data container
                    Vector< real > & tData = tField->data() ;

                    for( mesh::Block * tBlock : mMesh->blocks() )
                    {
                        for( mesh::Element * tElement : tBlock->elements() )
                        {
                            // convert value of node into float
                            tFChar = vtk::swap_byte_endian( tData( tElement->index() ) );
                            mFile.write( ( char * ) &tFChar, sizeof( float ));
                        }
                    }

                    for( mesh::SideSet * tCut : mMesh->cuts() )
                    {
                        // convert value of node into float
                        tFChar = vtk::swap_byte_endian( 0.0 );

                        for( index_t f=0 ; f<tCut->number_of_facets(); ++f )
                        {
                            mFile.write( ( char * ) &tFChar, sizeof( float ));
                        }
                    }

                    // create new line
                    mFile << std::endl;
                }
            }

        }

//------------------------------------------------------------------------------

        void
        VtkWriter::write_node_fields()
        {
            mFile << "POINT_DATA " << mNumberOfNodes << std::endl;

            mFile << "SCALARS NODE_ID int" << std::endl;
            mFile << "LOOKUP_TABLE default" << std::endl;

            // integer container
            int tIChar;
            float tFChar;

            // loop over all nodes
            for ( mesh::Node * tNode : mMesh->nodes() )
            {
                // check if node is used
                if ( tNode->is_flagged() )
                {
                    // remember index of this node
                    mNodeMap[ tNode->id() ] = mNumberOfNodes++;

                    // write node coordinates
                    tIChar = vtk::swap_byte_endian(( int ) tNode->id() );
                    mFile.write( ( char * ) &tIChar, sizeof( int ));
                }
            }

            // create new line
            mFile << std::endl;

            mFile << "SCALARS NODE_OWNER int" << std::endl;
            mFile << "LOOKUP_TABLE default" << std::endl;

            // loop over all nodes
            for ( mesh::Node * tNode : mMesh->nodes() )
            {
                // check if node is used
                if ( tNode->is_flagged() )
                {
                    // remember index of this node
                    mNodeMap[ tNode->id() ] = mNumberOfNodes++;

                    // write node coordinates
                    tIChar = vtk::swap_byte_endian(( int ) tNode->owner() );
                    mFile.write( ( char * ) &tIChar, sizeof( int ));
                }
            }

            // create new line
            mFile << std::endl;

            // get number of fields
            uint tNumFields = mMesh->number_of_fields() ;

            for( uint k=0; k<tNumFields; ++k )
            {
                // get field
                mesh::Field * tField = mMesh->field( k ) ;

                if( tField->entity_type() == EntityType::NODE )
                {
                    string tLabel = search_and_replace( mMesh->field( k )->label()," ", "_" );

                    mFile << "SCALARS " << tLabel << " float" << std::endl;
                    mFile << "LOOKUP_TABLE default" << std::endl;

                    // get data container
                    Vector< real > & tData = tField->data() ;

                    for( mesh::Node * tNode : mMesh->nodes() )
                    {
                        // check if node is used
                        if( tNode->is_flagged() )
                        {
                            // convert value of node into float
                            tFChar = vtk::swap_byte_endian( tData( tNode->index() ) );
                            mFile.write( ( char * ) &tFChar, sizeof( float ));
                        }
                    }

                    // create new line
                    mFile << std::endl;
                }
            }

        }

//------------------------------------------------------------------------------

    }
}