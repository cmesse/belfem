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

#include "assert.hpp"
#include "stringtools.hpp"
#include "cl_Mesh_VtkWriter.hpp"
#include "cl_Vector.hpp"
#include "cl_Node.hpp"
#include "meshtools.hpp"

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
            // write_time() is not called because its time record causes a bug
            // in recent ParaView versions.
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

            mFile << "# vtk DataFile Version 3.0" << std::endl;
            mFile << "GO BUFFS!" << std::endl;
            mFile << "BINARY" << std::endl;

        }

//------------------------------------------------------------------------------

        void
        VtkWriter::write_time()
        {

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

            Cell< mesh::Element * > & tElements = mMesh->elements();

            for ( mesh::Element * tElement : tElements )
            {
                tElement->flag_nodes();
            }

            Cell< mesh::Node * > & tNodes = mMesh->nodes();

            mNumberOfNodes = 0;

            for ( mesh::Node * tNode : tNodes )
            {
                if ( tNode->is_flagged() )
                {
                    ++mNumberOfNodes;
                }
            }

            mFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

            mFile << "POINTS " << mNumberOfNodes << " float" << std::endl;

            float tFChar;

            mNumberOfNodes = 0;

            for ( mesh::Node * tNode : mMesh->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    mNodeMap[ tNode->id() ] = mNumberOfNodes++;

                    tFChar = vtk::swap_byte_endian(( float ) tNode->x());
                    mFile.write(( char * ) &tFChar, sizeof( float ));
                    tFChar = vtk::swap_byte_endian(( float ) tNode->y());
                    mFile.write(( char * ) &tFChar, sizeof( float ));
                    tFChar = vtk::swap_byte_endian(( float ) tNode->z());
                    mFile.write(( char * ) &tFChar, sizeof( float ));
                }
            }

            mFile << std::endl;

        }

//------------------------------------------------------------------------------

        void
        VtkWriter::write_elements()
        {
            Cell< mesh::Block * > & tBlocks = mMesh->blocks();

            int tCount = 0;

            mNumberOfElements = 0 ;

            for ( mesh::Block * tBlock : tBlocks )
            {
                index_t tNumNodesPerElement = mesh::number_of_nodes( tBlock->element_type() );

                mNumberOfElements += tBlock->number_of_elements();

                tCount += tNumNodesPerElement * tBlock->number_of_elements() + tBlock->number_of_elements();
            }

            mFile << "CELLS " << mNumberOfElements << " " << tCount << std::endl;

            int tIChar;

            for ( mesh::Block * tBlock : tBlocks )
            {
                ElementType tType = tBlock->element_type();

                uint tNumNodesPerElement = mesh::number_of_nodes( tType );

                Cell< mesh::Element * > & tElements = tBlock->elements();

                Vector< id_t > tNodeIDs( tNumNodesPerElement );
                Vector< uint > tNodeIndices( tNumNodesPerElement );

                for ( mesh::Element * tElement : tElements )
                {
                    tIChar = vtk::swap_byte_endian(( int ) tNumNodesPerElement );
                    mFile.write(( char * ) &tIChar, sizeof( int ));

                    vtk::get_node_ids( tElement, tNodeIDs );

                    for ( uint k = 0; k < tNumNodesPerElement; ++k )
                    {
                        tIChar = vtk::swap_byte_endian(( int ) mNodeMap( tNodeIDs( k ) ) );
                        mFile.write( ( char * ) &tIChar, sizeof( int ));
                    }
                }
            }

            mFile << std::endl;

        }

//------------------------------------------------------------------------------

        void
        VtkWriter::write_element_fields()
        {
            Cell< mesh::Block * > & tBlocks = mMesh->blocks();

            mFile << "CELL_TYPES " << mNumberOfElements << std::endl;

            int tIChar ;
            float tFChar ;

            for ( mesh::Block * tBlock : tBlocks )
            {
                tIChar = vtk::swap_byte_endian(( int ) vtk::vtk_type( tBlock->element_type()));

                for ( index_t k = 0; k < tBlock->number_of_elements(); ++k )
                {
                    mFile.write(( char * ) &tIChar, sizeof( int ));
                }
            }

            mFile << std::endl;

            mFile << "CELL_DATA " << mNumberOfElements << std::endl;

            mFile << "SCALARS ELEMENT_ID int" << std::endl;
            mFile << "LOOKUP_TABLE default" << std::endl;

            for ( mesh::Block * tBlock : tBlocks )
            {
                Cell< mesh::Element * > & tElements = tBlock->elements();

                for ( mesh::Element * tElement : tElements )
                {
                    tIChar =  vtk::swap_byte_endian( ( int ) tElement->id() );
                    mFile.write(( char * ) &tIChar, sizeof( int ));
                }
            }

            mFile << std::endl;

            mFile << "SCALARS ELEMENT_OWNER int" << std::endl;
            mFile << "LOOKUP_TABLE default" << std::endl;

            for ( mesh::Block * tBlock : tBlocks )
            {
                Cell< mesh::Element * > & tElements = tBlock->elements();

                for ( mesh::Element * tElement : tElements )
                {
                    tIChar =  vtk::swap_byte_endian( ( int ) tElement->owner() );
                    mFile.write(( char * ) &tIChar, sizeof( int ));
                }
            }

            mFile << std::endl;

            mFile << "SCALARS BLOCK_ID int" << std::endl;
            mFile << "LOOKUP_TABLE default" << std::endl;

            for ( mesh::Block * tBlock : tBlocks )
            {
                tIChar =  vtk::swap_byte_endian( ( int ) tBlock->id() );

                for ( index_t e=0; e<tBlock->number_of_elements(); ++e )
                {
                    mFile.write(( char * ) &tIChar, sizeof( int ));
                }
            }

            mFile << std::endl;

            uint tNumFields = mMesh->number_of_fields() ;

            for( uint k=0; k<tNumFields; ++k )
            {
                mesh::Field * tField = mMesh->field( k ) ;

                if( tField->entity_type() == EntityType::ELEMENT )
                {
                    string tLabel = search_and_replace( mMesh->field( k )->label()," ", "_" );

                    mFile << "SCALARS " << tLabel << " float" << std::endl;
                    mFile << "LOOKUP_TABLE default" << std::endl;

                    Vector< real > & tData = tField->data() ;

                    for( mesh::Block * tBlock : mMesh->blocks() )
                    {
                        for( mesh::Element * tElement : tBlock->elements() )
                        {
                            tFChar = vtk::swap_byte_endian( static_cast< float > ( tData( tElement->index() ) ) );
                            mFile.write( ( char * ) &tFChar, sizeof( float ));
                        }
                    }

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

            int tIChar;
            float tFChar;

            for ( mesh::Node * tNode : mMesh->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    tIChar = vtk::swap_byte_endian(( int ) tNode->id() );
                    mFile.write( ( char * ) &tIChar, sizeof( int ));
                }
            }

            mFile << std::endl;

            mFile << "SCALARS NODE_OWNER int" << std::endl;
            mFile << "LOOKUP_TABLE default" << std::endl;

            for ( mesh::Node * tNode : mMesh->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    tIChar = vtk::swap_byte_endian(( int ) tNode->owner() );
                    mFile.write( ( char * ) &tIChar, sizeof( int ));
                }
            }

            mFile << std::endl;

            uint tNumFields = mMesh->number_of_fields() ;

            for( uint k=0; k<tNumFields; ++k )
            {
                mesh::Field * tField = mMesh->field( k ) ;

                if( tField->entity_type() == EntityType::NODE )
                {
                    string tLabel = search_and_replace( mMesh->field( k )->label()," ", "_" );

                    mFile << "SCALARS " << tLabel << " float" << std::endl;
                    mFile << "LOOKUP_TABLE default" << std::endl;

                    Vector< real > & tData = tField->data() ;

                    for( mesh::Node * tNode : mMesh->nodes() )
                    {
                        if( tNode->is_flagged() )
                        {
                            tFChar = vtk::swap_byte_endian( static_cast< float > ( tData( tNode->index() ) ) );
                            mFile.write( ( char * ) &tFChar, sizeof( float ));
                        }
                    }

                    mFile << std::endl;
                }
            }

        }

//------------------------------------------------------------------------------

    }
}