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
#include "filetools.hpp"
#include "cl_Matrix.hpp"

#include "commtools.hpp"

#include <cstdlib>
#include "cl_Mesh_GmshReader.hpp"
#include "cl_Element_Factory.hpp"
#include "fn_min.hpp"
#include "fn_max.hpp"
#include "fn_unique.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        GmshReader::GmshReader(
            const string & aPath,
            Mesh * aMesh,
            const bool aComputeConnectivities,
            const real aMeshScale ) :
            Ascii( aPath, FileMode::OPEN_RDONLY ),
            mMeshScale( aMeshScale )
        {
            mFilename = basename( aPath );

            this->tidy_up_buffer();
            this->read_version();

            if( mVersion == 2.2 )
            {
                this->read_tags_v22();
                this->check_tag_existence();
                this->read_nodes_v22();
                this->create_node_map();
                this->read_elements_v22();
            }
            else if ( mVersion == 4.1 )
            {
                this->read_mesh_v41();
                this->check_tag_existence();
            }
            else
            {
                BELFEM_ERROR( false,
                        "The version number of mesh %s is %f, however, only 2.2 and 4.1 are supported",
                             mFilename.c_str(),
                             ( float ) mVersion );
            }

            this->read_mesh_dimension();
            this->create_element_containers_per_dimension();

            this->create_group_ids();

            if( aMesh == NULL )
            {
                mMesh = new Mesh( mNumberOfDimensions, comm_rank() );
                mOwnMesh = true;
            }
            else
            {
                mMesh = aMesh;
            }

            // for higher order elements, the order is different in exodus
            this->convert_orders_for_quadratic_volume_elements_to_exo();
            this->create_blocks();
            this->create_sidesets();
            this->create_mesh();

            // only keep vertices whose node is connected to an element.
            // unconnected construction points would dangle: they carry no dof,
            // trigger a singular matrix if adopted as orphaned nodes, and break
            // the mesh save/load round-trip
            this->create_vertices();

            // free what the mesh did not adopt. Must stay after
            // create_vertices() ( which flags the kept vertices ) and before
            // finalize() ( which unflags the block elements )
            this->delete_unused_nodes_and_elements();

            if ( ! aComputeConnectivities )
            {
                mMesh->reset_connectivity( Connectivity::Compute );
            }
            mMesh->finalize();
        }

//------------------------------------------------------------------------------

        GmshReader::~GmshReader()
        {
            if( mOwnMesh )
            {
                delete mMesh;
            }
        }

//------------------------------------------------------------------------------

        void
        GmshReader::tidy_up_buffer()
        {
            for( auto tLine: mBuffer )
            {
                tLine = clean_string( tLine );
            }
        }

//------------------------------------------------------------------------------

        void
        GmshReader::check_tag_existence()
        {
            BELFEM_ERROR( mNodesTag > 0,
                    "$Nodes tag not found in file %s",
                    mFilename.c_str() );

            BELFEM_ERROR( mElementsTag > 0,
                         "$Elements tag not found in file %s",
                         mFilename.c_str() );
        }

//------------------------------------------------------------------------------

        void
        GmshReader::read_version()
        {
            size_t tCount = 0;

            for( auto tLine: mBuffer )
            {
                if ( tLine == "$MeshFormat")
                {
                    break;
                }
                ++tCount;
            }

            BELFEM_ERROR( tCount < mBuffer.size(),
                    "The file %s does not seem to be a GMSH mesh file",
                         mFilename.c_str() );

            mVersion = std::stod( first_word( mBuffer( ++tCount ) ) );
        }

//------------------------------------------------------------------------------

        void
        GmshReader::read_tags_v22()
        {
            size_t tCount = 0;

            size_t tN = mBuffer.size();

            while( tCount < tN )
            {
                if ( mBuffer( tCount ) == "$Nodes")
                {
                    this->read_nodes_tag_v22( tCount );
                }
                else if( mBuffer( tCount ) == "$Elements" )
                {
                    this->read_elements_tag_v22( tCount );
                }
                else if( mBuffer( tCount ) == "$PhysicalNames" )
                {
                    this->read_physical_tag( tCount );
                }

                ++tCount;
            }
        }

//------------------------------------------------------------------------------

        void
        GmshReader::read_nodes_tag_v22( size_t & aCount )
        {
            mNodesTag = aCount;

            mNumberOfNodes = std::stoi( mBuffer( ++aCount ) );

            aCount += mNumberOfNodes + 1;
        }

//------------------------------------------------------------------------------

        void
        GmshReader::read_elements_tag_v22( size_t & aCount )
        {
            mElementsTag = aCount;

            mNumberOfElements = std::stoi( mBuffer( ++aCount ) );

            aCount += mNumberOfElements + 1;
        }

//------------------------------------------------------------------------------

        void
        GmshReader::read_physical_tag( size_t & aCount )
        {
            mPhysicalTag = aCount;

            mNumberOfPhysicalGroups = std::stoi( mBuffer( ++aCount ) );

            aCount += mNumberOfPhysicalGroups + 1;
        }

//------------------------------------------------------------------------------

        void
        GmshReader::read_nodes_v22()
        {
            mNodes.set_size( mNumberOfNodes, nullptr );

            size_t tCount = 0;

            uint tStart = mNodesTag + 2;

            uint tEnd = tStart + mNumberOfNodes;

            for ( size_t k = tStart; k < tEnd; ++k )
            {
                Cell< string > tWords = string_to_words( mBuffer( k ));

                id_t tID = std::stoi( tWords( 0 ));

                real tX = std::stod( tWords( 1 )) * mMeshScale ;

                real tY = std::stod( tWords( 2 )) * mMeshScale ;

                real tZ = std::stod( tWords( 3 )) * mMeshScale ;

                mNodes( tCount++ ) = new Node( tID, tX, tY, tZ );
            }

        }

//------------------------------------------------------------------------------

        void
        GmshReader::read_elements_v22()
        {
            mElements.set_size( mNumberOfElements, nullptr );

            size_t tCount = 0;

            uint tStart = mElementsTag + 2;

            uint tEnd = tStart + mNumberOfElements;

            ElementFactory tFactory;

            for( uint k=tStart; k<tEnd; ++k )
            {
                Cell< string > tWords = string_to_words( mBuffer( k ) );

                id_t tID = std::stoi( tWords( 0 ) );

                ElementType tType = element_type_from_gmsh(  std::stoi( tWords( 1 ) ) );

                index_t tNumNodes = number_of_nodes( tType );

                uint tNumTags = std::stoi( tWords( 2 ) );

                BELFEM_ERROR( tWords.size()-3-tNumTags == tNumNodes,
                             "Number of nodes for element %u does not match. Is: %u, Expect %u",
                             ( unsigned int ) tID,
                             ( unsigned int ) tWords.size()-3-tNumTags,
                             ( unsigned int ) tNumNodes );

                if ( tType == ElementType::VERTEX && tNumTags == 2 )
                {
                    tID = std::stoi( tWords( 4 ) );
                }

                Element * tElement = tFactory.create_element( tType, tID );

                tElement->set_physical_tag( std::stoi( tWords( 3 ) ) );

                tElement->set_geometry_tag( std::stoi( tWords( 4 ) ) );

                index_t tOff = tNumTags + 3;

                for ( index_t i=0; i<tNumNodes; ++i )
                {
                    tElement->insert_node( mNodeMap( std::stoi( tWords( i + tOff ) )  ) , i );
                }

                mElements( tCount++ ) = tElement;
            }
        }

 //------------------------------------------------------------------------------

        void
        GmshReader::read_mesh_v41()
        {
            size_t tCount = 0;

            size_t tN = mBuffer.size();

            while( tCount < tN )
            {
                if ( mBuffer( tCount ) == "$PhysicalNames" )
                {
                    this->read_physical_tag( tCount );
                }
                else if ( mBuffer( tCount ) == "$Entities" )
                {
                    this->read_entity_tag_v41( tCount );
                }
                if ( mBuffer( tCount ) == "$Nodes" )
                {

                    this->read_nodes_v41( tCount );
                    this->create_node_map();
                }
                else if ( mBuffer( tCount ) == "$Elements" )
                {
                    this->read_elements_v41( tCount );
                }

                ++tCount;
            }
        }

//------------------------------------------------------------------------------

        void
        GmshReader::read_entity_tag_v41( size_t & aCount )
        {
            // header: numPoints numCurves numSurfaces numVolumes
            Cell< string > tWords = string_to_words( mBuffer( ++aCount ) );

            BELFEM_ERROR( tWords.size() == 4,
                "Malformed $Entities header in file %s", mFilename.c_str() );

            for( uint tDim=0; tDim<4; ++tDim )
            {
                size_t tNumEntities = std::stoi( tWords( tDim ) );

                // a point line reads: tag x y z numPhysicalTags [ tags ]
                // every other:        tag minX minY minZ maxX maxY maxZ
                //                     numPhysicalTags [ tags ] numBounding [ .. ]
                uint tPos = tDim == 0 ? 4 : 7 ;

                for( size_t k=0; k<tNumEntities; ++k )
                {
                    Cell< string > tLine = string_to_words( mBuffer( ++aCount ) );

                    BELFEM_ERROR( tLine.size() > tPos,
                        "Malformed $Entities line in file %s", mFilename.c_str() );

                    // an entity may carry several physical tags; an element
                    // holds one, so the first is kept, as the v2.2 path does.
                    // the tag is stored unsigned; a negative value would wrap
                    if( std::stoi( tLine( tPos ) ) > 0 )
                    {
                        BELFEM_ERROR( tLine.size() > tPos + 1,
                            "Malformed $Entities line in file %s", mFilename.c_str() );

                        mPhysicalTagOfEntity[ tDim ][ std::stoi( tLine( 0 ) ) ]
                            = std::abs( std::stoi( tLine( tPos + 1 ) ) );
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        GmshReader::read_nodes_v41( size_t & aCount )
        {
            mNodesTag = aCount++;

            Cell< string > tWords = string_to_words( mBuffer( aCount++ ) );

            mNumberOfNodes = std::stoi( tWords( 1 ) );

            uint tNodeCount = 0;

            Vector< id_t > mNodeIDs( mNumberOfNodes );
            Matrix< real > mNodeCoords( 3, mNumberOfNodes );

            while( tNodeCount < mNumberOfNodes )
            {
                tWords = string_to_words( mBuffer( aCount++ ) );

                uint tParametric = std::stoi( tWords( 2 ) );

                BELFEM_ERROR( tParametric == 0,
                    "The mesh %s contains parameteric nodes, which are not supported by BELFEM",
                    mFilename.c_str() );

                uint tN = std::stoi( tWords( 3 ) );

                for( uint k=0; k<tN; ++k )
                {
                    mNodeIDs( tNodeCount + k ) = std::stoi( mBuffer( aCount++ ) );
                }

                for( uint k=0; k<tN; ++k )
                {
                    tWords = string_to_words( mBuffer( aCount++ ) );
                    mNodeCoords( 0, tNodeCount + k ) = std::stod( tWords( 0 ) ) * mMeshScale ;
                    mNodeCoords( 1, tNodeCount + k ) = std::stod( tWords( 1 ) ) * mMeshScale ;
                    mNodeCoords( 2, tNodeCount + k ) = std::stod( tWords( 2 ) ) * mMeshScale ;
                }

                tNodeCount += tN;
            }

            mNodes.set_size( mNumberOfNodes, nullptr );

            for ( uint k=0; k<mNumberOfNodes; ++k )
            {
                mNodes( k ) = new Node(
                        mNodeIDs( k ),
                        mNodeCoords( 0, k ),
                        mNodeCoords( 1, k ),
                        mNodeCoords( 2, k ) );
            }

        }

//------------------------------------------------------------------------------

        void
        GmshReader::convert_orders_for_quadratic_volume_elements_to_exo()
        {
            Cell< Node * > tNodes;
            tNodes.set_size( 27, nullptr );

            // gmsh to exodus conversion table
            uint tHexIndex[ 27 ] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 9, 10,
                                     12, 14, 15, 16, 18, 19, 17, 26, 20, 25,
                                     22, 23, 21, 24 };

            uint tPentaIndex[ 18 ] ={ 0, 1, 2, 3, 4, 5, 6, 9, 7, 8, 10, 11,
                                      12, 14, 13, 15, 17, 16 };

            uint tPyraIndex[ 14 ] ={ 0, 1, 2, 3, 4, 5, 8, 10, 6, 7, 9, 11, 12, 13 };

            for ( Element * tElement : mElements )
            {
                switch( tElement->type() )
                {
                    case( ElementType::TET10 ):
                    {
                        // only 8 and 9 are swapped
                        Node * tNode = tElement->node( 8 );
                        tElement->insert_node( tElement->node( 9 ), 8 );
                        tElement->insert_node( tNode, 9 );
                        break;
                    }
                    case( ElementType::PENTA15 ):
                    case( ElementType::PENTA18 ):
                    {
                        uint tNumNodes = tElement->number_of_nodes();

                        for ( uint k=0; k<tNumNodes; ++k )
                        {
                            tNodes( k ) = tElement->node( k );
                        }

                        for ( uint k=0; k<tNumNodes; ++k )
                        {
                            tElement->insert_node( tNodes( tPentaIndex[ k ] ), k );
                        }

                        break;
                    }
                    case( ElementType::PYRA13 ):
                    case( ElementType::PYRA14 ):
                    {
                        uint tNumNodes = tElement->number_of_nodes();

                        for ( uint k=0; k<tNumNodes; ++k )
                        {
                            tNodes( k ) = tElement->node( k );
                        }

                        for ( uint k=0; k<tNumNodes; ++k )
                        {
                            tElement->insert_node( tNodes( tPyraIndex[ k ] ), k );
                        }

                        break;
                    }
                    case( ElementType::HEX20 ):
                    case( ElementType::HEX27 ):
                    {
                        uint tNumNodes = tElement->number_of_nodes();

                        for ( uint k=0; k<tNumNodes; ++k )
                        {
                            tNodes( k ) = tElement->node( k );
                        }

                        for ( uint k=0; k<tNumNodes; ++k )
                        {
                            tElement->insert_node( tNodes( tHexIndex[ k ] ), k );
                        }

                        break;
                    }
                    default:
                    {
                        break;
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        GmshReader::read_elements_v41( size_t & aCount )
        {
            // this routine assumes that a node map was already created.
            // therefore, check for that
            BELFEM_ERROR( mNodesTag > 0,
                "$Elements tag must be after $Nodes tag in file %s",
                mFilename.c_str() );

            mElementsTag = aCount++;

            Cell<string> tWords = string_to_words( mBuffer( aCount++ ));

            mNumberOfElements = std::stoi( tWords( 1 ));

            mElements.set_size( mNumberOfElements, nullptr );

            uint tElementCount = 0;

            ElementFactory tFactory;

            while ( tElementCount < mNumberOfElements )
            {
                tWords = string_to_words( mBuffer( aCount++ ));

                uint tDim = std::stoi( tWords( 0 ));
                uint tGeometryTag = std::stoi( tWords( 1 ));

                BELFEM_ERROR( tDim < 4,
                    "Invalid entity dimension %u in file %s", tDim, mFilename.c_str() );

                ElementType tType = element_type_from_gmsh( std::stoi(  tWords( 2 ) ));

                uint tNumElements = std::stoi( tWords( 3 ) );

                size_t tNumNodes = number_of_nodes( tType );

                for ( uint k = 0; k < tNumElements; ++k )
                {
                    tWords = string_to_words( mBuffer( aCount++ ) );

                    id_t tID = std::stoi( tWords( 0 ) );

                    BELFEM_ERROR( tWords.size() == tNumNodes + 1,
                                 "Number of nodes for element %ui does not match. Is: %u, Expect %u",
                                 ( unsigned int ) tID,
                                 ( unsigned int ) tWords.size()-1,
                                 ( unsigned int ) tNumNodes );

                    if( tType == ElementType::VERTEX )
                    {
                        BELFEM_ERROR( tNumElements == 1, "Don't understand input file." );

                        // we want to use the GeometryTag as Entity ID
                        tID = tGeometryTag ;
                    }

                    Element * tElement = tFactory.create_element( tType, tID );

                    for( uint i=0; i < tNumNodes; ++i )
                    {
                        tElement->insert_node( mNodeMap( std::stoi( tWords( i + 1 ) ) ), i );
                    }

                    tElement->set_geometry_tag( tGeometryTag );

                    if( mPhysicalTagOfEntity[ tDim ].key_exists( tGeometryTag ) )
                    {
                        tElement->set_physical_tag( mPhysicalTagOfEntity[ tDim ]( tGeometryTag ) );
                    }

                    mElements( tElementCount++ ) = tElement;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        GmshReader::create_node_map()
        {
            mNodeMap.clear();

            for ( Node* tNode: mNodes )
            {
                mNodeMap[ tNode->id() ] = tNode;
            }
        }

//------------------------------------------------------------------------------

        void
        GmshReader::read_mesh_dimension()
        {
            Vector< real > tZ( mNumberOfNodes );

            for ( index_t k=0; k<mNumberOfNodes; ++k )
            {
                tZ( k ) = mNodes( k )->z();

            }

            if( min( tZ ) == max( tZ ) )
            {
                mNumberOfDimensions = 2;
            }
            else
            {
                mNumberOfDimensions = 3;
            }
        }

//------------------------------------------------------------------------------

        void
        GmshReader::create_element_containers_per_dimension()
        {
            mNumberOfElementsPerDimension.set_size( 4, 0 );

            for( Element * tElement: mElements )
            {
                ++mNumberOfElementsPerDimension( dimension( tElement->type() ) );
            }

            mVertices.set_size( mNumberOfElementsPerDimension( 0 ), nullptr );
            mEdges.set_size( mNumberOfElementsPerDimension( 1 ), nullptr );
            mFaces.set_size( mNumberOfElementsPerDimension( 2 ), nullptr );
            mVolumes.set_size( mNumberOfElementsPerDimension( 3 ), nullptr );

            mNumberOfElementsPerDimension.fill( 0 );

            for( Element * tElement: mElements )
            {
                switch(  dimension( tElement->type() ) )
                {
                    case( 0 ):
                    {
                        mVertices( mNumberOfElementsPerDimension( 0 )++ ) = tElement;
                        break;
                    }
                    case( 1 ):
                    {
                        mEdges( mNumberOfElementsPerDimension( 1 )++ ) = tElement;

                        break;
                    }
                    case( 2 ):
                    {
                        mFaces( mNumberOfElementsPerDimension( 2 )++ ) = tElement;
                        break;
                    }
                    case( 3 ):
                    {
                        mVolumes( mNumberOfElementsPerDimension( 3 )++ ) = tElement;
                        break;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Invalid element dimension." );
                    }
                }
            }

        }

//------------------------------------------------------------------------------

        void
        GmshReader::create_group_ids()
        {
            // copy numbers for better readability
            index_t tN0 = mNumberOfElementsPerDimension( 0 );
            index_t tN1 = mNumberOfElementsPerDimension( 1 );
            index_t tN2 = mNumberOfElementsPerDimension( 2 );
            index_t tN3 = mNumberOfElementsPerDimension( 3 );

            if( tN0 > 0 )
            {
                mNodeGroupIDs.set_size( tN0, mVertices( 0 )->geometry_tag() );

                for( uint k=1; k<tN0; ++k )
                {
                    mNodeGroupIDs( k ) = mVertices( k )->geometry_tag();
                }
                unique( mNodeGroupIDs );
            }

            if( tN1 > 0 )
            {
                mEdgeGroupIDs.set_size( tN1, mEdges( 0 )->geometry_tag() );

                for( uint k=1; k<tN1; ++k )
                {
                    mEdgeGroupIDs( k ) = mEdges( k )->geometry_tag();
                }
                unique( mEdgeGroupIDs );
            }

            if( tN2 > 0 )
            {
                mFaceGroupIDs.set_size( tN2, mFaces( 0 )->geometry_tag() );

                for( uint k=1; k<tN2; ++k )
                {
                    mFaceGroupIDs( k ) = mFaces( k )->geometry_tag();
                }
                unique( mFaceGroupIDs );
            }

            if( tN3 > 0 )
            {
                mVolumeGroupIDs.set_size( tN3, mVolumes( 0 )->geometry_tag() );

                for( uint k=1; k<tN3; ++k )
                {
                    mVolumeGroupIDs( k ) = mVolumes( k )->geometry_tag();
                }
                unique( mVolumeGroupIDs );
            }

        }

//------------------------------------------------------------------------------

        void
        GmshReader::create_blocks()
        {
            if( mNumberOfDimensions == 2 )
            {
                uint tNumBlocks = mFaceGroupIDs.length();

                Map< uint, uint > tMap;
                for( uint k=0; k<tNumBlocks; ++k )
                {
                    tMap[ mFaceGroupIDs( k ) ] = k;
                }

                Vector< index_t > tElementsPerBlock( tNumBlocks, 0 );
                for( Element * tElement: mFaces )
                {
                    ++tElementsPerBlock( tMap( tElement->geometry_tag() ) );
                }

                mMesh->mBlocks.set_size( tNumBlocks, nullptr );
                for( uint k=0; k<tNumBlocks; ++k )
                {
                    mMesh->mBlocks( k ) = new Block( mFaceGroupIDs( k ), tElementsPerBlock( k ) );
                }

                for( Element * tElement: mFaces )
                {
                    mMesh->mBlocks( tMap( tElement->geometry_tag() ) )->insert_element( tElement );
                }
                
            }
            else if( mNumberOfDimensions == 3 )
            {
                uint tNumBlocks = mVolumeGroupIDs.length();

                Map< uint, uint > tMap;
                for( uint k=0; k<tNumBlocks; ++k )
                {
                    tMap[ mVolumeGroupIDs( k ) ] = k;
                }

                Vector< index_t > tElementsPerBlock( tNumBlocks, 0 );
                for( Element * tElement: mVolumes )
                {
                    ++tElementsPerBlock( tMap( tElement->geometry_tag() ) );
                }

                mMesh->mBlocks.set_size( tNumBlocks, nullptr );

                for( uint k=0; k<tNumBlocks; ++k )
                {
                    mMesh->mBlocks( k ) = new Block( mVolumeGroupIDs( k ), tElementsPerBlock( k ) );
                }

                for( Element * tElement: mVolumes )
                {
                    mMesh->mBlocks( tMap( tElement->geometry_tag() ) )->insert_element( tElement );
                }

                mMesh->mBoundaryEdges.set_size( mEdges.size(), nullptr );

                index_t tCount = 0;

                for( mesh::Element * tEdge : mEdges )
                {
                    // flag edge and its nodes: create_mesh() and create_vertices()
                    // keep only flagged entities
                    tEdge->flag();

                    for( uint k=0; k<tEdge->number_of_nodes(); ++k )
                    {
                        tEdge->node( k )->flag() ;
                    }

                    mMesh->mBoundaryEdges( tCount++ ) = tEdge;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        GmshReader::create_sidesets()
        {
            Vector< uint > tGroupIDs;

            if( mNumberOfDimensions == 2 )
            {
                mMesh->mFacets.set_size( mEdges.size(), nullptr );
                for( size_t k=0; k<mEdges.size(); ++k )
                {
                    mMesh->mFacets( k ) = new Facet( mEdges( k ) );
                }

                tGroupIDs = mEdgeGroupIDs;
            }
            else if( mNumberOfDimensions == 3 )
            {
                mMesh->mFacets.set_size( mFaces.size(), nullptr );
                for( size_t k=0; k<mFaces.size(); ++k )
                {
                    mMesh->mFacets( k ) = new Facet( mFaces( k ) );
                }

                tGroupIDs = mFaceGroupIDs;
            }
            mMesh->mNumberOfDimensions = mNumberOfDimensions;

            uint tNumSideSets = tGroupIDs.length();

            Map< uint, uint > tMap;
            for( uint k=0; k<tNumSideSets; ++k )
            {
                tMap[ tGroupIDs( k ) ] = k;
            }

            Vector< uint > mFacesPerSet( tNumSideSets, 0 );

            for( Facet * tFacet: mMesh->mFacets )
            {
                ++mFacesPerSet( tMap( tFacet->element()->geometry_tag() ) );
            }

            mMesh->mSideSets.set_size( tNumSideSets, nullptr );

            if( mNumberOfDimensions == 2 )
            {
                for( uint k=0; k<tNumSideSets; ++k )
                {
                    mMesh->mSideSets( k ) = new SideSet( mEdgeGroupIDs( k ), mFacesPerSet( k ) );
                }
            }
            else if ( mNumberOfDimensions == 3 )
            {
                for( uint k=0; k<tNumSideSets; ++k )
                {
                    mMesh->mSideSets( k ) = new SideSet( mFaceGroupIDs( k ), mFacesPerSet( k ) );
                }
            }

            for( Facet * tFacet: mMesh->mFacets )
            {
                mMesh->mSideSets( tMap( tFacet->element()->geometry_tag() ) )->insert_facet( tFacet );
            }
        }

//------------------------------------------------------------------------------

        void
        GmshReader::delete_unused_nodes_and_elements()
        {
            for ( Node * tNode: mNodes )
            {
                if ( tNode != nullptr && ! tNode->is_flagged() )
                {
                    delete tNode;
                }
            }

            for ( Element * tElement: mElements )
            {
                if ( tElement != nullptr && ! tElement->is_flagged() )
                {
                    delete tElement;
                }
            }

            // the reader's containers now hold dangling pointers next to
            // pointers the mesh owns: drop them all
            mNodes.clear();
            mNodeMap.clear();
            mElements.clear();
            mVertices.clear();
            mEdges.clear();
            mFaces.clear();
            mVolumes.clear();
        }

//------------------------------------------------------------------------------

        void
        GmshReader::create_mesh()
        {
            size_t tCount=0;

            for( Block * tBlock: mMesh->mBlocks )
            {
                for( Element * tElement: tBlock->mElements )
                {
                    tElement->flag();
                    tElement->flag_nodes();
                    ++tCount;
                }
            }

            mMesh->mElements.set_size( tCount, nullptr );
            tCount = 0;
            for( Block * tBlock: mMesh->mBlocks )
            {
                for( Element * tElement: tBlock->mElements )
                {
                    mMesh->mElements( tCount++ ) = tElement;
                }
            }

            tCount = 0;
            for( SideSet * tSideSet: mMesh->mSideSets )
            {
                for( Facet * tFacet: tSideSet->mFacets )
                {
                    tFacet->element()->flag();
                    tFacet->element()->flag_nodes();
                    ++tCount;
                }
            }

            BELFEM_ERROR( tCount == mMesh->mFacets.size(),
                "Some Faces have not been assigned in file %s",
                    mFilename.c_str() );

            mMesh->mFacets.set_size( tCount, nullptr );

            tCount = 0;
            for( SideSet * tSideSet: mMesh->mSideSets )
            {
                for( Facet * tFacet: tSideSet->mFacets )
                {
                    mMesh->mFacets( tCount++ ) = tFacet;
                }
            }

            tCount = 0;
            for( Node * tNode: mNodes )
            {
                if( tNode->is_flagged() )
                {
                    ++tCount;
                }
            }

            mMesh->mNodes.set_size( tCount, nullptr );
            tCount = 0;
            for( Node * tNode: mNodes )
            {
                if( tNode->is_flagged() )
                {
                    mMesh->mNodes( tCount ++ ) = tNode;
                }
            }
        }

 //------------------------------------------------------------------------------

        void
        GmshReader::create_vertices()
        {

            index_t tCount = 0;

            for ( Element * tElement : mVertices )
            {
                if( tElement->node( 0 )->is_flagged() )
                {
                    ++tCount;

                    tElement->flag();
                }
            }

            Cell< Element * > & tVertices = mMesh->mVertices;

            tVertices.set_size( tCount, nullptr );

            tCount = 0;

            for ( Element * tElement : mVertices )
            {
                 if( tElement->is_flagged() )
                 {
                     tVertices( tCount++ ) = tElement;
                 }
            }
        }

//------------------------------------------------------------------------------
    }
}
