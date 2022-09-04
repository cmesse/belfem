//
// Created by Christian Messe on 2019-07-25.
//

#include "assert.hpp"
#include "stringtools.hpp"
#include "filetools.hpp"
#include "cl_Matrix.hpp"

#include "commtools.hpp"

#include "cl_Mesh_GmshReader.hpp"
#include "cl_Element_Factory.hpp"


#include "fn_min.hpp"
#include "fn_max.hpp"
#include "fn_unique.hpp"
#include "cl_Mesh_OrientationChecker.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        GmshReader::GmshReader( const string & aPath, Mesh * aMesh, const real aMeshScale ) :
            Ascii( aPath, FileMode::OPEN_RDONLY ),
            mMeshScale( aMeshScale )
        {
            // remember name
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
                mMesh = new Mesh( comm_rank() );
                mOwnMesh = true;
            }
            else
            {
                mMesh = aMesh;
            }

            // for higher order elements, the order is different in exodus
            this->convert_orders_for_quadratic_volume_elements_to_exo();
            this->create_blocks();
            this->check_element_orientation();
            this->create_sidesets();
            this->create_mesh();

            this->create_vertices();

            // tidy up memory
            this->delete_unused_nodes_and_elements();


            // finish mesh
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
            // loop over all lines
            for( auto tLine: mBuffer )
            {
                // remove tabs, double spaces etc
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
            // counter
            size_t tCount = 0;

            // loop over all possible lines
            for( auto tLine: mBuffer )
            {
                if ( tLine == "$MeshFormat")
                {
                    break;
                }
                ++tCount;
            }

            // make sure that a version was found
            BELFEM_ERROR( tCount < mBuffer.size(),
                    "The file %s does not seem to be a GMSH mesh file",
                         mFilename.c_str() );

            // read version string
            mVersion = std::stod( first_word( mBuffer( ++tCount ) ) );
        }

//------------------------------------------------------------------------------

        void
        GmshReader::read_tags_v22()
        {
            // line counter
            size_t tCount = 0;

            // number of lines
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

                // increment counter
                ++tCount;
            }
        }

//------------------------------------------------------------------------------

        void
        GmshReader::read_nodes_tag_v22( size_t & aCount )
        {
            mNodesTag = aCount;

            // read number of nodes
            mNumberOfNodes = std::stoi( mBuffer( ++aCount ) );

            // increment counter
            aCount += mNumberOfNodes + 1;
        }

//------------------------------------------------------------------------------

        void
        GmshReader::read_elements_tag_v22( size_t & aCount )
        {
            mElementsTag = aCount;

            // read number of elements
            mNumberOfElements = std::stoi( mBuffer( ++aCount ) );

            // increment counter
            aCount += mNumberOfElements + 1;
        }

//------------------------------------------------------------------------------

        void
        GmshReader::read_physical_tag( size_t & aCount )
        {
            // remember tag
            mPhysicalTag = aCount;

            // read number of physical groups
            mNumberOfPhysicalGroups = std::stoi( mBuffer( ++aCount ) );

            // increment counter
            aCount += mNumberOfPhysicalGroups + 1;
        }

//------------------------------------------------------------------------------

        void
        GmshReader::read_nodes_v22()
        {
            // reserve memory for nodes
            mNodes.set_size( mNumberOfNodes, nullptr );

            // Node Counter
            size_t tCount = 0;

            // sart line
            uint tStart = mNodesTag + 2;

            // end line
            uint tEnd = tStart + mNumberOfNodes;

            for ( size_t k = tStart; k < tEnd; ++k )
            {
                Cell< string > tWords = string_to_words( mBuffer( k ));

                // read id of node
                id_t tID = std::stoi( tWords( 0 ));

                // read x-coordinate
                real tX = std::stod( tWords( 1 )) * mMeshScale ;

                // read y-coordinate
                real tY = std::stod( tWords( 2 )) * mMeshScale ;

                // read z-coordinate
                real tZ = std::stod( tWords( 3 )) * mMeshScale ;

                // create the node object
                mNodes( tCount++ ) = new Node( tID, tX, tY, tZ );
            }

        }

//------------------------------------------------------------------------------

        void
        GmshReader::read_elements_v22()
        {
            // reserve memory for elements
            mElements.set_size( mNumberOfElements, nullptr );

            // Element Counter
            size_t tCount = 0;

            // sart line
            uint tStart = mElementsTag + 2;

            // end line
            uint tEnd = tStart + mNumberOfElements;

            // create the factory
            ElementFactory tFactory;

            for( uint k=tStart; k<tEnd; ++k )
            {
                // create words
                Cell< string > tWords = string_to_words( mBuffer( k ) );

                // get the element ID
                id_t tID = std::stoi( tWords( 0 ) );

                // get the gmsh type
                ElementType tType = element_type_from_gmsh(  std::stoi( tWords( 1 ) ) );

                // number of nodes for this element
                index_t tNumNodes = number_of_nodes( tType );

                // get the number of tags
                uint tNumTags = std::stoi( tWords( 2 ) );

                // make sure that number of elements is correct
                BELFEM_ERROR( tWords.size()-3-tNumTags == tNumNodes,
                             "Number of nodes for element %u does not match. Is: %u, Expect %u",
                             ( unsigned int ) tID,
                             ( unsigned int ) tWords.size()-3-tNumTags,
                             ( unsigned int ) tNumNodes );


                // create an element
                Element * tElement = tFactory.create_lagrange_element( tType, tID );

                // set the physical tag of the element
                tElement->set_physical_tag( std::stoi( tWords( 3 ) ) );

                // set the geometry tag of the element
                tElement->set_geometry_tag( std::stoi( tWords( 4 ) ) );

                // offset in gmsh string
                index_t tOff = tNumTags + 3;

                // insert nodes into element
                for ( index_t i=0; i<tNumNodes; ++i )
                {
                    tElement->insert_node( mNodeMap( std::stoi( tWords( i + tOff ) )  ) , i );
                }

                // add element to list
                mElements( tCount++ ) = tElement;
            }
        }

 //------------------------------------------------------------------------------

        void
        GmshReader::read_mesh_v41()
        {
            // line counter
            size_t tCount = 0;

            // number of lines
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

                // increment counter
                ++tCount;
            }
        }

//------------------------------------------------------------------------------

        void
        GmshReader::read_entity_tag_v41( size_t & aCount )
        {
            // convert next line to words
            Cell< string > tWords = string_to_words( mBuffer( ++aCount ) );

            // increment counter
            for( uint k=0; k<tWords.size(); ++k )
            {
                aCount += std::stoi( tWords( k ) );
            }
        }

//------------------------------------------------------------------------------

        void
        GmshReader::read_nodes_v41( size_t & aCount )
        {
            // remember the tag
            mNodesTag = aCount++;

            Cell< string > tWords = string_to_words( mBuffer( aCount++ ) );

            // get number of elements
            mNumberOfNodes = std::stoi( tWords( 1 ) );

            // node counter
            uint tNodeCount = 0;

            // vector containing node IDs
            Vector< id_t > mNodeIDs( mNumberOfNodes );
            Matrix< real > mNodeCoords( 3, mNumberOfNodes );

            while( tNodeCount < mNumberOfNodes )
            {
                tWords = string_to_words( mBuffer( aCount++ ) );

                // read parametric flag
                uint tParametric = std::stoi( tWords( 2 ) );

                // make sure that there are no parametric nodes
                BELFEM_ERROR( tParametric == 0,
                    "The mesh %s contains parameteric nodes, which are not supported by BELFEM",
                    mFilename.c_str() );

                // read number of nodes
                uint tN = std::stoi( tWords( 3 ) );

                // read node ids
                for( uint k=0; k<tN; ++k )
                {
                    mNodeIDs( tNodeCount + k ) = std::stoi( mBuffer( aCount++ ) );
                }

                // read node coordinates
                for( uint k=0; k<tN; ++k )
                {
                    tWords = string_to_words( mBuffer( aCount++ ) );
                    mNodeCoords( 0, tNodeCount + k ) = std::stod( tWords( 0 ) ) * mMeshScale ;
                    mNodeCoords( 1, tNodeCount + k ) = std::stod( tWords( 1 ) ) * mMeshScale ;
                    mNodeCoords( 2, tNodeCount + k ) = std::stod( tWords( 2 ) ) * mMeshScale ;
                }

                tNodeCount += tN;
            }

            // create nodes
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
            // Allocate pointer with Nodes
            Cell< Node * > tNodes;
            tNodes.set_size( 27, nullptr );

            // gmsh to exodus conversion table
            uint tHexIndex[ 27 ] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 9, 10,
                                     12, 14, 15, 16, 18, 19, 17, 26, 20, 25,
                                     22, 23, 21, 24 };

            uint tPentaIndex[ 18 ] ={ 0, 1, 2, 3, 4, 5, 6, 9, 7, 8, 10, 11,
                                      12, 14, 13, 15, 17, 16 };

            // loop over all Elements
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

                        // copy pointers into temporary buffer
                        for ( uint k=0; k<tNumNodes; ++k )
                        {
                            tNodes( k ) = tElement->node( k );
                        }

                        // copy pointers back into Element
                        for ( uint k=0; k<tNumNodes; ++k )
                        {
                            tElement->insert_node( tNodes( tPentaIndex[ k ] ), k );
                        }

                        break;
                    }
                    case( ElementType::HEX20 ):
                    case( ElementType::HEX27 ):
                    {
                        uint tNumNodes = tElement->number_of_nodes();

                        // copy pointers into temporary buffer
                        for ( uint k=0; k<tNumNodes; ++k )
                        {
                            tNodes( k ) = tElement->node( k );
                        }

                        // copy pointers back into Element
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

            // remember the element tag
            mElementsTag = aCount++;

            Cell<string> tWords = string_to_words( mBuffer( aCount++ ));

            // get number of elements
            mNumberOfElements = std::stoi( tWords( 1 ));

            // allocate  container
            mElements.set_size( mNumberOfElements, nullptr );

            uint tElementCount = 0;

            // create the factory
            ElementFactory tFactory;

            while ( tElementCount < mNumberOfElements )
            {
                // get next line
                tWords = string_to_words( mBuffer( aCount++ ));

                // get geometry tag
                uint tGeometryTag = std::stoi( tWords( 1 ));

                // get element type
                ElementType tType = element_type_from_gmsh( std::stoi(  tWords( 2 ) ));

                // number of elements in this block
                uint tNumElements = std::stoi( tWords( 3 ) );

                // get number of nodes
                size_t tNumNodes = number_of_nodes( tType );

                for ( uint k = 0; k < tNumElements; ++k )
                {
                    // read buffer
                    tWords = string_to_words( mBuffer( aCount++ ) );

                    // read element ID
                    id_t tID = std::stoi( tWords( 0 ) );

                    // make sure that number of elements is correct
                    BELFEM_ERROR( tWords.size() == tNumNodes + 1,
                                 "Number of nodes for element %ui does not match. Is: %u, Expect %u",
                                 ( unsigned int ) tID,
                                 ( unsigned int ) tWords.size()-1,
                                 ( unsigned int ) tNumNodes );

                    // create an element
                    Element * tElement = tFactory.create_lagrange_element( tType, tID );

                    // add nodes
                    for( uint i=0; i < tNumNodes; ++i )
                    {
                        tElement->insert_node( mNodeMap( std::stoi( tWords( i + 1 ) ) ), i );
                    }

                    // set the geometry tag
                    tElement->set_geometry_tag( tGeometryTag );

                    // add element to list
                    mElements( tElementCount++ ) = tElement;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        GmshReader::create_node_map()
        {
            // reset map
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

            // loop over all nodes
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
            // init vector
            mNumberOfElementsPerDimension.set_size( 4, 0 );

            // loop over all elements and count elements per dimension
            for( Element * tElement: mElements )
            {
                // get dimension of element and increment counter
                ++mNumberOfElementsPerDimension( dimension( tElement->type() ) );
            }

            // perform some small sanity checks
            if( mNumberOfDimensions == 3 )
            {
                BELFEM_ERROR( mNumberOfElementsPerDimension( 3 ) > 0,
                        "The file %s has 3D node coordinates, but no 3D elements were found",
                        mFilename.c_str() );
            }
            else if( mNumberOfDimensions == 2 )
            {
                BELFEM_ERROR( mNumberOfElementsPerDimension( 2 ) > 0,
                             "The file %s has 2D node coordinates, but no 2D elements were found",
                             mFilename.c_str() );

                BELFEM_ERROR( mNumberOfElementsPerDimension( 3 ) == 0,
                             "The file %s has 2D node coordinates, but contains 3D element types",
                             mFilename.c_str() );
            }

            // allocate containers
            mVertices.set_size( mNumberOfElementsPerDimension( 0 ), nullptr );
            mEdges.set_size( mNumberOfElementsPerDimension( 1 ), nullptr );
            mFaces.set_size( mNumberOfElementsPerDimension( 2 ), nullptr );
            mVolumes.set_size( mNumberOfElementsPerDimension( 3 ), nullptr );

            // reset container
            mNumberOfElementsPerDimension.fill( 0 );

            // loop over all elements and add elements per dimension
            for( Element * tElement: mElements )
            {
                // get dimension of element
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
            // set initial values

            // get tags of 0D elements
            if( tN0 > 0 )
            {
                // initialize vector
                mNodeGroupIDs.set_size( tN0, mVertices( 0 )->geometry_tag() );

                for( uint k=1; k<tN0; ++k )
                {
                    mNodeGroupIDs( k ) = mVertices( k )->geometry_tag();
                }
                unique( mNodeGroupIDs );
            }

            // get tags of 1D elements
            if( tN1 > 0 )
            {
                // initialize vector
                mEdgeGroupIDs.set_size( tN1, mEdges( 0 )->geometry_tag() );

                for( uint k=1; k<tN1; ++k )
                {
                    mEdgeGroupIDs( k ) = mEdges( k )->geometry_tag();
                }
                unique( mEdgeGroupIDs );
            }

            // get tags of 2D elements
            if( tN2 > 0 )
            {
                // initialize vector
                mFaceGroupIDs.set_size( tN2, mFaces( 0 )->geometry_tag() );

                for( uint k=1; k<tN2; ++k )
                {
                    mFaceGroupIDs( k ) = mFaces( k )->geometry_tag();
                }
                unique( mFaceGroupIDs );
            }

            // get tags of 3D elements
            if( tN3 > 0 )
            {
                // initialize vector
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
                // number of blocks
                uint tNumBlocks = mFaceGroupIDs.length();

                // create id map
                Map< uint, uint > tMap;
                for( uint k=0; k<tNumBlocks; ++k )
                {
                    tMap[ mFaceGroupIDs( k ) ] = k;
                }

                // count elements per block
                Vector< index_t > tElementsPerBlock( tNumBlocks, 0 );
                for( Element * tElement: mFaces )
                {
                    ++tElementsPerBlock( tMap( tElement->geometry_tag() ) );
                }

                // allocate blocks
                mMesh->mBlocks.set_size( tNumBlocks, nullptr );
                for( uint k=0; k<tNumBlocks; ++k )
                {
                    mMesh->mBlocks( k ) = new Block( mFaceGroupIDs( k ), tElementsPerBlock( k ) );
                }

                // add elements to block
                for( Element * tElement: mFaces )
                {
                    mMesh->mBlocks( tMap( tElement->geometry_tag() ) )->insert_element( tElement );
                }
                
            }
            else if( mNumberOfDimensions == 3 )
            {
                // number of blocks
                uint tNumBlocks = mVolumeGroupIDs.length();

                // create id map
                Map< uint, uint > tMap;
                for( uint k=0; k<tNumBlocks; ++k )
                {
                    tMap[ mVolumeGroupIDs( k ) ] = k;
                }

                // count elements per block
                Vector< index_t > tElementsPerBlock( tNumBlocks, 0 );
                for( Element * tElement: mVolumes )
                {
                    ++tElementsPerBlock( tMap( tElement->geometry_tag() ) );
                }

                // allocate blocks
                mMesh->mBlocks.set_size( tNumBlocks, nullptr );

                for( uint k=0; k<tNumBlocks; ++k )
                {
                    mMesh->mBlocks( k ) = new Block( mVolumeGroupIDs( k ), tElementsPerBlock( k ) );
                }

                /// insert elements into Blocks
                for( Element * tElement: mVolumes )
                {
                    mMesh->mBlocks( tMap( tElement->geometry_tag() ) )->insert_element( tElement );
                }


                // allocate edges container
                mMesh->mBoundaryEdges.set_size( mEdges.size(), nullptr );

                // reset counter
                index_t tCount = 0;

                // copy edges in container
                for( mesh::Element * tEdge : mEdges )
                {
                    // flag edge so that it is not deleted during destruction of reader
                    tEdge->flag();

                    // make sure that all nodes on the edge are flagged
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
        GmshReader::check_element_orientation()
        {
            OrientationChecker tChecker ;
            for( Block * tBlock :  mMesh->mBlocks )
            {
                tChecker.set_element_type( tBlock->element_type() );
                for( Element * tElement : tBlock->elements() )
                {
                    tChecker.process_element( tElement );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        GmshReader::create_sidesets()
        {
            // container for ID numbers
            Vector< uint > tGroupIDs;

            // create facets
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

            // create the ID map
            Map< uint, uint > tMap;
            for( uint k=0; k<tNumSideSets; ++k )
            {
                tMap[ tGroupIDs( k ) ] = k;
            }

            // count faces per set
            Vector< uint > mFacesPerSet( tNumSideSets, 0 );

            for( Facet * tFacet: mMesh->mFacets )
            {
                ++mFacesPerSet( tMap( tFacet->element()->geometry_tag() ) );
            }

            // create sidesets
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

            // insert Facets into Sidesets

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
                if ( ! tNode->is_flagged() )
                {
                    delete tNode;
                }
            }

            for ( Element * tElement: mElements )
            {
                if ( ! tElement->is_flagged() )
                {
                    delete tElement;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        GmshReader::create_mesh()
        {
            // count elements
            size_t tCount=0;

            // flag elements and nodes blockwise
            for( Block * tBlock: mMesh->mBlocks )
            {
                for( Element * tElement: tBlock->mElements )
                {
                    tElement->flag();
                    tElement->flag_nodes();
                    ++tCount;
                }
            }

            // set memory for mesh
            mMesh->mElements.set_size( tCount, nullptr );
            tCount = 0;
            for( Block * tBlock: mMesh->mBlocks )
            {
                for( Element * tElement: tBlock->mElements )
                {
                    mMesh->mElements( tCount++ ) = tElement;
                }
            }

            // flag elements and nodes sideset wise
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

            // sanity check
            BELFEM_ERROR( tCount == mMesh->mFacets.size(),
                "Some Faces have not been assigned in file %s",
                    mFilename.c_str() );

            // set memory for mesh object
            mMesh->mFacets.set_size( tCount, nullptr );

            // add facets
            tCount = 0;
            for( SideSet * tSideSet: mMesh->mSideSets )
            {
                for( Facet * tFacet: tSideSet->mFacets )
                {
                    mMesh->mFacets( tCount++ ) = tFacet;
                }
            }

            // count used nodes
            tCount = 0;
            for( Node * tNode: mNodes )
            {
                if( tNode->is_flagged() )
                {
                    ++tCount;
                }
            }

            // copy nodes into mesh
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
            // count bearings that are connected to elements

            index_t tCount = 0;

            // loop over all bearings
            for ( Element * tElement : mVertices )
            {
                // test if the node on this elememnt is used by any other element
                if( tElement->node( 0 )->is_flagged() )
                {
                    // increment counter
                    ++tCount;

                    // flag this element
                    tElement->flag();
                }
            }

            // allocate bearing container
            // get vertex container
            Cell< Element * > & tVertices = mMesh->mVertices;

            tVertices.set_size( tCount, nullptr );

            // reset counter
            tCount = 0;

            // loop over all bearings
            for ( Element * tElement : mVertices )
            {
                 // test if element is flagged
                 if( tElement->is_flagged() )
                 {
                     // add this element to the container
                     tVertices( tCount++ ) = tElement;
                 }
            }
        }

//------------------------------------------------------------------------------
    }
}