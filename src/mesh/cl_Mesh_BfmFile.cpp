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
#include <cstdio>

#include "belfem_version.hpp"
#include "banner.hpp"
#include "cl_Logger.hpp"
#include "cl_Mesh.hpp"
#include "cl_ProtoMesh.hpp"
#include "cl_Mesh_BfmFile.hpp"

#include "cl_HDF5_Dataset.hpp"
#include "fn_sum.hpp"
#include "cl_Element_Factory.hpp"


namespace belfem
{
    namespace mesh
    {
        BfmFile::BfmFile( const string & aFilePath, Mesh * aMesh ) :
            mFilePath( aFilePath ),
            mMesh( aMesh ),
            mOwnMesh( aMesh == nullptr )
        {

        }

        BfmFile::~BfmFile()
        {
            if( mOwnMesh ) delete mMesh ;
            if ( mProto != nullptr ) delete mProto ;
            if ( mFile != nullptr ) delete mFile ;
        }

        size_t
        BfmFile::checksum() const
        {
            HDF5 * tFile = new HDF5( mFilePath, FileMode::OPEN_RDONLY ) ;

            size_t aChecksum ;
            tFile->select_group( "meta" );
            tFile->load_data( "checksum", aChecksum );
            tFile->close_active_group() ;
            tFile->close() ;
            delete tFile ;

            return aChecksum ;
        }

//------------------------------------------------------------------------------

        uint64_t
        BfmFile::config_tag() const
        {
            uint64_t tTag = 0 ;

#ifdef BELFEM_HDF5
            HDF5 * tFile = new HDF5( mFilePath, FileMode::OPEN_RDONLY ) ;

            tFile->select_group( "meta" );

            if ( hdf5::dataset_exists( tFile->active_group(), "config" ) )
            {
                tFile->load_data( "config", tTag );
            }

            tFile->close_active_group() ;
            tFile->close() ;
            delete tFile ;
#endif
            return tTag ;
        }

//------------------------------------------------------------------------------

        string
        BfmFile::config_text() const
        {
            string tText ;

#ifdef BELFEM_HDF5
            HDF5 * tFile = new HDF5( mFilePath, FileMode::OPEN_RDONLY ) ;

            tFile->select_group( "meta" );

            if ( hdf5::dataset_exists( tFile->active_group(), "config_text" ) )
            {
                tFile->load_data( "config_text", tText );
            }

            tFile->close_active_group() ;
            tFile->close() ;
            delete tFile ;
#endif
            return tText ;
        }

//------------------------------------------------------------------------------

        void
        BfmFile::save()
        {
            mMesh->update_node_indices();
            mMesh->update_element_indices();

            BELFEM_ASSERT( mFile == nullptr, "File already allocated" ) ;

            mFile = new HDF5( mFilePath, FileMode::NEW ) ;

            mProto = new ProtoMesh( mMesh ) ;


            this->save_meta_data();
            this->save_group_data( GroupType::BLOCK );
            this->save_group_data( GroupType::SIDESET );


            this->save_node_data() ;
            this->save_element_data();
            this->save_facet_data( false );
            this->save_edge_data();

            this->save_face_data();

            this->save_control_point_data();

            this->save_node_duplicate_data();

            this->save_hanging_entities();

            this->save_periodicity_data();

            this->save_thinshell_data();

            this->save_vertex_data();
            this->save_curve_data();

            delete mProto ;
            mProto = nullptr;
        }

        void
        BfmFile::load()
        {
#ifdef BELFEM_HDF5
            BELFEM_ASSERT( mFile == nullptr, "File already allocated" ) ;
            mFile = new HDF5( mFilePath, FileMode::OPEN_RDONLY );

            mFile->select_group( "meta" );
            uint tNumDimensions ;
            mFile->load_data( "dimensions" , tNumDimensions );
            mFile->close_active_group() ;

            if ( mMesh == nullptr )
            {
                mMesh = new Mesh( tNumDimensions );
            }

            mProto = new ProtoMesh( mMesh ) ;

            this->load_meta_data();

            mMesh->force_checksum( mProto->meta_data()->mChecksum );

            // carry the settings stamp back onto the mesh. Without this a mesh
            // that is loaded and then saved again drops its tag, and the file it
            // writes is treated as untagged for ever after -- it would rebuild
            // on every run and never become valid
            mMesh->set_config_tag( mProto->meta_data()->mConfigTag,
                                   mProto->meta_data()->mConfigText );

            this->load_group_data( GroupType::BLOCK );
            this->load_group_data( GroupType::SIDESET );

            this->load_node_data() ;
            this->load_element_data();
            this->load_facet_data();
            this->load_edge_data();
            this->load_face_data();

            this->load_control_point_data();
            this->load_node_duplicate_data();
            this->load_hanging_entities();

            // the hangs exist only now: re-point 1:1 edge-on-edge sources to
            // the twin favored by the slot reconstruction ( B8 ruling )
            mProto->normalize_edge_hangs();

            this->load_periodicity_data();
            this->load_thinshell_data();

            this->load_vertex_data();
            this->load_curve_data();

            mFile->close();

            delete mProto ;
            mProto = nullptr;


            mMesh->finalize();

            // we trust BFM meshes: the checker ran before the file was
            // written, so the stored element orientation is correct by
            // construction and must never be re-checked ( the mesh is
            // enriched — a rerun would abort on the edges-exist guard )
            mMesh->set_mesh_checker_flag();
#endif
        }

        void
        BfmFile::save_meta_data()
        {
#ifdef BELFEM_HDF5
            mProto->populate_meta_data();
            mFile->create_group( "meta" );
            mFile->save_data( "dimensions" , mProto->meta_data()->mNumberOfDimensions );
            mFile->save_data( "entities" , mProto->meta_data()->mNumberOfEntities );
            mFile->save_data( "groups" , mProto->meta_data()->mNumberOfGroups );
            mFile->save_data( "checksum", mMesh->checksum() );

            // we save the version of BELFEM for future reference
            mFile->save_data( "belfem", version() );

            // git provenance of the writing build: the checksum only covers
            // the base mesh, so this stamp is what identifies a file written
            // by an older construction state
            if ( is_built_from_git() )
            {
                mFile->save_data( "git", git_is_dirty() ?
                    git_commit_hash() + "-dirty" : git_commit_hash() );
                mFile->save_data( "branch", git_branch() );
            }

            // fingerprint of the settings this mesh was enriched with. The
            // checksum above is the base mesh only, so this is what catches an
            // edited deck on an unchanged .msh. The text is stored too, so a
            // mismatch can name the setting that changed
            if ( mMesh->config_tag() != 0 )
            {
                mFile->save_data( "config", mMesh->config_tag() );
                mFile->save_data( "config_text", mMesh->config_text() );
            }
            mFile->close_active_group();
#endif
        }

        void
        BfmFile::load_meta_data()
        {
#ifdef BELFEM_HDF5
            mFile->select_group( "meta" );
            mFile->load_data( "dimensions" , mProto->meta_data()->mNumberOfDimensions );
            mFile->load_data( "entities" , mProto->meta_data()->mNumberOfEntities );
            mFile->load_data( "groups" , mProto->meta_data()->mNumberOfGroups );
            mFile->load_data( "checksum", mProto->meta_data()->mChecksum );

            // provenance stamps ( absent on files written before 0.9.0 )
            if ( hdf5::dataset_exists( mFile->active_group(), "belfem" ) )
            {
                mFile->load_data( "belfem", mProto->meta_data()->mBelfemVersion );
            }
            if ( hdf5::dataset_exists( mFile->active_group(), "git" ) )
            {
                mFile->load_data( "git", mProto->meta_data()->mGitHash );
            }

            // settings fingerprint ( absent on files written before it existed;
            // the reuse gate treats absent as a mismatch and rebuilds )
            if ( hdf5::dataset_exists( mFile->active_group(), "config" ) )
            {
                mFile->load_data( "config", mProto->meta_data()->mConfigTag );
            }
            if ( hdf5::dataset_exists( mFile->active_group(), "config_text" ) )
            {
                mFile->load_data( "config_text", mProto->meta_data()->mConfigText );
            }
            mFile->close_active_group();

            // surface the provenance: the first question on any reload
            // oddity is which build wrote the file
            if ( mProto->meta_data()->mBelfemVersion.size() > 0 )
            {
                // the hash goes on its own line: 40 hex plus "-dirty" does
                // not fit beside the version inside 80 columns
                message( InfoLevel::Default,
                    "    Mesh file written by BELFEM %s\n    ( git %s )",
                    mProto->meta_data()->mBelfemVersion.c_str(),
                    mProto->meta_data()->mGitHash.size() > 0 ?
                        mProto->meta_data()->mGitHash.c_str() : "unknown" );

                // a file from a NEWER code may carry schema this build does
                // not know: unknown datasets are skipped silently by the
                // dataset_exists guards, so the user must be told
                int tMajor = 0 ;
                int tMinor = 0 ;
                int tPatch = 0 ;
                if ( std::sscanf( mProto->meta_data()->mBelfemVersion.c_str(),
                        "%d.%d.%d", & tMajor, & tMinor, & tPatch ) == 3 )
                {
                    const bool tIsNewer =
                           tMajor > gVersionMajor
                        || ( tMajor == gVersionMajor && tMinor > gVersionMinor )
                        || ( tMajor == gVersionMajor && tMinor == gVersionMinor
                             && tPatch > gVersionPatch ) ;

                    if ( tIsNewer )
                    {
                        message( InfoLevel::Default,
                            "    Warning: this mesh file was written by BELFEM %s, but this build is %s.",
                            mProto->meta_data()->mBelfemVersion.c_str(),
                            gVersionString );
                        message( InfoLevel::Default,
                            "             Datasets unknown to this build are ignored silently." );
                    }
                }
            }
#endif
        }

        void
        BfmFile::save_node_data()
        {
#ifdef BELFEM_HDF5
            mProto->populate_node_data();
            mFile->create_group( "nodes" );

            mFile->save_data( "ids" ,    mProto->node_data()->mIDs );
            mFile->save_data( "coords" , mProto->node_data()->mCoords, true );

            if ( mMesh->abstract_nodes().size() > 0 )
            {
                Cell< id_t > tData( mMesh->abstract_nodes().size() );

                for ( Node * tNode : mMesh->abstract_nodes() )
                {
                    tData.push( tNode->id());
                }
                mFile->save_data( "abstract" , tData );
            }

            if ( mMesh->autopins().size() > 0 )
            {
                Cell< id_t > tData( mMesh->autopins().size() );

                for ( Node * tNode : mMesh->autopins() )
                {
                    tData.push( tNode->id());
                }
                mFile->save_data( "pinned" , tData );
            }

            if ( mMesh->orphaned_nodes().size() > 0 )
            {
                Cell< id_t > tData( mMesh->orphaned_nodes().size() );

                for ( Node * tNode : mMesh->orphaned_nodes() )
                {
                    tData.push( tNode->id() );
                }
                mFile->save_data( "orphaned" , tData );
            }

            mFile->close_active_group();
            mProto->reset_node_data();
#endif
        }

        void
        BfmFile::load_node_data()
        {
#ifdef BELFEM_HDF5
            mFile->select_group( "nodes" );

            mFile->load_data( "ids" , mProto->node_data()->mIDs );
            mFile->load_data( "coords" , mProto->node_data()->mCoords, true );



            mProto->create_nodes();

            if ( hdf5::group_exists( mFile->active_group(), "abstract" ) )
            {
                Cell< id_t > tData ;
                mFile->load_data( "abstract", tData );
                Cell< Node * > & tNodes = mMesh->abstract_nodes();
                tNodes.reserve( tData.size() );
                for ( id_t tID : tData )
                {
                    tNodes.push( mProto->node( tID ) );
                }
            }

            if ( hdf5::group_exists( mFile->active_group(), "pinned" ) )
            {
                Cell< id_t > tData ;
                mFile->load_data( "pinned", tData );
                Cell< Node * > & tNodes = mMesh->autopins();
                tNodes.reserve( tData.size() );
                for ( id_t tID : tData )
                {
                    tNodes.push( mProto->node( tID ) );
                }
            }

            if ( hdf5::group_exists( mFile->active_group(), "orphaned" ) )
            {
                Cell< id_t > tData ;
                mFile->load_data( "orphaned", tData );
                Cell< Node * > & tNodes = mMesh->orphaned_nodes();
                tNodes.reserve( tData.size() );
                for ( id_t tID : tData )
                {
                    tNodes.push( mProto->node( tID ) );
                }
            }

            mFile->close_active_group();
#endif
        }

        void
        BfmFile::save_group_data( const GroupType aType )
        {
#ifdef BELFEM_HDF5
            Cell< proto::GroupData > & tGroups =
                aType == GroupType::BLOCK ? mProto->block_data() : mProto->sideset_data() ;

            if ( aType == GroupType::BLOCK )
            {
                mProto->populate_block_data();
            }
            else
            {
                mProto->populate_sideset_data();
            }
            uint tNumGroups = tGroups.size();

            Cell< id_t >   tIDs( tNumGroups );
            Cell< string > tLabels( tNumGroups );
            Cell< index_t > tNumElems( tNumGroups );
            Cell< uint > tElemTypes( tNumGroups );
            Cell< uint > tDomainTypes( tNumGroups );
            Cell< uchar > tTags( tNumGroups );

            for ( proto::GroupData & tData : tGroups )
            {
                BELFEM_ERROR( tData.mID != 0 && tData.mID != gNoID, "invalid group id" );

                tIDs.push( tData.mID );
                tLabels.push( tData.mLabel );
                tNumElems.push( tData.mNumElements );
                tElemTypes.push( static_cast< uint >( tData.mElementType ) );
                tDomainTypes.push( static_cast< uint >( tData.mDomainType ) );

                uchar tTag = 0 ;
                if ( tData.mHidden )   tTag += 1 ;
                if ( tData.mHasEdges ) tTag += 2 ;
                if ( tData.mHasFaces ) tTag += 4 ;
                tTags.push( tTag );
            }


            mFile->create_group( aType == GroupType::BLOCK ? "blocks" : "sidesets" );

            mFile->save_data( "ids" , tIDs );
            mFile->save_data( "labels" , tLabels );
            mFile->save_data( "elements" , tNumElems );
            mFile->save_data( "types" , tElemTypes );
            mFile->save_data( "domains" , tDomainTypes );
            mFile->save_data( "tags" , tTags );


            mFile->close_active_group();
#endif
        }

        void
        BfmFile::load_group_data( const GroupType aType )
        {
#ifdef BELFEM_HDF5
            Cell< id_t >   tIDs ;
            Cell< string > tLabels ;
            Cell< index_t > tNumElems ;
            Cell< uint > tElemTypes ;
            Cell< uint > tDomainTypes ;
            Cell< uchar > tTags ;

            mFile->select_group( aType == GroupType::BLOCK ? "blocks" : "sidesets" );

            mFile->load_data( "ids" , tIDs );
            mFile->load_data( "labels" , tLabels );
            mFile->load_data( "elements" , tNumElems );
            mFile->load_data( "types" , tElemTypes );
            mFile->load_data( "domains" , tDomainTypes );
            mFile->load_data( "tags" , tTags );

            uint tNumGroups = tIDs.size();

            Cell< proto::GroupData > & tGroups =
              aType == GroupType::BLOCK ? mProto->block_data() : mProto->sideset_data() ;

            tGroups.set_size( tNumGroups );

            uint g = 0 ;
            for ( proto::GroupData & tGroup : tGroups )
            {
                tGroup.mID = tIDs( g );
                tGroup.mLabel = tLabels( g );
                tGroup.mNumElements = tNumElems( g );
                tGroup.mElementType = static_cast< ElementType >( tElemTypes( g ) );
                tGroup.mDomainType  = static_cast< DomainType >( tDomainTypes( g ) );

                // decode the tag bits ( symmetric with save_group_data )
                uchar tTag = tTags( g );
                tGroup.mHidden   = tTag & ( 1u << 0 );
                tGroup.mHasEdges = tTag & ( 1u << 1 );
                tGroup.mHasFaces = tTag & ( 1u << 2 );

                ++g ;
            }

            mFile->close_active_group();
#endif
        }

        void
        BfmFile::save_element_data()
        {
#ifdef BELFEM_HDF5
            mProto->populate_element_data( false, false );

            mFile->create_group( "elements" );

            mFile->save_data( "ids", mProto->element_data()->mIDs );

            // check if the physical tags are even set, otherwise we don't
            // need to save them

            if ( mProto->element_data()->mPhysicalTags.size() > 0 )
            {
                mFile->save_data( "physical", mProto->element_data()->mPhysicalTags );
            }

            // we store the topology unflattened. This has to be done manually
            hsize_t tNumElems = mProto->element_data()->mIDs.size();

            index_t tCount = 0 ;

            hdf5::Dataset< id_t > tData( mFile->active_group(), "topology", tNumElems );

            for ( Block * tBlock : mMesh->blocks() )
            {
                Cell< Element * > & tElements = tBlock->elements();
                hsize_t n = number_of_nodes( tBlock->element_type() );
                for ( Element * tElement : tElements )
                {
                    id_t * tNodes =  tData.set_size( tCount++, n );

                    for ( hsize_t k=0; k<n; ++k )
                    {
                        tNodes[ k ] = tElement->node( k )->id();
                    }
                }
            }

            tData.save();

            mFile->close_active_group() ;
            mProto->reset_element_data();
#endif
        }

        void
        BfmFile::load_element_data()
        {
#ifdef BELFEM_HDF5
            mFile->select_group( "elements" );

            mFile->load_data( "ids", mProto->element_data()->mIDs );

            index_t tNumElems = mProto->element_data()->mIDs.size();

            // check if physical dataset exists
            if ( hdf5::dataset_exists( mFile->active_group(), "physical" ) )
            {
                mFile->load_data( "physical", mProto->element_data()->mPhysicalTags );
            }
            else
            {
                mProto->element_data()->mPhysicalTags.set_size( tNumElems, 0 );
            }

            // populate block data
            index_t tCount = 0 ;
            Cell< uint > & tBlockIDs = mProto->element_data()->mGeometryTags ;
            tBlockIDs.set_size( tNumElems, 0 );

            for ( proto::GroupData & tData : mProto->block_data() )
            {
                for ( index_t k = 0 ; k<tData.mNumElements; ++k )
                {
                    tBlockIDs( tCount++ ) = tData.mID ;
                }
            }


            hdf5::Dataset< id_t > tData( mFile->active_group(), "topology" );

            tCount = 0 ;

            for ( proto::GroupData & tBlockData : mProto->block_data() )
            {
                tCount += tBlockData.mNumElements * number_of_nodes( tBlockData.mElementType );
            }

            Cell< id_t > & tTopology = mProto->element_data()->mTopology;
            tTopology.set_size( tCount );
            tCount = 0 ;

            index_t tElemCount = 0 ;
            for ( proto::GroupData & tBlockData : mProto->block_data() )
            {
                for ( index_t e=0; e<tBlockData.mNumElements; ++e )
                {
                    uint n = tData.length( tElemCount );
                    const id_t * tNodes = tData[ tElemCount++ ];

                    for ( uint k=0; k<n; ++k )
                    {
                        tTopology( tCount++ ) = tNodes[ k ];
                    }
                }
            }


            tData.close();

            mFile->close_active_group() ;

            mProto->create_elements();
            mProto->edge();

            mMesh->set_connectivity( Connectivity::ElementToNode );
#endif
        }

        void
        BfmFile::save_facet_data( const bool aSaveTopology )
        {
#ifdef BELFEM_HDF5

            mProto->populate_facet_data( false, false, false );

            mFile->create_group( "facets" );

            mFile->save_data( "ids", mProto->facet_data()->mIDs );

            if ( mProto->facet_data()->mPhysicalTags.size() > 0 )
            {
                mFile->save_data( "physical", mProto->facet_data()->mPhysicalTags );
            }

            hsize_t tNumFacets = mProto->facet_data()->mIDs.size() ;

            if ( aSaveTopology )
            {
                index_t tCount = 0 ;

                hdf5::Dataset< id_t > tData( mFile->active_group(), "topology", tNumFacets );
                for ( SideSet * tSideSet : mMesh->sidesets() )
                {
                    Cell< Facet * > & tFacets= tSideSet->facets() ;

                    hsize_t n = number_of_nodes( tSideSet->element_type() );
                    for ( Facet * tFacet : tFacets )
                    {
                        id_t * tNodes = tData.set_size( tCount++, n );

                        for ( hsize_t k=0; k<n; ++k )
                        {
                            tNodes[ k ] = tFacet->node( k )->id();
                        }
                    }
                }

                // write data
                tData.save();
            }

            hdf5::Dataset<id_t>  tElements( mFile->active_group(), "elements", tNumFacets );
            hdf5::Dataset<uchar> tIndices( mFile->active_group(), "indices", tNumFacets );

            index_t tCount = 0 ;

            for ( SideSet * tSideSet : mMesh->sidesets() )
            {
                Cell< Facet * > & tFacets= tSideSet->facets() ;

                for ( Facet * tFacet : tFacets )
                {
                    // every facet must have a master: the no-topology load path
                    // re-derives the facet nodes from the master element
                    BELFEM_ERROR( tFacet->has_master() || aSaveTopology,
                        "Facet %lu has no master element; cannot save without node topology",
                        ( long unsigned int ) tFacet->id() );

                    uchar tCase = 0 ;
                    tCase +=     tFacet->has_master() ;
                    tCase += 2 * tFacet->has_slave() ;

                    switch ( tCase )
                    {
                        case 1 : // only master
                        {

                            id_t * tIDs  = tElements.set_size( tCount, 1 );
                            tIDs[ 0 ]    = tFacet->master()->id() ;
                            uchar * tIDX = tIndices.set_size( tCount++, 1 );
                            tIDX[ 0 ] = tFacet->index_on_master() ;
                            break ;
                        }
                        case 2 : // only slave
                        {

                            id_t * tIDs  = tElements.set_size( tCount, 1 );
                            tIDs[ 0 ]    = tFacet->slave()->id() ;
                            uchar * tIDX = tIndices.set_size( tCount++, 2 );
                            tIDX[ 0 ] = tFacet->index_on_slave() ;
                            tIDX[ 1 ] = tFacet->orientation_on_slave();
                            break ;
                        }
                        case 3 : // master and slave
                        {
                            id_t * tIDs  = tElements.set_size( tCount, 2 );
                            tIDs[ 0 ]    = tFacet->master()->id() ;
                            tIDs[ 1 ]    = tFacet->slave()->id() ;

                            uchar * tIDX = tIndices.set_size( tCount++, 3 );
                            tIDX[ 0 ] = tFacet->index_on_master() ;
                            tIDX[ 1 ] = tFacet->index_on_slave() ;
                            tIDX[ 2 ] = tFacet->orientation_on_slave();
                            break ;
                        }
                        default:
                        {
                            BELFEM_ERROR(  false, "Facet %lu has neither master nor slave" ,
                                ( long unsigned int ) tFacet->id() );
                        }
                    }
                }
            }

            tElements.save();
            tIndices.save();

            mFile->close_active_group() ;

            mProto->reset_facet_data();
#endif
        }

        void
        BfmFile::load_facet_data()
        {
#ifdef BELFEM_HDF5

            if ( ! hdf5::group_exists( mFile->active_group(), "facets" ) ) return ;

            mFile->select_group( "facets" );
            mFile->load_data( "ids", mProto->facet_data()->mIDs );

            // populate sideset data
            hsize_t tNumFacets = mProto->facet_data()->mIDs.size();
            Cell< uint > & tSideSetIDs = mProto->facet_data()->mGeometryTags ;
            tSideSetIDs.set_size( tNumFacets, 0 );
            // facet element types are uniform per sideset and not stored on
            // file; re-derive them from the sideset element type
            Cell< uchar > & tTypes = mProto->facet_data()->mTypes ;
            tTypes.set_size( tNumFacets, 0 );

            index_t tCount = 0 ;
            for ( proto::GroupData & tData : mProto->sideset_data() )
            {
                for ( index_t k = 0 ; k<tData.mNumElements; ++k )
                {
                    tSideSetIDs( tCount )   = tData.mID ;
                    tTypes( tCount++ )      = static_cast< uchar >( tData.mElementType );
                }
            }

            if ( hdf5::dataset_exists( mFile->active_group(), "physical" ) )
            {
                mFile->load_data( "physical", mProto->facet_data()->mPhysicalTags );
            }

            if ( hdf5::dataset_exists( mFile->active_group(), "topology" ) )
            {
                                tCount = 0 ;

                hdf5::Dataset< id_t > tData( mFile->active_group(), "topology" );

                for ( proto::GroupData & tSideSetData : mProto->sideset_data() )
                {
                    tCount += tSideSetData.mNumElements * number_of_nodes( tSideSetData.mElementType );
                }

                Cell< id_t > & tTopology = mProto->facet_data()->mTopology;
                tTopology.set_size( tCount );
                tCount = 0 ;

                index_t tFacetCount = 0 ;
                for ( proto::GroupData & tSideSetData : mProto->sideset_data() )
                {
                    for ( index_t e=0; e<tSideSetData.mNumElements; ++e )
                    {
                        const id_t * tNodes = tData[ tFacetCount ];

                        uint n = tData.length( tFacetCount++ );
                        for ( uint k=0; k<n; ++k )
                        {
                            tTopology( tCount++ ) = tNodes[ k ];
                        }
                    }
                }

                tData.close();
            }


            if ( hdf5::dataset_exists( mFile->active_group(), "elements" ) &&
                 hdf5::dataset_exists( mFile->active_group(), "indices" ) )
            {
                hdf5::Dataset< id_t > tElements( mFile->active_group(), "elements" );
                hdf5::Dataset< uchar > tIndices( mFile->active_group(), "indices" );

                Cell< id_t > & tMasterIDs      = mProto->facet_data()->mMasterIDs ;
                tMasterIDs.set_size( tNumFacets, gNoID );

                Cell< uchar > & tMasterIndices = mProto->facet_data()->mIndicesOnMaster ;
                tMasterIndices.set_size( tNumFacets, BELFEM_UCHAR_MAX );

                Cell< id_t >  & tSlaveIDs      = mProto->facet_data()->mSlaveIDs ;
                tSlaveIDs.set_size( tNumFacets, gNoID );

                Cell< uchar > & tSlaveIndices = mProto->facet_data()->mIndicesOnSlave ;
                tSlaveIndices.set_size( tNumFacets, BELFEM_UCHAR_MAX );

                Cell< uchar > & tOrientations  = mProto->facet_data()->mOrientationsOnSlave ;
                tOrientations.set_size( tNumFacets, BELFEM_UCHAR_MAX );

                tCount = 0 ;

                for ( proto::GroupData & tData : mProto->sideset_data() )
                {
                    for ( index_t f=0; f<tData.mNumElements; ++f )
                    {
                        const id_t  * tIDs = tElements[ tCount ] ;
                        const uchar * tIDX = tIndices[ tCount ] ;

                        switch ( tIndices.length( tCount ) )
                        {
                            case 1 :
                            {
                                tMasterIDs( tCount ) = tIDs[ 0 ];
                                tMasterIndices( tCount ) = tIDX[ 0 ];
                                break ;
                            }
                            case 2 :
                            {
                                tSlaveIDs( tCount )     = tIDs[ 0 ];
                                tSlaveIndices( tCount ) = tIDX[ 0 ];
                                tOrientations( tCount ) = tIDX[ 1 ];
                                break ;
                            }
                            case 3 :
                            {
                                tMasterIDs( tCount ) = tIDs[ 0 ];
                                tMasterIndices( tCount ) = tIDX[ 0 ];
                                tSlaveIDs( tCount ) = tIDs[ 1 ];
                                tSlaveIndices( tCount ) = tIDX[ 1 ];
                                tOrientations( tCount ) = tIDX[ 2 ];
                                break ;
                            }
                            default:
                            {
                                BELFEM_ERROR( false, "invalid facet topology" );
                            }
                        }

                        ++tCount ;
                    }
                }

                tElements.close();
                tIndices.close();
            }

            mFile->close_active_group() ;

            mProto->create_facets();

            // keep zero-facet sidesets: thin-shell source sidesets are empty
            // by design, but read_domain_types looks them up by id on reload
            mProto->create_sidesets( true );

            mMesh->set_connectivity( Connectivity::FacetToNode );
            mMesh->set_connectivity( Connectivity::FacetToElement );
#endif
        }

        void
        BfmFile::save_edge_data()
        {
#ifdef BELFEM_HDF5

            if ( ! mMesh->edges_exist() ) return ;

            mFile->create_group( "edges" );
            mProto->populate_edge_data( false );

            mFile->save_data( "ids", mProto->edge_data()->mIDs );

            hsize_t tNumEdges = mProto->edge_data()->mIDs.size() ;

            index_t tCount = 0 ;

            hdf5::Dataset< id_t > tData( mFile->active_group(), "topology", tNumEdges );
            // NOTE: this branch re-packs an already-flattened proto topology into
            // the variable-length container. It is currently NOT in use, since
            // save_edge_data() always calls populate_edge_data( false ); it is
            // kept for discipline / symmetry with the live-mesh branch below.
            if ( mProto->edge_data()->mTopology.size() > 0 )
            {
                Cell< id_t > & tTopology = mProto->edge_data()->mTopology;

                for ( hsize_t e=0; e<tNumEdges; ++e )
                {
                    hsize_t n = tTopology( tCount++ );

                    id_t * tNodes = tData.set_size( e , n );

                    for ( hsize_t k=0; k<n; ++k )
                    {
                        tNodes[ k ] = tTopology( tCount++ );
                    }
                }
            }
            else
            {
                Cell< Edge * > & tEdges = mMesh->edges() ;

                for ( hsize_t e=0; e<tNumEdges; ++e )
                {
                    Edge * tEdge = tEdges( e );
                    hsize_t n = tEdge->number_of_nodes();

                    id_t * tNodes = tData.set_size( e , n );

                    for ( hsize_t k=0; k<n; ++k )
                    {
                        tNodes[ k ] = tEdge->node( k )->id();
                    }
                }
            }

            tData.save();

            mFile->close_active_group() ;

            mProto->reset_edge_data();
#endif
        }

        void
        BfmFile::load_edge_data()
        {
#ifdef BELFEM_HDF5
            if ( ! hdf5::group_exists( mFile->active_group(), "edges" ) ) return ;

            mFile->select_group( "edges" );

            mFile->load_data( "ids", mProto->edge_data()->mIDs );

            hsize_t tNumEdges = mProto->edge_data()->mIDs.size();

            hdf5::Dataset< id_t > tData( mFile->active_group(), "topology" );

            // count memory
            index_t tCount = tNumEdges ;

            for ( hsize_t e=0; e<tNumEdges; ++e )
            {
                tCount += tData.length( e );
            }

            Cell< id_t > & tTopology = mProto->edge_data()->mTopology;
            tTopology.set_size( tCount );

            tCount = 0 ;

            for ( hsize_t e=0; e<tNumEdges; ++e )
            {
                uint n = tData.length( e );
                tTopology( tCount++ ) = n ;

                id_t * tNodes = tData[ e ];

                for ( uint k=0; k<n; ++k )
                {
                    tTopology( tCount++ ) = tNodes[ k ];
                }
            }

            tData.close();

            mFile->close_active_group() ;

            mProto->create_edges();
            mProto->reconstruct_edge_connectivity();
#endif
        }

        void
        BfmFile::save_face_data()
        {
#ifdef BELFEM_HDF5

            if ( ! mMesh->faces_exist() ) return ;

            mProto->populate_face_data( false );

            mFile->create_group( "faces" );


            mFile->save_data( "ids", mProto->face_data()->mIDs );

            Cell< Face * > & tFaces = mMesh->faces();

            hsize_t tNumFaces = tFaces.size();

            index_t tCount = 0 ;

            hdf5::Dataset<id_t> tElements( mFile->active_group(), "elements", tNumFaces );
            hdf5::Dataset<uchar> tIndices( mFile->active_group(), "indices", tNumFaces );

            for ( Face * tFace : tFaces )
            {
                // the case of a face without master is not forseen in the code
                // at this time. We keep the slave-olny case for symmetry reasons with the Facet storing scheme
                BELFEM_ERROR( tFace->master() != nullptr,
                    "Face %lu has no master element; cannot save",
                    ( long unsigned int ) tFace->id() );

                uchar tCase = 0 ;
                tCase +=       tFace->master() != nullptr ;
                tCase += 2 * ( tFace->slave()  != nullptr );

                switch ( tCase )
                {
                    case 1 : // only master
                    {

                        id_t * tIDs  = tElements.set_size( tCount, 1 );
                        tIDs[ 0 ]    = tFace->master()->id() ;

                        uchar * tIDX = tIndices.set_size( tCount++, 1 );
                        tIDX[ 0 ] = tFace->index_on_master() ;
                        break ;
                    }
                    case 2 : // only slave (should not happen, written for completeness)
                    {
                        id_t * tIDs  = tElements.set_size( tCount, 1 );
                        tIDs[ 0 ]    = tFace->slave()->id() ;

                        uchar * tIDX = tIndices.set_size( tCount++, 2 );
                        tIDX[ 0 ] = tFace->index_on_slave() ;
                        tIDX[ 1 ] = tFace->orientation_on_slave();

                        break ;
                    }
                    case 3 : // master and slave
                    {
                        id_t * tIDs  = tElements.set_size( tCount, 2 );
                        tIDs[ 0 ]    = tFace->master()->id() ;
                        tIDs[ 1 ]    = tFace->slave()->id() ;

                        uchar * tIDX = tIndices.set_size( tCount++, 3 );
                        tIDX[ 0 ] = tFace->index_on_master() ;
                        tIDX[ 1 ] = tFace->index_on_slave() ;
                        tIDX[ 2 ] = tFace->orientation_on_slave();

                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR(  false, "Face %lu has neither master nor slave" ,
                            ( long unsigned int ) tFace->id() );
                    }
                }
            }

            tElements.save();
            tIndices.save();

            mFile->close_active_group() ;

            mProto->reset_face_data();
#endif
        }

        void
        BfmFile::load_face_data()
        {
#ifdef BELFEM_HDF5

            if ( ! hdf5::group_exists( mFile->active_group(), "faces" ) ) return ;

            mFile->select_group( "faces" );
            mFile->load_data( "ids", mProto->face_data()->mIDs );
            hsize_t tNumFaces = mProto->face_data()->mIDs.size();

            hdf5::Dataset<id_t> tElements( mFile->active_group(), "elements" );
            hdf5::Dataset<uchar> tIndices( mFile->active_group(), "indices" );

            Cell< id_t > & tMasterIDs      = mProto->face_data()->mMasterIDs ;
            tMasterIDs.set_size( tNumFaces, gNoID );

            Cell< uchar > & tMasterIndices = mProto->face_data()->mIndicesOnMaster ;
            tMasterIndices.set_size( tNumFaces, BELFEM_UCHAR_MAX );

            Cell< id_t >  & tSlaveIDs      = mProto->face_data()->mSlaveIDs ;
            tSlaveIDs.set_size( tNumFaces, gNoID );

            Cell< uchar > & tSlaveIndices = mProto->face_data()->mIndicesOnSlave ;
            tSlaveIndices.set_size( tNumFaces, BELFEM_UCHAR_MAX );

            Cell< uchar > & tOrientations  = mProto->face_data()->mOrientationsOnSlave ;
            tOrientations.set_size( tNumFaces, BELFEM_UCHAR_MAX );

            index_t tCount = 0 ;

            for ( size_t f=0; f<tNumFaces; ++f )
            {
                id_t  * tIDs = tElements[ tCount ];
                uchar * tIDX = tIndices[ tCount ];

                switch ( tIndices.length( tCount ) )
                {
                    case 1 :
                    {
                        tMasterIDs( tCount ) = tIDs[ 0 ];
                        tMasterIndices( tCount ) = tIDX[ 0 ];
                        break ;
                    }
                    case 2 :
                    {
                        tSlaveIDs( tCount )     = tIDs[ 0 ];
                        tSlaveIndices( tCount ) = tIDX[ 0 ];
                        tOrientations( tCount ) = tIDX[ 1 ];
                        break ;
                    }
                    case 3 :
                    {
                        tMasterIDs( tCount ) = tIDs[ 0 ];
                        tMasterIndices( tCount ) = tIDX[ 0 ];
                        tSlaveIDs( tCount ) = tIDs[ 1 ];
                        tSlaveIndices( tCount ) = tIDX[ 1 ];
                        tOrientations( tCount ) = tIDX[ 2 ];
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "invalid face topology" );
                    }
                }

                ++tCount ;
            }

            tElements.close();
            tIndices.close();

            mFile->close_active_group() ;

            mProto->create_faces();
            mProto->reconstruct_face_connectivity();
#endif
        }

        void
        BfmFile::save_control_point_data()
        {
#ifdef BELFEM_HDF5
            if ( mMesh->number_of_control_points() == 0 ) return ;

            mProto->populate_control_point_data( false  );

            mFile->create_group( "control_points" );
            mFile->save_data( "ids", mProto->control_point_data()->mIDs );
            mFile->save_data( "coords", mProto->control_point_data()->mCoords , true );

            hsize_t tNumElements = mMesh->number_of_elements();

            hdf5::Dataset<id_t> tTopo( mFile->active_group(), "topology", tNumElements );

            index_t tCount = 0 ;

            for ( Block * tBlock : mMesh->blocks() )
            {
                Cell< Element * > & tElements = tBlock->elements();
                for ( Element * tElement : tElements )
                {
                    uint n = tElement->number_of_control_points();
                    if ( n > 0 )
                    {
                        id_t * tIDs = tTopo.set_size( tCount,  n );
                        for ( uint k=0; k<n; ++k )
                        {
                            tIDs[ k ] = tElement->control_point( k )->id();
                        }
                    }
                    ++tCount ;
                }
            }

            tTopo.save();
            mFile->close_active_group();
            mProto->reset_control_point_data();
#endif
        }

        void
        BfmFile::load_control_point_data()
        {
#ifdef BELFEM_HDF5

            if ( ! hdf5::group_exists( mFile->active_group(), "control_points" ) ) return;

            mFile->select_group( "control_points" );
            mFile->load_data( "ids", mProto->control_point_data()->mIDs );
            mFile->load_data( "coords", mProto->control_point_data()->mCoords , true );

            // create the control point entities ( + mControlPointMap ) from the
            // loaded ids/coords. The element->control_point incidence is restored
            // below from the vlen "topology" dataset, so create_control_points'
            // own ( mElementTopology ) incidence restore no-ops here.
            mProto->create_control_points();

            hdf5::Dataset<id_t> tTopo( mFile->active_group(), "topology" );

            index_t tCount = 0 ;
            for ( proto::GroupData & tBlockData : mProto->block_data() )
            {
                Block * tBlock = mProto->block( tBlockData.mID );

                Cell< Element * > & tElements = tBlock->elements() ;

                for ( Element * tElement : tElements )
                {
                    uint n = tTopo.length( tCount );
                    if ( n == 0 ) { ++tCount; continue; };

                    const id_t * tIDs = tTopo[ tCount++ ];

                    tElement->allocate_control_points_container( n );
                    for ( uint k=0; k<n; ++k )
                    {
                        tElement->insert_control_point( mProto->control_point( tIDs[ k ] ) , k );
                    }
                }
            }
            tTopo.close();

            mFile->close_active_group();
#endif
        }

        void
        BfmFile::save_node_duplicate_data()
        {
#ifdef BELFEM_HDF5

            hsize_t tCount = 0 ;
            Cell< Node * > & tNodes = mMesh->nodes();
            for ( Node * tNode : tNodes )
            {
                if ( tNode->is_duplicate() ) continue;
                if ( tNode->number_of_duplicates() > 0 ) ++tCount ;
            }

            if ( tCount == 0 ) return ;

            mFile->select_group( "nodes" );

            hdf5::Dataset< id_t > tData( mFile->active_group(), "duplicates", tCount );

            // reset: tCount now becomes the row index for the write loop
            tCount = 0 ;

            for ( Node * tNode : tNodes )
            {
                if ( tNode->is_duplicate() ) continue;
                hsize_t n = tNode->number_of_duplicates();

                if ( n == 0 ) continue;


                id_t * tIDs  = tData.set_size( tCount, n + 1 );

                tIDs[ 0 ] = tNode->id();

                for ( hsize_t d=0; d<n; ++d )
                {
                    tIDs[ d + 1 ] = tNode->duplicate( d )->id();
                }


                ++tCount;
            }

            tData.save();

            mFile->close_active_group();

#endif
        }

        void
        BfmFile::load_node_duplicate_data()
        {
#ifdef BELFEM_HDF5
            mFile->select_group( "nodes" );

            if ( ! hdf5::dataset_exists( mFile->active_group(), "duplicates" ) )
            {
                mFile->close_active_group();
                return;
            }

            hdf5::Dataset< id_t > tData( mFile->active_group(), "duplicates" );

            hsize_t tNumNodes = tData.size();

            for ( hsize_t k = 0; k<tNumNodes; ++k )
            {
                BELFEM_ERROR( tData.length( k ) >= 1,
                    "Malformed duplicate row %lu in node duplicate data",
                    ( long unsigned int ) k );

                const id_t * tIDs = tData[ k ];

                Node * tOrg = mProto->node( tIDs[ 0 ] );

                uint n = tData.length( k ) - 1 ;
                tOrg->allocate_duplicate_container( n );

                for ( uint d=0; d<n; ++d )
                {
                    Node * tDup = mProto->node( tIDs[ d + 1 ] );
                    tOrg->add_duplicate( tDup );
                    tDup->set_original( tOrg );
                }
            }

            tData.close();
            mFile->close_active_group();
#endif
        }

        void
        BfmFile::save_hanging_entities()
        {
#ifdef BELFEM_HDF5

            Cell< Node * > & tNodes = mMesh->nodes();
            Cell< Edge * > & tEdges = mMesh->edges();
            Cell< Face * > & tFaces = mMesh->faces();
            Cell< Facet * > & tFacets = mMesh->facets();
            Cell< ControlPoint * > & tControlPoints = mMesh->control_points();

            index_t tNumNodes = 0 ;
            for ( Node * tNode : tNodes )
            {
                if ( tNode->is_hanging() ) ++ tNumNodes ;
            }

            index_t tNumEdges = 0 ;
            for ( Edge * tEdge : tEdges )
            {
                if ( tEdge->is_hanging() ) ++ tNumEdges ;
            }

            index_t tNumFaces = 0 ;

            for ( Face * tFace : tFaces )
            {
                if ( tFace->is_hanging() ) ++ tNumFaces ;
            }

            index_t tNumFacets = 0 ;
            for ( Facet * tFacet : tFacets )
            {
                if ( tFacet->is_hanging() ) ++ tNumFacets;
            }

            index_t tNumControlPoints = 0 ;
            for ( ControlPoint * tPoint : tControlPoints )
            {
                if ( tPoint->is_hanging() ) ++ tNumControlPoints;
            }

            if ( tNumNodes == 0 && tNumEdges == 0 && tNumFaces == 0 && tNumFacets == 0 && tNumControlPoints == 0 ) return;

            mFile->create_group( "hanging" );

            if ( tNumNodes > 0 ) this->save_hanging_nodes( tNumNodes );
            if ( tNumEdges > 0 ) this->save_hanging_edges( tNumEdges );
            if ( tNumFaces > 0 ) this->save_hanging_faces( tNumFaces );
            if ( tNumFacets > 0 ) this->save_hanging_facets( tNumFacets );
            if ( tNumControlPoints > 0 ) this->save_hanging_control_points( tNumControlPoints );

            mFile->close_active_group();
#endif
        }

        void
        BfmFile::load_hanging_entities()
        {
#ifdef BELFEM_HDF5
            if ( ! hdf5::group_exists( mFile->active_group(), "hanging" ) ) return;

            mFile->select_group( "hanging" );
            if ( hdf5::group_exists( mFile->active_group(), "nodes" ) ) this->load_hanging_nodes() ;
            if ( hdf5::group_exists( mFile->active_group(), "edges" ) ) this->load_hanging_edges() ;
            if ( hdf5::group_exists( mFile->active_group(), "faces" ) ) this->load_hanging_faces() ;
            if ( hdf5::group_exists( mFile->active_group(), "facets" ) ) this->load_hanging_facets() ;
            if ( hdf5::group_exists( mFile->active_group(), "control_points" ) ) this->load_hanging_control_points() ;

            mFile->close_active_group();
#endif
        }

        void
        BfmFile::save_hanging_nodes( const hsize_t aNumNodes )
        {
#ifdef BELFEM_HDF5
            if ( aNumNodes == 0 ) return;

            Cell< id_t > tIDs( aNumNodes, gNoID );

            Cell< Node * > & tNodes = mMesh->nodes();

            mFile->create_group( "nodes" );

            index_t tCount = 0 ;
            hid_t tGroup = mFile->active_group();
            hdf5::Dataset< id_t >  tTopo( tGroup, "topology",  aNumNodes );
            hdf5::Dataset< suint > tTypes( tGroup, "types", aNumNodes );
            hdf5::Dataset< real >  tWeights( tGroup, "weights", aNumNodes );

            for ( Node * tNode : tNodes )
            {
                if ( ! tNode->is_hanging() ) continue ;



                uint n = tNode->number_of_sources();
                id_t   * tSrc  = tTopo.set_size( tCount, n );
                suint  * tType = tTypes.set_size( tCount, n );
                double * tWeight = tWeights.set_size( tCount, n );

                tIDs( tCount++ ) = tNode->id();

                for ( uint s=0; s<n; ++s )
                {
                    tSrc[ s ]    = tNode->source( s )->id();
                    tType[ s ]   = static_cast< suint >( tNode->source( s )->entity_type() );
                    tWeight[ s ] = tNode->weight( s );
                }
            }

            mFile->save_data( "ids", tIDs );
            tTopo.save();
            tTypes.save();
            tWeights.save();
            mFile->close_active_group();

#endif
        }
        void
        BfmFile::load_hanging_nodes()
        {
#ifdef BELFEM_HDF5
            if ( ! hdf5::group_exists( mFile->active_group(), "nodes" ) ) return;

            mFile->select_group( "nodes" );

            // capture the group handle AFTER selecting "nodes" so the datasets
            // are opened in hanging/nodes, not the parent
            hid_t tGroup = mFile->active_group();

            Cell< id_t > tIDs ;

            mFile->load_data( "ids", tIDs );
            hdf5::Dataset< id_t >  tTopo( tGroup, "topology" );
            hdf5::Dataset< suint > tTypes( tGroup, "types" );
            hdf5::Dataset< real >  tWeights( tGroup, "weights" );

            hsize_t tNumNodes = tIDs.size();

            for ( hsize_t k=0; k<tNumNodes; ++k )
            {
                Node * tNode = mProto->node( tIDs( k ) );

                uint n = tTopo.length( k );
                tNode->allocate_source_container( n );

                const id_t   * tID     = tTopo[ k ];
                const suint  * tType   = tTypes[ k ];
                const real   * tWeight = tWeights[ k ];
                for ( uint s=0; s<n; ++s )
                {
                    tNode->add_source( mProto->basis( static_cast< EntityType >( tType[ s ] ), tID[ s ] ), tWeight[ s ] );
                }
            }

            tTopo.close();
            tTypes.close();
            tWeights.close();
            mFile->close_active_group();
#endif
        }

        void
        BfmFile::save_hanging_edges( const hsize_t aNumEdges )
        {
#ifdef BELFEM_HDF5
            if ( aNumEdges == 0 ) return;

            Cell< id_t > tIDs( aNumEdges, gNoID );

            Cell< Edge * > & tEdges = mMesh->edges();

            mFile->create_group( "edges" );

            index_t tCount = 0 ;
            hid_t tGroup = mFile->active_group();
            hdf5::Dataset< id_t >  tTopo( tGroup, "topology",  aNumEdges );
            hdf5::Dataset< suint > tTypes( tGroup, "types", aNumEdges );
            hdf5::Dataset< real >  tWeights( tGroup, "weights", aNumEdges );

            for ( Edge * tEdge : tEdges )
            {
                if ( ! tEdge->is_hanging() ) continue ;



                uint n = tEdge->number_of_sources();
                id_t   * tSrc  = tTopo.set_size( tCount, n );
                suint  * tType = tTypes.set_size( tCount, n );
                double * tWeight = tWeights.set_size( tCount, n );

                tIDs( tCount++ ) = tEdge->id();

                for ( uint s=0; s<n; ++s )
                {
                    tSrc[ s ]    = tEdge->source( s )->id();
                    tType[ s ]   = static_cast< suint >( tEdge->source( s )->entity_type() );
                    tWeight[ s ] = tEdge->weight( s );
                }
            }

            mFile->save_data( "ids", tIDs );
            tTopo.save();
            tTypes.save();
            tWeights.save();
            mFile->close_active_group();

#endif
        }

        void
        BfmFile::load_hanging_edges()
        {
#ifdef BELFEM_HDF5
            if ( ! hdf5::group_exists( mFile->active_group(), "edges" ) ) return;

            mFile->select_group( "edges" );

            // capture the group handle AFTER selecting "edges" so the datasets
            // are opened in hanging/edges, not the parent
            hid_t tGroup = mFile->active_group();

            Cell< id_t > tIDs ;

            mFile->load_data( "ids", tIDs );
            hdf5::Dataset< id_t >  tTopo( tGroup, "topology" );
            hdf5::Dataset< suint > tTypes( tGroup, "types" );
            hdf5::Dataset< real >  tWeights( tGroup, "weights" );

            hsize_t tNumEdges = tIDs.size();

            for ( hsize_t k=0; k<tNumEdges; ++k )
            {
                Edge * tEdge = mProto->edge( tIDs( k ) );

                uint n = tTopo.length( k );
                tEdge->allocate_source_container( n );

                const id_t   * tID     = tTopo[ k ];
                const suint  * tType   = tTypes[ k ];
                const real   * tWeight = tWeights[ k ];
                for ( uint s=0; s<n; ++s )
                {
                    tEdge->add_source( mProto->basis( static_cast< EntityType >( tType[ s ] ), tID[ s ] ), tWeight[ s ] );
                }
            }

            tTopo.close();
            tTypes.close();
            tWeights.close();
            mFile->close_active_group();
#endif
        }

        void
        BfmFile::save_hanging_faces( const hsize_t aNumFaces )
        {
#ifdef BELFEM_HDF5
            if ( aNumFaces == 0 ) return;

            Cell< id_t > tIDs( aNumFaces, gNoID );

            Cell< Face * > & tFaces = mMesh->faces();

            mFile->create_group( "faces" );

            index_t tCount = 0 ;
            hid_t tGroup = mFile->active_group();
            hdf5::Dataset< id_t >  tTopo( tGroup, "topology",  aNumFaces );
            hdf5::Dataset< suint > tTypes( tGroup, "types", aNumFaces );
            hdf5::Dataset< real >  tWeights( tGroup, "weights", aNumFaces );

            for ( Face * tFace : tFaces )
            {
                if ( ! tFace->is_hanging() ) continue ;

                uint n = tFace->number_of_sources();
                id_t   * tSrc  = tTopo.set_size( tCount, n );
                suint  * tType = tTypes.set_size( tCount, n );
                double * tWeight = tWeights.set_size( tCount, n );

                tIDs( tCount++ ) = tFace->id();

                for ( uint s=0; s<n; ++s )
                {
                    tSrc[ s ]    = tFace->source( s )->id();
                    tType[ s ]   = static_cast< suint >( tFace->source( s )->entity_type() );
                    tWeight[ s ] = tFace->weight( s );
                }
            }

            mFile->save_data( "ids", tIDs );
            tTopo.save();
            tTypes.save();
            tWeights.save();
            mFile->close_active_group();

#endif
        }
        void
        BfmFile::load_hanging_faces()
        {
#ifdef BELFEM_HDF5
            if ( ! hdf5::group_exists( mFile->active_group(), "faces" ) ) return;

            mFile->select_group( "faces" );

            // capture the group handle AFTER selecting "faces" so the datasets
            // are opened in hanging/faces, not the parent
            hid_t tGroup = mFile->active_group();

            Cell< id_t > tIDs ;

            mFile->load_data( "ids", tIDs );
            hdf5::Dataset< id_t >  tTopo( tGroup, "topology" );
            hdf5::Dataset< suint > tTypes( tGroup, "types" );
            hdf5::Dataset< real >  tWeights( tGroup, "weights" );

            hsize_t tNumFaces = tIDs.size();

            for ( hsize_t k=0; k<tNumFaces; ++k )
            {
                Face * tFace = mProto->face( tIDs( k ) );

                uint n = tTopo.length( k );
                tFace->allocate_source_container( n );

                const id_t   * tID     = tTopo[ k ];
                const suint  * tType   = tTypes[ k ];
                const real   * tWeight = tWeights[ k ];
                for ( uint s=0; s<n; ++s )
                {
                    tFace->add_source( mProto->basis( static_cast< EntityType >( tType[ s ] ), tID[ s ] ), tWeight[ s ] );
                }
            }

            tTopo.close();
            tTypes.close();
            tWeights.close();
            mFile->close_active_group();
#endif
        }

        void
        BfmFile::save_hanging_facets( const hsize_t aNumFacets )
        {
#ifdef BELFEM_HDF5
            if ( aNumFacets == 0 ) return;

            Cell< id_t > tIDs( aNumFacets, gNoID );

            Cell< Facet * > & tFacets = mMesh->facets();

            mFile->create_group( "facets" );

            index_t tCount = 0 ;
            hid_t tGroup = mFile->active_group();
            hdf5::Dataset< id_t >  tTopo( tGroup, "topology",  aNumFacets );
            hdf5::Dataset< suint > tTypes( tGroup, "types", aNumFacets );
            hdf5::Dataset< real >  tWeights( tGroup, "weights", aNumFacets );

            for ( Facet * tFacet : tFacets )
            {
                if ( ! tFacet->is_hanging() ) continue ;



                uint n = tFacet->number_of_sources();
                id_t   * tSrc  = tTopo.set_size( tCount, n );
                suint  * tType = tTypes.set_size( tCount, n );
                double * tWeight = tWeights.set_size( tCount, n );

                tIDs( tCount++ ) = tFacet->id();

                for ( uint s=0; s<n; ++s )
                {
                    tSrc[ s ]    = tFacet->source( s )->id();
                    tType[ s ]   = static_cast< suint >( tFacet->source( s )->entity_type() );
                    tWeight[ s ] = tFacet->weight( s );
                }
            }

            mFile->save_data( "ids", tIDs );
            tTopo.save();
            tTypes.save();
            tWeights.save();
            mFile->close_active_group();

#endif
        }
        void
        BfmFile::load_hanging_facets()
        {
#ifdef BELFEM_HDF5
            if ( ! hdf5::group_exists( mFile->active_group(), "facets" ) ) return;

            mFile->select_group( "facets" );

            // capture the group handle AFTER selecting "facets" so the datasets
            // are opened in hanging/facets, not the parent
            hid_t tGroup = mFile->active_group();

            Cell< id_t > tIDs ;

            mFile->load_data( "ids", tIDs );
            hdf5::Dataset< id_t >  tTopo( tGroup, "topology" );
            hdf5::Dataset< suint > tTypes( tGroup, "types" );
            hdf5::Dataset< real >  tWeights( tGroup, "weights" );

            hsize_t tNumFacets = tIDs.size();

            for ( hsize_t k=0; k<tNumFacets; ++k )
            {
                Facet * tFacet = mProto->facet( tIDs( k ) );

                uint n = tTopo.length( k );
                tFacet->allocate_source_container( n );

                const id_t   * tID     = tTopo[ k ];
                const suint  * tType   = tTypes[ k ];
                const real   * tWeight = tWeights[ k ];
                for ( uint s=0; s<n; ++s )
                {
                    tFacet->add_source( mProto->basis( static_cast< EntityType >( tType[ s ] ), tID[ s ] ), tWeight[ s ] );
                }
            }

            tTopo.close();
            tTypes.close();
            tWeights.close();
            mFile->close_active_group();
#endif
        }

        void
        BfmFile::save_hanging_control_points( const hsize_t aNumControlPoints )
        {
#ifdef BELFEM_HDF5
            if ( aNumControlPoints == 0 ) return;

            Cell< id_t > tIDs( aNumControlPoints );

            Cell< ControlPoint * > & tControlPoints = mMesh->control_points();

            mFile->create_group( "control_points" );

            index_t tCount = 0 ;
            hid_t tGroup = mFile->active_group();
            hdf5::Dataset< id_t >  tTopo( tGroup, "topology",  aNumControlPoints );
            hdf5::Dataset< suint > tTypes( tGroup, "types", aNumControlPoints );
            hdf5::Dataset< real >  tWeights( tGroup, "weights", aNumControlPoints );

            for ( ControlPoint * tControlPoint : tControlPoints )
            {
                if ( ! tControlPoint->is_hanging() ) continue ;



                uint n = tControlPoint->number_of_sources();
                id_t   * tSrc  = tTopo.set_size( tCount, n );
                suint  * tType = tTypes.set_size( tCount, n );
                double * tWeight = tWeights.set_size( tCount, n );

                tIDs( tCount++ ) = tControlPoint->id();

                for ( uint s=0; s<n; ++s )
                {
                    tSrc[ s ]    = tControlPoint->source( s )->id();
                    tType[ s ]   = static_cast< suint >( tControlPoint->source( s )->entity_type() );
                    tWeight[ s ] = tControlPoint->weight( s );
                }
            }

            mFile->save_data( "ids", tIDs );
            tTopo.save();
            tTypes.save();
            tWeights.save();
            mFile->close_active_group();

#endif
        }
        void
        BfmFile::load_hanging_control_points()
        {
#ifdef BELFEM_HDF5
            if ( ! hdf5::group_exists( mFile->active_group(), "control_points" ) ) return;

            mFile->select_group( "control_points" );

            // capture the group handle AFTER selecting "control_points" so the datasets
            // are opened in hanging/control_points, not the parent
            hid_t tGroup = mFile->active_group();

            Cell< id_t > tIDs ;

            mFile->load_data( "ids", tIDs );
            hdf5::Dataset< id_t >  tTopo( tGroup, "topology" );
            hdf5::Dataset< suint > tTypes( tGroup, "types" );
            hdf5::Dataset< real >  tWeights( tGroup, "weights" );

            hsize_t tNumControlPoints = tIDs.size();

            for ( hsize_t k=0; k<tNumControlPoints; ++k )
            {
                ControlPoint * tControlPoint = mProto->control_point( tIDs( k ) );

                uint n = tTopo.length( k );
                tControlPoint->allocate_source_container( n );

                const id_t   * tID     = tTopo[ k ];
                const suint  * tType   = tTypes[ k ];
                const real   * tWeight = tWeights[ k ];
                for ( uint s=0; s<n; ++s )
                {
                    tControlPoint->add_source( mProto->basis( static_cast< EntityType >( tType[ s ] ), tID[ s ] ), tWeight[ s ] );
                }
            }

            tTopo.close();
            tTypes.close();
            tWeights.close();
            mFile->close_active_group();
#endif
        }

        void
        BfmFile::save_periodicity_data()
        {
#ifdef BELFEM_HDF5
            if ( ! mMesh->has_periodicity() ) return ;

            mFile->create_group( "periodic" );

            this->save_periodic_planes();

            this->save_periodic_nodes();
            this->save_periodic_edges();
            this->save_periodic_faces();
            this->save_periodic_facets();

            mFile->close_active_group();
#endif
        }

        void
        BfmFile::load_periodicity_data()
        {
#ifdef BELFEM_HDF5
            if ( ! hdf5::group_exists( mFile->active_group(), "periodic" ) ) return ;

            mFile->select_group( "periodic" );

            this->load_periodic_planes();
            this->load_periodic_nodes();
            this->load_periodic_edges();
            this->load_periodic_faces();
            this->load_periodic_facets();

            mFile->close_active_group();

            mProto->create_periodicitiy( true );
#endif
        }
        void
        BfmFile::save_periodic_planes()
        {
#ifdef BELFEM_HDF5
            Matrix< id_t > tPeriodicPlane( 2, 3 );

            Cell< Node * > & tMaster = mMesh->periodicity()->master_plane();
            Cell< Node * > & tSlave  = mMesh->periodicity()->slave_plane();
            for ( uint k=0; k<3; ++k )
            {
                tPeriodicPlane( 0, k ) = tMaster( k )->id();
                tPeriodicPlane( 1, k ) = tSlave( k )->id();
            }

            mFile->save_data( "planes", tPeriodicPlane );
#endif
        }


        void
        BfmFile::load_periodic_planes()
        {
#ifdef BELFEM_HDF5

            Matrix< id_t > tPeriodicPlane( 2, 3 );
            mFile->load_data( "planes", tPeriodicPlane );


            Cell< id_t > & tMaster = mProto->periodicity_data()->mMasterPlane ;
            Cell< id_t > & tSlave  = mProto->periodicity_data()->mSlavePlane ;

            tMaster.set_size( 3 );
            tSlave.set_size( 3 );
            for ( uint k=0; k<3; ++k )
            {
                tMaster( k ) = tPeriodicPlane( 0, k );
                tSlave( k )  = tPeriodicPlane( 1, k );
            }
#endif
        }

        void
        BfmFile::save_periodic_nodes()
        {
#ifdef BELFEM_HDF5

            if ( mMesh->periodicity()->master_nodes().size() == 0 ) return ;

            Cell< Node * > & tNodes = mMesh->periodicity()->master_nodes();
            Matrix< id_t > tData( 2, tNodes.size() );

            index_t tCount = 0 ;
            for ( Node * tNode : tNodes )
            {
                tData( 0, tCount ) = tNode->id();
                tData( 1, tCount ) = tNode->periodic()->id();
                ++tCount ;
            }
            mFile->save_data( "nodes", tData, true );
#endif
        }

        void
        BfmFile::load_periodic_nodes()
        {
#ifdef BELFEM_HDF5
            if ( ! hdf5::group_exists( mFile->active_group(), "nodes" ) ) return ;

            Matrix< id_t > tData ;
            mFile->load_data( "nodes", tData, true );

            Cell< id_t > & tMaster = mProto->periodicity_data()->mMasterNodes ;
            Cell< id_t > & tSlave  = mProto->periodicity_data()->mSlaveNodes ;

            index_t tNumNodes = tData.n_cols();

            tMaster.set_size( tNumNodes );
            tSlave.set_size( tNumNodes );
            for ( index_t k=0; k<tNumNodes; ++k )
            {
                tMaster( k ) = tData( 0, k );
                tSlave( k )  = tData( 1, k );
            }
#endif
        }

        void
        BfmFile::save_periodic_edges()
        {
#ifdef BELFEM_HDF5

            if ( mMesh->periodicity()->master_edges().size() == 0 ) return ;

            Cell< Edge * > & tEdges = mMesh->periodicity()->master_edges();


            Matrix< id_t > tData( 2, tEdges.size() );

            index_t tCount = 0 ;
            for ( Edge * tEdge : tEdges )
            {
                tData( 0, tCount ) = tEdge->id();
                tData( 1, tCount ) = tEdge->periodic()->id();
                ++tCount ;
            }

            mFile->save_data( "edges", tData, true );

#endif
        }
        void
        BfmFile::load_periodic_edges()
        {
#ifdef BELFEM_HDF5
            if ( ! hdf5::group_exists( mFile->active_group(), "edges" ) ) return ;

            Matrix< id_t > tData ;
            mFile->load_data( "edges", tData, true );

            Cell< id_t > & tMaster = mProto->periodicity_data()->mMasterEdges ;
            Cell< id_t > & tSlave  = mProto->periodicity_data()->mSlaveEdges ;

            index_t tNumEdges = tData.n_cols();

            tMaster.set_size( tNumEdges );
            tSlave.set_size( tNumEdges );
            for ( index_t k=0; k<tNumEdges; ++k )
            {
                tMaster( k ) = tData( 0, k );
                tSlave( k )  = tData( 1, k );
            }
#endif
        }

        void
        BfmFile::save_periodic_faces()
        {
#ifdef BELFEM_HDF5

            if ( mMesh->periodicity()->master_faces().size() == 0 ) return ;

            Cell< Face * > & tFaces = mMesh->periodicity()->master_faces();

            Matrix< id_t > tData( 2, tFaces.size() );
            index_t tCount = 0 ;
            for ( Face * tFace : tFaces )
            {
                tData( 0, tCount ) = tFace->id();
                tData( 1, tCount ) = tFace->periodic()->id();
                ++tCount ;
            }

            mFile->save_data( "faces", tData, true );
#endif
        }

        void
        BfmFile::load_periodic_faces()
        {
#ifdef BELFEM_HDF5

            if ( ! hdf5::group_exists( mFile->active_group(), "faces" ) ) return ;

            Matrix< id_t > tData ;
            mFile->load_data( "faces", tData, true );

            Cell< id_t > & tMaster = mProto->periodicity_data()->mMasterFaces ;
            Cell< id_t > & tSlave  = mProto->periodicity_data()->mSlaveFaces ;

            index_t tNumFaces = tData.n_cols();

            tMaster.set_size( tNumFaces );
            tSlave.set_size( tNumFaces );
            for ( index_t k=0; k<tNumFaces; ++k )
            {
                tMaster( k ) = tData( 0, k );
                tSlave( k )  = tData( 1, k );
            }
# endif
        }

        void
        BfmFile::save_periodic_facets()
        {
#ifdef BELFEM_HDF5
            if ( mMesh->periodicity()->master_facets().size() == 0 ) return ;

            Cell< Facet * > & tFacets = mMesh->periodicity()->master_facets();

            Matrix< id_t > tData( 2, tFacets.size() );

            index_t tCount = 0 ;
            for ( Facet * tFacet : tFacets )
            {
                tData( 0, tCount ) = tFacet->id();
                tData( 1, tCount ) = tFacet->periodic()->id();
                ++tCount ;
            }

            mFile->save_data( "facets", tData, true );
#endif
        }

        void
        BfmFile::load_periodic_facets()
        {
#ifdef BELFEM_HDF5

            if ( ! hdf5::group_exists( mFile->active_group(), "facets" ) ) return ;

            Matrix< id_t > tData ;
            mFile->load_data( "facets", tData, true );

            Cell< id_t > & tMaster = mProto->periodicity_data()->mMasterFacets ;
            Cell< id_t > & tSlave  = mProto->periodicity_data()->mSlaveFacets ;

            index_t tNumFacets = tData.n_cols();

            tMaster.set_size( tNumFacets );
            tSlave.set_size( tNumFacets );
            for ( index_t k=0; k<tNumFacets; ++k )
            {
                tMaster( k ) = tData( 0, k );
                tSlave( k )  = tData( 1, k );
            }
#endif
        }

        void
        BfmFile::save_thinshell_data()
        {
#ifdef BELFEM_HDF5

            Cell< ThinShell * > & tShells = mMesh->thin_shells();

            uint n = tShells.size();

            if ( n == 0 ) return ;

            mFile->create_group( "thinshells" );

            Vector< id_t > tShellIDs( n, 0 );
            Vector< id_t > tGhostIDs( n, 0 );

            hdf5::Dataset< id_t > tBlockIDs( mFile->active_group(), "layers" , n );
            hdf5::Dataset< real > tThicknesses( mFile->active_group(), "thicknesses", n );

            Cell< string > tMaterials ;

            bool tHaveCoatings = false ;

            for ( uint s=0; s<n; ++s )
            {
                ThinShell * tShell = tShells( s );
                tShellIDs( s ) = tShell->id();
                if ( tShell->ghost_id() != gNoID )
                {
                    tGhostIDs( s ) = tShell->ghost_id();
                }

                uint m = tShell->blocks().size();

                id_t * tIDs   = tBlockIDs.set_size( s, m );
                real * tThick = tThicknesses.set_size( s, m );

                for ( uint i=0; i<m; ++i )
                {
                    tIDs[ i ]   = tShell->blocks()( i )->id();
                    tThick[ i ] = tShell->thicknesses()( i );
                    tMaterials.push( tShell->materials()( i ) );
                }

                if ( tShell->side_connector_blocks().size() > 0 ) tHaveCoatings = true ;
            }



            mFile->save_data( "ids", tShellIDs );
            mFile->save_data( "ghost", tGhostIDs );
            mFile->save_data( "materials", tMaterials );

            tBlockIDs.save();
            tThicknesses.save();

            if ( tHaveCoatings )
            {
                hdf5::Dataset< id_t > tCoatingIDs( mFile->active_group(), "coatings" , n );
                hdf5::Dataset< real > tWidths( mFile->active_group(), "widths", n );
                hdf5::Dataset< id_t > tSeamIDs( mFile->active_group(), "seams" , n );

                for ( uint s=0; s<n; ++s )
                {
                    ThinShell * tShell = tShells( s );

                    uint m = tShell->side_connector_blocks().size() ;

                    if ( m == 0 ) continue ;


                    id_t * tIDs   = tCoatingIDs.set_size( s, m );
                    id_t * tSeam  = tSeamIDs.set_size( s, m );
                    real * tWidth = tWidths.set_size( s, m );

                    for ( uint k=0; k<m; ++k )
                    {
                        Block * tBlock =  tShell->side_connector_blocks()( k ) ;

                        tIDs[ k ] = tBlock->id() ;
                        tSeam[ k ] = tShell->side_connector_sidesets()( k )->id() ;
                        tWidth[ k ] = tBlock->thickness() ;
                    }
                }

                tCoatingIDs.save();
                tSeamIDs.save();
                tWidths.save();
            }
            mFile->close_active_group();

#endif
        }

        void
        BfmFile::load_thinshell_data()
        {
#ifdef BELFEM_HDF5

            if ( ! hdf5::group_exists( mFile->active_group(), "thinshells" ) ) return;

            mFile->select_group( "thinshells" );

            Vector< id_t > tShellIDs ;
            mFile->load_data( "ids", tShellIDs );

            Vector< id_t > tGhostIDs ;
            mFile->load_data( "ghost", tGhostIDs );

            Cell< string > tMaterials ;
            mFile->load_data( "materials", tMaterials );

            hdf5::Dataset< id_t > tBlockIDs( mFile->active_group(), "layers" );
            hdf5::Dataset< real > tThicknesses( mFile->active_group(), "thicknesses" );

            uint n = tShellIDs.length();

            DynamicBitset tBitset( mMesh->number_of_nodes() );

            mMesh->thin_shells().reserve( n );

            index_t tCount = 0 ;

            bool tHaveCoatings = hdf5::dataset_exists( mFile->active_group(), "coatings" );

            for ( uint s=0; s<n; ++s )
            {
                SideSet * tSideSet = mProto->sideset( tShellIDs( s ) );
                SideSet * tGhostSideset = tGhostIDs( s ) == 0 ? nullptr : mProto->sideset( tGhostIDs( s ) );

                ThinShell * tThinShell
                    = new ThinShell( tSideSet, tGhostSideset );

                for ( Facet * tFacet : tSideSet->facets() )
                {
                    for ( uint k=0; k<tFacet->number_of_nodes(); ++k )
                    {
                        tBitset.set( tFacet->node( k )->original()->index() );
                    }
                }

                Cell< index_t > tIndices ;
                tBitset.where( tIndices );
                tThinShell->move_node_indices( tIndices );

                hsize_t m = tBlockIDs.length( s );
                Vector< real > tThick( m );

                const id_t * tIDs = tBlockIDs[ s ];
                tThinShell->blocks().reserve( m );
                Cell< string > tMats ;
                tMats.reserve( m );

                for ( hsize_t b=0; b<m; ++b )
                {
                    Block * tBlock = mProto->block( tIDs[ b ] ) ;
                    tThinShell->blocks().push( tBlock );
                    tThick( b ) = tThicknesses[ s ][ b ];
                    // no set_thickness here: ThinShell::set_thicknesses below
                    // propagates the vector onto the layer blocks
                    //tBlock->set_thickness( tThick( b ) );
                    tMats.push( tMaterials( tCount++) );
                }
                tThinShell->set_thicknesses( tThick );
                tThinShell->set_materials( tMats );
                mMesh->thin_shells().push( tThinShell );
                if ( s+1 < n ) tBitset.reset();
            }

            if ( tHaveCoatings )
            {
                hdf5::Dataset< id_t > tCoatingIDs( mFile->active_group(), "coatings" );
                hdf5::Dataset< id_t > tSeamIDs( mFile->active_group(), "seams" );
                hdf5::Dataset< real > tWidths( mFile->active_group(), "widths" );

                Cell< Edge * > tEdges ;

                for ( uint s=0; s<n; ++s )
                {
                    uint m = tCoatingIDs.length( s );
                    if ( m == 0 ) continue ;

                    const id_t * tIDs     = tCoatingIDs[ s ];
                    const id_t * tIDs2   = tSeamIDs[ s ];
                    const real * tWidth  = tWidths[ s ];

                    ThinShell * tShell = mMesh->thin_shells()( s );

                    tShell->side_connector_blocks().reserve( m );
                    tShell->side_connector_sidesets().reserve( m );

                    for ( uint k=0; k<m; ++k )
                    {
                        Block * tBlock = mProto->block( tIDs[ k ] ) ;
                        tBlock->set_thickness( tWidth[ k ] );
                        tShell->side_connector_blocks().push( tBlock ) ;

                        SideSet * tSeam = mProto->sideset( tIDs2[ k ] ) ;
                        tShell->side_connector_sidesets().push( tSeam ) ;

                        // the recovery facets carry the two longitudinal dof
                        // edges of the wall's inner face. The file does not
                        // store them: the slave face determines them
                        // completely, and the wall's edge slots were already
                        // rebuilt by reconstruct_edge_connectivity
                        for ( Facet * tFacet : tSeam->facets() )
                        {
                            tFacet->slave()->get_edges_of_facet(
                                tFacet->index_on_slave(), tEdges );

                            tFacet->element()->allocate_edge_container();
                            tFacet->element()->insert_edge( tEdges( 0 ), 0 );
                            tFacet->element()->insert_edge( tEdges( 1 ), 1 );
                        }
                    }
                }
            }

            mFile->close_active_group();
#endif
        }

        void
        BfmFile::save_vertex_data()
        {
#ifdef BELFEM_HDF5
            Cell< Element * > & tVertices = mMesh->vertices();

            index_t n = tVertices.size();

            if ( n == 0 ) return ;
            Vector< id_t > tIDs( n );
            Vector< id_t > tNodes( n );

            index_t tCount = 0 ;
            for ( Element * tVertex : tVertices )
            {
                tIDs( tCount )   = tVertex->id();
                tNodes( tCount )   = tVertex->node( 0 )->id();
                ++tCount;
            }

            mFile->create_group( "vertices" );
            mFile->save_data( "ids", tIDs );
            mFile->save_data( "nodes", tNodes );
            mFile->close_active_group();

#endif
        }

        void
        BfmFile::load_vertex_data()
        {
#ifdef BELFEM_HDF5
            if ( ! hdf5::group_exists( mFile->active_group(), "vertices" ) ) return ;

            mFile->select_group( "vertices" );

            Cell< id_t > tIDs ;
            Cell< id_t > tNodes ;

            mFile->load_data( "ids", tIDs );
            mFile->load_data( "nodes", tNodes );
            BELFEM_ERROR(  tIDs.size() == tNodes.size() , "Vertex data size mismatch" );

            mFile->close_active_group();

            ElementFactory tFactory ;

            index_t n = tIDs.size();
            Cell< Element * > & tVertices = mMesh->vertices();

            tVertices.set_size(  n );
            for ( index_t k=0; k<n; ++k )
            {
                Element * tVertex = tFactory.create_element( ElementType::VERTEX, tIDs( k ) );
                tVertex->insert_node( mProto->node( tNodes( k ) ), 0 );
                tVertices( k ) = tVertex;
            }

#endif
        }

        void
        BfmFile::save_curve_data()
        {
#ifdef BELFEM_HDF5
            Cell< Curve * > & tCurves = mMesh->curves();

            index_t n = tCurves.size();

            if ( n == 0 ) return ;

            mFile->create_group( "curves" );

            Cell< id_t > tIDs( n, 0 );
            Cell< string > tLabels( n );
            Matrix< id_t > tSideSets( 2, n, 0);
            Cell< uchar > tTypes( n );
            Cell< uchar > tClosed( n );

            //hdf5::Dataset< id_t > tEdgeIDs( mFile->active_group(), "edges", n );
            hdf5::Dataset< id_t > tSegmentIDs( mFile->active_group(), "segments", n );
            hdf5::Dataset< real > tLengths( mFile->active_group(), "lengths", n );

            index_t tCount = 0 ;
            for ( Curve * tCurve : tCurves )
            {
                tCount += tCurve->segments().size();
            }

            hdf5::Dataset< id_t > tTopology( mFile->active_group(), "topology", tCount );

            uint c = 0;
            tCount = 0 ;

            for ( Curve * tCurve : tCurves )
            {
                tIDs( c ) = tCurve->id();
                tLabels.push( tCurve->label() );

                if ( tCurve->sideset_a() != nullptr )
                {
                    tSideSets( 0, c ) = tCurve->sideset_a()->id();
                }

                if ( tCurve->sideset_b() != nullptr )
                {
                    tSideSets( 1, c ) =tCurve->sideset_b()->id();
                }

                tTypes.push( static_cast< uchar >( tCurve->element_type() ) );
                tClosed.push( tCurve->is_closed() );

                Cell< Segment * > & tSegments = tCurve->segments();
                hsize_t p = tSegments.size();

                id_t * tSeg  = tSegmentIDs.set_size( c, p  );
                //id_t * tEdge = tEdgeIDs.set_size( c, p  );

                index_t q = 0 ;
                for ( Segment * tSegment : tSegments )
                {
                    //tEdge[ q ] = tSegment->edge()->id();
                    tSeg[ q++ ] = tSegment->id();


                    id_t * tTopo = tTopology.set_size( tCount++, tSegment->number_of_nodes() );
                    for ( uint k=0; k<tSegment->number_of_nodes(); ++k )
                    {
                        tTopo[ k ] = tSegment->node( k )->id();
                    }
                }

                const Vector< real > & tArcLengh = tCurve->arclength();
                hsize_t s = tArcLengh.length();
                real * tS = tLengths.set_size( c, s );
                for ( index_t i = 0 ; i < s ; ++i )
                {
                    tS[ i ] = tArcLengh( i );
                }
                ++c ;
            }

            mFile->save_data( "ids", tIDs );
            mFile->save_data( "labels", tLabels );
            mFile->save_data( "sidesets", tSideSets, true );
            mFile->save_data( "types", tTypes );
            mFile->save_data( "closed", tClosed );
            //tEdgeIDs.save();
            tSegmentIDs.save();
            tLengths.save();
            tTopology.save();

            mFile->close_active_group();
#endif
        }

        void
        BfmFile::load_curve_data()
        {
#ifdef BELFEM_HDF5

            if ( ! hdf5::group_exists( mFile->active_group(), "curves" ) ) return;

            mFile->select_group( "curves" );

            Cell< id_t > tIDs ;
            mFile->load_data( "ids", tIDs );

            Cell< string > tLabels ;
            mFile->load_data( "labels", tLabels );

            Matrix< id_t > tSideSets ;
            mFile->load_data( "sidesets", tSideSets, true );

            Cell< uchar > tTypes ;
            mFile->load_data( "types", tTypes );

            Cell< uchar > tClosed ;
            mFile->load_data( "closed", tClosed );

            //hdf5::Dataset< id_t > tEdgeIDs( mFile->active_group(), "edges" );
            hdf5::Dataset< id_t > tSegmentIDs( mFile->active_group(), "segments" );
            hdf5::Dataset< real > tLengths( mFile->active_group(), "lengths" );
            hdf5::Dataset< id_t > tTopology( mFile->active_group(), "topology" );

            index_t n = tIDs.size();

            Cell< Curve * > & tCurves = mMesh->curves();
            tCurves.set_size( n, nullptr );

            ElementFactory tFactory ;

            index_t tCount = 0;

            for ( index_t c=0; c<n; ++c )
            {
                ElementType tType =  static_cast< ElementType >( tTypes( c ) );
                Curve * tCurve = new Curve( tIDs( c ), tType );

                tCurve->label() = tLabels( c );

                if ( tSideSets( 0, c ) != 0 ) tCurve->sideset_a( mProto->sideset( tSideSets( 0, c ) ));
                if ( tSideSets( 1, c ) != 0 ) tCurve->sideset_b( mProto->sideset( tSideSets( 1, c ) ));

                hsize_t m = tSegmentIDs.length( c );

                Cell< Segment * > & tSegments = tCurve->segments();

                tSegments.set_size( m, nullptr );
                const id_t * tSids = tSegmentIDs[ c ];
                //const id_t * tEids = tEdgeIDs[ c ];
                const real * tLen  = tLengths[ c ];
                for ( index_t s=0; s<m; ++s )
                {
                    Element * tElement = tFactory.create_element( tType , tSids[ s ] );

                    const id_t * tNids = tTopology[ tCount++ ];
                    for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                    {
                        tElement->insert_node( mProto->node( tNids[ k ] ) , k );
                    }

                    Segment * tSegment = new Segment( tElement );

                    tSegments( s ) = tSegment;
                }

                tCurve->set_closed_flag( tClosed( c ) != 0 );
                //tCurve->assign_edges(); only the homology generator needs edges on curves.

                m = tLengths.length( c );
                Vector< real > & tS = tCurve->arclength() ;
                tS.set_size( m );
                for ( hsize_t s=0; s<m; ++s )
                {
                    tS( s ) = tLen[ s ];
                }

                tCurves( c ) = tCurve;
            }

            tSegmentIDs.close();
            tLengths.close();

            mFile->close_active_group();

#endif
        }

    }
}