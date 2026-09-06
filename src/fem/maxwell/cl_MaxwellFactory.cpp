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

#include <cmath>
#include <filesystem>
#include <memory>
#include <utility>


#include "cl_MaxwellFactory.hpp"
#include "fn_check_facet_orientation.hpp"
#include "fn_check_unit.hpp"
#include "meshtools.hpp"

#include "cl_CutFactory.hpp"
#include "cl_Mesh_PeriodicityFactory.hpp"
#include "cl_Mesh_Periodicity.hpp"
#include "cl_CutProcessorManual.hpp"
#include "cl_FEM_SideSet.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_FEM_Calculator.hpp"
#include "cl_Curve.hpp"
#include "cl_CurveFactory.hpp"
#include "cl_IWG_MaxwellPostproc.hpp"
#include "fn_to_master_orientation.hpp"
#include "cl_FEM_Postprocessor.hpp"
#include "cl_FEM_Controller.hpp"

#include "fn_unique.hpp"
#include "fn_max.hpp"
#include "fn_combine.hpp"
#include "cl_MaterialFactory.hpp"
#include "cl_MaxwellPostprocessor.hpp"
#include "cl_Logger.hpp"
#include "cl_ThinShellFactory.hpp"
#include "fn_FEM_ghost_switch.hpp"
#include "cl_Maxwell_TMatrix.hpp"
#include "cl_MeshChecker.hpp"

#include "cl_Mesh_BfmFile.hpp"
#include "cl_Mesh_ConnectivityCalculator.hpp"
#include "cl_Queue.hpp"
#include "en_DomainType.hpp"
#include "constants.hpp"
#include "fn_Graph_find_connected_partitions.hpp"
#include "fn_Graph_symrcm.hpp"
#include "fn_Mesh_symrcm_nodes.hpp"
#include "fn_mesh_config_tag.hpp"
#include "fn_Graph_multibfs.hpp"
#include "op_Graph_Vertex_ID.hpp"
#include "op_Graph_Vertex_Index.hpp"

namespace belfem
{
    namespace fem {
    //------------------------------------------------------------------------------

    MaxwellFactory::MaxwellFactory( const string &aInputFile ) :
        mCommRank( comm_rank() ),
        mInputFile( new InputFile( aInputFile ) )
    {
        this->create_materials();

        // synchronize data
        Vector< uint > tIData( 3 );
        if ( mCommRank == 0 )
        {
            // read mesh, either a gmsh mesh or a belfem mesh
            mMesh = this->read_mesh();

            //mMesh->save( "mesh.exo" );

            mTopology = new mesh::Topology( mMesh );

            // create the commtable
            uint n = comm_size() - 1;
            if ( n > 0 )
            {
                mCommTable.set_size( n );
                for ( uint k = 0; k < n; ++k )
                {
                    mCommTable( k ) = k + 1;
                }
            }
            mMeshDimension = mMesh->number_of_dimensions();
            mElementOrder  = mMesh->max_element_order();
            tIData( 0 ) = mMeshDimension ;
            tIData( 1 ) = mElementOrder ;
            tIData( 2 ) = mComputeCohomologies == 1;

            comm_barrier();

            broadcast( tIData );
        }
        else
        {
            mTopology = new mesh::Topology( nullptr );

            comm_barrier();
            broadcast( tIData );

            mMeshDimension = tIData( 0 );
            mElementOrder = tIData( 1 );
            mComputeCohomologies = tIData( 2 ) == 1;

            mMesh = new Mesh( mMeshDimension, 0 );
        }
        if ( mCommRank == 0 )
        {
            this->read_domain_types();
            if ( mMesh->curves().size() == 0 )
            {
                // domain_type() accepts "curve" as well as "curves", so
                // honour both spellings here -- looking up only the plural
                // left a singular "curve { }" section silently unread
                const input::Section * tTopo = mInputFile->section( "topology" ) ;

                if ( tTopo->section_exists( "curves" ) )
                {
                    this->create_curves( tTopo->section( "curves" ) ) ;
                }
                else if ( tTopo->section_exists( "curve" ) )
                {
                    this->create_curves( tTopo->section( "curve" ) ) ;
                }
            }

            if ( ! mMesh->has_periodicity() )
            {
                if ( mInputFile->section( "topology" )->section_exists( "periodic" ) )
                {
                    this->create_periodic( mInputFile->section( "topology" )->section( "periodic" ) ) ;
                }
            }
        }

        //Create the Maxwell BC Factory
        BELFEM_ERROR(mInputFile->section_exists("boundary conditions" ),"No boundary conditions defined") ;
        if ( mInputFile->section( "boundary conditions" )->section_exists( "maxwell" ) )
        {
            mBoundaryConditionFactory = new MaxwellBoundaryConditionFactory( mInputFile->section( "boundary conditions" )->section( "maxwell" ),
                                                                                mMesh->number_of_dimensions() ) ;
        }
        else
        {
            mBoundaryConditionFactory = new MaxwellBoundaryConditionFactory( mInputFile->section( "boundary conditions" ),
                                                                                mMesh->number_of_dimensions() ) ;
        }

        BELFEM_ERROR(mInputFile->section_exists("boundary conditions" ),"No boundary conditions defined") ;

        // Read temperature value if it exists
        if ( mInputFile->section_exists("initial conditions" ) )
        {
            const input::Section * tSection = mInputFile->section( "initial conditions" ) ;
            for ( uint i = 0 ; i < tSection->num_keys() ; ++i )
            {
                string tKey = tSection->key( i ) ;
                if (tKey == "t" || tKey == "temp" || tKey == "temperature" )
                {
                    string tUnitDef = tSection->get_units( tKey ) ;
                    value tValue = unit_to_si( tUnitDef );
                    BELFEM_ERROR(check_unit( tValue, "K" ),
                                "required unit for the temperature is: K");

                    gTbulk = tSection->get_value( tKey,"K" ).first ;
                }
            }
        }


        //Set the sideset domain type for the appropriate boundary conditions
        if ( mCommRank == 0 )
        {
            for ( PhysicalBoundaryCondition * tBC : mBoundaryConditionFactory->boundary_conditions() )
            {
                if ( tBC->type() == BoundaryConditionType::Background )
                {
                    for (id_t tID : tBC->domains())
                    {
                        mMesh->sideset( tID )->set_domain_type( DomainType::BackgroundField );

                        // A deck that states the background as a flux density
                        // is read as H = B/mu0, which is exact only where the
                        // permeability IS mu0. Of the block domain types only
                        // Ferro carries a B-H curve -- air and coil are
                        // free-space, a conductor is a non-magnetic HTS, and a
                        // buffer is magnesia -- so ferro adjacency is the one
                        // case where the conversion misstates the field, by the
                        // relative permeability. Warn and continue: the deck is
                        // still runnable and the imposed H is well defined; it
                        // is the B the user THINKS they asked for that is not
                        if ( tBC->amplitude_is_flux_density() )
                        {
                            mesh::SideSet * tSideSet = mMesh->sideset( tID );

                            for ( mesh::Facet * tFacet : tSideSet->facets() )
                            {
                                id_t tFerroBlock = gNoID ;

                                // Both guards are load-bearing, and neither is
                                // redundant with the other. This walk runs in the
                                // CONSTRUCTOR: mMaster is nullable and a .bfm may
                                // carry an only-slave facet, and an element that
                                // belongs to no registered block keeps block id 0,
                                // for which Mesh::block() asserts and then throws
                                // out of the map. Topology::detect_sideset_types
                                // dereferences master() unguarded, but it runs
                                // much later on a fully wired mesh.
                                //
                                // This is a DIAGNOSTIC. It must never abort the
                                // run it is only commenting on, so it asks before
                                // it looks -- and stays silent when it cannot tell
                                if ( tFacet->has_master()
                                    && mMesh->block_exists( tFacet->master()->block_id() )
                                    && mMesh->block( tFacet->master()->block_id() )
                                        ->domain_type() == DomainType::Ferro )
                                {
                                    tFerroBlock = tFacet->master()->block_id() ;
                                }
                                else if ( tFacet->has_slave()
                                    && mMesh->block_exists( tFacet->slave()->block_id() )
                                    && mMesh->block( tFacet->slave()->block_id() )
                                        ->domain_type() == DomainType::Ferro )
                                {
                                    tFerroBlock = tFacet->slave()->block_id() ;
                                }

                                if ( tFerroBlock != gNoID )
                                {
                                    message( InfoLevel::Default,
                                        "\n    Warning: the background condition on sideset %lu states its\n"
                                        "             amplitude as a flux density, but that sideset borders\n"
                                        "             ferro block %lu. B = mu0*H holds only outside a magnetic\n"
                                        "             material, so the imposed field strength H = B/mu0 does\n"
                                        "             NOT correspond to the requested B there. State the\n"
                                        "             amplitude in A/m if that is not what you meant.",
                                        ( long unsigned int ) tID,
                                        ( long unsigned int ) tFerroBlock );

                                    // one warning per sideset, not per facet
                                    break ;
                                }
                            }
                        }
                    }
                }
            }
        }


        // create a parameter object
        mKernelParameters = new KernelParameters( mMesh );
    }

    const string &
    MaxwellFactory::mesh_path() const
    {
        return mInputFile->section( "mesh" )->get_string( "file" );
    }

    //------------------------------------------------------------------------------

    MaxwellFactory::~MaxwellFactory()
    {
        delete mInputFile;

        delete mBoundaryConditionFactory;

        delete mPeriodicFactory ;

        if ( mTopology != nullptr )
        {
            delete mTopology;
        }
        for ( mesh::SideSet *tSideSet : mTapes )
        {
            for ( mesh::Facet *tFacet : tSideSet->facets() )
            {
                delete tFacet;
            }
            delete tSideSet;
        }

        for ( Protoshell *tShell : mProtoshells )
        {
            delete tShell;
        }

        this->delete_unused_materials();

        for ( Domain *tDomain : mDomains )
        {
            delete tDomain;
        }
        if ( mOwnKernelParameters && mKernelParameters != nullptr )
        {
            delete mKernelParameters;
        }
        if ( mOwnMagneticEquation && mMagneticEquation != nullptr )
        {
            delete mMagneticEquation;
        }
    }

    //------------------------------------------------------------------------------

    Mesh *
    MaxwellFactory::read_mesh()
    {
        // get name of mesh
        string tMeshPath = mInputFile->section( "mesh" )->get_string(
            "file" );

        mLabel = filename( tMeshPath );
        mLabel = mLabel.substr( 0, mLabel.find_last_of( "." ) );

        // check if file is configured as bfm
        if ( tMeshPath.substr( tMeshPath.find_last_of( "." ) + 1 ) == "bfm" )
        {
            // load bfm mesh
            mesh::BfmFile tFile( tMeshPath );

            // the deck names a processed mesh directly, so there is no .msh
            // to validate against and nothing to rebuild from -- this is the
            // "I was handed a prepared mesh" path. Warn if the settings the
            // file was built with disagree with this deck, then continue:
            // the user asked for this file explicitly
            const uint64_t tFileTag = tFile.config_tag();

            if ( tFileTag != 0
                 && tFileTag != maxwell::mesh_config_tag( *mInputFile ) )
            {
                message( InfoLevel::Default,
                         "\n    Warning: %s was built with a different mesh configuration\n"
                         "             than this input file describes; using the file as given",
                         filename( tMeshPath ).c_str() );
            }

            // the one setting the tag carries that changes the DISCRETIZATION
            // of the thin shells, not just the run: a file built with ghost
            // facets ( duplicate interface dofs ) under a deck that asks for
            // none, or the reverse, is refused rather than warned about --
            // the dof layout would silently differ from what the deck says
            {
                // a file written before the switch existed ( no line in its
                // config text ) was built by the old factory, whose flag was
                // hardcoded ON -- so "no line" means ON, not unknown, and a
                // legacy file cannot slip past this check ( audit finding )
                const string tText = tFile.config_text();
                const bool tHaveLine  = tText.find( "thinshell.ghost = " ) != string::npos ;
                const bool tFileGhost = tHaveLine ?
                    tText.find( "thinshell.ghost = on" ) != string::npos : true ;
                const bool tDeckGhost = fem::ghost_facets_requested( *mInputFile );

                // and a ghost-ON layout cannot be loaded at all -- see the
                // cache branch below for why the reload degrades it
                BELFEM_ERROR( ! tDeckGhost,
                    "%s cannot be used with the thin-shell ghost on: a reloaded mesh loses the "
                    "duplicate interface edges. Name the .msh in the deck instead.",
                    filename( tMeshPath ).c_str() );

                BELFEM_ERROR( tFileGhost == tDeckGhost,
                    "%s was built with the thin-shell ghost facets %s%s, but this deck asks for them "
                    "%s ( nonlinear magnetic { nitsche ghost penalty { eta } } ). Rebuild the file "
                    "from the .msh, or match the deck to the file.",
                    filename( tMeshPath ).c_str(),
                    tFileGhost ? "ON" : "OFF",
                    tHaveLine ? "" : " ( a file from before 2026-09-01, when ON was the only layout )",
                    tDeckGhost ? "ON" : "OFF" );
            }

            tFile.load();

            // no need to compute cohomologies for bfm meshes
            mComputeCohomologies = false;

            return tFile.get();
        }

        // next we check if there is a corresponding Bfm file
        std::filesystem::path tBfmFilePath( tMeshPath );
        tBfmFilePath.replace_extension( ".bfm" );
        mMeshPath = tBfmFilePath.string();

        // load original mesh
        Mesh * aMesh =  new Mesh( tMeshPath, 0, false, false );


        // we need to scale the original mesh using the given unit
        string tUnit = mInputFile->section( "mesh" )->get_string(
               "unit" );

        // test if the unit is valid
        value tValue = unit_to_si( tUnit );
        BELFEM_ERROR(
            check_unit( tValue, "m" ),
            "need a length unit for the mesh, is: %s",
            tUnit.c_str() );

        // make sure the mesh is scaled to SI
        aMesh->scale_mesh( tValue.first );

        // test the checksum
        size_t tChecksum = aMesh->checksum();

        // ... and the settings the mesh would be enriched with. The checksum
        // identifies the BASE mesh only, so it cannot see an edited deck on
        // an unchanged .msh; the tag cannot see a changed .msh. Cache reuse
        // needs both to match
        const string   tConfigText = maxwell::mesh_config_text( *mInputFile );
        const uint64_t tConfigTag  = maxwell::mesh_config_hash( tConfigText );

        // carry them along so they land in the file we are about to write
        aMesh->set_config_tag( tConfigTag, tConfigText );

        // check if bfm file exists
        if ( std::filesystem::exists( tBfmFilePath ) )
        {
            // OK, the mesh exists, now we need to load the original mesh
            // and test the checksums
            // get the unit

            mesh::BfmFile tBfmFile( tBfmFilePath );

            // if the checksum is identical
            // we can skip the cohomology computation
            // and load the processed mesh directly
            if ( tChecksum == tBfmFile.checksum() )
            {
                const uint64_t tFileTag = tBfmFile.config_tag();

                // a cache is never reused while the ghost is ON: the reload
                // relinks elements to edges by node pair
                // ( ProtoMesh::reconstruct_edge_connectivity ), and the
                // duplicate interface edges of a DG stack share their node
                // pairs with the originals -- they come back as orphans and
                // the stack silently degrades to shared edges ( measured
                // 2026-09-01: a fresh ghost-ON build carries 5 x 24908
                // duplicate edge dofs on the tape deck, every reload of it
                // carried none ). Rebuilding is the only correct answer until
                // the reload preserves twin edges
                if ( tFileTag == tConfigTag && fem::ghost_facets_requested( *mInputFile ) )
                {
                    message( InfoLevel::Default,
                             "\n    %s is not reused: the thin-shell ghost is on, and a reloaded mesh\n"
                             "    loses the duplicate interface edges ( rebuilding from the .msh )",
                             filename( tBfmFilePath.string() ).c_str() );
                }
                else if ( tFileTag == tConfigTag )
                {
                    delete aMesh ;
                    tBfmFile.load();
                    mComputeCohomologies = false ;
                    aMesh = tBfmFile.get();

                    message( InfoLevel::Verbose,
                             "    Reusing processed mesh %s ( base mesh and settings unchanged )",
                             filename( tBfmFilePath.string() ).c_str() );

                    // BfmFile::load() has set the mesh checker flag:
                    // we trust BFM meshes
                    return aMesh;
                }

                // the base mesh matches but the settings do not: rebuilding
                // is the only correct answer, and the user has to be told,
                // because the run is about to take much longer than usual
                if ( tFileTag == 0 )
                {
                    message( InfoLevel::Default,
                             "\n    %s predates the mesh-configuration check and is being rebuilt",
                             filename( tBfmFilePath.string() ).c_str() );
                }
                else
                {
                    message( InfoLevel::Default,
                             "\n    %s was built with a different mesh configuration and is being rebuilt",
                             filename( tBfmFilePath.string() ).c_str() );

                    this->report_config_difference(
                            tBfmFile.config_text(), tConfigText );
                }
            }
        }

        // if the BFM file doesn't exist or doesn't match
        // we need to finalize the mesh and run the cohomology computation

        if ( ! aMesh->mesh_checker_flag() )
        {
            // make sure that the mesh is oriented correctly
            MeshChecker tCheck( aMesh );
        }

        aMesh->unfinalize();
        aMesh->set_connectivity( Connectivity::Compute );
        aMesh->finalize();

        mComputeCohomologies = true;
        return aMesh;

    }

    void
    MaxwellFactory::report_config_difference(
            const string & aStored,
            const string & aCurrent )
    {
        // both texts are sorted, one setting per line, so a plain set
        // difference names exactly what moved
        Cell< string > tStored  = string_to_words( aStored, '\n' );
        Cell< string > tCurrent = string_to_words( aCurrent, '\n' );

        // cap the report: a wholesale topology rewrite must not bury the
        // console in a diff nobody reads
        const uint tMaxLines = 8 ;
        uint tCount = 0 ;

        for ( const string & tLine : tCurrent )
        {
            bool tFound = false ;
            for ( const string & tOther : tStored )
            {
                if ( tLine == tOther )
                {
                    tFound = true ;
                    break ;
                }
            }

            if ( ! tFound )
            {
                if ( tCount++ < tMaxLines )
                {
                    message( InfoLevel::Default, "        now : %s", tLine.c_str() );
                }
            }
        }

        for ( const string & tLine : tStored )
        {
            bool tFound = false ;
            for ( const string & tOther : tCurrent )
            {
                if ( tLine == tOther )
                {
                    tFound = true ;
                    break ;
                }
            }

            if ( ! tFound )
            {
                if ( tCount++ < tMaxLines )
                {
                    message( InfoLevel::Default, "        was : %s", tLine.c_str() );
                }
            }
        }

        if ( tCount > tMaxLines )
        {
            message( InfoLevel::Default, "        ... and %u more",
                     ( unsigned int ) ( tCount - tMaxLines ) );
        }
    }

    //------------------------------------------------------------------------------

    void
    MaxwellFactory::read_domain_types()
    {
        uint n = mInputFile->section( "topology" )->num_sections();

        mDomains.set_size( n, nullptr );
        for ( uint d = 0; d < n; ++d )
        {
            mDomains( d ) = new Domain(
                mInputFile->section( "topology" )->section( d ) );
        }

        for ( Domain * tDomain : mDomains )
        {
            if ( tDomain->is_block() )
            {
                for ( id_t tID : tDomain->group_ids() )
                {
                    mMesh->block( tID )->set_domain_type( tDomain->type() );
                }
            }
            else if (tDomain->is_sideset() )
            {
                for ( id_t tID : tDomain->group_ids() )
                {
                    mMesh->sideset( tID )->set_domain_type( tDomain->type() );
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    MaxwellFactory::create_curves( const input::Section * aSection )
    {
        uint n = aSection->num_keys() ;
        mesh::CurveFactory tFactory = mesh::CurveFactory( mMesh ) ;

        for ( uint i = 0; i < n; ++i )
        {
            mesh::Curve * tCurve ;
            std::string tKey = aSection->key( i ) ;
            id_t tID = std::stoi( tKey ) ;
            if (mMesh->number_of_dimensions() == 2)
            {
                Vector < id_t > aIDs ;
                aSection->get_ids( tKey, aIDs ) ;
                tCurve = tFactory.from_2d_sidesets( aIDs , tID ) ;
            }
            else
            {
                std::pair<id_t, id_t> tPair = aSection->get_intersection_ids( tKey ) ;
                tCurve = tFactory.intersect( tPair.first, tPair.second, tID ) ;
            }

            // Add these curves to the list and the mesh
            mCurves.push(tCurve) ;
            mMesh->curves().push(tCurve) ;
        }
        mMesh->curves().shrink_to_fit();
    }

    //------------------------------------------------------------------------------

    void
    MaxwellFactory::create_periodic( const input::Section * aSection )
    {

        //Read source and target nodes
        BELFEM_ERROR( aSection->key_exists( "source" ), "Undefined source nodes for periodic" ) ;
        BELFEM_ERROR( aSection->key_exists( "target" ), "Undefined target nodes for periodic" ) ;

        //Get source nodes
        Vector<id_t> tSourceIDs;
        aSection->get_ids( "source", tSourceIDs) ;

        Vector<id_t> tTargetIDs;
        aSection->get_ids( "target", tTargetIDs) ;

        BELFEM_ERROR( tSourceIDs.length() == 3, "Exactly 3 nodes must be defined on source" ) ;
        BELFEM_ERROR( tTargetIDs.length() == 3, "Exactly 3 nodes must be defined on target" ) ;


        // Create the periodicity
        mPeriodicFactory = new mesh::PeriodicityFactory( mMesh );

        mPeriodicFactory->set_master_plane( tSourceIDs(0), tSourceIDs(1), tSourceIDs(2) );
        mPeriodicFactory->set_slave_plane( tTargetIDs(0), tTargetIDs(1), tTargetIDs(2) );

        mPeriodicFactory->tag_periodic_sidesets();
    }

    //------------------------------------------------------------------------------

    IWG_Maxwell *
    MaxwellFactory::create_equation(
        const maxwell::Formulation aFormulation )
    {
        mFormulation = aFormulation;
        mDimensionality = mMeshDimension == 2
            ? ModelDimensionality::TwoD
            : ModelDimensionality::ThreeD;

        return new IWG_Maxwell( mFormulation, mDimensionality,
                                mElementOrder == 2, mUseEnrichment );

    }

    //------------------------------------------------------------------------------

    std::shared_ptr< Kernel >
    MaxwellFactory::create_magnetic_kernel()
    {
        mMagneticEquation = this->create_equation( mFormulation );

        if ( mComputeCohomologies )
        {
            this->create_cuts();
        }
        else
        {
            // the fresh path fills the type map inside create_cuts(), BEFORE
            // the enrichment factories run. On the reloaded (enriched) mesh
            // we must reproduce that pre-enrichment view: skip the thin-shell
            // layer/buffer blocks and tape/ghost sidesets, and do not
            // re-detect sideset types. synchronize_maps() is the COLLECTIVE
            // broadcast of the map — without it the non-root ranks have an
            // empty type map and groups(Conductor) fails.
            mTopology->run_on_enriched_mesh() ;
            mTopology->synchronize_maps() ;

            if ( mCommRank == 0 )
            {
                this->create_block_to_material_map();

                mMagneticEquation->
                    set_abstract_nodes( mMesh->abstract_nodes() );
                mMagneticEquation->
                    set_orphaned_nodes( mMesh->orphaned_nodes() );
            }
        }

        this->create_edges_and_faces_on_mesh();
        this->create_thinshells();

        this->synch_material_map();

        //if ( mCommRank == 0 ) mMesh->save( "debug.exo" );

        comm_barrier();


        this->set_block_types_in_magnetic_equation();

        comm_barrier();

        if ( mCommRank == 0 && mMesh->has_periodicity() )
        {
            if ( ! mMesh->periodicity()->node_pairs_restored() )
            {
                mMesh->periodicity()->update() ;
                mMesh->periodicity()->set_entity_dependencies();
            }
        }

        if ( mCommRank == 0 && mComputeCohomologies )
        {
            mMesh->flag_curved_elements();
            this->create_hanging_edges_and_facets();

            mesh::symrcm( mMesh->nodes() );

            // rebuild indices
            mMesh->update_node_indices() ;
            mMesh->update_edge_indices() ;
            mMesh->update_face_indices() ;
            mMesh->update_element_indices() ;
            mMesh->update_facet_indices() ;

            // cut enrichment is currently unsupported: the bubble
            // functions in BELFEM are not the right space for it, and
            // the .bfm topology snapshot does not represent enriched
            // meshes. Refuse to cache rather than reload wrongly
            BELFEM_ERROR( ! mUseEnrichment,
                "storing a .bfm of an enriched mesh is not supported" );

            mMesh->save( mMeshPath );

        }

        comm_barrier();

        // selecting the blocks
        Cell< id_t > tBlocks ;
        for ( mesh::Block * tBlock : mMesh->blocks() )
        {
            switch ( tBlock->domain_type() )
            {
            case DomainType::Air :
            case DomainType::Buffer :
            case DomainType::Conductor :
            case DomainType::Ferro :
            case DomainType::ThinShell :
            {
                tBlocks.push( tBlock->id() ) ;
                break ;
            }
            default:
            {
                break ;
            }
            }
        }
        mKernelParameters->select_blocks( tBlocks );

        mMagneticKernel = std::make_shared< Kernel >( mKernelParameters );

        Mesh * tMesh = mMagneticKernel->mesh();

        mMagneticKernel->claim_parameter_ownership( true );
        mOwnKernelParameters = false ;

        mMagneticKernel->mesh()->flag_curved_elements();

        // link equation to kernel and surrender ownership
        mMagneticKernel->add_equation( mMagneticEquation );
        mOwnMagneticEquation = false ;

        mMagneticField = mMagneticKernel->create_field( mMagneticEquation );

        mMagneticField->solver_data()->use_full_force( true );

        // get the solver section from the input file
        const input::Section *tSolverSection = mInputFile->section(
            "solver" );

        // when we have thermal, we will be able to chose different settings for the thermal and the magnetic part
        const input::Section *tLinearSection = tSolverSection->
            section_exists( "linear magnetic" )
            ? tSolverSection->section( "linear magnetic" )
            : tSolverSection->section( "linear" );

        //SolverParameters tSolverParameters( tLinearSection );

        this->configure_solver( tLinearSection, mMagneticField );

        //extract the abstract node data to access it later in the controller
        mMagneticField->extract_abstract_dofs_from_mesh() ;

        // update the block types
        for( Block * tGroup : mMagneticKernel->dofmgr()->blocks() )
        {
            switch( tMesh->block( tGroup->id() )->domain_type() )
            {
            case DomainType::Air :
            case DomainType::Buffer :
            case DomainType::Conductor :
            case DomainType::Ferro :
            case DomainType::ThinShell :
            case DomainType::LeftCoating :
            case DomainType::RightCoating :
            {
                tGroup->set_domain_type( tMesh->block( tGroup->id() )->domain_type() );
                break ;
            }
            default:
            {
                // activation mode is set automatically based on domain type
                tGroup->set_activation_mode( GroupActivationMode::Inactive ) ;
                break ;
            }
            }

        }

        for ( Group *tGroup : mMagneticKernel->dofmgr()->blocks() )
        {
            tGroup->set_domain_type( tMesh->block( tGroup->id() )->domain_type() );
        }

        for ( Group *tGroup : mMagneticKernel->dofmgr()->sidesets() )
        {
            tGroup->set_domain_type( tMesh->sideset( tGroup->id() )->domain_type() );
        }

        // get enriched order
        uint tOrder = mUseEnrichment ? ( tMesh->max_element_order() == 1 ? 5 : 13 ) : 0 ;

        if ( mUseEnrichment )
        {
            for ( Block *tBlock : mMagneticKernel->dofmgr()->blocks() )
            {
                tBlock->set_integration_order( tOrder );
            }
        }


        // deactivate all antisymmetry BCs and unwanted interface
        for( SideSet * tGroup : mMagneticKernel->dofmgr()->sidesets())
        {
            if (   tGroup->domain_type() == DomainType::AirAntiSymmetry
                || tGroup->domain_type() == DomainType::BufferAntiSymmetry
                || tGroup->domain_type() == DomainType::FerroAntiSymmetry
                || tGroup->domain_type() == DomainType::ConductorAntiSymmetry )
            {
                tGroup->set_activation_mode( GroupActivationMode::Inactive );
            }
        }


        for( SideSet * tGroup : mMagneticKernel->dofmgr()->sidesets())
        {
            switch ( tMesh->sideset( tGroup->id() )->domain_type() )
            {
            case DomainType::ConductorSymmetry :
            case DomainType::BackgroundField :
            case DomainType::EnrichedInterface :
            {
                tGroup->calculator()->set_integration_order( tOrder ) ;
                tGroup->set_domain_type( tMesh->sideset( tGroup->id() )->domain_type() );
                break;
            }
            case DomainType::AirSymmetry :
            case DomainType::BufferSymmetry :
            case DomainType::FerroSymmetry :
            {
                //Phi symmetry corresponds to a constant phi on the symmetry plane
                tGroup->impose_dirichlet( 0.0 ) ;
                break;
            }
            default:
            {
                // activation mode is set automatically based on domain type
                break;
            }
            }
        }

        //Send the dofmngr to the boundary conditions
        mBoundaryConditionFactory->set_fields( mMagneticField );

        //Send boundary conditions to the kernel
        for ( PhysicalBoundaryCondition *tBC : mBoundaryConditionFactory->
              boundary_conditions() )
        {
            mMagneticKernel->add_boundary_condition( tBC );
        }

        //Impose the bearing directly here
        bool tHaveBearing = false ;

        for ( PhysicalBoundaryCondition * tBC : mBoundaryConditionFactory->boundary_conditions() )
        {
            if ( tBC->type() == BoundaryConditionType::Bearing )
            {
                for (id_t tID : tBC->domains())
                {
                    mMagneticKernel->dofmgr()->bearing( tID )->impose_dirichlet( 0.0 );
                }
                tHaveBearing = true ;
            }
            //Initialize the Dirichlet conditions, to set the fixed Dofs
            else if ( tBC->type() == BoundaryConditionType::BackgroundDirichlet ||
                      tBC->type() == BoundaryConditionType::Dirichlet)
            {
                for (id_t tID: tBC->domains()) {
                    mMagneticKernel->dofmgr()->sideset( tID )->set_domain_type( DomainType::BackgroundField );
                    for (mesh::Node * tNode: mMagneticKernel->dofmgr()->sideset( tID )->nodes())
                    {
                        if (!tNode->is_duplicate())
                        {
                            Dof * tDof = mMagneticKernel->dofmgr()->dof(
                                    mMagneticKernel->dofmgr()->calculate_dof_id( tNode, 0 ) );

                            // a single-source unit-weight condensation
                            // must pin its SOURCE, or the source stays
                            // FREE when initialize() freezes the split
                            // and the per-timestep imposition can no
                            // longer reroute safely. The
                            // shared helper applies the SAME shape test
                            // as SideSet::impose_dirichlet, so the set
                            // pre-pinned here is exactly the set the
                            // per-timestep path reroutes
                            index_t tNumFirstFlips = 0 ;
                            if ( ! pin_dirichlet_dof(
                                    mMagneticKernel->dofmgr(), tDof, 0.0,
                                    tID, tNode->id(), tNumFirstFlips ) )
                            {
                                tDof->fix( 0.0 );
                            }
                        }
                    }
                }
            }
        }

        if ( ! tHaveBearing )
        {
            for ( mesh::Node * tNode : mMesh->autopins() )
            {
                Dof * tDof = reinterpret_cast< Dof * >( tNode->dof( 0 ) ) ;

                if ( tDof->is_hanging() )
                {
                    bool tFoundNodeSource = false ;
                    for ( uint k=0; k<tDof->number_of_sources() ; ++k )
                    {
                        Dof * tOther = tDof->source( k ) ;
                        if ( tOther->entity_type() == EntityType::NODE )
                        {
                            if ( ! tOther->is_fixed() )
                            {
                                tOther->fix( 0.0 ) ;
                            }
                            tFoundNodeSource = true ;
                            break ;
                        }
                    }
                    BELFEM_ERROR( tFoundNodeSource, "Node %lu seems to be hanging, but no source seems to be a node." ,
                        ( long unsigned int ) tNode->id() );
                }
                else if ( ! tDof->is_fixed() )
                {
                    // the regular case: an interior, free phi dof.
                    // ( an already-fixed dof keeps its value — e.g. a
                    // Dirichlet sideset got there first )
                    tDof->fix( 0.0 ) ;
                }
            }

            // backwards compatibility: a .bfm written before the autopin
            // feature carries no "pinned" data, and without a bearing in
            // the deck the potential would float undetected —
            // say so instead of running silently
            if ( mMesh->autopins().size() == 0 && mCommRank == 0 )
            {
                message( InfoLevel::Default,
                    "    Warning: no bearing in the input file and no automatic gauge\n"
                    "             pins on the mesh ( .bfm written before the autopin\n"
                    "             feature? ). The magnetic potential is unpinned -\n"
                    "             delete the .bfm to recompute the pins, or set a\n"
                    "             bearing in the input file." );
            }
        }

        comm_barrier();

        mMagneticField->create_fields( mMagneticEquation );
        this->init_fields();

        // hide internal scratch fields from the output
        if ( mMagneticKernel->mesh()->field_exists( "element_rho" ) )
        {
            mMagneticKernel->mesh()->field( "element_rho" )
                ->set_write_to_file_flag( false );
        }

        this->assign_materials();

        this->create_postprocessors();

        mMagneticKernel->dofmgr()->disconnect_dofs_from_mesh();

        this->delete_unused_materials();

        // create one mesh global variable per labeled boundary condition
        // whose type publishes one ( has_block_global ). The terminal
        // conditions — including the circuit terminal pairs, which the
        // circuit factory pushed into the same container — carry a label
        // too, but it is consumed by Controller::save_IV, which writes
        // their I/U pairs. RANK 0 ONLY, by design ( 2026-08-15 ): the
        // globals' consumers — the Exodus writer, the memdump, ParaView —
        // all read rank 0's mesh, no worker code reads them, and keeping
        // worker copies meant replicating label strings over MPI for
        // nothing. update_global() no-ops on ranks where the variable does
        // not exist. A pre-existing variable at this point is a duplicate
        // BC label, never a restart leftover: the executables load the
        // memdump AFTER factory construction
        if ( comm_rank() == 0 )
        {
            for ( PhysicalBoundaryCondition * tBC :
                  mBoundaryConditionFactory->boundary_conditions() )
            {
                if ( tBC->label().size() > 0 && has_block_global( tBC->type() ) )
                {
                    Mesh * tMesh = mMagneticKernel->dofmgr()->mesh() ;

                    BELFEM_ERROR( ! tMesh->global_variable_exists( tBC->label() ),
                        "duplicate boundary condition global '%s' - label the sections to disambiguate",
                        tBC->label().c_str() );

                    // the temperature global is created by Controller::save
                    // under one of these two names; a condition taking either
                    // would be overwritten at the first save
                    BELFEM_ERROR( tBC->label() != "T_max" && tBC->label() != "T_bulk",
                        "a boundary condition is labelled '%s' - rename it, "
                        "the temperature publishes under that name",
                        tBC->label().c_str() );

                    tMesh->create_global_variable( tBC->label(), tBC->value() );
                }
            }
        }

        BELFEM_ERROR( ! mMagneticKernel->dofmgr()->mesh()->global_variable_exists( "dotQ" ),
                   "duplicate mesh global 'dotQ' - label the sections to disambiguate" );

        mMagneticKernel->dofmgr()->mesh()->create_global_variable( "dotQ", 0 );

        return mMagneticKernel;

    }


//------------------------------------------------------------------------------

        std::shared_ptr< Controller >
        MaxwellFactory::create_controller()
        {
            std::shared_ptr< Controller > aControl = std::make_shared<
                Controller >( mMagneticKernel.get() );
            aControl->set_params( mInputFile->section( "solver" ) );

            mMagneticEquation->set_timestepping_method( aControl->euler_method() );

            // register the kernels for the maxwell data helpers; if a thermal
            // kernel is attached later, Controller::set_thermal_kernel relinks
            for ( Block * tBlock : mMagneticKernel->dofmgr()->blocks() )
            {
                if ( tBlock->calculator() != nullptr )
                {
                    tBlock->calculator()->link_maxwell( mMagneticKernel.get() );
                }
            }

            return aControl;
        }

//------------------------------------------------------------------------------

        Cell< PhysicalBoundaryCondition *> &
        MaxwellFactory::boundary_conditions()
        {
            return mBoundaryConditionFactory->boundary_conditions();
        }

//------------------------------------------------------------------------------

        Cell< PhysicalBoundaryCondition *>
        MaxwellFactory::current_BCs()
        {
            return mBoundaryConditionFactory->current_boundary_conditions();
        }

//------------------------------------------------------------------------------

        Mesh *
        MaxwellFactory::mesh()
        {
            return mMesh ;
        }

//------------------------------------------------------------------------------
        string
        MaxwellFactory::label() const
        {
            return mLabel ;
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::create_cuts()
        {
            if ( mCommRank == 0 )
            {
                // creating the protoshells
                this->read_thin_shell_data();

                // set fields in topology class
                mTopology->run() ;

                if ( mComputeCohomologies )
                {
                    mTopology->synchronize_maps();

                    BELFEM_ERROR( mInputFile->section_exists( "homology" ),
                                  "No homology section found in input file" );
                    const input::Section *tSection = mInputFile->section(
                        "homology" );

                    if ( tSection->key_exists( "algorithm" ) )
                    {
                        mCutAlgorithm = mesh::to_cut_algorithm(
                            tSection->get_string( "algorithm" ) );
                    }

                    BELFEM_ERROR(
                        mCutAlgorithm != mesh::CutAlgorithm::UNDEFINED,
                        "Invalid cut algorithm" );

                    this->fix_facet_masters();

                    // deck-signed sideset flips must run AFTER the master
                    // normalization ( which would otherwise undo them ) and
                    // BEFORE the cut / thin-shell pipeline consumes the facet
                    // windings for layer normals
                    this->flip_thin_shell_sidesets();

                    this->create_cuts_sub_master();
                    this->create_block_to_material_map();

                    mMesh->unflag_all_nodes() ;
                }

                mMagneticEquation->
                    set_abstract_nodes( mMesh->abstract_nodes() );
                mMagneticEquation->
                    set_orphaned_nodes( mMesh->orphaned_nodes() );
            }
            else
            {
                if ( mComputeCohomologies )
                {
                    mTopology->synchronize_maps();
                    this->create_cuts_sub_slave();
                }
            }
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::create_cuts_sub_master()
        {
            BELFEM_ASSERT( mCommRank == 0, "This function should only be called on rank 0" );

            mesh::CutFactory tFactory( mMesh, mTopology, mProtoshells, mCutAlgorithm,
                                                mUseEnrichment );

            tFactory.set_periodicity( mPeriodicFactory ) ;

            //Open the thin shells for cohomology purposes
            if (tFactory.create_thin_shell_cuts())
            {
                // for debugging
                //tFactory.save_curve_debug_meshes();

                mThinShellMasterNodes = std::move(tFactory.thin_shell_master_nodes()) ;
                mThinShellSlaveNodes  = std::move( tFactory.thin_shell_slave_nodes() );
            }

            //Create the terminal information to send to the cut factory
            this->create_terminal_list() ;

            tFactory.set_terminals( mTerminals, mThinShellTerminalIndices ) ;

            comm_barrier();
            tFactory.run();

            for ( mesh::SideSet *tCut : tFactory.cuts() )
            {
                tCut->set_domain_type( DomainType::Cut );
            }
            this->set_block_and_sideset_names();


            comm_barrier();
        }

        void
        MaxwellFactory::set_block_and_sideset_names()
        {
            // create the format
            string tFormat = "%s_" + format_with_leading_zeros(
                mMesh->max_block_and_sideset_id() );

            for ( mesh::Block *tBlock : mMesh->blocks() )
            {
                if ( tBlock->domain_type() != DomainType::Default || tBlock->domain_type() == DomainType::UNDEFINED )
                {
                    string tName = string_to_lower(
                        to_string( tBlock->domain_type() ) );

                    tBlock->label() = sprint( tFormat.c_str(), tName.c_str(),
                                              ( uint ) tBlock->id() );
                }
            }
            for ( mesh::SideSet *tSideSet : mMesh->sidesets() )
            {
                if ( tSideSet->domain_type() != DomainType::Default || tSideSet->domain_type() == DomainType::UNDEFINED )
                {
                    string tName = string_to_lower(
                        to_string( tSideSet->domain_type() ) );

                    tSideSet->label() = sprint( tFormat.c_str(), tName.c_str(),
                                                ( uint ) tSideSet->id() );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::create_cuts_sub_slave()
        {
            BELFEM_ASSERT( mCommRank != 0, "This function should not be called on rank 0" );

            // in this section, we only need to create and run the factory
            // on the non-root ranks so they participate in the collective
            // mesh distribution.

            mesh::CutFactory tFactory( mMesh, mTopology, mProtoshells,
                                       mesh::CutAlgorithm::UNDEFINED,
                                       mUseEnrichment );

            comm_barrier();
            tFactory.run();

            comm_barrier();
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::create_thinshells()
        {
            if ( mCommRank == 0 )
            {
                if ( mMesh->thin_shells().size() != 0 )
                {
                    // shells came from the .bfm: don't recreate them, but still
                    // register their layer block->material so assign_materials
                    // can resolve and label them ( mirrors the create loop below )
                    for ( mesh::ThinShell * tShell : mMesh->thin_shells() )
                    {
                        const Cell< string > & tMaterials = tShell->materials();
                        uint b = 0;
                        for ( mesh::Block * tBlock : tShell->blocks() )
                        {
                            if ( ! mMaterialBlockAssignment.key_exists( tBlock->id() ) )
                            {
                                mMaterialBlockAssignment[ tBlock->id() ] = tMaterials( b );
                            }
                            ++b;
                        }

                        // edge-coating walls inherit the outer layer material
                        // ( mirrors the create loop below )
                        for ( mesh::Block * tBlock : tShell->side_connector_blocks() )
                        {
                            if ( ! mMaterialBlockAssignment.key_exists( tBlock->id() ) )
                            {
                                mMaterialBlockAssignment[ tBlock->id() ] = tMaterials.first();
                            }
                        }
                    }

                    comm_barrier();
                    return ;
                }
                // whether the layer interfaces get duplicate dofs and ghost
                // facets is the deck's call ( nitsche ghost penalty { eta } ),
                // read through the same helper the controller and the mesh
                // cache tag use
                const bool tGhostFacets = fem::ghost_facets_requested( *mInputFile );

                mesh::ThinShellFactory tFactory(
                    mMesh,
                    mThinShellMasterNodes,
                    mThinShellSlaveNodes,
                       & mMaterialMap,
                    tGhostFacets );

                if ( mProtoshells.size() > 0 && mMesh->has_periodicity() )
                {
                    tFactory.flag_periodic_nodes();
                }
                for ( Protoshell * tProtoshell : mProtoshells )
                {
                    mesh::ThinShell * tShell = tFactory.create( tProtoshell );

                    // add shell to map
                    for ( id_t tID : tProtoshell->sidesets() )
                    {
                        mThinShellMap[ tID ] = tShell;
                    }

                    tShell->set_label( tProtoshell->label() == ""
                        ? to_string( DomainType::ThinShell )
                        : tProtoshell->label() );

                    uint b = 0;
                    const Cell< string > &tMaterials = tShell->materials();
                    for ( mesh::Block *tBlock : tShell->blocks() )
                    {
                        // b tracks the layer index, so it must advance for every
                        // block even when the assignment is already set
                        if ( ! mMaterialBlockAssignment.key_exists( tBlock->id() ) )
                        {
                            mMaterialBlockAssignment[ tBlock->id() ] = tMaterials( b );
                        }
                        ++b;
                    }

                    // edge-coating walls inherit the outer layer material:
                    // the side coating deposits in the same plating step
                    // ( top == bottom enforced in the ThinShellFactory )
                    for ( mesh::Block * tBlock : tShell->side_connector_blocks() )
                    {
                        if ( ! mMaterialBlockAssignment.key_exists( tBlock->id() ) )
                        {
                            mMaterialBlockAssignment[ tBlock->id() ] = tMaterials.first();
                        }
                    }
                }
                if ( mProtoshells.size() > 0 )
                {
                    mMesh->unfinalize();
                    mMesh->finalize();

                    if ( mMesh->has_periodicity() )
                    {
                        tFactory.unflag_periodic_nodes();
                    }
                }
                this->find_autopins( mMesh->autopins() ) ;


            }
            comm_barrier();
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::create_terminal_list()
        {
            uint tNumConditions = 0;

            //Current boundary conditions are first
            for (PhysicalBoundaryCondition * tBC : mBoundaryConditionFactory->boundary_conditions())
            {
                if (tBC->type() == BoundaryConditionType::Current )
                {
                    mTerminals.push(Cell<id_t>()) ;
                    for (id_t tID : tBC->domains())
                    {
                        mTerminals(tNumConditions).push(tID);
                    }

                    if (tBC->is_thinshell())
                    {
                        //Keeping track of the indices that correspond to thin shell in mTerminals
                        mThinShellTerminalIndices.push(tNumConditions) ;
                    }
                    tNumConditions++;
                }
            }

            //Current conditions associated with circuit then
            for (PhysicalBoundaryCondition * tBC : mBoundaryConditionFactory->boundary_conditions())
            {
                if (tBC->type() == BoundaryConditionType::CircuitCurrent )
                {
                    mTerminals.push(Cell<id_t>()) ;
                    for (id_t tID : tBC->domains())
                    {
                        mTerminals(tNumConditions).push(tID);
                    }

                    if (tBC->is_thinshell())
                    {
                        //Keeping track of the indices that correspond to thin shell in mTerminals
                        mThinShellTerminalIndices.push(tNumConditions) ;
                    }
                    tNumConditions++;
                }
            }

            //Finally, voltage conditions
            for (PhysicalBoundaryCondition * tBC : mBoundaryConditionFactory->boundary_conditions())
            {
                if (tBC->type() == BoundaryConditionType::Voltage )
                {
                    mTerminals.push(Cell<id_t>()) ;
                    for (id_t tID : tBC->domains())
                    {
                        mTerminals(tNumConditions).push(tID);
                    }

                    if (tBC->is_thinshell())
                    {
                        //Keeping track of the indices that correspond to thin shell in mTerminals
                        mThinShellTerminalIndices.push(tNumConditions) ;
                    }
                    tNumConditions++;
                }
            }

            for ( PhysicalBoundaryCondition * tBC : mBoundaryConditionFactory->boundary_conditions() )
            {
                if (tBC->type() == BoundaryConditionType::CircuitVoltage )
                {
                    mTerminals.push(Cell<id_t>()) ;
                    for (id_t tID : tBC->domains())
                    {
                        mTerminals(tNumConditions).push(tID);
                    }

                    if (tBC->is_thinshell())
                    {
                        //Keeping track of the indices that correspond to thin shell in mTerminals
                        mThinShellTerminalIndices.push(tNumConditions) ;
                    }
                    tNumConditions++;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::fix_facet_masters()
        {
            BELFEM_ASSERT( mCommRank == 0, "This function should only be called on rank 0" );

            // same-type propagation is 3D-only: in 2D the sole consumer of
            // same-type facet orientation is the thin-shell cut pipeline,
            // which needs a UNIFORM master side per tape ( inherited from
            // the per-block element numbering ) — selective pairwise flips
            // can only break it. Cross-type facets are still normalized by
            // domain type below
            const bool tPropagate = mMesh->number_of_dimensions() == 3 ;

            bool tHaveFacetToFacet = mMesh->test_connectivity( Connectivity::FacetToFacet );

            if ( tPropagate && ! tHaveFacetToFacet )
            {
                mesh::ConnectivityCalculator tCalc( mMesh );
                tCalc.connect_facets_to_facets();
            }

            mMesh->unflag_all_facets();

            Queue< mesh::Facet * > tQueue;
            Cell< mesh::Facet * > & tFacets = mMesh->facets();

            for ( mesh::Facet * tFacet : tFacets )
            {
                if ( ! tFacet->has_slave() ) continue ;
                uint tM = static_cast< uint >( mMesh->block( tFacet->master()->block_id() )->domain_type() );
                uint tS = static_cast< uint >( mMesh->block( tFacet->slave()->block_id() )->domain_type() );

                if ( tM < tS )
                {
                    tQueue.push( tFacet );
                    tFacet->flip();
                }
                else if ( tM > tS )
                {
                    tQueue.push( tFacet );
                }
                else if ( tPropagate )
                {
                    tFacet->flag();
                }
            }

            while ( tPropagate && ! tQueue.empty() )
            {
                mesh::Facet * tFacet = tQueue.pop();

                for ( uint f=0; f<tFacet->number_of_facets(); ++f )
                {
                    mesh::Facet * tOther = tFacet->facet( f );

                    // skip facets that are not of interest
                    if ( ! tOther->is_flagged() ) continue;

                    // add other to Queue
                    tQueue.push( tOther );

                    tOther->unflag();

                    // check if facets have aligned normals
                    if ( mesh::check_facet_orientation( tFacet, tOther ) ) continue;
                    tOther->flip();
                }
            }

            // the flags are worklist scratch: facets with equal domain types on
            // both sides that were never reached by the queue are still flagged
            mMesh->unflag_all_facets();

            if ( tPropagate && ! tHaveFacetToFacet )
            {
                for ( mesh::Facet * tFacet : tFacets )
                {
                    tFacet->reset_facet_container();
                }
                mMesh->reset_connectivity( Connectivity::FacetToFacet );
            }
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::create_edges_and_faces_on_mesh()
        {
            if ( !mMesh->edges_exist() || !mMesh->faces_exist() )
            {
                if ( !mMesh->edges_exist() )
                {
                    mMesh->create_edges( false,
                        mTopology->groups( DomainType::Conductor ),
                                         Vector<id_t>(),
                                         false );
                }

                if ( !mMesh->faces_exist() && mMesh->max_element_order() > 1 )
                {
                    mMesh->create_faces( false,
                        mTopology->groups( DomainType::Conductor ) ,
                                         Vector<id_t>() );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::create_hanging_edges_and_facets()
        {
            // Initialize: clear any existing flags on edges and facets
            mMesh->unflag_all_edges();
            mMesh->unflag_all_facets( 0 );
            mMesh->unflag_all_facets( 1 );
            mMesh->unflag_all_facets( 2 );
            // Node containers for interface node matching
            Cell< mesh::Node * > tMasterNodes;   // Nodes from master side (typically conducting)
            Cell< mesh::Node * > tSlaveNodes;    // Nodes from slave side (reordered to match master)
            Cell< mesh::Node * > & tSlaveNodesTemp = tMasterNodes; // Temporary reuse of tMasterNodes buffer

            // Only linear elements supported currently
            BELFEM_ERROR( mMesh->max_element_order() == 1, "Not implemented for higher order" );

            // =========================================================================
            // PART 1: Process conductor-air and conductor-ferro interfaces (sidesets)
            // =========================================================================
            for ( mesh::SideSet * tSideSet : mMesh->sidesets() )
            {
                const DomainType tType = tSideSet->domain_type();

                // Handle InterfaceCondAir or InterfaceCondFerro (HPhi only)
                // These are interfaces where H-field DOFs exist on conductor side,
                // and we need to couple them to air/ferro nodes
                if ( tType == DomainType::InterfaceCondAir ||
                    ( tType == DomainType::InterfaceCondFerro &&
                        mFormulation == maxwell::Formulation::HPhi ) )
                {
                    Cell< mesh::Facet * > &tFacets = tSideSet->facets();
                    Cell< mesh::Edge * > tEdges;

                    for ( mesh::Facet *tFacet : tFacets )
                    {
                        // Ensure facet knows which element is master/slave
                        tFacet->compute_orientation();

                        // --- Node retrieval with proper orientation matching ---

                        // Step 1: Get slave nodes in their local element ordering
                        tFacet->slave()->get_nodes_of_facet(
                            tFacet->index_on_slave(), tSlaveNodesTemp );

                        // Step 2: Reorder slave nodes to match master facet's orientation
                        // This ensures node correspondence: tSlaveNodes[k] corresponds to tMasterNodes[k]
                        mesh::to_master_orientation(
                            tFacet, tSlaveNodesTemp, tSlaveNodes );

                        // Step 3: Get the master nodes (conducting side)
                        tFacet->master()->get_nodes_of_facet(
                            tFacet->index_on_master(), tMasterNodes );

                        // Get edges from the master side (these will become "hanging")
                        tFacet->master()->get_edges_of_facet(
                            tFacet->index_on_master(), tEdges );

                        // --- Build correspondence map: master node ID → slave node pointer ---
                        Map< id_t, mesh::Node * > tNodeMap;
                        for ( uint k = 0; k < tMasterNodes.size(); ++k )
                        {
                            // Use original() to handle periodic boundary conditions
                            tNodeMap[ tMasterNodes( k )->original()->id() ] =
                                tSlaveNodes( k );
                        }

                        // --- Create hanging edges ---
                        // Each edge on the master side gets "source" nodes from the slave side
                        for ( mesh::Edge *tEdge : tEdges )
                        {
                            if ( !tEdge->is_flagged() )  // Process each edge only once
                            {
                                // a periodic slave edge already hangs on its
                                // master-plane partner; keep that tie and let the
                                // t-matrix cascade resolve the chain to the
                                // interface nodes ( periodic wins on the slave rim )
                                if ( tEdge->number_of_sources() > 0 &&
                                     tEdge->source( 0 )->entity_type() == EntityType::EDGE )
                                {
                                    tEdge->flag();
                                    continue;
                                }

                                // Allocate storage for source nodes and their weights
                                tEdge->allocate_source_container(
                                    tEdge->number_of_nodes() );

                                // For each node on this master edge, find corresponding slave node
                                for ( uint k = 0; k < tEdge->number_of_nodes(); ++k )
                                {
                                    mesh::Node * tNode = tNodeMap(
                                        tEdge->node( k )->original()->id() ); // yes, original is correct here

                                    // note: we determine the weights later in
                                    // DofData::create_dofwise_t_matrices_master()
                                    tEdge->add_source( tNode );
                                }

                                tEdge->flag();  // Mark as processed
                            }
                        }

                        tFacet->flag();  // Mark facet as processed
                    }

                    // the h-side is coupled to the phi-side by the condensation
                    // above, not by a weak form, so this sideset must not become
                    // an active group. InterfaceCondAir never reaches one because
                    // it is not admitted to the magnetic equation, but
                    // InterfaceCondFerro is, so it needs the same explicit
                    // deactivation as InterfaceFerroAir below. HPhi is implied
                    // by the branch condition
                    if ( tType == DomainType::InterfaceCondFerro )
                    {
                        tSideSet->set_domain_type( DomainType::Inactive );
                    }
                }
                // Deactivate ThinShell sidesets (handled separately below)
                else if ( tType == DomainType::ThinShell )
                {
                    tSideSet->set_domain_type( DomainType::Inactive );
                }
                // Deactivate InterfaceFerroAir in HPhi formulation
                else if ( tType == DomainType::InterfaceFerroAir &&
                    mFormulation == maxwell::Formulation::HPhi )
                {
                    tSideSet->set_domain_type( DomainType::Inactive );
                }
            }

            // =========================================================================
            // PART 2: Process thin-shell elements
            // =========================================================================
#ifdef DEBUG
            if ( mMesh->thin_shells().size() > 0 )
            {
                for ( mesh::Node * tNode : mMesh->nodes() )
                {
                    tNode->set_index( gNoIndex );
                }
                for ( mesh::Edge * tEdge : mMesh->edges() )
                {
                    tEdge->set_index( gNoIndex );
                }
            }
#endif
            for ( mesh::ThinShell * tShell : mMesh->thin_shells() )
            {
                Cell< mesh::Facet * > & tFacets = tShell->facets();
                if ( tFacets.size() == 0 ) break;

                // Verify linear elements (TRI3 in 3D, LINE2 in 2D)
                BELFEM_ERROR( tFacets( 0 )->element()->type() == ElementType::TRI3 ||
                              tFacets( 0 )->element()->type() == ElementType::LINE2,
                              "Only linear triangular elements are supported" );

                DynamicBitset tBitset( tFacets.size() );

                // we begin with the upper ones
                index_t tCount = 0 ;
                for ( mesh::Facet * tFacet : tFacets )
                {
                    tFacet->compute_orientation();
                    if ( mMesh->block( tFacet->master()->block_id() )->domain_type() == DomainType::Conductor )
                    {
                        tBitset.set( tCount );
                    }
                    ++tCount ;
                }

                EdgeWorkData tWork ;

                Cell< index_t > tIndices ;
                Cell< mesh::Element * > & tElementsBottom = tShell->blocks().first()->elements() ;
                Cell< mesh::Element * > & tElementsTop = tShell->blocks().last()->elements() ;

                // Process air-master facets FIRST (node-based coupling).
                // At conductor/air transition boundaries, shared edges must
                // resolve to node DOFs (phi), not stay as independent H-edge
                // DOFs.  Processing air before conductor ensures transition
                // edges get edge-to-node coupling; pure conductor-interior
                // edges are never touched by the air path.
                tBitset.flip();
                tBitset.where( tIndices );

                for ( index_t k : tIndices )
                {
                    this->hang_thinshell_edges_on_nodes_bottom( tWork, tFacets( k ), tElementsBottom( k ) );
                }

                // Then process conductor-master facets (edge-based coupling).
                tBitset.flip();
                tBitset.where( tIndices );

                for ( index_t k : tIndices )
                {
                    this->hang_thinshell_edges_on_edges_bottom( tWork, tFacets( k ), tElementsBottom( k ) );
                }

                tBitset.reset();
                tCount = 0 ;
                for ( mesh::Facet * tFacet : tFacets )
                {
                    if ( mMesh->block( tFacet->slave()->block_id() )->domain_type() == DomainType::Conductor )
                    {
                        tBitset.set( tCount );
                    }
                    ++tCount ;
                }

                // Process air-slave facets FIRST (same reasoning as bottom)
                tBitset.flip();
                tBitset.where( tIndices );

                for ( index_t k : tIndices )
                {
                    this->hang_thinshell_edges_on_nodes_top( tWork, tFacets( k ), tElementsTop( k ) );
                }

                // Then conductor-slave facets
                tBitset.flip();
                tBitset.where( tIndices );

                for ( index_t k : tIndices )
                {
                    this->hang_thinshell_edges_on_edges_top( tWork, tFacets( k ), tElementsTop( k ) );
                }

            }

            // =========================================================================
            // PART 3: Higher-order facet processing (3D quadratic elements only)
            // =========================================================================
            if ( mMesh->number_of_dimensions() == 3 &&
                 mMesh->max_element_order() == 2 )
            {
                // Use TMatrix to compute proper interpolation weights for facets
                maxwell::TMatrix tTmatrix( mMesh );

                Cell< mesh::Node * > tNodes( 12, nullptr );

                for ( mesh::Facet *tFacet : mMesh->facets() )
                {
                    if ( tFacet->is_flagged() )  // Only process flagged facets from above
                    {
                        // Get slave nodes in local ordering (reusing tMasterNodes temporarily)
                        tFacet->slave()->get_nodes_of_facet(
                            tFacet->index_on_slave(), tMasterNodes );

                        // Reorder slave node      s to match master orientation
                        mesh::to_master_orientation(
                            tFacet, tSlaveNodes, tMasterNodes );

                        // Get actual master nodes
                        tFacet->master()->get_nodes_of_facet(
                            tFacet->index_on_master(), tMasterNodes );

                        // Collect all nodes (master + slave)
                        uint k = 0;
                        for ( mesh::Node *tNode : tMasterNodes )
                        {
                            tNodes( k++ ) = tNode;
                        }
                        for ( mesh::Node *tNode : tSlaveNodes )
                        {
                            tNodes( k++ ) = tNode;
                        }

                        // Compute interpolation weights using TMatrix
                        const Vector< real > &tWeights = tTmatrix.process( tFacet );

                        // Set weighted sources on facet
                        tFacet->set_sources( tNodes, tWeights );
                    }
                }
            }

            // Cleanup: unflag all edges and facets
            mMesh->unflag_all_edges( );
            mMesh->unflag_all_facets( 0 );
            mMesh->unflag_all_facets( 1 );
            mMesh->unflag_all_facets( 2 );
            mMesh->update_node_indices() ;
            mMesh->update_edge_indices() ;

        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::hang_thinshell_edges_on_nodes_bottom(
            EdgeWorkData  & aWork,
            mesh::Facet   * aFacet,
            mesh::Element * aElement )
        {
            Cell< mesh::Node * > & tNodesOnThinShell = aWork.NodesOnThinShell;
            Cell< mesh::Node * > & tNodesOnVolume    = aWork.NodesOnVolume;
            Cell< mesh::Edge * > & tEdgesOnThinShell = aWork.EdgesOnThinShell;

            // Shell-side nodes and edges on the bottom layer. The helpers
            // return them in the natural label order of the shell's bottom
            // nodes (0, 1, ...), which — under the thin-shell extrusion
            // convention — is position-wise aligned with the air-below's
            // top-facet ordering. This hides the CW/CCW asymmetry baked
            // into get_nodes_of_facet(bottom) (stored CW for the outward-
            // down normal in normal_penta's Jacobian formula).
            mesh::get_bottom_nodes( aElement, tNodesOnThinShell );

            // get volume nodes (index_on_master is correct here because
            // aFacet->master() is the volume element)
            aFacet->master()->get_nodes_of_facet( aFacet->index_on_master(), tNodesOnVolume );

            uint n = aFacet->number_of_nodes();
            for ( uint k=0; k<n; ++k )
            {
                tNodesOnThinShell( k )->set_index( k );
            }

            mesh::get_bottom_edges( aElement, tEdgesOnThinShell );

            for ( mesh::Edge * tEdge : tEdgesOnThinShell )
            {
                // check if edge has already been processed
                if ( tEdge->is_flagged() ) continue ;

                // a periodic slave edge already hangs on its master-plane
                // partner ( including the pure half-cut ties from
                // match_edges ); keep that tie and let the t-matrix cascade
                // resolve it — same guard as the InterfaceCondAir pass
                if ( tEdge->number_of_sources() > 0 &&
                     tEdge->source( 0 )->entity_type() == EntityType::EDGE )
                {
                    tEdge->flag();
                    continue;
                }

                tEdge->allocate_source_container( tEdge->number_of_nodes() );

                for ( uint k=0; k<tEdge->number_of_nodes(); ++k )
                {
                    mesh::Node * tNode = tNodesOnVolume( tEdge->node( k )->index() );

                    // note: we determine the weights later in
                    // DofData::create_dofwise_t_matrices_master()
                    tEdge->add_source( tNode );
                }

                // tag edge as processed
                tEdge->flag();
            }

#ifdef DEBUG
            // this will provoke an error if we mess things up
            for ( uint k=0; k<n; ++k )
            {
                tNodesOnThinShell( k )->original()->set_index( gNoIndex );
            }
#endif
        }

        void
        MaxwellFactory::hang_thinshell_edges_on_nodes_top(
            EdgeWorkData  & aWork,
            mesh::Facet   * aFacet,
            mesh::Element * aElement )
        {
            Cell< mesh::Node * > & tNodesOnThinShell = aWork.NodesOnThinShell;
            Cell< mesh::Node * > & tNodesOnVolume    = aWork.NodesOnVolume;
            Cell< mesh::Node * > & tTemporaryNodes   = aWork.NodesOnThinShell;
            Cell< mesh::Edge * > & tEdgesOnThinShell = aWork.EdgesOnThinShell;

            // Step 1: Get slave nodes in local ordering
            aFacet->slave()->get_nodes_of_facet(
                aFacet->index_on_slave(), tTemporaryNodes );

            // Step 2: Reorder slave nodes to match master orientation
            mesh::to_master_orientation(
                aFacet, tTemporaryNodes, tNodesOnVolume );

            // Step 3: Get shell nodes/edges on the top layer. The helpers
            // return them in the position-wise alignment that matches
            // tNodesOnVolume (master orientation, CCW top-facet convention).
            mesh::get_top_nodes( aElement, tNodesOnThinShell );

            uint n = aFacet->number_of_nodes();
            for ( uint k=0; k<n; ++k )
            {
                tNodesOnThinShell( k )->set_index( k );
            }

            mesh::get_top_edges( aElement, tEdgesOnThinShell );

            for ( mesh::Edge * tEdge : tEdgesOnThinShell )
            {
                // check if edge has already been processed
                if ( tEdge->is_flagged() ) continue ;

                // a periodic slave edge already hangs on its master-plane
                // partner ( including the pure half-cut ties from
                // match_edges ); keep that tie and let the t-matrix cascade
                // resolve it — same guard as the InterfaceCondAir pass
                if ( tEdge->number_of_sources() > 0 &&
                     tEdge->source( 0 )->entity_type() == EntityType::EDGE )
                {
                    tEdge->flag();
                    continue;
                }

                tEdge->allocate_source_container( tEdge->number_of_nodes() );

                for ( uint k=0; k<tEdge->number_of_nodes(); ++k )
                {
                    // note: we determine the weights later in
                    // DofData::create_dofwise_t_matrices_master()
                    mesh::Node * tNode = tNodesOnVolume( tEdge->node( k )->index() );
                    tEdge->add_source( tNode );
                }

                // tag edge as processed
                tEdge->flag();
            }

#ifdef DEBUG
            // this will provoke an error if we mess things up
            for ( uint k=0; k<n; ++k )
            {
                tNodesOnThinShell( k )->set_index( gNoIndex );
            }
#endif
        }

        void
        MaxwellFactory::hang_thinshell_edges_on_edges_bottom(
            EdgeWorkData  & aWork,
            mesh::Facet   * aFacet,
            mesh::Element * aElement )
        {
            Cell< mesh::Node * > & tNodesOnThinShell = aWork.NodesOnThinShell;
            Cell< mesh::Node * > & tNodesOnVolume    = aWork.NodesOnVolume;
            Cell< mesh::Edge * > & tEdgesOnThinShell = aWork.EdgesOnThinShell;
            Cell< mesh::Edge * > & tEdgesOnVolume    = aWork.EdgesOnVolume;

            // Shell-side nodes/edges on the bottom layer, in the position-wise
            // alignment that matches the air-below master's top-facet ordering
            // (hides the CW/CCW asymmetry baked into get_nodes_of_facet for
            // the outward-down normal).
            mesh::get_bottom_nodes( aElement, tNodesOnThinShell );

            // get volume nodes (index_on_master is correct here because
            // aFacet->master() is the volume element)
            aFacet->master()->get_nodes_of_facet( aFacet->index_on_master(), tNodesOnVolume );

            uint n = aFacet->number_of_nodes();
            uint m = aFacet->number_of_edges();

            for ( uint k=0; k<n; ++k )
            {
                tNodesOnThinShell( k )->set_index( k );
                tNodesOnVolume( k )->set_index( k );
            }

            mesh::get_bottom_edges( aElement, tEdgesOnThinShell );

            // get volume edges
            aFacet->master()->get_edges_of_facet( aFacet->index_on_master(), tEdgesOnVolume );

            for ( uint e=0; e<m; ++e )
            {
                tEdgesOnThinShell( e )->set_index( e );
            }

            real tWeight ;

            for ( mesh::Edge * tEdge : tEdgesOnThinShell )
            {


                // check if edge has already been processed
                if ( tEdge->is_flagged() ) continue ;

                // a periodic slave edge already hangs on its master-plane
                // partner ( including the pure half-cut ties from
                // match_edges ); keep that tie and let the t-matrix cascade
                // resolve it — same guard as the InterfaceCondAir pass
                if ( tEdge->number_of_sources() > 0 &&
                     tEdge->source( 0 )->entity_type() == EntityType::EDGE )
                {
                    tEdge->flag();
                    continue;
                }

                tEdge->allocate_source_container( 1 );

                mesh::Edge * tOther = tEdgesOnVolume( tEdge->index() );

                if ( tEdge->node( 0 )->index() == tOther->node( 0 )->index() &&
                     tEdge->node( 1 )->index() == tOther->node( 1 )->index() )
                {
                    tWeight = 1.0 ;
                }
                else if ( tEdge->node( 0 )->index() == tOther->node( 1 )->index() &&
                     tEdge->node( 1 )->index() == tOther->node( 0 )->index() )
                {
                    tWeight = -1.0 ;
                }
                else
                {
                    BELFEM_ERROR( false, "Could not determine edge orientation" );
                    tWeight = BELFEM_QUIET_NAN ;
                }
                tEdge->add_source( tOther, tWeight );

                tEdge->flag();
            }
#ifdef DEBUG
            // this will provoke an error if we mess things up
            for ( uint k=0; k<n; ++k )
            {
                tNodesOnThinShell( k )->set_index( gNoIndex );
                tNodesOnVolume( k )->set_index( gNoIndex );
            }

            // this will provoke an error if we mess things up
            for ( uint e=0; e<m; ++e )
            {
                tEdgesOnThinShell( e )->set_index( gNoIndex );
            }
#endif
        }


        void
        MaxwellFactory::hang_thinshell_edges_on_edges_top(
            EdgeWorkData  & aWork,
            mesh::Facet   * aFacet,
            mesh::Element * aElement )
        {
            Cell< mesh::Node * > & tNodesOnThinShell = aWork.NodesOnThinShell;
            Cell< mesh::Node * > & tNodesOnVolume    = aWork.NodesOnVolume;
            Cell< mesh::Node * > & tTemporaryNodes   = aWork.NodesOnThinShell;
            Cell< mesh::Edge * > & tEdgesOnThinShell = aWork.EdgesOnThinShell;
            Cell< mesh::Edge * > & tEdgesOnVolume    = aWork.EdgesOnVolume;
            Cell< mesh::Edge * > & tTemporaryEdges   = aWork.EdgesOnThinShell;

            uint n = aFacet->number_of_nodes();
            uint m = aFacet->number_of_edges();

            // Step 1: Get slave nodes in local ordering
            aFacet->slave()->get_nodes_of_facet(
                aFacet->index_on_slave(), tTemporaryNodes );

            // Step 2: Reorder slave nodes to match master orientation
            mesh::to_master_orientation(
                aFacet, tTemporaryNodes, tNodesOnVolume );

            // Step 3: Get shell nodes on the top layer (position-wise aligned
            // with tNodesOnVolume in master orientation).
            mesh::get_top_nodes( aElement, tNodesOnThinShell );

            for ( uint k=0; k<n; ++k )
            {
                tNodesOnThinShell( k )->set_index( k );
                tNodesOnVolume( k )->set_index( k );
            }

            // get volume edges
            aFacet->slave()->get_edges_of_facet( aFacet->index_on_slave(), tTemporaryEdges );
            to_master_orientation( aFacet, tTemporaryEdges, tEdgesOnVolume );

            // get shell edges on the top layer
            mesh::get_top_edges( aElement, tEdgesOnThinShell );

            for ( uint e=0; e<m; ++e )
            {
                tEdgesOnThinShell( e )->set_index( e );
            }

            real tWeight = 0.0 ;

            for ( mesh::Edge * tEdge : tEdgesOnThinShell )
            {
                //
                // check if edge has already been processed
                if ( tEdge->is_flagged() ) continue ;

                // a periodic slave edge already hangs on its master-plane
                // partner ( including the pure half-cut ties from
                // match_edges ); keep that tie and let the t-matrix cascade
                // resolve it — same guard as the InterfaceCondAir pass
                if ( tEdge->number_of_sources() > 0 &&
                     tEdge->source( 0 )->entity_type() == EntityType::EDGE )
                {
                    tEdge->flag();
                    continue;
                }

                tEdge->allocate_source_container( 1 );


                mesh::Edge * tOther = tEdgesOnVolume( tEdge->index() );

                if ( tEdge->node( 0 )->index() == tOther->node( 0 )->index() &&
                     tEdge->node( 1 )->index() == tOther->node( 1 )->index() )
                {
                    tWeight = 1.0 ;
                }
                else if ( tEdge->node( 0 )->index() == tOther->node( 1 )->index() &&
                          tEdge->node( 1 )->index() == tOther->node( 0 )->index() )
                {
                    tWeight = -1.0 ;
                }
                else
                {
                    BELFEM_ERROR( false, "Could not determine edge orientation" );
                    tWeight = BELFEM_QUIET_NAN ;
                }
                tEdge->add_source( tOther, tWeight );
                tEdge->flag();

            }

#ifdef DEBUG
            // this will provoke an error if we mess things up
            for ( uint k=0; k<n; ++k )
            {
                tNodesOnThinShell( k )->set_index( gNoIndex );
                tNodesOnVolume( k )->set_index( gNoIndex );
            }

            // this will provoke an error if we mess things up
            for ( uint e=0; e<m; ++e )
            {
                tEdgesOnThinShell( e )->set_index( gNoIndex );
            }
#endif
        }

        void
        MaxwellFactory::set_block_types_in_magnetic_equation()
        {
            // coil blocks are excluded: they are magnetically inert by design
            // ( no material, and the FieldList's Coil dof list is empty, so
            // they carry no dofs ). Handing them to the equation anyway made
            // the dof manager build blocks and calculators for them, and the
            // MaxwellData constructor then dereferenced their null material.
            // The kernel-parameter selection in create_magnetic_kernel already
            // excludes them; this closes the same hole on the equation side.
            // Only Coil is skipped -- the two lists are deliberately NOT
            // identical: the coating blocks, for instance, belong here even
            // though the kernel-parameter filter omits them
            index_t tNumBlocks = 0;
            for ( mesh::Block *tBlock : mMesh->blocks() )
            {
                if ( tBlock->domain_type() != DomainType::Coil )
                {
                    ++tNumBlocks;
                }
            }

            Vector< id_t > tIDs( tNumBlocks );
            Cell< DomainType > tTypes( tNumBlocks,
                                       DomainType::UNDEFINED );

            index_t tCount = 0;

            for ( mesh::Block *tBlock : mMesh->blocks() )
            {
                if ( tBlock->domain_type() == DomainType::Coil )
                {
                    continue;
                }
                tIDs( tCount ) = tBlock->id();
                tTypes( tCount++ ) = tBlock->domain_type();
            }
            mMagneticEquation->set_blocks( tIDs, tTypes );

            for ( uint s = 0; s < 2; ++s )
            {
                tCount = 0;
                for ( mesh::SideSet * tSideset : mMesh->sidesets() )
                {
                    switch ( tSideset->domain_type() )
                    {
                        case DomainType::AirPeriodic :
                        case DomainType::AirSymmetry :
                        case DomainType::BackgroundField :
                        case DomainType::AirAntiSymmetry :
                        case DomainType::BufferPeriodic :
                        case DomainType::BufferSymmetry :
                        case DomainType::BufferAntiSymmetry :
                        case DomainType::ConductorPeriodic :
                        case DomainType::ConductorSymmetry :
                        case DomainType::ConductorAntiSymmetry :
                        case DomainType::FerroPeriodic :
                        case DomainType::FerroSymmetry :
                        case DomainType::FerroAntiSymmetry :
                        case DomainType::InterfaceCondFerro :
                        case DomainType::InterfaceFerroAir :
                        case DomainType::ThinShell :
                        case DomainType::GeometryOnly :
                        case DomainType::Ghost :
                        {
                            if( s == 1 )
                            {
                                tIDs( tCount ) = tSideset->id();
                                tTypes( tCount ) = tSideset->domain_type();
                            }
                            ++tCount ;
                            break;
                        }
                        default:
                        {
                            // pass
                        }
                    }
                }

                if( s == 0 )
                {
                    tIDs.set_size( tCount );
                    tTypes.set_size( tCount, DomainType::UNDEFINED );
                }
            }

            mMagneticEquation->set_sidesets( tIDs, tTypes );
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::create_postprocessors()
        {
            if ( mCommRank == 0 )
            {
                // we can now disconnect the abstract nodes from the nodes
                mMesh->unflag_all_nodes();

                const Vector< id_t > & tAirBlocks = mTopology->groups( DomainType::Air );

                for ( id_t tID : tAirBlocks )
                {
                    mMesh->block( tID )->flag_nodes();
                }

                // Buffer blocks use the same scalar-phi DOFs as Air and
                // need their nodes flagged for hanging-node detection.
                const Vector< id_t > & tBufferBlocks = mTopology->groups( DomainType::Buffer );
                for ( id_t tID : tBufferBlocks )
                {
                    mMesh->block( tID )->flag_nodes();
                }

                if ( mFormulation == maxwell::Formulation::HPhi )
                {
                    const Vector< id_t > & tFerroBlocks = mTopology->groups( DomainType::Ferro );

                    for ( id_t tID : tFerroBlocks )
                    {
                        mMesh->block( tID )->flag_nodes();
                    }
                }
                Vector< real > tCoeffs;
                Cell< mesh::Node * > tSources;

                for ( mesh::Node *tNode : mMesh->nodes() )
                {
                    if ( tNode->is_hanging() )
                    {
                        uint tCount = 0;

                        // check for abstract nodes
                        for ( uint s = 0; s < tNode->number_of_sources(); ++s )
                        {
                            if ( tNode->source( s )->is_flagged() && tNode->
                                entity_type() == EntityType::NODE )
                            {
                                ++tCount;
                            }
                        }
                        tCoeffs.set_size( tCount );
                        tSources.set_size( tCount, nullptr );
                        tCount = 0;

                        // relink nodes without abstract nodes
                        for ( uint s = 0; s < tNode->number_of_sources(); ++s )
                        {
                            if ( tNode->source( s )->is_flagged() && tNode->
                                entity_type() == EntityType::NODE )
                            {
                                tSources( tCount ) = reinterpret_cast<
                                    mesh::Node * >( tNode->source( s ) );
                                tCoeffs( tCount++ ) = tNode->weight( s );
                            }
                        }

                        tNode->reset_source_container();
                        tNode->set_sources( tSources, tCoeffs );
                    }
                }
            }

            comm_barrier();

            const Vector< id_t > & tFerroBlocks = mTopology->groups( DomainType::Ferro );
            const Vector< id_t > & tCondBlocks = mTopology->groups( DomainType::Conductor );

            // removes connections at interface
            if ( tFerroBlocks.length() > 0 )
            {
                mMesh->unflag_all_nodes();
                for ( id_t tID : tFerroBlocks )
                {
                    if ( mMesh->block_exists( tID ) )
                    {
                        mMesh->block( tID )->flag_nodes();
                    }
                }
                for ( mesh::Node *tNode : mMesh->nodes() )
                {
                    if ( tNode->is_flagged() && tNode->is_hanging() )
                    {
                        tNode->reset_source_container();
                    }
                }
            }

            // the postprocessor constructor requires the fields existing

            Mesh * tMesh = mMagneticKernel->mesh();

            if ( !tMesh->field_exists( "Hx" ) )
                tMesh->create_field( "Hx" );
            if ( !tMesh->field_exists( "Hy" ) )
                tMesh->create_field( "Hy" );
            if ( tMesh->number_of_dimensions() == 3 && !tMesh->field_exists(
                "Hz" ) )
                tMesh->create_field( "Hz" );
            if ( !tMesh->field_exists( "Bx" ) )
                tMesh->create_field( "Bx" );
            if ( !tMesh->field_exists( "By" ) )
                tMesh->create_field( "By" );
            if ( tMesh->number_of_dimensions() == 3 && !tMesh->field_exists(
                "Bz" ) )
                tMesh->create_field( "Bz" );

            // note, this class also works, but the one above should be faster
            mMagneticField->postprocessors().push(
                       new MaxwellPostprocessor(
                           mMagneticKernel.get(),
                           mTopology->block_map(),
                           mMaterialBlockAssignment,
                           MaxwellPostprocessorType::Air ) );

            if ( tFerroBlocks.length() > 0 )
            {
                if ( mFormulation == maxwell::Formulation::HPhi )
                {
                    mMagneticField->postprocessors().push(
                        new MaxwellPostprocessor(
                            mMagneticKernel.get(),
                            mTopology->block_map(),
                            mMaterialBlockAssignment,
                            MaxwellPostprocessorType::Ferro ) );

                }
            }

            // the thin shell blocks are not part of the topology.
            // we must collect them manually
            Vector< id_t > tThinShellBlocks;

            // the edge coating walls ride along with their TRUE domain types
            // ( Left/RightCoating ): the types must travel explicitly since
            // non-root ranks rebuild ThinShells without the connector record
            Vector< id_t > tCoatingBlocks;
            Vector< id_t > tCoatingTypes;

            if ( mCommRank == 0 )
            {
                index_t tCount = 0 ;
                index_t tCountC = 0 ;

                for ( mesh::ThinShell * tShell : mMesh->thin_shells() )
                {
                    tCount += tShell->blocks().size() ;
                    tCountC += tShell->side_connector_blocks().size() ;
                }
                tThinShellBlocks.set_size( tCount );
                tCoatingBlocks.set_size( tCountC );
                tCoatingTypes.set_size( tCountC );
                tCount = 0 ;
                tCountC = 0 ;
                for ( mesh::ThinShell * tShell : mMesh->thin_shells() )
                {
                    for ( mesh::Block * tBlock : tShell->blocks() )
                    {
                        tThinShellBlocks( tCount++ ) = tBlock->id() ;
                    }
                    for ( mesh::Block * tBlock : tShell->side_connector_blocks() )
                    {
                        tCoatingBlocks( tCountC ) = tBlock->id() ;
                        tCoatingTypes( tCountC++ ) =
                            static_cast< id_t >( tBlock->domain_type() );
                    }
                }
            }

            // synch data with other procs
            broadcast( tThinShellBlocks );
            broadcast( tCoatingBlocks );
            broadcast( tCoatingTypes );

            // we need to add the thin shell types to the block map here.
            // The walls get their OWN map and their own postprocessor pass:
            // a Postprocessor instance is single-element-type by construction
            // ( compute_node_matrices sizes the patch matrices once from
            // mElementType ), so PENTA6TS layers and HEX8TB walls must never
            // share an instance
            Map< id_t, DomainType > tThinShellMap ;
             for ( id_t tID : tThinShellBlocks )
            {
                tThinShellMap[ tID ] = DomainType::ThinShell ;
            }

            Map< id_t, DomainType > tCoatingMap ;
            for ( index_t k=0; k<tCoatingBlocks.length(); ++k )
            {
                tCoatingMap[ tCoatingBlocks( k ) ] =
                    static_cast< DomainType >( tCoatingTypes( k ) );
            }

            // s = 0 : volume conductors, s = 1 : thin shell layers,
            // s = 2 : edge coating walls ( separate instance, see above )
            for ( uint s = 0; s < 3; ++s )
            {
                const Vector < id_t > & tIDs = s==0 ? tCondBlocks :
                                               s==1 ? tThinShellBlocks : tCoatingBlocks ;
                MaxwellPostprocessorType tC = s==0 ? MaxwellPostprocessorType::Conductor :
                                              s==1 ? MaxwellPostprocessorType::ThinShellConductor :
                                                     MaxwellPostprocessorType::SideConnector ;
                // walls are pure metal by the factory gate, so the s==2
                // superconductor count is always zero and tS never fires there
                MaxwellPostprocessorType tS = s==0 ? MaxwellPostprocessorType::SuperConductor : MaxwellPostprocessorType::ThinShellSuperConductor ;
                const Map< id_t, DomainType > & tMap = s==0 ? mTopology->block_map() :
                                                       s==1 ? tThinShellMap : tCoatingMap ;
                if ( tIDs.length() > 0 )
                {
                    // count conductor types
                    uint tNumNormalConductors = 0;
                    uint tNumSuperConductors = 0;

                    for ( id_t tID : tIDs )
                    {

                        if ( !mMagneticField->block_exists( tID ) )
                            continue;

                        if ( mMagneticKernel->material(
                            mMaterialBlockAssignment( tID ) )->have(MaterialProperty::jc) )
                        {
                            ++tNumSuperConductors;
                        }
                        else
                        {
                            ++tNumNormalConductors;
                        }
                    }

                    if ( comm_size() > 1 )
                    {
                        if ( mCommRank == 0 )
                        {
                            Vector< uint > tAllNumConductors;
                            Vector< uint > tAllNumSuperConductors;
                            collect( tAllNumConductors,
                                     tNumNormalConductors );
                            collect( tAllNumSuperConductors,
                                     tNumSuperConductors );

                            tNumNormalConductors = max( tAllNumConductors );
                            tNumSuperConductors = max( tAllNumSuperConductors );
                            tAllNumConductors.fill( tNumNormalConductors );
                            tAllNumSuperConductors.fill( tNumSuperConductors );

                            comm_barrier();
                            distribute( tAllNumConductors );
                            distribute( tAllNumSuperConductors );
                        }
                        else
                        {
                            send( tNumNormalConductors );
                            send( tNumSuperConductors );
                            comm_barrier();
                            receive( tNumNormalConductors );
                            receive( tNumSuperConductors );
                        }
                    }
                    if ( tNumNormalConductors > 0 )
                    {
                        mMagneticField->postprocessors().push(
                            new MaxwellPostprocessor(
                                mMagneticKernel.get(),
                                tMap,
                                mMaterialBlockAssignment,
                                tC ) );
                    }
                    if ( tNumSuperConductors > 0 )
                    {
                        mMagneticField->postprocessors().push(
                            new MaxwellPostprocessor(
                                mMagneticKernel.get(),
                                tMap,
                                mMaterialBlockAssignment,
                                tS ) );
                    }
                }
            }

            comm_barrier();
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::init_fields()
        {
            // initialize phi values
            Mesh * tMesh = mMagneticKernel->mesh() ;

            Vector< real > &tPhi = tMesh->field_exists( "phi" )
                ? tMesh->field_data( "phi" )
                : tMesh->create_field( "phi" );

            tMesh->unflag_all_nodes();

            const Vector< id_t > & tAirBlocks = mTopology->groups( DomainType::Air );
            const Vector< id_t > & tBufferBlocks = mTopology->groups( DomainType::Buffer );
            const Vector< id_t > & tFerroBlocks = mTopology->groups( DomainType::Ferro );

            for ( id_t tID : tAirBlocks )
            {
                if ( tMesh->block_exists( tID ) )
                {
                    tMesh->block( tID )->flag_nodes();
                }
            }

            // Buffer blocks carry scalar-phi DOFs and need 0.0 init,
            // otherwise their nodes start at BELFEM_QUIET_NAN.
            for ( id_t tID : tBufferBlocks )
            {
                if ( tMesh->block_exists( tID ) )
                {
                    tMesh->block( tID )->flag_nodes();
                }
            }

            if ( mFormulation == maxwell::Formulation::HPhi )
            {
                for ( id_t tID : tFerroBlocks )
                {
                    if ( tMesh->block_exists( tID ) )
                    {
                        tMesh->block( tID )->flag_nodes();
                    }
                }
            }

            for ( mesh::Node * tNode : tMesh->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    tPhi( tNode->index() ) = 0.0 ;
                }
                else
                {
                    tPhi( tNode->index() ) = BELFEM_QUIET_NAN ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::configure_solver( const input::Section *aSection,
                                          DofManager *aField )
        {
            // read solver settings
            SolverParameters tParams( aSection );

            // set the solver of the field
            aField->set_solver( tParams );
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::create_block_to_material_map()
        {
            Vector< id_t > tBlockIDs;
            Cell< string > tMaterials;

            this->collect_material_labels_from_domains( tBlockIDs, tMaterials );

            index_t tCount = 0;
            mMaterialBlockAssignment.clear();

            for ( id_t tID : tBlockIDs )
            {
                mMaterialBlockAssignment[ tID ] = tMaterials( tCount );
                ++tCount;
            }
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::set_physical_tags_for_elements()
        {
            Cell< string > tMaterials;
            for ( auto tPair : mMaterialBlockAssignment )
            {
                tMaterials.push(  tPair.second );
            }
            unique( tMaterials );

            id_t tCount = 1 ;

            for ( const string & tMaterial : tMaterials )
            {
                mMaterialIDs[ tMaterial ] = tCount++ ;
            }

            for ( auto tPair : mMaterialBlockAssignment )
            {
                Cell< mesh::Element * > & tElements = mMesh->block( tPair.first )->elements() ;

                for ( mesh::Element * tElement : tElements )
                {
                    tElement->set_physical_tag( mMaterialIDs( tPair.second) );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::collect_material_labels_from_domains(
            Vector< id_t > &aBlockIDs,
            Cell< string > &aMaterialLabels )
        {
            // count memory needs
            index_t tCount = 0;
            for ( Domain *tDomain : mDomains )
            {
                if ( tDomain->type() == DomainType::ThinShell ) continue ;
                tCount += tDomain->group_ids().length();
            }

            aBlockIDs.set_size( tCount );
            aMaterialLabels.set_size( tCount, "" );
            tCount = 0;
            for ( Domain *tDomain : mDomains )
            {
                switch ( tDomain->type() )
                {
                    case DomainType::Air:
                    {
                        for ( uint k = 0; k < tDomain->group_ids().length(); ++k )
                        {
                            aBlockIDs( tCount ) = tDomain->group_ids()( k );
                            aMaterialLabels( tCount++ ) = "air";
                        }
                        break ;
                    }
                    case DomainType::Conductor:
                    case DomainType::Ferro:
                    case DomainType::Buffer :
                    {
                        for ( uint k = 0; k < tDomain->group_ids().length(); ++k )
                        {
                            aBlockIDs( tCount ) = tDomain->group_ids()( k );
                            aMaterialLabels( tCount++ ) = tDomain->material();
                        }
                        break ;
                    }
                    case DomainType::ThinShell :
                    {
                        // pass
                        break ;
                    }
                    default:
                    {
                        for ( uint k = 0; k < tDomain->group_ids().length(); ++k )
                        {
                            aBlockIDs( tCount ) = tDomain->group_ids()( k );
                            aMaterialLabels( tCount++ ) = "inactive";
                        }
                        break ;
                    }
                }
            }

        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::synch_material_map()
        {
            if ( mCommRank == 0 )
            {
                Vector< id_t > tBlockIDs( mMaterialBlockAssignment.size() );
                Cell< string > tMaterialLabels( mMaterialBlockAssignment.size(),
                                                "" );

                index_t tCount = 0;

                for ( auto tPair : mMaterialBlockAssignment )
                {
                    tBlockIDs( tCount ) = tPair.first;
                    tMaterialLabels( tCount++ ) = tPair.second;
                }
                comm_barrier();
                broadcast( tBlockIDs );
                broadcast( tMaterialLabels );
            }
            else
            {
                comm_barrier();
                Vector< id_t > tBlockIDs;
                Cell< string > tMaterialLabels;
                broadcast( tBlockIDs );
                broadcast( tMaterialLabels );

                index_t tCount = 0;
                mMaterialBlockAssignment.clear();
                for ( id_t tID : tBlockIDs )
                {
                    mMaterialBlockAssignment[ tID ] = tMaterialLabels(
                        tCount++ );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::create_materials()
        {
            if ( mInputFile->section_exists( "materials" ) )
            {
                MaterialFactory tFactory( mInputFile->section( "materials" ) );
                mMaterialMap = tFactory.materials() ;
            }
        }

        void
        MaxwellFactory::assign_materials()
        {
            // first we collect all materials that we need
            Cell< string > tMaterialLabels ;

            Mesh * tMesh = mMagneticKernel->mesh() ;

            for ( auto tPair : mMaterialBlockAssignment )
            {
                if (tPair.second == "inactive" || ! tMesh->block_exists( tPair.first ) ) continue ;
                switch ( tMesh->block( tPair.first )->domain_type() )
                {
                    case DomainType::Conductor :
                    case DomainType::ThinShell :
                    case DomainType::Ferro :
                    case DomainType::Buffer :
                    {
                        // Buffer blocks use the phi formulation in the
                        // magnetic kernel and never query mat->rho, but the
                        // material (Magnesia by default) still needs to be
                        // registered so the thermal kernel can look it up
                        // by label later.
                        tMaterialLabels.push( tPair.second );
                        break;
                    }
                    case DomainType::LeftCoating :
                    case DomainType::RightCoating :
                    {
                        // edge-coating wall: the kernel evaluates the generic
                        // rho( T, B, beta ) law, so the material must be a
                        // pure metal. Gate on the TYPE - YBCO subclasses
                        // Metal but is MaterialType::HTS
                        Material * tMat = mMaterialMap( tPair.second );

                        BELFEM_ERROR( tMat->type() == MaterialType::PureMetal,
                            "edge coating material %s on block %lu must be a pure metal",
                            tPair.second.c_str(),
                            ( long unsigned int ) tPair.first );

                        if (    string_to_lower( tPair.second ) != "copper"
                             && string_to_lower( tPair.second ) != "cu" )
                        {
                            message( InfoLevel::Default,
                                "    Warning: edge coating material %s is not copper, but surround plating is usually Cu",
                                tPair.second.c_str() );
                        }

                        tMaterialLabels.push( tPair.second );
                        break;
                    }
                    default:
                        break;
                }
            }
            unique( tMaterialLabels );

            for ( const string &tLabel : tMaterialLabels )
            {
                Material * tMat = mMaterialMap( tLabel ) ;
                tMat->flag( 7 );
                // tagging the material so that the
                // destructor doesn't delete it
                mMagneticKernel->add_material( tLabel, tMat );
            }

            for ( uint d = 0; d < mMagneticKernel->number_of_dof_managers(); ++
                  d )
            {
                Cell< Block * > &tBlocks = mMagneticKernel->dofmgr( d )->
                    blocks();

                for ( Block *tBlock : tBlocks )
                {
                    // get material label
                    const string &tMaterialLabel = mMaterialBlockAssignment(
                        tBlock->id() );
                    if ( tMaterialLabel != "air" && tMaterialLabel !=
                        "inactive" && tMesh->block_exists( tBlock->id() ) )
                    {
                        tBlock->set_material(
                            mMaterialBlockAssignment( tBlock->id() ) );
                    }
                }
            }

            // A thin shell next to an h-conductor recovers its normal field
            // from the conductor's own Nedelec trace ( compute_h_trace ),
            // which equals the shell's normal field only where [ n . B ] = 0
            // reduces to n . H continuity, i.e. for a nonmagnetic conductor.
            // Same predicate as MaxwellData's constant-mu0 branch. Setup
            // code, runs once: always-active check. The thin-shell records
            // live on the master mesh.
            if ( mCommRank == 0 )
            {
                for ( mesh::ThinShell * tShell : tMesh->thin_shells() )
                {
                    for ( mesh::Facet * tFacet : tShell->facets() )
                    {
                        mesh::Element * tSides[ 2 ] = { tFacet->master(), tFacet->slave() };

                        for ( mesh::Element * tVolume : tSides )
                        {
                            if ( tVolume == nullptr ) continue ;

                            const id_t tBlockID = tVolume->block_id() ;

                            if ( tMesh->block( tBlockID )->domain_type() != DomainType::Conductor ) continue ;

                            /* commenting these errors because they thow on Lawrencium without cause
                            // The recovery gathers one edge dof per edge of the
                            // conductor element ( nedelec_data_master_h /
                            // _slave_h ) into a vector sized for ALL of its
                            // Nedelec dofs. A higher-order conductor ( TET10,
                            // TRI6: face and second edge dofs ) would leave the
                            // tail of that vector stale, and E * q would read
                            // it. The length check inside the gather is a
                            // debug-only assert; this is the always-active
                            // rejection of the configuration ( no thin-shell
                            // element above linear exists today )
                            BELFEM_ERROR( mesh::interpolation_order( tVolume->type() ) == InterpolationOrder::LINEAR,
                                "conductor block %lu ( %s ) is a volume neighbor of thin shell %s but is not linear ( %s ): the shell's normal-field recovery reads one edge dof per edge of the conductor element",
                                ( long unsigned int ) tBlockID,
                                mMaterialBlockAssignment( tBlockID ).c_str(),
                                tShell->label().c_str(),
                                to_string( tVolume->type() ).c_str() );

                            const Material * tMat = mMagneticField->block( tBlockID )->material() ;

                            BELFEM_ERROR( tMat != nullptr
                                && tMat->is_constant( MaterialProperty::mu )
                                && std::abs( tMat->constant_property( MaterialProperty::mu ) - constant::mu0 ) < BELFEM_EPSILON,
                                "conductor block %lu ( %s ) is a volume neighbor of thin shell %s but is not nonmagnetic: the shell's normal-field recovery needs mu = mu0 on an h-conductor side",
                                ( long unsigned int ) tBlockID,
                                mMaterialBlockAssignment( tBlockID ).c_str(),
                                tShell->label().c_str() ); */
                        }
                    }
                }
            }

            // rename mesh blocks based on material assignment.
            // iterate over all material assignments (which includes all
            // thin shell layer blocks) so that thin shells from every
            // tape get labeled, not just the ones present in the dof
            // manager block list.
            for ( auto tPair : mMaterialBlockAssignment )
            {
                const id_t tID = tPair.first ;
                const string & tMaterialLabel = tPair.second ;

                if ( tMaterialLabel == "air" || tMaterialLabel == "inactive" )
                    continue ;

                if ( ! tMesh->block_exists( tID ) )
                    continue ;

                mesh::Block * tMeshBlock = tMesh->block( tID );
                if ( tMeshBlock == nullptr )
                    continue ;

                string tDomain = string_to_lower(
                        to_string( tMeshBlock->domain_type() ) );
                tMeshBlock->label() = sprint(
                        "%s_%02u_%s",
                        tDomain.c_str(),
                        ( unsigned int ) tMeshBlock->id(),
                        tMaterialLabel.c_str() );
            }
        }

        void
        MaxwellFactory::delete_unused_materials()
        {
            for ( auto tPair : mMaterialMap )
            {
                if ( ! tPair.second->is_flagged( 7 ) )
                {
                    delete tPair.second;
                }
            }
            mMaterialMap.clear();
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::read_signed_sidesets(
                const input::Section * aSection,
                const string         & aKey,
                Protoshell           * aShell )
        {
            // gmsh-style signed sideset list, e.g.
            //     sidesets : -5, -6, 7:20 ;
            // A NEGATIVE id requests that the orientation of every facet of
            // that sideset be flipped after the master normalization ( see
            // flip_thin_shell_sidesets ). This exists because the layer stack
            // is laid along the facet normal, and the domain-type master rule
            // is mirror-symmetric: a tape stack bounded by air on both faces
            // cannot come out uniform — one outer tape always flips. Which
            // side the layers face is user intent that no geometry rule can
            // derive ( a corc wrap has no meaningful mean normal at all ), so
            // the sign carries it per sideset, and an unsigned deck behaves
            // exactly as before.
            Cell< string > tWords = string_to_words( search_and_replace(
                aSection->get_string( aKey ), " ", "" ), ',' );

            Cell< id_t > tIDs ;
            Cell< id_t > tFlipped ;

            Vector< id_t > tExpanded ;

            for ( string & tWord : tWords )
            {
                if ( tWord.length() > 0 && tWord[ 0 ] == '-' )
                {
                    const string tNumber = tWord.substr( 1 );

                    // inside a range the scope of the sign is ambiguous
                    // ( -5:8 ? -5:-8 ? ), so it is refused rather than guessed
                    BELFEM_ERROR( tNumber.find( ":" ) == string::npos,
                        "thin shell %s, key %s: the sign must bind to a single id, not to the range %s — list the flipped ids individually",
                        aShell->label().c_str(),
                        aKey.c_str(),
                        tWord.c_str() );

                    BELFEM_ERROR( is_integer( tNumber ),
                        "thin shell %s, key %s: %s is not a valid signed sideset id",
                        aShell->label().c_str(),
                        aKey.c_str(),
                        tWord.c_str() );

                    const id_t tID = ( id_t ) std::stoi( tNumber );
                    tIDs.push( tID );
                    tFlipped.push( tID );
                }
                else
                {
                    // plain id or a:b range: the standard expansion applies
                    aSection->ids_from_string( tWord, tExpanded );
                    for ( index_t k = 0; k < tExpanded.length(); ++k )
                    {
                        tIDs.push( tExpanded( k ) );
                    }
                }
            }

            // hand over as vectors; sidesets() carries the absolute values,
            // so every consumer downstream stays sign-agnostic
            aShell->sidesets().set_size( tIDs.size() );
            for ( index_t k = 0; k < tIDs.size(); ++k )
            {
                aShell->sidesets()( k ) = tIDs( k );
            }

            aShell->flipped_sidesets().set_size( tFlipped.size() );
            for ( index_t k = 0; k < tFlipped.size(); ++k )
            {
                aShell->flipped_sidesets()( k ) = tFlipped( k );
            }
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::flip_thin_shell_sidesets()
        {
            BELFEM_ASSERT( mCommRank == 0,
                "This function should only be called on rank 0" );

            for ( Protoshell * tShell : mProtoshells )
            {
                const Vector< id_t > & tFlipped = tShell->flipped_sidesets() ;

                for ( index_t k = 0; k < tFlipped.length(); ++k )
                {
                    mesh::SideSet * tSideSet = mMesh->sideset( tFlipped( k ) );

                    for ( mesh::Facet * tFacet : tSideSet->facets() )
                    {
                        // a thin-shell facet always sits between two volume
                        // elements; without a slave there is nothing to flip
                        // to and the model is broken anyway
                        BELFEM_ERROR( tFacet->has_slave(),
                            "thin shell %s: cannot flip facet %lu of sideset %lu, it has no slave element",
                            tShell->label().c_str(),
                            ( long unsigned int ) tFacet->id(),
                            ( long unsigned int ) tSideSet->id() );

                        tFacet->flip() ;
                    }

                    message( InfoLevel::Default,
                        "    thin shell %s: flipped orientation of sideset %lu ( %lu facets )",
                        tShell->label().c_str(),
                        ( long unsigned int ) tSideSet->id(),
                        ( long unsigned int ) tSideSet->number_of_facets() );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        MaxwellFactory::read_thin_shell_data()
        {
            // get the topology section in the input file
            const input::Section *tTopo = mInputFile->section( "topology" );

            Map< string, string > tTapeMap;

            // count number of thin shells
            uint tCount = 0;
            for ( uint k = 0; k < tTopo->num_sections(); ++k )
            {
                const input::Section *tSection = tTopo->section( k );
                // gate on the domain type, not the spelling: domain_type()
                // also accepts "tape" and "shell", and comparing the raw
                // string here made those sections build no shell at all
                if ( domain_type( tSection->type() ) == DomainType::ThinShell )
                {
                    tTapeMap[ tSection->label() ] = tSection->
                        key_exists( "sideset" )
                        ? tSection->get_string( "sideset" )
                        : tSection->get_string( "sidesets" );
                    ++tCount;
                }
            }

            if ( tCount == 0 )
                return;

            tCount = 0;

            // get sections that describe layers
            Map< string, const input::Section * > tLayers;
            for ( uint k = 0; k < mInputFile->num_sections(); ++k )
            {
                const input::Section *tSection = mInputFile->section( k );

                if ( tSection->type() == "layers" )
                {
                    tLayers[ tSection->label() ] = tSection;
                }
            }

            // read the topology data
            id_t tID = 0;
            for ( uint j = 0; j < tTopo->num_sections(); ++j )
            {
                const input::Section *tSection = tTopo->section( j );
                // see the note on the counting loop above
                if ( domain_type( tSection->type() ) == DomainType::ThinShell )
                {
                    Protoshell *tShell = new Protoshell( ++tID );
                    tShell->label() = tSection->label();

                    if ( tSection->key_exists( "sideset" ) xor tSection->
                        key_exists( "sidesets" ) )
                    {
                        if ( tSection->key_exists( "sideset" ) )
                        {
                            this->read_signed_sidesets( tSection, "sideset", tShell );
                        }
                        else if ( tSection->key_exists( "sidesets" ) )
                        {
                            this->read_signed_sidesets( tSection, "sidesets", tShell );
                        }
                        else
                        {
                            BELFEM_ERROR(
                                false,
                                "need sideset or sidesets defined for thin shell %s, not both.",
                                tSection->label().c_str() );
                        }

                        // read the layer data
                        BELFEM_ERROR( tLayers.key_exists( tSection->label() ),
                                      "Layers for %s not defined.",
                                      tShell->label().c_str() );

                        const input::Section *tLayer = tLayers[ tShell->
                            label() ];

                        uint tNumLayers = tLayer->num_keys();
                        tShell->thicknesses().set_size( tNumLayers, 0.0 );
                        tShell->materials().set_size( tNumLayers, "" );

                        //Read the material and thicknesses
                        const Cell< string > & tBuffer = tLayer->buffer();
                        tCount = 0 ;
                        for ( index_t k = tLayer->start(); k<tLayer->end(); ++k )
                        {
                            Cell< string > tWords = string_to_words(  tBuffer( k ) );

                            const string & tStrMat = tWords( 0 );
                            real tThickness = to_real( tWords( 2 ) );

                            value tUnit = unit_to_si( tWords( 3 ) );

                            BELFEM_ERROR( check_unit( tUnit, "m" ),
                            "Invalid thickness unit for %s, expect %s or same dimension",
                            tStrMat.c_str(), "m" );

                            tThickness *= tUnit.first ;

                            // A layer thickness must be positive, and this is the
                            // place to say so: the deck line is still in hand. The
                            // downstream guard is Kernel::compute_element_volumes,
                            // a whole setup phase later, and all it can report is
                            // "elements with negative volume on mesh" -- true, but
                            // it names neither the layer nor the value. A negative
                            // thickness inverts the layer stack and flips the sign
                            // of the QUAD4TS Jacobian.
                            //
                            // Phrased as a positive test on purpose: to_real hands
                            // back NaN for a malformed number, and NaN > 0 is false,
                            // so that is refused here too instead of propagating
                            BELFEM_ERROR( tThickness > 0.0,
                                "Layer thickness for %s in layer stack %s must be positive, but is %s %s",
                                tStrMat.c_str(),
                                tShell->label().c_str(),
                                tWords( 2 ).c_str(),
                                tWords( 3 ).c_str() );

                            tShell->thicknesses()( tCount ) = tThickness ;
                            tShell->materials()( tCount++ ) = tStrMat ;

                        }

                        // edge coating: opt-in surround plating walls on the
                        // tape slit edges ( 3D only ); wall material and width
                        // derive from the outer stabilizer layer unless the
                        // width key overrides
                        if ( tSection->key_exists( "edge coating" ) )
                        {
                            tShell->edge_coating() =
                                tSection->get_bool( "edge coating" );

                            BELFEM_ERROR( ! tShell->edge_coating()
                                || mMesh->number_of_dimensions() == 3,
                                "edge coating on thin shell %s requires a 3D mesh",
                                tShell->label().c_str() );
                        }
                        if ( tSection->key_exists( "edge coating width" ) )
                        {
                            BELFEM_ERROR( tShell->edge_coating(),
                                "edge coating width is set on thin shell %s, but edge coating is not switched on",
                                tShell->label().c_str() );

                            tShell->edge_coating_width() =
                                tSection->get_value( "edge coating width", "m" ).first ;

                            BELFEM_ERROR( tShell->edge_coating_width() > 0.0,
                                "edge coating width on thin shell %s must be positive",
                                tShell->label().c_str() );
                        }

                        //Associate the open curves to the terminals of the thin shell
                        for ( id_t tSideSetID : tShell->sidesets() )
                        {
                            for ( mesh::Curve * tCurve : mMesh->curves() )
                            {
                                if ( mMesh->number_of_dimensions() == 2)
                                {
                                    if ( tSideSetID == tCurve->sideset_a()->id() )
                                    {
                                        tShell->terminal_curves().push(tCurve) ;
                                    }
                                }
                                else
                                {
                                    if (tCurve->sideset_a()->id() == tSideSetID )
                                    {
                                        //Reverse SideSetA and SideSetB such that SideSetB is always the thin shell
                                        id_t tIDa =  tCurve->sideset_a()->id() ;
                                        id_t tIDb =  tCurve->sideset_b()->id() ;
                                        tCurve->sideset_a( mMesh->sideset( tIDb ) );
                                        tCurve->sideset_b( mMesh->sideset( tIDa ) );

                                        tShell->terminal_curves().push(tCurve) ;
                                    }
                                    else if (tCurve->sideset_b()->id() == tSideSetID)
                                    {
                                        tShell->terminal_curves().push(tCurve) ;
                                    }
                                }
                            }
                        }
                        mProtoshells.push( tShell );
                    }
                    else
                    {
                        BELFEM_ERROR(
                            false, "no sidesets defined for thin shell %s.",
                            tSection->label().c_str() );
                    }
                }
            }

            // set sideset types in mesh
            for ( auto tPair : mProtoshells )
            {
                for ( id_t s : tPair->sidesets() )
                {
                    mMesh->sideset( s )->set_domain_type( DomainType::ThinShell );
                }
            }

        }


//------------------------------------------------------------------------------

        void
        MaxwellFactory::find_autopins( Cell< mesh::Node * > & aPins )
        {
            mMesh->unflag_all_nodes();
            mMesh->unflag_all_elements();
            mMesh->update_node_indices();

            Cell< graph::Vertex * > tGraph ;

            // tag all nodes that sit on air domains or buffer layers
            for ( auto tBlock : mMesh->blocks() )
            {
                if (    tBlock->domain_type() == DomainType::Buffer
                    ||  tBlock->domain_type() == DomainType::Air
                    ||  tBlock->domain_type() == DomainType::Ferro )
                {
                    {
                        Cell< mesh::Element * > & tElements = tBlock->elements() ;

                        for ( auto tElement : tElements )
                        {
                            tElement->flag();
                            for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                            {
                                mesh::Node * tNode = tElement->node( k )->original() ;
                                if ( ! tNode->is_flagged() )
                                {
                                    tNode->flag() ;
                                    tGraph.push( tNode );
                                }
                            }
                        }
                    }
                }
            }

            // cancel if no points were found ( a pure h-conductor deck has
            // no phi region; the tail reads tGraph.last() )
            if ( tGraph.size() == 0 ) { return; }

            tGraph.shrink_to_fit();
            sort( tGraph, opVertexIndex );

            for ( auto tNode : tGraph )
            {
                tNode->unflag();
            }

            Graph tAdj ;

            for ( auto tVertex : tGraph )
            {
                mesh::Node * tNode = reinterpret_cast< mesh::Node * >( tVertex );

                tAdj.clear();
                tNode->flag();

                for ( uint e=0; e<tNode->number_of_elements(); ++e )
                {
                    mesh::Element * tElement = tNode->element( e );

                    if ( ! tElement->is_flagged() ) continue;

                    for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                    {
                        mesh::Node * tOther = tElement->node( k )->original() ;

                        if ( ! tOther->is_flagged() )
                        {
                            tOther->flag();
                            tAdj.push( tOther );
                        }
                    }
                }

                for ( uint d=0; d<tNode->number_of_duplicates(); ++d )
                {
                    mesh::Node * tDup = tNode->duplicate( d );
                    for ( uint e=0; e<tDup->number_of_elements(); ++e )
                    {
                        mesh::Element * tElement = tDup->element( e );

                        if ( ! tElement->is_flagged() ) continue;

                        for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                        {
                            mesh::Node * tOther = tElement->node( k )->original() ;

                            if ( ! tOther->is_flagged() )
                            {
                                tOther->flag();
                                tAdj.push( tOther );
                            }
                        }
                    }
                }

                if ( tNode->is_periodic() )
                {
                    mesh::Node * tPeriodic = tNode->periodic()->original();

                    for ( uint e=0; e<tPeriodic->number_of_elements(); ++e )
                    {
                        mesh::Element * tElement = tPeriodic->element( e );

                        if ( ! tElement->is_flagged() ) continue;

                        for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                        {
                            mesh::Node * tOther = tElement->node( k )->original() ;

                            if ( ! tOther->is_flagged() )
                            {
                                tOther->flag();
                                tAdj.push( tOther );
                            }
                        }
                    }

                    for ( uint d=0; d<tPeriodic->number_of_duplicates(); ++d )
                    {
                        mesh::Node * tDup = tPeriodic->duplicate( d );
                        for ( uint e=0; e<tDup->number_of_elements(); ++e )
                        {
                            mesh::Element * tElement = tDup->element( e );

                            if ( ! tElement->is_flagged() ) continue;

                            for ( uint k=0; k<tElement->number_of_nodes(); ++k )
                            {
                                mesh::Node * tOther = tElement->node( k )->original() ;

                                if ( ! tOther->is_flagged() )
                                {
                                    tOther->flag();
                                    tAdj.push( tOther );
                                }
                            }
                        }
                    }
                }

                sort( tAdj, opVertexIndex );

                tVertex->init_vertex_container( tAdj.size() );

                for ( auto tOther : tAdj )
                {
                    tOther->unflag();
                    tVertex->insert_vertex( tOther );
                }
                tVertex->unflag();
            }

            // element flags were scratch for the adjacency
            mMesh->unflag_all_elements();

            sort( tGraph , opVertexID );

            // backup ownerships and levels, just in case.
            // two-arg ctor: the one-arg Cell ctor only RESERVES ( size
            // stays 0 and the indexed writes below would be out of range )
            Cell< proc_t >  tOwners( tGraph.size(), 0 );
            Cell< index_t > tLevels( tGraph.size(), 0 );
            index_t tCount = 0 ;
            for ( auto tVertex : tGraph )
            {
                tOwners( tCount )   = tVertex->owner();
                tLevels( tCount++ ) = tVertex->level();
            }

            graph::find_connected_partitions( tGraph );
            index_t tNumPins = tGraph.last()->owner() + 1 ;
            aPins.set_size( tNumPins, nullptr );

            tCount = 0 ;
            proc_t tLastOwner = gNoOwner ;
            for ( auto tVertex : tGraph )
            {
                if ( tVertex->owner() != tLastOwner )
                {
                    aPins( tCount++ ) = reinterpret_cast< mesh::Node * >( tVertex );
                    if ( tCount == tNumPins )
                    {
                        break ;
                    }
                    tLastOwner = tVertex->owner();
                }
            }

            mMesh->unflag_all_nodes();
            for ( auto tSideSet : mMesh->sidesets() )
            {
                tSideSet->flag_all_nodes();
            }

            tCount = 0 ;

            for ( mesh::Node * tStart : aPins )
            {
                graph::Vertex * tPin = graph::multibfs( tGraph, tStart->owner() );

                BELFEM_ERROR( tPin != nullptr,
                    "autopin: phi component %u vanished from the graph",
                    ( unsigned int ) tStart->owner() );

                // winner-level semantics from multibfs: 0 = the component
                // touches no sideset at all ( enclosed — any node is safe,
                // the tie-break returned its smallest id ); 1 = EVERY node
                // sits on a sideset, nothing is safely pinnable; >= 2 = an
                // interior node, the regular case
                BELFEM_ERROR( tPin->level() != 1,
                    "autopin: every node of phi component %u lies on a sideset - cannot select a safe gauge pin",
                    ( unsigned int ) tStart->owner() );

                aPins( tCount++ ) = reinterpret_cast< mesh::Node * >( tPin );
            }

            // tidy up
            sort( tGraph , opVertexID );
            tCount = 0 ;
            for ( auto tNode : tGraph )
            {
                tNode->reset_vertex_container() ;
                tNode->set_level( tLevels( tCount  ) );
                tNode->set_owner( tOwners( tCount++ ) );
            }

            // the seed pass flagged EVERY sideset node mesh-wide, including
            // nodes outside the graph — clear both scratch planes globally
            // ( house rule: flags are cleared before AND after )
            mMesh->unflag_all_nodes( 0 );
            mMesh->unflag_all_nodes( 1 );

            mMesh->update_node_indices() ;

        }

    }
}
