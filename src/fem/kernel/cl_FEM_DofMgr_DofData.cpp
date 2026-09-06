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

#include "cl_FEM_DofMgr_DofData.hpp"
#include "cl_FEM_DofMgr_SolverData.hpp"

#include "cl_FEM_DofManager.hpp"
#include "cl_Mesh.hpp"
#include "cl_IWG.hpp"
#include "cl_FEM_Dof.hpp"
#include "cl_Vertex.hpp"

#include "fn_entity_type.hpp"
#include "assert.hpp"
#include "commtools.hpp"
#include "cl_FEM_Kernel.hpp"
#include "fn_max.hpp"
#include "fn_unique.hpp"

#include "fn_to_master_orientation.hpp"

#include "fn_Graph_METIS.hpp"
#include "fn_Graph_ParMETIS.hpp"
#include "fn_Graph_SCOTCH.hpp"
#include "fn_Graph_PTSCOTCH.hpp"


#include "fn_Graph_symrcm.hpp"
#include "op_Graph_Vertex_Index.hpp"
#include "fn_sort.hpp"

namespace belfem
{
    namespace fem
    {

        namespace dofmgr
        {
//------------------------------------------------------------------------------

            DofData::DofData(  DofManager * aParent,
                               Parameters * aParams ) :
                mParent( aParent ),
                mKernel( aParent->parent() ),
                mMesh( aParent->parent()->mesh() ),
                mParams( aParams ),
                mCommRank( comm_rank() ),
                mCommSize( comm_size() )
            {
                
            }

//------------------------------------------------------------------------------

            DofData::~DofData()
            {
                // restore factory settings
                this->reset() ;
            }

//------------------------------------------------------------------------------

            void
            DofData::reset()
            {
                // delete the maps
                mDofMap.clear() ;
                mDofTypeToField.clear() ;

                // delete the tables
                mDofIndexTables.clear() ;

                // delete the dofs
                for ( Dof * tDof : mDOFs )
                {
                    delete tDof ;
                }
                mDOFs.clear() ;

                // delete the hanging
                for ( Dof * tDof : mHangingDOFs )
                {
                    delete tDof ;
                }
                mHangingDOFs.clear() ;


                mAbstractDOFs.clear() ;

                // reset the offsets
                mEdgeDofOffset   = gNoID ;
                mFaceDofOffset   = gNoID ;
                mCellDofOffset   = gNoID ;
                mLambdaDofOffset = gNoID ;
            }

//------------------------------------------------------------------------------

            void
            DofData::create_dofs( IWG * aIWG )
            {
                BELFEM_ERROR( comm_size() < BELFEM_MAX_NUMPROCS,
                             "This program is not supposed to run on more that %u procs. \nRedefine BELFEM_MAX_NUMPROCS in cl_IWG.hpp and compile again",
                             ( unsigned int ) BELFEM_MAX_NUMPROCS );


                Vector< id_t >    tNodeDofIDs ;
                Vector< id_t >    tNodeDofEntityIDs;
                Vector< index_t > tNodeDofTypes;

                Vector< id_t >    tEdgeDofIDs ;
                Vector< id_t >    tEdgeDofEntityIDs;
                Vector< index_t > tEdgeDofTypes;

                Vector< id_t >    tFaceDofIDs ;
                Vector< id_t >    tFaceDofEntityIDs;
                Vector< index_t > tFaceDofTypes;

                Vector< id_t >    tCellDofIDs ;
                Vector< id_t >    tCellDofEntityIDs;
                Vector< index_t > tCellDofTypes;

                Vector< id_t >    tLambdaDofIDs ;
                Vector< id_t >    tLambdaDofEntityIDs;
                Vector< index_t > tLambdaDofTypes;

                // compute offsets
                this->compute_dof_offsets( aIWG );

                if( mCommRank == 0 )
                {
                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // Identify Node-Based DOFs
                    // - - - - - - - - - - - - - - - - - - - - - - - - - -


                    Cell< Bitset< BELFEM_MAX_NUMPROCS > > tNodeProcs;
                    index_t tNumNodeDofs =
                            this->count_node_dofs( aIWG,
                                                   tNodeDofEntityIDs,
                                                   tNodeDofTypes,
                                                   tNodeProcs );

                    // make list with NODE IDs
                    tNodeDofIDs.set_size( tNumNodeDofs );
                    for ( index_t k = 0; k < tNumNodeDofs; ++k )
                    {
                        tNodeDofIDs( k ) = this->node_dof_id(
                                tNodeDofEntityIDs( k ),
                                tNodeDofTypes( k ) );
                    }

                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // Identify Edge-Based DOFs
                    // - - - - - - - - - - - - - - - - - - - - - - - - - -);

                    Cell< Bitset< BELFEM_MAX_NUMPROCS > > tEdgeProcs;
                    index_t tNumEdgeDofs =
                            this->count_edge_dofs( aIWG,
                                                   tEdgeDofEntityIDs,
                                                   tEdgeDofTypes,
                                                   tEdgeProcs );


                    tEdgeDofIDs.set_size( tNumEdgeDofs );

                    // make list with EDGE IDs
                    for ( index_t k = 0; k < tNumEdgeDofs; ++k )
                    {
                        tEdgeDofIDs( k ) = this->edge_dof_id(
                                tEdgeDofEntityIDs( k ),
                                tEdgeDofTypes( k ) );

                    }

                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // Identify Face-Based DOFs
                    // - - - - - - - - - - - - - - - - - - - - - - - - - -


                    Cell< Bitset< BELFEM_MAX_NUMPROCS > > tFaceProcs;
                    index_t tNumFaceDofs =
                            this->count_face_dofs( aIWG,
                                                   tFaceDofEntityIDs,
                                                   tFaceDofTypes,
                                                   tFaceProcs );

                    tFaceDofIDs.set_size( tNumFaceDofs );

                    // make list with FACE IDs
                    for ( index_t k = 0; k < tNumFaceDofs; ++k )
                    {
                        tFaceDofIDs( k ) = this->face_dof_id(
                                tFaceDofEntityIDs( k ),
                                tFaceDofTypes( k ) );
                    }

                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // Identify Cell-Based DOFs
                    // - - - - - - - - - - - - - - - - - - - - - - - - - -

                    Cell< Bitset< BELFEM_MAX_NUMPROCS > > tCellProcs;
                    index_t tNumCellDofs =
                            this->count_cell_dofs( aIWG,
                                                   tCellDofEntityIDs,
                                                   tCellDofTypes,
                                                   tCellProcs );

                    tCellDofIDs.set_size( tNumCellDofs );

                    // make list with FACE IDs
                    for ( index_t k = 0; k < tNumCellDofs; ++k )
                    {
                        tCellDofIDs( k ) = this->cell_dof_id(
                                tCellDofEntityIDs( k ),
                                tCellDofTypes( k ) );
                    }

                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // Identify Lambda DOFs
                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Cell< Bitset< BELFEM_MAX_NUMPROCS > > tLambdaProcs;
                    index_t tNumLambdaDofs =
                            this->count_lambda_dofs( aIWG,
                                                     tLambdaDofEntityIDs,
                                                     tLambdaDofTypes,
                                                     tLambdaProcs );

                    // make list with Lambda IDs
                    tLambdaDofIDs.set_size( tNumLambdaDofs );
                    for ( index_t k = 0; k < tNumLambdaDofs; ++k )
                    {
                        tLambdaDofIDs( k ) = this->lambda_dof_id(
                                tLambdaDofEntityIDs( k ),
                                tLambdaDofTypes( k ) );
                    }

                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // Allocate Containers for Data
                    // - - - - - - - - - - - - - - - - - - - - - - - - - -

                    proc_t tNumberOfProcs = mKernel->number_of_procs() ;

                    Cell< Vector< id_t > > tAllDofIDs( tNumberOfProcs, Vector< id_t >() );
                    Cell< Vector< id_t > > tAllEntityIDs( tNumberOfProcs, Vector< id_t >() );
                    Cell< Vector< index_t > > tAllDofTypes( tNumberOfProcs, Vector< index_t >() );

                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // Send Node-Based DOFs
                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    for( proc_t p=0; p<tNumberOfProcs; ++p )
                    {

                        // only do something if this is not the master proc
                        if ( p != mCommRank )
                        {

                            // collect dof entities
                            this->count_dofs_for_proc( p,
                                                       tNodeDofIDs,
                                                       tNodeDofEntityIDs,
                                                       tNodeDofTypes,
                                                       tNodeProcs,
                                                       tAllDofIDs( p ),
                                                       tAllEntityIDs( p ),
                                                       tAllDofTypes( p ) );
                        }
                    }

                    // wait for other procs
                    comm_barrier();

                    // send data to other procs
                    distribute( tAllDofIDs );
                    distribute( tAllEntityIDs );
                    distribute( tAllDofTypes );

                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // Send Edge-Based DOFs
                    // - - - - - - - - - - - - - - - - - - - - - - - - - -

                    for( proc_t p=0; p<tNumberOfProcs; ++p )
                    {

                        // only do something if this is not the master proc
                        if ( p != mCommRank )
                        {

                            // collect dof entities
                            this->count_dofs_for_proc( p,
                                                       tEdgeDofIDs,
                                                       tEdgeDofEntityIDs,
                                                       tEdgeDofTypes,
                                                       tEdgeProcs,
                                                       tAllDofIDs( p ),
                                                       tAllEntityIDs( p ),
                                                       tAllDofTypes( p ) );
                        }
                    }

                    // wait for other procs
                    comm_barrier();

                    // send data to other procs
                    distribute( tAllDofIDs );
                    distribute( tAllEntityIDs );
                    distribute( tAllDofTypes );

                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // Send Face-Based DOFs
                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    for( proc_t p=0; p<tNumberOfProcs; ++p )
                    {

                        // only do something if this is not the master proc
                        if ( p != mCommRank )
                        {

                            // collect dof entities
                            this->count_dofs_for_proc( p,
                                                       tFaceDofIDs,
                                                       tFaceDofEntityIDs,
                                                       tFaceDofTypes,
                                                       tFaceProcs,
                                                       tAllDofIDs( p ),
                                                       tAllEntityIDs( p ),
                                                       tAllDofTypes( p ) );
                        }
                    }

                    // wait for other procs
                    comm_barrier();

                    // send data to other procs
                    distribute( tAllDofIDs );
                    distribute( tAllEntityIDs );
                    distribute( tAllDofTypes );

                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // Send Cell-Based DOFs
                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    for( proc_t p=0; p<tNumberOfProcs; ++p )
                    {
                        // only do something if this is not the master proc
                        if ( p != mCommRank )
                        {
                            // collect dof entities
                            this->count_dofs_for_proc( p,
                                                       tCellDofIDs,
                                                       tCellDofEntityIDs,
                                                       tCellDofTypes,
                                                       tCellProcs,
                                                       tAllDofIDs( p ),
                                                       tAllEntityIDs( p ),
                                                       tAllDofTypes( p ) );
                        }
                    }

                    // wait for other procs
                    comm_barrier();

                    // send data to other procs
                    distribute( tAllDofIDs );
                    distribute( tAllEntityIDs );
                    distribute( tAllDofTypes );

                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // Send LambdaDOFs
                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    for( proc_t p=0; p<tNumberOfProcs; ++p )
                    {
                        // only do something if this is not the master proc
                        if ( p != mCommRank )
                        {
                            // collect dof entities
                            this->count_dofs_for_proc( p,
                                                       tLambdaDofIDs,
                                                       tLambdaDofEntityIDs,
                                                       tLambdaDofTypes,
                                                       tLambdaProcs,
                                                       tAllDofIDs( p ),
                                                       tAllEntityIDs( p ),
                                                       tAllDofTypes( p ) );
                        }
                    }

                    // wait for other procs
                    comm_barrier();

                    // send data to other procs
                    distribute( tAllDofIDs );
                    distribute( tAllEntityIDs );
                    distribute( tAllDofTypes );


                }
                else
                {
                    // wait for other procs
                    comm_barrier();

                    // receive node info from master
                    receive( tNodeDofIDs );
                    receive( tNodeDofEntityIDs );
                    receive( tNodeDofTypes );

                    // wait for other procs
                    comm_barrier();

                    // receive edge info from master
                    receive( tEdgeDofIDs );
                    receive( tEdgeDofEntityIDs );
                    receive( tEdgeDofTypes );

                    // wait for other procs
                    comm_barrier();

                    // receive face info from master
                    receive( tFaceDofIDs );
                    receive( tFaceDofEntityIDs );
                    receive( tFaceDofTypes );

                    // wait for other procs
                    comm_barrier();

                    // receive cell info from master
                    receive( tCellDofIDs );
                    receive( tCellDofEntityIDs );
                    receive( tCellDofTypes );

                    // wait for other procs
                    comm_barrier();

                    receive( tLambdaDofIDs );
                    receive( tLambdaDofEntityIDs );
                    receive( tLambdaDofTypes );
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Crete the dofs based on the given IDs
                // - - - - - - - - - - - - - - - - - - - - - - - - - -

                // total number of dofs ( master: all, other: on this proc )
                index_t tNumNodeDofs = tNodeDofIDs.length() ;
                index_t tNumEdgeDofs = tEdgeDofIDs.length() ;
                index_t tNumFaceDofs = tFaceDofIDs.length() ;
                index_t tNumCellDofs = tCellDofIDs.length() ;

                index_t tNumLambdaDofs = tLambdaDofIDs.length() ;

                // total number of dofs
                index_t tNumDofs =
                          tNumNodeDofs
                        + tNumEdgeDofs
                        + tNumFaceDofs
                        + tNumCellDofs
                        + tNumLambdaDofs ;

                BELFEM_ASSERT( tNumDofs > 0 || ! mKernel->is_master(), "No dofs exist. What now?" );

                index_t tCount = 0;


                mDOFs.set_size( tNumDofs, nullptr );

                // create the node DOFs
                for( index_t k=0; k<tNumNodeDofs; ++k )
                {
                    // create a new dof
                    Dof * tDof = new Dof( tNodeDofIDs( k ),
                                          tNodeDofTypes( k ),
                                          mMesh->node( tNodeDofEntityIDs( k ) ) );

                    // set the dof index. This will be overwritten during Jacobi initialization
                    tDof->set_index( tCount );

                    // add the dof to the container
                    mDOFs( tCount++ ) = tDof ;
                }

                // the ID of the old entity
                id_t tID = 0 ;

                // counter for multiplicity
                index_t i = 0 ;

                // multiplicity counter
                index_t tMult = aIWG->edge_multiplicity() ;

                // create the edge dofs
                for( index_t k=0; k<tNumEdgeDofs; ++k )
                {
                    // get edge
                    mesh::Edge * tEdge = mMesh->edge( tEdgeDofEntityIDs( k ) );

                    // check if this is a new edge
                    if( tID != tEdge->id() )
                    {
                        // remember new id
                        tID = tEdge->id() ;

                        // reset counter
                        i = 0 ;
                    }

                    // compute the field index
                    index_t tIndexOnField =  tEdge->index() * tMult + i;

                    // create a new dof
                    Dof * tDof = new Dof( tEdgeDofIDs( k ),
                                          tEdgeDofTypes( k ),
                                          tEdge ,
                                          i++,
                                          tIndexOnField );

                    // set the dof index. This will be overwritten during Jacobi initialization
                    tDof->set_index( tCount );

                    // add the dof to the container
                    mDOFs( tCount++ ) = tDof ;
                }

                // reset counters
                tID = 0 ;
                i = 0 ;
                tMult = aIWG->face_multiplicity() ;

                // create the facedofs
                for( index_t k=0; k<tNumFaceDofs; ++k )
                {
                    mesh::Face * tFace = mMesh->face( tFaceDofEntityIDs( k ) );

                    // check if this is a new face
                    if( tID != tFace->id() )
                    {
                        // remember new id
                        tID = tFace->id() ;

                        // reset counter
                        i = 0 ;
                    }

                    // compute new index
                    index_t tIndexOnField =  tFace->index()    * tMult + i;

                    // create a new dof
                    Dof * tDof = new Dof( tFaceDofIDs( k ),
                                          tFaceDofTypes( k ),
                                          tFace,
                                          i++,
                                          tIndexOnField );

                    // set the dof index. This will be overwritten during Jacobi initialization
                    tDof->set_index( tCount );

                    // add the dof to the container
                    mDOFs( tCount++ ) = tDof ;
                }

                // create cell dofs
                // reset counters
                tID = 0 ;
                i = 0 ;
                tMult = aIWG->cell_multiplicity() ;

                // create the cell dofs
                for( index_t k=0; k<tNumCellDofs; ++k )
                {
                    mesh::Element * tCell = mMesh->element( tCellDofEntityIDs( k ) );

                    // check if this is a new face
                    if( tID != tCell->id() )
                    {
                        // remember new id
                        tID = tCell->id() ;

                        // reset counter
                        i = 0 ;
                    }

                    // compute new index
                    index_t tIndexOnField =  tCell->index() * tMult + i;

                    // create a new dof
                    Dof * tDof = new Dof( tCellDofIDs( k ),
                                          tCellDofTypes( k ),
                                          tCell,
                                          i++,
                                          tIndexOnField );

                    // set the dof index. This will be overwritten during Jacobi initialization
                    tDof->set_index( tCount );

                    // add the dof to the container
                    mDOFs( tCount++ ) = tDof ;
                }

                tMult = aIWG->lambda_multiplicity() ;
                tID = 0 ;
                i = 0 ;

                // create the lambda DOFs
                for( index_t k=0; k<tNumLambdaDofs; ++k )
                {
                    mesh::Facet * tFacet = mMesh->facet( tLambdaDofEntityIDs( k ) ) ;

                    // check if this is a new facet
                    if( tID != tFacet->id() )
                    {
                        // remember new id
                        tID = tFacet->id() ;

                        // reset counter
                        i = 0 ;
                    }

                    index_t tIndexOnField =  tFacet->index() * tMult + i;

                    // create a new dof
                    Dof * tDof = new Dof( tLambdaDofIDs( k ),
                                          tLambdaDofTypes( k ),
                                          tFacet,
                                          i++,
                                          tIndexOnField );

                    // set the dof index. This will be overwritten during Jacobi initialization
                    tDof->set_index( tCount );

                    // add the dof to the container
                    mDOFs( tCount++ ) = tDof ;
                }

                // create the dof map
                this->create_dof_map();
            }

//------------------------------------------------------------------------------

            void
            DofData::create_field_map( IWG  * aIwg )
            {
                // get the fields that belong to the IWG
                const Cell< string > & tLabels = aIwg->dof_fields();

                // number of fields
                uint tNumFields = tLabels.size() ;

                // field counter
                index_t tCount = 0;

                // reset the map
                mDofTypeToField.clear() ;

                for( uint f=0; f<tNumFields; ++f )
                {
                    // get the label
                    const string & tLabel = tLabels( f );

                    // get the entity type of the field
                    EntityType tType = entity_type( tLabel );

                    // get field index
                    index_t tFieldIndex = mMesh->field( tLabel )->index() ;


                    switch( tType )
                    {
                        case( EntityType::NODE ) :
                        {
                            // add index to map
                            mDofTypeToField[ tCount++ ] = tFieldIndex ;

                            break ;
                        }
                        case( EntityType::EDGE ) :
                        {
                            for( uint k=0; k<aIwg->edge_multiplicity(); ++k )
                            {
                                mDofTypeToField[ tCount++ ] = tFieldIndex ;
                            }
                            break ;
                        }
                        case( EntityType::FACE ) :
                        {
                            for( uint k=0; k<aIwg->face_multiplicity(); ++k )
                            {
                                mDofTypeToField[ tCount++ ] = tFieldIndex ;
                            }
                            break ;
                        }
                        case( EntityType::CELL ) :
                        {
                            for( uint k=0; k<aIwg->cell_multiplicity(); ++k )
                            {
                                mDofTypeToField[ tCount++ ] = tFieldIndex ;
                            }
                            break ;
                        }
                        case( EntityType::FACET ) : // for lambda dofs
                        {
                            for( uint k=0; k<aIwg->lambda_multiplicity(); ++k )
                            {
                                mDofTypeToField[ tCount++ ] = tFieldIndex ;
                            }
                            break ;
                        }
                        case( EntityType::ELEMENT ) :
                        {
                            // element is not a dof. Use Cell instead
                            break ;
                        }
                        default:
                        {
                            BELFEM_ERROR( false, "Invalid entity type");
                        }
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            DofData::create_dof_map()
            {
                // reset the map
                mDofMap.clear() ;

                // reindex entities
                for( Dof * tDOF : mDOFs )
                {
                    mDofMap[ tDOF->id() ] = tDOF ;
                }
            }

//------------------------------------------------------------------------------

            void
            DofData::connect_dofs_to_mesh()
            {
                // count dofs per basis
                for( Dof * tDof : mDOFs )
                {
                    tDof->mesh_basis()->increment_dof_counter() ;
                }

                // allocate memory
                allocate_dof_containers( mMesh->nodes() );
                allocate_dof_containers( mMesh->edges() );
                allocate_dof_containers( mMesh->faces() );
                allocate_dof_containers( mMesh->facets() );
                allocate_dof_containers( mMesh->elements() );

                for( Dof * tDof : mDOFs )
                {
                    tDof->mesh_basis()->insert_dof( tDof );
                }
            }

//------------------------------------------------------------------------------

            void
            DofData::disconnect_dofs_from_mesh()
            {
                for( Dof * tDof : mDOFs )
                {
                    tDof->mesh_basis()->reset_dof_container();
                }
            }

//------------------------------------------------------------------------------

            void
            DofData::compute_dof_offsets( IWG  * aIWG )
            {
                // wait
                comm_barrier() ;

                BELFEM_ASSERT( mParent->iwg()->is_initialized(),
                     "Equation object has not been initialized");

                if( comm_rank() == 0 )
                {
                    Vector< index_t > tMaxIDs ;
                    this->compute_max_ids( tMaxIDs );


                    mNumDofTypes = aIWG->dof_entity_types().length() ;

                    // compute the offsets
                    mEdgeDofOffset = ( tMaxIDs( 0 ) + 1 ) * mNumDofTypes ;

                    mFaceDofOffset = mEdgeDofOffset+
                                     ( tMaxIDs( 1 ) + 1 ) * mNumDofTypes ;

                    mCellDofOffset = mFaceDofOffset +
                                     ( tMaxIDs( 2 ) + 1 ) * mNumDofTypes ;


                    mLambdaDofOffset = mCellDofOffset +
                                       ( tMaxIDs( 3 ) + 1 ) * mNumDofTypes ;


                    Vector< index_t > tOffsets( 5 ) ;

                    tOffsets( 0 ) = mNumDofTypes ;
                    tOffsets( 1 ) = mEdgeDofOffset ;
                    tOffsets( 2 ) = mFaceDofOffset ;
                    tOffsets( 3 ) = mCellDofOffset ;
                    tOffsets( 4 ) = mLambdaDofOffset  ;

                    // send info to other procs
                    broadcast( tOffsets );
                }
                else
                {
                    // receive data from master
                    Vector< index_t > tOffsets( 5 ) ;
                    broadcast( tOffsets );

                    mNumDofTypes  = tOffsets( 0 ) ;
                    mEdgeDofOffset = tOffsets( 1 ) ;
                    mFaceDofOffset = tOffsets( 2 ) ;
                    mCellDofOffset = tOffsets( 3 ) ;
                    mLambdaDofOffset  = tOffsets( 4 ) ;
                }

                // wait
                comm_barrier() ;
            }

//------------------------------------------------------------------

            void
            DofData::compute_max_ids( Vector< id_t > & aMaxEntityIDs )
            {
                aMaxEntityIDs.set_size( 4, 0 );

                // compute max node id
                id_t & tMaxNodeID = aMaxEntityIDs( 0 );

                // loop over all nodes in mesh
                for( mesh::Node * tNode : mMesh->nodes() )
                {
                    if( tNode->id() > tMaxNodeID )
                    {
                        tMaxNodeID = tNode->id() ;
                    }
                }

                // check if edges have been generated
                if ( mMesh->edges_exist() )
                {
                    // compute max edge ID
                    id_t & tMaxEdgeID = aMaxEntityIDs( 1 );

                    for( mesh::Edge * tEdge : mMesh->edges() )
                    {
                        if( tEdge->id() > tMaxEdgeID )
                        {
                            tMaxEdgeID = tEdge->id() ;
                        }
                    }
                }

                // check if faces have been generated
                if( mMesh->faces_exist() )
                {
                    // compute max face ID
                    id_t & tMaxFaceID = aMaxEntityIDs( 2 );

                    for( mesh::Face * tFace : mMesh->faces() )
                    {
                        if( tFace->id() > tMaxFaceID )
                        {
                            tMaxFaceID = tFace->id() ;
                        }
                    }
                }


                // compute max element id
                id_t & tMaxElementID = aMaxEntityIDs( 3 );
                for( mesh::Element * tElement : mMesh->elements() )
                {
                    if( tElement->id() > tMaxElementID )
                    {
                        tMaxElementID = tElement->id() ;
                    }
                }
            }

//-----------------------------------------------------------------------------

            index_t
            DofData::count_node_dofs( IWG  * aIWG,
                                      Vector< id_t >    & aEntityIDs,
                                      Vector< index_t > & aDofTypes,
                                      Cell< Bitset<BELFEM_MAX_NUMPROCS> > & aProcFlags )
            {
                /* depending on which function this is, "entity" refers to
                   a node, an edge, face or element

                   this routine performs xxx steps:

                   step 1 : flag all mesh entities that sit on the selected blocks

                   step 2:  count the selected entities and create a lookup table

                   step 3:  create the dof table: in order to find out which dofs exist, we create an
                            array of bitsets and flip the bitsets for each mesh entity

                   step 4:  count how many dofs have to be created

                   step 5:  polpulate the containers for the entity ids and dof types

                   step 6:  determine which dofs are visible on which proc
                   */
                // unflag all entities on the mesh
                mMesh->unflag_all_nodes() ;

                // we remember the ID of reach entity
                Vector< id_t > tEntityIDs ;

                // first, we need to figure out if entities dofs exist at all
                // if so, relevant entities are flagged
                bool tHaveEntityDofs = false ;

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 1: flag all mesh entities the dofs are based on
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                for ( id_t tID : aIWG->selected_blocks() )
                {
                    if ( aIWG->number_of_dofs_per_node( tID ) > 0 )
                    {
                        tHaveEntityDofs = true;
                        mMesh->block( tID )->flag_nodes();
                    }
                }

                // also flag abstract nodes if they exist (special for Maxwell)
                for( mesh::Node * tNode : aIWG->abstract_nodes() )
                {
                    tNode->flag() ;
                }

                // also flag orphaned nodes if they exist  (special for Maxwell)
                for( mesh::Node * tNode : aIWG->orphaned_nodes() )
                {
                    tNode->flag() ;
                }

                // exit the routine if no entity dofs exist
                if( ! tHaveEntityDofs )
                {
                    return 0 ;
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 2: count the selected entities and create tables
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // initialize counter
                index_t tEntityCount = 0 ;

                // the entity maps connects the ids to the counter
                Map< id_t, index_t > tEntityMap ;

                // now we count the number of flagged entities
                for( mesh::Node * tNode : mMesh->nodes() )
                {
                    if( tNode->is_flagged() )
                    {
                        // write counter into map and increment it
                        tEntityMap[ tNode->id() ] = tEntityCount++ ;
                    }
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 3: flip entity-wise bitsets for used dofs
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // the bitset array
                Cell< Bitset<BELFEM_MAX_DOFTYPES> > tDofFlags( tEntityCount, Bitset<BELFEM_MAX_DOFTYPES>() );

                // we loop over all blocks and get the block-wise dofs
                for( id_t tBlockID: aIWG->selected_blocks() )
                {
                    // grab block pointer on mesh
                    mesh::Block * tBlock = mMesh->block( tBlockID );

                    // ask IWG about selected dofs
                    const Vector< index_t > & tSelectedDofs = aIWG->dofs_per_node( tBlockID );

                    // check if any dofs are selected
                    if( tSelectedDofs.length() > 0 )
                    {
                        // number of nodes per element on this block ( all nodes carry dofs )
                        uint tNumEntities = mesh::number_of_nodes( tBlock->element_type() );

                        // now we loop over all elements on the block and the selected number of nodes
                        for( mesh::Element * tElement : tBlock->elements() )
                        {
                            // loop over all mesh entities on this element
                            for( uint k=0; k<tNumEntities; ++k )
                            {
                                // get the temporary index of this entity
                                index_t tIndex = tEntityMap( tElement->node( k )->id() );

                                // grab the corresponding bitset
                                Bitset< BELFEM_MAX_DOFTYPES > & tBitset = tDofFlags( tIndex );

                                // flip the bits that represent the selected dofs
                                for( index_t d : tSelectedDofs )
                                {
                                    tBitset.set( d );
                                }
                            }
                        }
                    }
                }

                for( id_t tSideSetID: aIWG->selected_sidesets() )
                {
                    // grab sideset pointer on mesh
                    mesh::SideSet * tSideSet = mMesh->sideset( tSideSetID );

                    // ask IWG about selected dofs
                    const Vector< index_t > & tSelectedDofs = aIWG->dofs_per_node_on_sideset( tSideSetID, true );

                    // check if any dofs are selected
                    if( tSelectedDofs.length() > 0 )
                    {
                        // get the number of edges per element on this sideset
                        uint tNumEntities = mesh::number_of_nodes( tSideSet->element_type() );

                        // now we loop over all elements on the block and the selected number oedges
                        for( mesh::Facet * tFacet : tSideSet->facets() )
                        {
                            // loop over all mesh entities on this element
                            for( uint k=0; k<tNumEntities; ++k )
                            {
                                // get the temporary index of this entity
                                index_t tIndex = tEntityMap( tFacet->node( k )->id() );

                                // grab the corresponding bitset
                                Bitset< BELFEM_MAX_DOFTYPES > & tBitset = tDofFlags( tIndex );

                                // flip the bits that represent the selected dofs
                                for( index_t d : tSelectedDofs )
                                {
                                    tBitset.set( d );
                                }
                            }
                        }
                    }
                }

                // dealing with abstract and orphaned nodes
                if( aIWG->abstract_dof_type() != gNoIndex )
                {
                    index_t d = aIWG->abstract_dof_type() ;

                    for( mesh::Node * tNode : aIWG->abstract_nodes() )
                    {

                        // get the temporary index of this entity
                        index_t tIndex = tEntityMap( tNode->id() );

                        // grab the corresponding bitset
                        Bitset< BELFEM_MAX_DOFTYPES > & tBitset = tDofFlags( tIndex );

                        tBitset.set( d );
                    }

                    for( mesh::Node * tNode : aIWG->orphaned_nodes() )
                    {
                        // get the temporary index of this entity
                        index_t tIndex = tEntityMap( tNode->id() );

                        // grab the corresponding bitset
                        Bitset< BELFEM_MAX_DOFTYPES > & tBitset = tDofFlags( tIndex );

                        tBitset.set( d );
                    }
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 4: count how many dofs have to be created
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                index_t aDofCount = 0 ;

                // we loop over all dof bitsets
                for( index_t k=0; k<tEntityCount; ++k )
                {
                    // grab the current bitset
                    Bitset< BELFEM_MAX_DOFTYPES > & tBitset = tDofFlags( k );

                    // count the dofs
                    aDofCount += tBitset.count() ;
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 5: populate the containers for the entity ids,
                //         the dof types and the used procs
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // allocate memory
                aEntityIDs.set_size( aDofCount );
                aDofTypes.set_size( aDofCount );


                // this is a temporary map which is needed to assign procs and dofs
                Map< luint, index_t > tDofMap ;


                // reset counters
                aDofCount = 0 ;

                index_t tIndex ;
                index_t tBitCountA ;
                index_t tBitCountB ;

                // loop over all flagged entities
                for( mesh::Node * tNode : mMesh->nodes() )
                {
                    // check if node is flagged
                    if( tNode->is_flagged() )
                    {
                        // get index in array
                        tIndex = tEntityMap( tNode->id() );

                        // get the corresponding bitset
                        Bitset< BELFEM_MAX_DOFTYPES > & tBitset = tDofFlags( tIndex );

                        tBitCountA = 0 ;
                        tBitCountB = tBitset.count() ;

                        // loop over all bits
                        for( uint d=0; d<BELFEM_MAX_DOFTYPES; ++d )
                        {
                            if( tBitset.test( d ) )
                            {
                                // remember ID
                                aEntityIDs( aDofCount ) = tNode->id() ;

                                // remember dof type
                                aDofTypes( aDofCount )  = d ;

                                // create a unique and map it
                                tDofMap[ tNode->id() * BELFEM_MAX_DOFTYPES + d ] = aDofCount++ ;

                                // cancel loop if we are done
                                if( ++tBitCountA == tBitCountB )
                                {
                                    break ;
                                }
                            }
                        }
                    }
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 6: determine which dofs are visible on which proc
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // allocate memory
                aProcFlags.set_size( aDofCount, Bitset< BELFEM_MAX_NUMPROCS >() );

                // loop over all blocks
                // define which dof is used by which proc
                for( id_t tBlockID : aIWG->selected_blocks() )
                {
                    // grab block on mesh
                    mesh::Block * tBlock = mMesh->block( tBlockID );

                    // grab the DOFs for this entity
                    const Vector< index_t > & tSelectedDofs = aIWG->dofs_per_node( tBlockID );

                    if ( tSelectedDofs.length() > 0 )
                    {
                        // number of nodes per element on this block ( all nodes carry dofs )
                        uint tNumEntities = mesh::number_of_nodes( tBlock->element_type() );


                        // loop over all elements on this block
                        for( mesh::Element * tElement : tBlock->elements() )
                        {
                            // loop over all entities on this element
                            for( uint k=0; k<tNumEntities; ++k )
                            {
                                // loop over all selected dofs
                                for( index_t d : tSelectedDofs )
                                {
                                    // compute the unique dof ID
                                    luint tID = tElement->node( k )->id() * BELFEM_MAX_DOFTYPES + d ;

                                    // grab the corresponding bitset
                                    Bitset< BELFEM_MAX_NUMPROCS > & tBitset = aProcFlags( tDofMap( tID ) );

                                    // set the proc bitset
                                    tBitset.set( tElement->owner() );
                                }
                            }
                        }
                    }
                }

                for( id_t tSideSetID : aIWG->selected_sidesets() )
                {
                    // grab sideset on mesh
                    mesh::SideSet * tSideSet = mMesh->sideset( tSideSetID );

                    // grab the DOFs for this entity
                    const Vector< index_t > & tSelectedDofs = aIWG->dofs_per_node_on_sideset( tSideSetID, true );


                    if ( tSelectedDofs.length() > 0 )
                    {
                        // number of nodes per element on this block ( all nodes carry dofs )
                        uint tNumEntities = mesh::number_of_nodes( tSideSet->element_type());


                        // loop over all elements on this block
                        for ( mesh::Facet * tFacet: tSideSet->facets())
                        {
                            // loop over all entities on this element
                            for ( uint k = 0; k < tNumEntities; ++k )
                            {
                                // loop over all selected dofs
                                for ( index_t d: tSelectedDofs )
                                {
                                    // compute the unique dof ID
                                    luint tID = tFacet->node( k )->id() * BELFEM_MAX_DOFTYPES + d;

                                    // grab the corresponding bitset
                                    Bitset< BELFEM_MAX_NUMPROCS > & tBitset = aProcFlags( tDofMap( tID ));

                                    // set the proc bitset
                                    tBitset.set( tFacet->owner());
                                }
                            }
                        }
                    }
                }

                // return the number of DOFs
                return aDofCount ;
            }

//------------------------------------------------------------------------

            index_t
            DofData::count_edge_dofs( IWG  * aIWG,
                                      Vector< id_t >    & aEntityIDs,
                                      Vector< index_t > & aDofTypes,
                                      Cell< Bitset<BELFEM_MAX_NUMPROCS> > & aProcFlags )
            {

                /* depending on which function this is, "entity" refers to
                   a node, an edge, face or element

                   this routine performs xxx steps:

                   step 1 : flag all mesh entities that sit on the selected blocks

                   step 2:  count the selected entities and create a lookup table

                   step 3:  create the dof table: in order to find out which dofs exist, we create an
                            array of bitsets and flip the bitsets for each mesh entity

                   step 4:  count how many dofs have to be created

                   step 5:  polpulate the containers for the entity ids and dof types

                   step 6:  determine which dofs are visible on which proc
                   */
                // unflag all entities on the mesh
                mMesh->unflag_all_edges() ;

                // we remember the ID of each entity
                Vector< id_t > tEntityIDs ;

                // first, we need to figure out if entities dofs exist at all
                // if so, relevant entities are flagged
                bool tHaveEntityDofs = false ;

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 1: flag all mesh entities the dofs are based on
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                for( id_t tID : aIWG->selected_blocks() )
                {
                    if ( aIWG->number_of_dofs_per_edge( tID ) > 0 )
                    {
                        tHaveEntityDofs = true ;
                        mMesh->block( tID )->flag_edges() ;
                    }
                }

                // exit the routine if no entity dofs exist
                if( ! tHaveEntityDofs )
                {
                    return 0 ;
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 2: count the selected entities and create tables
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // initialize counter
                index_t tEntityCount = 0 ;

                // the entity maps connects the ids to the counter
                Map< id_t, index_t > tEntityMap ;

                // now we count the number of flagged entities
                for( mesh::Edge * tEdge : mMesh->edges() )
                {
                    if( tEdge->is_flagged() )
                    {
                        // write counter into map and increment it
                        tEntityMap[ tEdge->id() ] = tEntityCount++ ;
                    }
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 3: flip entity-wise bitsets for used dofs
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // the bitset array
                Cell< Bitset<BELFEM_MAX_DOFTYPES> > tDofFlags( tEntityCount, Bitset<BELFEM_MAX_DOFTYPES>() );

                index_t tMultiplicity = aIWG->edge_multiplicity();

                // we loop over all blocks and get the block-wise dofs
                for( id_t tBlockID: aIWG->selected_blocks() )
                {
                    // grab block pointer on mesh
                    mesh::Block * tBlock = mMesh->block( tBlockID );

                    // ask IWG about selected dofs
                    const Vector< index_t > & tSelectedDofs = aIWG->dofs_per_edge( tBlockID );

                    // check if any dofs are selected
                    if( tSelectedDofs.length() > 0 )
                    {
                        // get the number of edges per element on this block
                        uint tNumEntities = mesh::number_of_edges( tBlock->element_type() );

                        // now we loop over all elements on the block and the selected number oedges
                        for( mesh::Element * tElement : tBlock->elements() )
                        {
                            // loop over all mesh entities on this element
                            for( uint k=0; k<tNumEntities; ++k )
                            {
                                // get the temporary index of this entity
                                index_t tIndex = tEntityMap( tElement->edge( k )->id() );

                                // grab the corresponding bitset
                                Bitset< BELFEM_MAX_DOFTYPES > & tBitset = tDofFlags( tIndex );

                                // flip the bits that represent the selected dofs
                                for( index_t s : tSelectedDofs )
                                {
                                    for( index_t i=0; i<tMultiplicity; ++i )
                                    {
                                        tBitset.set( s + i );
                                    }
                                }

                            }
                        }
                    }
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 4: count how many dofs have to be created
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                index_t aDofCount = 0 ;

                // we loop over all dof bitsets
                for( index_t k=0; k<tEntityCount; ++k )
                {
                    // grab the current bitset
                    Bitset< BELFEM_MAX_DOFTYPES > & tBitset = tDofFlags( k );

                    // count the dofs
                    aDofCount += tBitset.count() ;
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 5: populate the containers for the entity ids,
                //         the dof types and the used procs
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // allocate memory
                aEntityIDs.set_size( aDofCount );
                aDofTypes.set_size( aDofCount );


                // this is a temporary map which is needed to assign procs and dofs
                Map< luint, index_t > tDofMap ;


                // reset counters
                aDofCount = 0 ;

                index_t tIndex ;
                index_t tBitCountA ;
                index_t tBitCountB ;

                // loop over all flagged entities
                for( mesh::Edge * tEdge : mMesh->edges() )
                {
                    // check if edge is flagged
                    if( tEdge->is_flagged() )
                    {

                        // get index in array
                        tIndex = tEntityMap( tEdge->id() );

                        // get the corresponding bitset
                        Bitset< BELFEM_MAX_DOFTYPES > & tBitset = tDofFlags( tIndex );

                        tBitCountA = 0 ;
                        tBitCountB = tBitset.count() ;

                        // loop over all bits
                        for( uint d=0; d<BELFEM_MAX_DOFTYPES; ++d )
                        {

                            if( tBitset.test( d ) )
                            {
                                // remember ID
                                aEntityIDs( aDofCount ) = tEdge->id() ;

                                // remember dof type
                                aDofTypes( aDofCount )  = d ;


                                // create a unique and map it
                                tDofMap[ tEdge->id() * BELFEM_MAX_DOFTYPES + d ] = aDofCount++ ;

                                // cancel loop if we are done
                                if( ++tBitCountA == tBitCountB )
                                {
                                    break ;
                                }
                            }
                        }
                    }
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 6: determine which dofs are visible on which proc
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // allocate memory
                aProcFlags.set_size( aDofCount, Bitset< BELFEM_MAX_NUMPROCS >() );

                // loop over all blocks
                // define which dof is used by which proc
                for( id_t tBlockID : aIWG->selected_blocks() )
                {
                    // grab block on mesh
                    mesh::Block * tBlock = mMesh->block( tBlockID );

                    // grab the DOFs for this entity
                    const Vector< index_t > & tSelectedDofs = aIWG->dofs_per_edge( tBlockID );

                    if ( tSelectedDofs.length() > 0 )
                    {
                        // get the number of edges per element
                        uint tNumEntities = mesh::number_of_edges( tBlock->element_type() );

                        // loop over all elements on this block
                        for( mesh::Element * tElement : tBlock->elements() )
                        {
                            // loop over all entities on this element
                            for( uint k=0; k<tNumEntities; ++k )
                            {
                                // loop over all selected dofs
                                for( index_t s : tSelectedDofs )
                                {
                                    for( index_t i=0; i<tMultiplicity; ++i )
                                    {
                                        // compute the unique dof ID
                                        luint tID = tElement->edge( k )->id() * BELFEM_MAX_DOFTYPES + s + i ;

                                        // grab the corresponding bitset
                                        Bitset< BELFEM_MAX_NUMPROCS > & tBitset = aProcFlags( tDofMap( tID ));

                                        // set the proc bitset
                                        tBitset.set( tElement->owner() );
                                    }
                                }
                            }
                        }
                    }
                }

                // return the number of DOFs
                return aDofCount ;
            }

//------------------------------------------------------------------------

            index_t
            DofData::count_face_dofs( IWG  * aIWG,
                                      Vector< id_t >    & aEntityIDs,
                                      Vector< index_t > & aDofTypes,
                                      Cell< Bitset<BELFEM_MAX_NUMPROCS> > & aProcFlags )
            {
                /* depending on which function this is, "entity" refers to
                                  a node, an face, face or element

                                  this routine performs xxx steps:

                                  step 1 : flag all mesh entities that sit on the selected blocks

                                  step 2:  count the selected entities and create a lookup table

                                  step 3:  create the dof table: in order to find out which dofs exist, we create an
                                           array of bitsets and flip the bitsets for each mesh entity

                                  step 4:  count how many dofs have to be created

                                  step 5:  polpulate the containers for the entity ids and dof types

                                  step 6:  determine which dofs are visible on which proc
                                  */
                // unflag all entities on the mesh
                mMesh->unflag_all_faces() ;

                // we remember the ID of reach entity
                Vector< id_t > tEntityIDs ;

                // first, we need to figure out if entities dofs exist at all
                // if so, relevant entities are flagged
                bool tHaveEntityDofs = false ;

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 1: flag all mesh entities the dofs are based on
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                for( id_t tID : aIWG->selected_blocks() )
                {
                    if ( aIWG->number_of_dofs_per_face( tID ) > 0 )
                    {
                        tHaveEntityDofs = true ;
                        mMesh->block( tID )->flag_faces();
                    }
                }

                // exit the routine if no entity dofs exist
                if( ! tHaveEntityDofs )
                {
                    return 0 ;
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 2: count the selected entities and create tables
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // initialize counter
                index_t tEntityCount = 0 ;

                // the entity maps connects the ids to the counter
                Map< id_t, index_t > tEntityMap ;

                // now we count the number of flagged entities
                for( mesh::Face * tFace : mMesh->faces() )
                {
                    if( tFace->is_flagged() )
                    {
                        // write counter into map and increment it
                        tEntityMap[ tFace->id() ] = tEntityCount++ ;
                    }
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 3: flip entity-wise bitsets for used dofs
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // the bitset array
                Cell< Bitset<BELFEM_MAX_DOFTYPES> > tDofFlags( tEntityCount, Bitset<BELFEM_MAX_DOFTYPES>() );

                index_t tMultiplicity = aIWG->face_multiplicity();

                // we loop over all blocks and get the block-wise dofs
                for( id_t tBlockID: aIWG->selected_blocks() )
                {
                    // grab block pointer on mesh
                    mesh::Block * tBlock = mMesh->block( tBlockID );

                    // ask IWG about selected dofs
                    const Vector< index_t > & tSelectedDofs = aIWG->dofs_per_face( tBlockID );

                    // check if any dofs are selected
                    if( tSelectedDofs.length() > 0 )
                    {
                        // get the number of faces per element on this block
                        uint tNumEntities = mesh::number_of_faces( tBlock->element_type() );


                        // now we loop over all elements on the block and the selected number ofaces
                        for( mesh::Element * tElement : tBlock->elements() )
                        {
                            // loop over all mesh entities on this element
                            for( uint k=0; k<tNumEntities; ++k )
                            {
                                // get the temporary index of this entity
                                index_t tIndex = tEntityMap( tElement->face( k )->id() );

                                // grab the corresponding bitset
                                Bitset< BELFEM_MAX_DOFTYPES > & tBitset = tDofFlags( tIndex );

                                // flip the bits that represent the selected dofs
                                for( index_t s : tSelectedDofs )
                                {
                                    for( index_t i=0; i<tMultiplicity; ++i )
                                    {
                                        tBitset.set( s + i );
                                    }
                                }
                            }
                        }
                    }
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 4: count how many dofs have to be created
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                index_t aDofCount = 0 ;

                // we loop over all dof bitsets
                for( index_t k=0; k<tEntityCount; ++k )
                {
                    // grab the current bitset
                    Bitset< BELFEM_MAX_DOFTYPES > & tBitset = tDofFlags( k );

                    // count the dofs
                    aDofCount += tBitset.count() ;
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 5: populate the containers for the entity ids,
                //         the dof types and the used procs
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // allocate memory
                aEntityIDs.set_size( aDofCount );
                aDofTypes.set_size( aDofCount );


                // this is a temporary map which is needed to assign procs and dofs
                Map< luint, index_t > tDofMap ;


                // reset counters
                aDofCount = 0 ;

                index_t tIndex ;
                index_t tBitCountA ;
                index_t tBitCountB ;

                // loop over all flagged entities
                for( mesh::Face * tFace : mMesh->faces() )
                {
                    // check if face is flagged
                    if( tFace->is_flagged() )
                    {
                        // get index in array
                        tIndex = tEntityMap( tFace->id() );

                        // get the corresponding bitset
                        Bitset< BELFEM_MAX_DOFTYPES > & tBitset = tDofFlags( tIndex );

                        tBitCountA = 0 ;
                        tBitCountB = tBitset.count() ;

                        // loop over all bits
                        for( uint d=0; d<BELFEM_MAX_DOFTYPES; ++d )
                        {
                            if( tBitset.test( d ) )
                            {
                                // remember ID
                                aEntityIDs( aDofCount ) = tFace->id() ;

                                // remember dof type
                                aDofTypes( aDofCount )  = d ;

                                // create a unique id and map it
                                tDofMap[ tFace->id() * BELFEM_MAX_DOFTYPES + d ] = aDofCount++ ;

                                // cancel loop if we are done
                                if( ++tBitCountA == tBitCountB )
                                {
                                    break ;
                                }
                            }
                        }
                    }
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 6: determine which dofs are visible on which proc
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // allocate memory
                aProcFlags.set_size( aDofCount, Bitset< BELFEM_MAX_NUMPROCS >() );

                // loop over all blocks
                // define which dof is used by which proc
                for( id_t tBlockID : aIWG->selected_blocks() )
                {
                    // grab block on mesh
                    mesh::Block * tBlock = mMesh->block( tBlockID );

                    // grab the DOFs for this entity
                    const Vector< index_t > & tSelectedDofs = aIWG->dofs_per_face( tBlockID );

                    if ( tSelectedDofs.length() > 0 )
                    {
                        // get the number of faces per element
                        uint tNumEntities = mesh::number_of_faces( tBlock->element_type() );

                        // loop over all elements on this block
                        for( mesh::Element * tElement : tBlock->elements() )
                        {
                            // loop over all entities on this element
                            for( uint k=0; k<tNumEntities; ++k )
                            {
                                // loop over all selected dofs
                                for( index_t s : tSelectedDofs )
                                {
                                    for( index_t i=0; i<tMultiplicity; ++i )
                                    {
                                        // compute the unique dof ID
                                        luint tID = tElement->face( k )->id() * BELFEM_MAX_DOFTYPES + s + i ;

                                        // grab the corresponding bitset
                                        Bitset< BELFEM_MAX_NUMPROCS > & tBitset = aProcFlags( tDofMap( tID ));

                                        // set the proc bitset
                                        tBitset.set( tElement->owner() );
                                    }
                                }
                            }
                        }
                    }
                }

                // return the number of DOFs
                return aDofCount ;
            }

//------------------------------------------------------------------------

            index_t
            DofData::count_cell_dofs( IWG  * aIWG,
                                      Vector< id_t >    & aEntityIDs,
                                      Vector< index_t > & aDofTypes,
                                      Cell< Bitset<BELFEM_MAX_NUMPROCS> > & aProcFlags )
            {
                /* depending on which function this is, "entity" refers to
                                  a node, an face, face or element

                                  this routine performs xxx steps:

                                  step 1 : flag all mesh entities that sit on the selected blocks

                                  step 2:  count the selected entities and create a lookup table

                                  step 3:  create the dof table: in order to find out which dofs exist, we create an
                                           array of bitsets and flip the bitsets for each mesh entity

                                  step 4:  count how many dofs have to be created

                                  step 5:  polpulate the containers for the entity ids and dof types

                                  step 6:  determine which dofs are visible on which proc */

                // unflag all entities on the mesh
                mMesh->unflag_all_elements() ;

                // we remember the ID of reach entity
                Vector< id_t > tEntityIDs ;

                // first, we need to figure out if entities dofs exist at all
                // if so, relevant entities are flagged
                bool tHaveEntityDofs = false ;

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 1: flag all mesh entities the dofs are based on
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                for( id_t tID : aIWG->selected_blocks() )
                {
                    if ( aIWG->number_of_dofs_per_cell( tID ) > 0 )
                    {
                        tHaveEntityDofs = true ;
                        mMesh->block( tID )->flag_elements();
                    }
                }

                // exit the routine if no entity dofs exist
                if( ! tHaveEntityDofs )
                {
                    return 0 ;
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 2: count the selected entities and create tables
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // initialize counter
                index_t tEntityCount = 0 ;

                // the entity maps connects the ids to the counter
                Map< id_t, index_t > tEntityMap ;

                // get element list
                Cell< mesh::Element * > & tElements
                        = mMesh->elements() ;

                // now we count the number of flagged entities
                for( mesh::Element * tElement : tElements )
                {
                    if( tElement->is_flagged() )
                    {
                        // write counter into map and increment it
                        tEntityMap[ tElement->id() ] = tEntityCount++ ;
                    }
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 3: flip entity-wise bitsets for used dofs
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // the bitset array
                Cell< Bitset<BELFEM_MAX_DOFTYPES> > tDofFlags( tEntityCount, Bitset<BELFEM_MAX_DOFTYPES>() );

                index_t tMultiplicity = aIWG->cell_multiplicity();

                // we loop over all blocks and get the block-wise dofs
                for( id_t tBlockID: aIWG->selected_blocks() )
                {
                    // grab block pointer on mesh
                    mesh::Block * tBlock = mMesh->block( tBlockID );

                    // ask IWG about selected dofs
                    const Vector< index_t > & tSelectedDofs = aIWG->dofs_per_cell( tBlockID );

                    // check if any dofs are selected
                    if( tSelectedDofs.length() > 0 )
                    {
                        // now we loop over all elements on the block and the selected number ofaces
                        for( mesh::Element * tElement : tBlock->elements() )
                        {
                            // get the temporary index of this entity
                            index_t tIndex = tEntityMap( tElement->id() );

                            // grab the corresponding bitset
                            Bitset< BELFEM_MAX_DOFTYPES > & tBitset = tDofFlags( tIndex );

                            // flip the bits that represent the selected dofs
                            for( index_t s : tSelectedDofs )
                            {
                                for( index_t i=0; i<tMultiplicity; ++i )
                                {
                                    tBitset.set( s + i );
                                }
                            }
                        }
                    }
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 4: count how many dofs have to be created
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                index_t aDofCount = 0 ;

                // we loop over all dof bitsets
                for( index_t k=0; k<tEntityCount; ++k )
                {
                    // grab the current bitset
                    Bitset< BELFEM_MAX_DOFTYPES > & tBitset = tDofFlags( k );

                    // count the dofs
                    aDofCount += tBitset.count() ;
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 5: populate the containers for the entity ids,
                //         the dof types and the used procs
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // allocate memory
                aEntityIDs.set_size( aDofCount );
                aDofTypes.set_size( aDofCount );

                // this is a temporary map which is needed to assign procs and dofs
                Map< luint, index_t > tDofMap ;

                // reset counters
                aDofCount = 0 ;

                index_t tIndex ;
                index_t tBitCountA ;
                index_t tBitCountB ;

                // loop over all flagged entities
                for( mesh::Element * tElement : tElements )
                {
                    // check if face is flagged
                    if( tElement->is_flagged() )
                    {
                        // get index in array
                        tIndex = tEntityMap( tElement->id() );

                        // get the corresponding bitset
                        Bitset< BELFEM_MAX_DOFTYPES > & tBitset = tDofFlags( tIndex );

                        tBitCountA = 0 ;
                        tBitCountB = tBitset.count() ;

                        // loop over all bits
                        for( uint d=0; d<BELFEM_MAX_DOFTYPES; ++d )
                        {
                            if( tBitset.test( d ) )
                            {
                                // remember ID
                                aEntityIDs( aDofCount ) = tElement->id() ;

                                // remember dof type
                                aDofTypes( aDofCount )  = d ;

                                // create a unique id and map it
                                tDofMap[ tElement->id() * BELFEM_MAX_DOFTYPES + d ] = aDofCount++ ;

                                // cancel loop if we are done
                                if( ++tBitCountA == tBitCountB )
                                {
                                    break ;
                                }
                            }
                        }
                    }
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 6: determine which dofs are visible on which proc
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // allocate memory
                aProcFlags.set_size( aDofCount, Bitset< BELFEM_MAX_NUMPROCS >() );

                // loop over all blocks
                // define which dof is used by which proc
                for( id_t tBlockID : aIWG->selected_blocks() )
                {
                    // grab block on mesh
                    mesh::Block * tBlock = mMesh->block( tBlockID );

                    // grab the DOFs for this entity
                    const Vector< index_t > & tSelectedDofs = aIWG->dofs_per_cell( tBlockID );

                    if ( tSelectedDofs.length() > 0 )
                    {
                        // loop over all elements on this block
                        for( mesh::Element * tElement : tBlock->elements() )
                        {
                            // loop over all selected dofs
                            for( index_t s : tSelectedDofs )
                            {
                                for( index_t i=0; i<tMultiplicity; ++i )
                                {
                                    // compute the unique dof ID
                                    luint tID = tElement->id() * BELFEM_MAX_DOFTYPES + s + i ;

                                    // grab the corresponding bitset
                                    Bitset< BELFEM_MAX_NUMPROCS > & tBitset = aProcFlags( tDofMap( tID ) );

                                    // set the proc bitset
                                    tBitset.set( tElement->owner() );

                                    // also make element visible of neighbor procs
                                    for( uint e=0; e<tElement->number_of_elements(); ++e )
                                    {
                                        // get neighbor
                                        mesh::Element * tNeighbor = tElement->element( e );

                                        if( tNeighbor->owner() != tElement->owner() && tNeighbor->is_flagged() )
                                        {
                                            tBitset.set( tNeighbor->owner() );
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // return the number of DOFs
                return aDofCount ;
            }

//------------------------------------------------------------------------

            index_t
            DofData::count_lambda_dofs( IWG  * aIWG,
                                      Vector< id_t >    & aEntityIDs,
                                      Vector< index_t > & aDofTypes,
                                      Cell< Bitset<BELFEM_MAX_NUMPROCS> > & aProcFlags )
            {
                /* depending on which function this is, "entity" refers to
                                   a node, an edge, face or element

                                   this routine performs xxx steps:

                                   step 1 : flag all mesh entities that sit on the selected blocks

                                   step 2:  count the selected entities and create a lookup table

                                   step 3:  create the dof table: in order to find out which dofs exist, we create an
                                            array of bitsets and flip the bitsets for each mesh entity

                                   step 4:  count how many dofs have to be created

                                   step 5:  polpulate the containers for the entity ids and dof types

                                   step 6:  determine which dofs are visible on which proc
                                   */

                // unflag all entities on the mesh
                mMesh->unflag_all_facets() ;

                // we remember the ID of reach entity
                Vector< id_t > tEntityIDs ;

                // first, we need to figure out if entities dofs exist at all
                // if so, relevant entities are flagged
                bool tHaveEntityDofs = false ;

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 1: flag all mesh entities the dofs are based on
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                for( id_t tID : aIWG->selected_sidesets() )
                {
                    // grab sideset on mesh
                    mesh::SideSet * tSideSet = mMesh->sideset( tID );
                    if( aIWG->number_of_lambda_dofs( tID ) > 0 )
                    {
                        tHaveEntityDofs = true ;
                        tSideSet->flag_all_facets() ;
                    }
                }

                // exit the routine if no entity dofs exist
                if( ! tHaveEntityDofs )
                {
                    return 0 ;
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 2: count the selected entities and create tables
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // initialize counter
                index_t tEntityCount = 0 ;

                // the entity maps connects the ids to the counter
                Map< id_t, index_t > tEntityMap ;

                // now we count the number of flagged entities
                for( mesh::Facet * tFacet : mMesh->facets() )
                {
                    if( tFacet->is_flagged() )
                    {
                        // write counter into map and increment it
                        tEntityMap[ tFacet->id() ] = tEntityCount++ ;
                    }
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 3: flip entity-wise bitsets for used dofs
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // the bitset array
                Cell< Bitset<BELFEM_MAX_DOFTYPES> > tDofFlags( tEntityCount, Bitset<BELFEM_MAX_DOFTYPES>() );

                // we loop over all blocks and get the block-wise dofs
                for( id_t tSideSetID : aIWG->selected_sidesets() )
                {
                    // grab block pointer on mesh
                    mesh::SideSet * tSideSet = mMesh->sideset( tSideSetID );

                    // ask IWG about selected dofs
                    const Vector< index_t > & tSelectedDofs = aIWG->lambda_dofs( tSideSetID );


                    // check if any dofs are selected
                    if( tSelectedDofs.length() > 0 )
                    {
                        for( mesh::Facet * tFacet : tSideSet->facets() )
                        {
                             // get the temporary index of this entity
                            index_t tIndex = tEntityMap( tFacet->id() );

                            // grab the corresponding bitset
                            Bitset< BELFEM_MAX_DOFTYPES > & tBitset = tDofFlags( tIndex );

                            // flip the bits that represent the selected dofs
                            for( index_t d : tSelectedDofs )
                            {
                                tBitset.set( d );
                            }
                        }
                    }
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 4: count how many dofs have to be created
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                index_t aDofCount = 0 ;

                // we loop over all dof bitsets
                for( index_t k=0; k<tEntityCount; ++k )
                {
                    // grab the current bitset
                    Bitset< BELFEM_MAX_DOFTYPES > & tBitset = tDofFlags( k );

                    // count the dofs
                    aDofCount += tBitset.count() ;
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 5: populate the containers for the entity ids,
                //         the dof types and the used procs
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // allocate memory
                aEntityIDs.set_size( aDofCount );
                aDofTypes.set_size( aDofCount );

                // this is a temporary map which is needed to assign procs and dofs
                Map< luint, index_t > tDofMap ;

                // reset counters
                aDofCount = 0 ;

                index_t tIndex ;
                index_t tBitCountA ;
                index_t tBitCountB ;

                // loop over all flagged entities
                Cell< mesh::Facet * > & tFacets = mMesh->facets();

                for ( mesh::Facet * tFacet: tFacets )
                {
                    // check if facet is flagged
                    if ( tFacet->is_flagged())
                    {
                        // get index in array
                        tIndex = tEntityMap( tFacet->id());

                        // get the corresponding bitset
                        Bitset< BELFEM_MAX_DOFTYPES > & tBitset = tDofFlags( tIndex );

                        tBitCountA = 0;
                        tBitCountB = tBitset.count();

                        // loop over all bits
                        for ( uint d = 0; d < BELFEM_MAX_DOFTYPES; ++d )
                        {
                            if ( tBitset.test( d ))
                            {
                                // remember ID
                                aEntityIDs( aDofCount ) = tFacet->id();

                                // remember dof type
                                aDofTypes( aDofCount ) = d;

                                // create a unique and map it
                                tDofMap[ tFacet->id() * BELFEM_MAX_DOFTYPES + d ] = aDofCount++;

                                // cancel loop if we are done
                                if ( ++tBitCountA == tBitCountB )
                                {
                                    break;
                                }
                            }
                        }
                    }
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // Step 6: determine which dofs are visible on which proc
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // allocate memory
                aProcFlags.set_size( aDofCount, Bitset< BELFEM_MAX_NUMPROCS >() );

                // loop over all sidesets
                for( id_t tSideSetId : aIWG->selected_sidesets() )
                {
                    // grab sideset on mesh
                    mesh::SideSet * tSideSet = mMesh->sideset( tSideSetId );

                    // grab the DOFs for this entity
                    const Vector< index_t > & tSelectedDofs = aIWG->lambda_dofs( tSideSetId );

                    if ( tSelectedDofs.length() > 0 )
                    {
                        // loop over all elements on this block
                        for( mesh::Facet * tFacet : tSideSet->facets() )
                        {
                            // loop over all selected dofs
                            for( index_t d : tSelectedDofs )
                            {
                                // compute the unique dof ID
                                luint tID = tFacet->id() * BELFEM_MAX_DOFTYPES + d ;

                                // grab the corresponding bitset
                                Bitset< BELFEM_MAX_NUMPROCS > & tBitset = aProcFlags( tDofMap( tID ) );

                                // set the proc bitset
                                tBitset.set( tFacet->owner() );
                            }
                        }
                    }
                }

                // return the number of DOFs
                return aDofCount ;
            }
//------------------------------------------------------------------------------

            index_t
            DofData::count_dofs_for_proc(
                    const index_t                               aProc,
                    const Vector< id_t >                      & aDofIDs,
                    const Vector< id_t >                      & aEntityIDs,
                    const Vector< index_t >                   & aDofTypes,
                    const Cell< Bitset<BELFEM_MAX_NUMPROCS> > & aProcFlags,
                    Vector< id_t >                            & aProcDofIDs,
                    Vector< id_t >                            & aProcEntityIDs,
                    Vector< index_t >                         & aProcDofTypes  )
            {
                // get the number of entitues
                index_t tNumEntities = aEntityIDs.length() ;

                // counter for entities
                index_t aCount = 0 ;

                // loop over all entities
                for( index_t k=0; k<tNumEntities; ++k )
                {
                    // grab bitset
                    const Bitset<BELFEM_MAX_NUMPROCS> & tBitset = aProcFlags( k );

                    // check if dof is used on proc
                    if ( tBitset.test( aProc ) )
                    {
                        // increment the counter
                        ++aCount ;
                    }
                }

                // allocate memory
                aProcDofIDs.set_size( aCount );
                aProcEntityIDs.set_size( aCount );
                aProcDofTypes.set_size( aCount );

                // cancel if there is nothing to do
                if( aCount == 0 )
                {
                    return 0 ;
                }

                // reset the counter
                aCount = 0;
                for( index_t k=0; k<tNumEntities; ++k )
                {
                    // grab bitset
                    const Bitset<BELFEM_MAX_NUMPROCS> & tBitset = aProcFlags( k );

                    // check if dof is used on proc
                    if ( tBitset.test( aProc ) )
                    {
                        // add data to list
                        aProcDofIDs( aCount )    = aDofIDs( k );
                        aProcEntityIDs( aCount ) = aEntityIDs( k );
                        aProcDofTypes( aCount )  = aDofTypes( k );

                        // increment the counter
                        ++aCount;
                    }
                }

                return aCount ;
            }

//------------------------------------------------------------------------------


            void
            DofData::collect_hanging_dofs()
            {
                for( Dof * tDof : mDOFs )
                {
                    tDof->unflag() ;
                }

                if( mCommRank == 0 )
                {

                    // special treatment for second order faces
                    if (   mMesh->faces_exist()
                        && mMesh->number_of_dimensions() == 3 )
                    {
                        mMesh->unflag_all_facets();

                        // flag faces that belong to interfaces
                        for ( mesh::SideSet * tSideSet : mMesh->sidesets() )
                        {
                            if ( tSideSet->number_of_facets() == 0 ) continue ;

                            // get first facet
                            mesh::Facet * tFacet0 = tSideSet->facet_by_index( 0 );

                            if ( ! tFacet0->has_slave() ) continue ;

                            // todo: might have to consider periodic BC condition here too
                            // todo: also would have to add thin shells here
                            Cell< mesh::Node * > tNodesInSlaveOrientation( tFacet0->number_of_nodes(), nullptr );
                            Cell< mesh::Node * > tNodesInMasterOrientation( tFacet0->number_of_nodes(), nullptr );
                            Cell< mesh::Node * > tSources( 2 * tFacet0->number_of_nodes(), nullptr );

                            Vector< real > tWeights( 2 * tFacet0->number_of_nodes(), 0.0 );
                            for ( uint k=0; k<tFacet0->number_of_nodes(); ++k )
                            {
                                tWeights( k ) = BELFEM_QUIET_NAN ;
                            }
                            if ( tFacet0->master()->has_faces() && ! tFacet0->slave()->has_faces() )
                            {
                                for ( mesh::Facet * tFacet : tSideSet->facets() )
                                {

                                    // first, we grab the nodes from the slave, we use tMasterNodes as temporary container
                                    tFacet->slave()->get_nodes_of_facet( tFacet->index_on_slave(), tNodesInSlaveOrientation );

                                    // now we align the nodes in the orientation of the master, after that, we can overwrite tMasterNodes
                                    mesh::to_master_orientation( tFacet, tNodesInSlaveOrientation, tNodesInMasterOrientation );

                                    // collect nodes
                                    for ( uint k=0; k<tFacet->number_of_sources(); ++k )
                                    {
                                        tWeights( k ) = tFacet->weight( k );
                                    }
                                    for ( uint k=0; k<tFacet->number_of_nodes(); ++k )
                                    {
                                        tSources( k ) = tNodesInMasterOrientation( k );
                                        tSources( k + tFacet->number_of_nodes() ) = tNodesInMasterOrientation( k );
                                    }

                                    tFacet->master()->face( tFacet->index_on_master() )->set_sources( tSources, tWeights );
                                }
                            }
                        }
                    }

                    // first, we count the hanging dofs.
                    // note that they are not yet linked to the other dofs
                    // so we need to identify the over their mesh bases

                    DynamicBitset * tBefore = new DynamicBitset( mDOFs.size() );

                    mNumberOfHangingDofs = 0 ;
                    mNumberOfFixedDofs = 0 ;
                    mNumberOfFreeDofs = 0 ;
                    index_t tCount = 0 ;
                    index_t tIndex = 0 ;
                    for ( Dof * tDof : mDOFs )
                    {
                        // check if basis is hanging
                        if ( tDof->basis_is_hanging() )
                        {
                            ++mNumberOfHangingDofs ;
                            tDof->set_index( gNoIndex );
                            tBefore->set( tIndex );
                        }
                        else
                        {
                            tDof->set_index( tCount++ );
                        }

                        ++tIndex ;
                    }
                    tBefore->lock();

                    Cell< index_t > tIndices ;
                    tBefore->where( tIndices );
                    mNumberOfHangingDofs = tIndices.size() ;
                    mHangingDOFs.set_size( mNumberOfHangingDofs, nullptr );

                    tCount = 0 ;
                    for ( index_t k : tIndices )
                    {
                        mHangingDOFs( tCount++ ) = mDOFs( k );
                    }


                    // creates the sources and weights of hanging dofs
                    // based on mesh information
                    this->create_dofwise_t_matrices_master() ;

                    // Remove DoFs from mHangingDOFs that ended up not being hanging
                    // (sources couldn't be set because of missing matching DoF types)

                    DynamicBitset * tAfter = new DynamicBitset( mDOFs.size() );

                    tIndex = 0 ;
                    tCount = 0 ;

                    for ( Dof * tDof : mDOFs )
                    {
                        // check if dof is hanging
                        if ( tDof->is_hanging() )
                        {
                            tDof->set_index( gNoIndex );
                            tAfter->set( tIndex );
                        }
                        else
                        {
                            tDof->set_index( tCount++ );
                        }

                        ++tIndex ;
                    }
                    tAfter->lock();
                    mNumberOfHangingDofs = mDOFs.size() - tCount ;


                    if ( *tBefore != *tAfter )
                    {
                        tAfter->where(tIndices  );

                        mHangingDOFs.set_size( mNumberOfHangingDofs, nullptr );

                        tCount = 0 ;
                        for ( index_t k : tIndices )
                        {
                            Dof * tDof = mDOFs( k ) ;
                            mHangingDOFs( tCount++ ) = tDof ;
                        }
                        mHangingDOFs.shrink_to_fit();
                    }

                    delete tBefore ;

                    // removing the hanging dofs from the dof container
                    if( comm_size() < 2 )
                    {
                        tAfter->unlock();
                        tAfter->flip();
                        tAfter->where(tIndices  );

                        Cell< Dof * > tDofs = std::move( mDOFs );
                        mDOFs.set_size( tIndices.size(), nullptr );
                        tCount = 0 ;
                        for ( index_t k : tIndices )
                        {
                            mDOFs( tCount++ ) = tDofs( k ) ;
                        }
                        mDOFs.shrink_to_fit();
                        delete tAfter ;
                        return ;
                    }

                    delete tAfter ;

                    // having the dofs, we can now create am id list for the other procsmKernel->comm_table()
                    tCount = 0 ;
                    Vector< id_t > tMasterHangingDofIDs( mNumberOfHangingDofs );
                    for( Dof * tDof : mHangingDOFs )
                    {
                        tMasterHangingDofIDs( tCount++ ) = tDof->id() ;
                    }

                    comm_barrier() ; // barrier 1

                    // Send the master hanging dof IDs to all other procs

                    proc_t tCommSize = comm_size() ;

                    for ( proc_t p=1; p<tCommSize; ++p )
                    {
                        send( tMasterHangingDofIDs, p );
                    }

                    comm_barrier() ; // barrier 2

                    uint n = comm_size() ;

                    Cell< Vector< id_t > > tAllHangingDofIDs( n, {} );
                    collect( tAllHangingDofIDs );

                    // now we crate a list with all source dofs that are needed per proc

                    Cell< Vector< id_t > > tAllSourceDofIDs( n, {} );

                    for ( uint p=0; p<n; ++p )
                    {
                        // get the hanging dof IDs
                        Vector< id_t > & tHangingDofIDs = tAllHangingDofIDs( p );

                        // get the container for the source dof IDs
                        Vector< id_t > & tSourceDofIDs = tAllSourceDofIDs( p );

                        // count memory needs
                        tCount = 0 ;
                        for ( id_t tID : tHangingDofIDs )
                        {
                            tCount += mDofMap( tID )->number_of_sources() ;
                        }

                        // populate the data
                        tSourceDofIDs.set_size( tCount );
                        tCount = 0 ;
                        for ( id_t tID : tHangingDofIDs )
                        {
                            // get the dof
                            Dof * tDof = mDofMap( tID ) ;
                            for ( uint k=0; k<tDof->number_of_sources(); ++k )
                            {
                                tSourceDofIDs( tCount++ ) = tDof->source( k )->id() ;
                            }
                        }

                        unique( tSourceDofIDs );
                    }

                    comm_barrier() ; // barrier 3
                    distribute( tAllSourceDofIDs );

                    // reset the tables
                    for ( uint p=0; p<n; ++p )
                    {
                        tAllSourceDofIDs( p ) = {} ;
                    }

                    comm_barrier() ; // barrier 4

                    collect( tAllSourceDofIDs );

                    Cell< Vector< id_t > > tAllIdata( n, {} );

                    for ( uint p=0; p<n; ++p )
                    {
                        // get the container for the source dof IDs
                        Vector< id_t > & tSourceDofIDs = tAllSourceDofIDs( p );

                        Vector< id_t > & tIdata = tAllIdata( p );
                        tIdata.set_size( tSourceDofIDs.length() * 5 );

                        // count memory needs
                        tCount = 0 ;
                        for ( id_t tID : tSourceDofIDs )
                        {
                            Dof * tDof = mDofMap( tID ) ;
                            tIdata( tCount++ ) = tDof->index() ;
                            tIdata( tCount++ ) = tDof->type_id() ;
                            tIdata( tCount++ ) = static_cast< id_t >( tDof->entity_type() );
                            tIdata( tCount++ ) = tDof->index_on_entity() ;
                            tIdata( tCount++ ) = tDof->mesh_basis()->id() ;
                        }
                    }

                    comm_barrier() ; // barrier 5
                    distribute( tAllIdata );

                    // reset containers
                    for ( uint p=0; p<n; ++p )
                    {
                        tAllIdata( p ) = {} ;
                    }

                    Cell< Vector< real > > tAllRdata( n, {} );

                    for ( uint p=0; p<n; ++p )
                    {
                        // get the list of hanging dofs for this proc
                        Vector< id_t > & tHangingDofIDs = tAllHangingDofIDs( p );
                        Vector< id_t > & tIdata = tAllIdata( p );
                        Vector< real > & tRdata = tAllRdata( p );

                        index_t tCountI = 0 ;
                        index_t tCountR = 0 ;

                        for ( id_t tID : tHangingDofIDs )
                        {
                            // get the dof
                            Dof * tDof = mDofMap( tID ) ;

                            tCountI += tDof->number_of_sources() + 2 ;
                            tCountR += tDof->number_of_sources() ;
                        }

                        tIdata.set_size( tCountI );
                        tRdata.set_size( tCountR );

                        tCountI = 0 ;
                        tCountR = 0 ;

                        for ( id_t tID : tHangingDofIDs )
                        {
                            // get the dof
                            Dof * tDof = mDofMap( tID ) ;

                            tIdata( tCountI++ ) = tDof->id();
                            tIdata( tCountI++ ) = tDof->number_of_sources() ;

                            for ( uint k=0; k<tDof->number_of_sources(); ++k )
                            {
                                tIdata( tCountI++ ) = tDof->source( k )->id() ;
                                tRdata( tCountR++ ) = tDof->weight( k );
                            }
                        }
                    }

                    comm_barrier() ; // barrier 6
                    distribute( tAllIdata );
                    distribute( tAllRdata );
                }
                else
                {
                    comm_barrier() ; // barrier 1

                    // receive list of all hanging dofs from the master
                    Vector< id_t > tMasterHangingDofIDs ;
                    receive( tMasterHangingDofIDs );
                    // now we check which dofs we have
                    mNumberOfHangingDofs = 0 ;

                    for ( id_t tID : tMasterHangingDofIDs )
                    {
                        if ( mDofMap.key_exists( tID ) )
                        {
                            // flag this dof
                            mDofMap( tID )->flag() ;

                            // increment hanging dof counter
                            ++mNumberOfHangingDofs ;
                        }
                    }

                    // next, we collect all hanging dofs
                    mHangingDOFs.set_size( mNumberOfHangingDofs, nullptr );

                    // container for ids
                    Vector< id_t > tMyHangingDofIDs( mNumberOfHangingDofs );
                    index_t tCount = 0 ;
                    for( Dof * tDof : mDOFs )
                    {
                        if( tDof->is_flagged() )
                        {
                            // unflag the dof and add it to the container
                            tDof->unflag() ;
                            tDof->set_index( gNoIndex );
                            tMyHangingDofIDs( tCount ) = tDof->id() ;
                            mHangingDOFs( tCount++ ) = tDof ;
                        }
                    }

                    comm_barrier() ; // barrier 2

                    send( tMyHangingDofIDs );

                    comm_barrier() ; // barrier 3

                    Vector< id_t > tMasterSourceDofIDs ;
                    receive( tMasterSourceDofIDs );

                    // having all sources from the master that are relevant, we can now
                    // check which ones we need

                    tCount = 0 ;
                    for ( id_t tID : tMasterSourceDofIDs )
                    {
                        if ( ! mDofMap.key_exists( tID ) )
                        {
                            ++tCount ;
                        }
                    }

                    // populate the list of required new source dofs
                    Vector< id_t > tMyNewSourceDofIDs( tCount );
                    tCount = 0 ;
                    for ( id_t tID : tMasterSourceDofIDs )
                    {
                        if ( ! mDofMap.key_exists( tID ) )
                        {
                            tMyNewSourceDofIDs( tCount++ ) = tID ;
                        }
                    }

                    comm_barrier() ; // barrier 4

                    send( tMyNewSourceDofIDs );

                    comm_barrier() ; // barrier 5

                    Vector< id_t > tIdata ;
                    receive( tIdata );
                    tCount = 0 ;
                    index_t tICount = 0 ;

                    // create the new dofs
                    Cell< Dof * > tNewDofs( tMyNewSourceDofIDs.length(), nullptr );

                    for ( id_t tID : tMyNewSourceDofIDs )
                    {
                        index_t tIndex      = tIdata( tICount++ );
                        uint tType          = tIdata( tICount++ );
                        EntityType tEntity  = static_cast< EntityType >( tIdata( tICount++ ));
                        uint tIndexOnEntity = tIdata( tICount++ );
                        id_t tMeshID        = tIdata( tICount++ );

                        switch ( tEntity )
                        {
                            case EntityType::NODE :
                            {
                                tNewDofs( tCount ) = new Dof( tID, tType, mMesh->node( tMeshID ) );
                                break ;
                            }
                            case EntityType::EDGE :
                            {
                                mesh::Edge * tEdge = mMesh->edge( tMeshID );
                                index_t tIndexOnMesh = tEdge->index() * mParent->iwg()->edge_multiplicity() + tIndexOnEntity ;
                                tNewDofs( tCount ) = new Dof( tID, tType, tEdge, tIndexOnEntity, tIndexOnMesh );
                                break ;
                            }
                            case EntityType::FACE :
                            {
                                mesh::Face * tFace = mMesh->face( tMeshID );
                                index_t tIndexOnMesh = tFace->index() * mParent->iwg()->face_multiplicity() + tIndexOnEntity ;
                                tNewDofs( tCount ) = new Dof( tID, tType, tFace, tIndexOnEntity, tIndexOnMesh );
                                break ;
                            }
                            case EntityType::CELL :
                            case EntityType::ELEMENT :
                            {
                                mesh::Element * tElement = mMesh->element( tMeshID );
                                index_t tIndexOnMesh = tElement->index() * mParent->iwg()->cell_multiplicity() + tIndexOnEntity ;
                                tNewDofs( tCount ) = new Dof( tID, tType, tElement, tIndexOnEntity, tIndexOnMesh );
                                break ;
                            }
                            case EntityType::FACET :
                            {
                                mesh::Facet * tFacet = mMesh->facet( tMeshID );
                                index_t tIndexOnMesh = tFacet->index() * mParent->iwg()->lambda_multiplicity() + tIndexOnEntity ;
                                tNewDofs( tCount ) = new Dof( tID, tType, tFacet, tIndexOnEntity, tIndexOnMesh );
                                break ;
                            }
                            default :
                            {
                                BELFEM_ERROR( false, "don't know how to create source dof");
                            }
                        }

                        tNewDofs( tCount )->set_index( tIndex );

                        // add dof to map
                        mDofMap[ tID ] = tNewDofs( tCount++ ) ;

                    }

                    // add new dofs to dof container
                    append( mDOFs, tNewDofs );

                    comm_barrier() ; // barrier 6

                    tIdata = {};
                    receive( tIdata );

                    Vector< real > tRdata ;
                    receive( tRdata );

                    tICount = 0 ;
                    index_t tRCount = 0 ;

                    Cell< Dof * > tSources ;
                    Vector< real > tCoeffs ;

                    // set coefficients for hanging dofs
                    while ( tICount < tIdata.length() )
                    {
                        Dof * tDof = mDofMap[ tIdata( tICount++ ) ] ;

                        uint tNumberOfSources = tIdata( tICount++ ) ;
                        tSources.set_size( tNumberOfSources, nullptr );
                        tCoeffs.set_size( tNumberOfSources, 0.0 );

                        for ( uint k=0; k<tNumberOfSources; ++k )
                        {
                            tSources( k ) = mDofMap[ tIdata( tICount++ ) ] ;
                            tCoeffs( k )  = tRdata( tRCount++ ) ;
                        }

                        tDof->set_sources( tSources, tCoeffs );
                    }

                }

                this->remove_hanging_dofs_from_container() ;

                // barrier 7
                comm_barrier() ;
            }

//------------------------------------------------------------------------------

            void
            DofData::init_dof_values( const bool aFreeDofsOnly )
            {
                // collect field list from IWG
                const Cell< string > & tLabels = mParent->iwg()->dof_fields();

                uint tNumFields = tLabels.size() ;

                // unflag all dofs
                for( Dof * tDof : mDOFs )
                {
                    tDof->unflag() ;
                }

                // loop over all number of dofs
                for( uint k=0; k<tNumFields; ++k )
                {
                    // grab field
                    mesh::Field * tField = mMesh->field( tLabels( k ) );

                    Vector< real > & tValues = tField->data();

                    // loop over all dofs
                    for ( Dof * tDof : mDOFs )
                    {
                        if ( !tDof->is_flagged() && tDof->entity_type() == tField->entity_type() )
                        {
                            BELFEM_ASSERT( tDof->dof_index_on_field() < tValues.length(),
                                          "Field index %lu of dof %lu ( %s : %lu ) of %s field %s is out of bounds (expect < %lu )",
                                          ( long unsigned int ) tDof->dof_index_on_field(),
                                          ( long unsigned int ) tDof->id(),
                                          to_string( tDof->entity_type()).c_str(),
                                          ( long unsigned int ) tDof->mesh_basis()->id(),
                                          to_string( tField->entity_type()).c_str(),
                                          tLabels( k ).c_str(),
                                          ( long unsigned int ) tValues.length() );

                            // check of type equals index
                            if ( tDof->type_id() == k )
                            {
                                if ( tDof->is_fixed() )
                                {
                                    // dof -> field: skipped in seeding mode,
                                    // where the field is the truth
                                    if ( ! aFreeDofsOnly )
                                    {
                                        tValues( tDof->dof_index_on_field() ) = tDof->value();
                                    }
                                }
                                else
                                {
                                    tDof->value() = tValues( tDof->dof_index_on_field() );
                                }
                                tDof->flag();
                            }
                        }
                    }
                } // end loop over all fields

                // unflag all dofs
                for( Dof * tDof : mDOFs )
                {
                    tDof->unflag() ;
                }
            }

//------------------------------------------------------------------------------
            /**
             * Temporarily splits the unified DOF container (mDOFs) into separate
             * free and fixed DOF graphs for independent reordering.
             *
             * Purpose:
             *   - Enables independent application of graph reordering (symrcm)
             *     to free and fixed DOFs
             *   - Simplifies sorting by separating DOF types
             *
             * Actions:
             *   1. Voids index and my_index of every regular DOF (gNoIndex) so that
             *      a stale index faults instead of aliasing; reorder_dofs() rebuilds them.
             *      Hanging DOFs get my_index = 0..h-1
             *   2. Populates the output containers in mDOFs order (fixed / free split)
             *   3. Clears mDOFs (restored by restore_dof_container())
             *
             * Postconditions:
             *   - aFreeDofs contains all free DOFs (size = n)
             *   - aFixedDofs contains all fixed DOFs (size = m)
             *   - mDOFs is empty
             *   - index and my_index of every regular DOF are gNoIndex
             *
             * @param aFreeDofs  Output: graph container for free DOFs
             * @param aFixedDofs Output: graph container for fixed DOFs
             */
            void
            DofData::split_dof_container(
                Cell< graph::Vertex * > & aFreeDofs,
                Cell< graph::Vertex * > & aFixedDofs )
            {
                // Count and assign initial SEPARATE mode indices
                mMyNumberOfFreeDofs = 0 ;
                mMyNumberOfFixedDofs = 0 ;

                for ( Dof * tDof : mDOFs )
                {
                    // we void all indices to provoke an error
                    // if we mess up. reorder_dofs has to rebuild them
                    tDof->set_index( gNoIndex );
                    tDof->set_my_index( gNoIndex ) ;
                    if ( tDof->is_fixed() )
                    {
                        mMyNumberOfFixedDofs++;
                    }
                    else
                    {
                        mMyNumberOfFreeDofs++;
                    }

                }

                // Handle hanging DOFs separately
                mMyNumberOfHangingDofs = 0 ;
                for (Dof * tDof : mHangingDOFs )
                {
                    tDof->set_my_index( mMyNumberOfHangingDofs++ );
                }

                // Allocate graph containers
                aFreeDofs.set_size( mMyNumberOfFreeDofs, nullptr );
                aFixedDofs.set_size( mMyNumberOfFixedDofs, nullptr );

                // Populate graphs using my_index as position
                mMyNumberOfFreeDofs = 0 ;
                mMyNumberOfFixedDofs = 0 ;
                for ( Dof * tDof : mDOFs )
                {
                    if ( tDof->is_fixed() )
                    {
                        aFixedDofs( mMyNumberOfFixedDofs++ ) = tDof ;
                    }
                    else
                    {
                        aFreeDofs( mMyNumberOfFreeDofs++ ) = tDof ;
                    }
                }

                // Clear mDOFs to save memory during reordering phase
                mDOFs.clear() ;
            }

            /**
             * Restores the unified DOF container (mDOFs) from separate free and
             * fixed DOF graphs after reordering is complete.
             *
             * Purpose:
             *   - Rebuilds mDOFs for use during assembly phase
             *   - Maintains consistent storage layout: free DOFs first, then fixed DOFs
             *
             * Preconditions:
             *   - aFreeDofs and aFixedDofs must be sorted by their reordered index
             *   - Each DOF's my_index must be set to its final position (0..n-1 or 0..m-1)
             *
             * Actions:
             *   1. Allocates mDOFs with size n+m
             *   2. Copies free DOFs to positions 0..n-1 (in iteration order)
             *   3. Copies fixed DOFs to positions n..n+m-1 (in iteration order)
             *   4. Clears input graph containers
             *
             * Postconditions:
             *   - mDOFs contains all DOFs in sorted order
             *   - mDOFs[i] for i<n: free DOF with my_index=i
             *   - mDOFs[n+j] for j<m: fixed DOF with my_index=j
             *   - aFreeDofs and aFixedDofs are empty
             *
             * Storage layout after restore:
             *   mDOFs[0..n-1]     : free DOFs in reordered sequence
             *   mDOFs[n..n+m-1]   : fixed DOFs in reordered sequence
             *
             * @param aFreeDofs  Input: sorted graph of free DOFs (will be cleared)
             * @param aFixedDofs Input: sorted graph of fixed DOFs (will be cleared)
             */
            void
            DofData::restore_dof_container(
                Cell< graph::Vertex * > & aFreeDofs,
                Cell< graph::Vertex * > & aFixedDofs )
            {
                // Allocate unified container with total size
                mDOFs.set_size( mMyNumberOfFreeDofs + mMyNumberOfFixedDofs, nullptr );

                // Copy free DOFs to positions 0..n-1
                index_t tCount = 0 ;
                for ( graph::Vertex * tVertex : aFreeDofs )
                {
                    mDOFs( tCount++ ) = reinterpret_cast< Dof * >( tVertex );
                }

                // Copy fixed DOFs to positions n..n+m-1
                for ( graph::Vertex * tVertex : aFixedDofs )
                {
                    mDOFs( tCount++ ) = reinterpret_cast< Dof * >( tVertex );
                }

                // Clear temporary graph containers
                aFreeDofs.clear() ;
                aFixedDofs.clear() ;
#ifdef DEBUG
                for ( Dof * tDof : mDOFs )
                {
                    BELFEM_ASSERT( tDof->index() != gNoIndex, "DOF index not set after restore" );
                    BELFEM_ASSERT( tDof->my_index() != gNoIndex, "DOF my_index not set after restore" );
                }
#endif
            }

//------------------------------------------------------------------------------

            void
            DofData::synchronize_dirichlet_bcs()
            {

                uint tNumProcs = mKernel->number_of_procs() ;

                if( tNumProcs > 1 )
                {
                    if ( mCommRank == 0 )
                    {
                        Cell< Vector< id_t > > tDofTables( tNumProcs, {} );

                        // cell with dof indices of fixed dofs
                        Cell< Vector< id_t > > tAllIDs( tNumProcs, {} );

                        // cell with dof values
                        Cell< Vector< real > > tAllValues( tNumProcs, {} );
                        comm_barrier() ;

                        collect( tAllIDs );
                        collect( tAllValues );
                        collect( tDofTables );

                        for ( uint p = 1; p < tNumProcs; ++p )
                        {

                            Vector< id_t > & tIDs = tAllIDs( p );
                            Vector< real > & tValues = tAllValues( p );

                            index_t tNumDofs = tIDs.length();

                            for ( index_t k = 0; k < tNumDofs; ++k )
                            {
                                this->dof( tIDs( k ) )->fix( tValues( k ) );
                            }
                        }

                        // now, we make sure that the data is consistent
                        for ( uint p = 1; p < tNumProcs; ++p )
                        {
                            // get ID vector
                            Vector< id_t > & tDofIDs = tDofTables( p );

                            // count fixed dofs
                            index_t tCount = 0 ;
                            for( id_t tID : tDofIDs )
                            {
                                if( this->dof( tID )->is_fixed() )
                                {
                                    ++tCount ;
                                }
                            }

                            // allocate containers
                            Vector< id_t > & tIDs = tAllIDs( p );
                            Vector< real > & tValues = tAllValues( p );

                            tIDs.set_size( tCount );
                            tValues.set_size( tCount );

                            // reset counter
                            tCount = 0 ;

                            // get values and IDs of fixed dofs
                            for( id_t tID : tDofIDs )
                            {
                                // get dof
                                Dof * tDof = this->dof( tID );

                                if( tDof->is_fixed() )
                                {
                                    tIDs( tCount ) = tDof->id() ;
                                    tValues( tCount++ ) = tDof->value() ;
                                }
                            }


                        }
                        comm_barrier() ;

                        // send containers to other procs
                        distribute( tAllIDs ) ;
                        distribute( tAllValues );
                    }
                    else
                    {
                        // count fixed dofs
                        index_t tCount = 0 ;
                        Vector< id_t > tAllIDs( mDOFs.size() );
                        for ( Dof * tDOF : mDOFs )
                        {
                            tAllIDs( tCount++ ) = tDOF->id() ;
                        }

                        // count fixed dofs
                        tCount = 0;
                        for ( Dof * tDOF : mDOFs )
                        {
                            if ( tDOF->is_fixed() )
                            {
                                ++tCount;
                            }
                        }

                        // container for IDs
                        Vector< id_t > tIDs( tCount );

                        // container for values
                        Vector< real > tValues( tCount );

                        // reset counter
                        tCount = 0;

                        for ( Dof * tDOF : mDOFs )
                        {
                            if ( tDOF->is_fixed() )
                            {
                                tIDs( tCount )    = tDOF->id();
                                tValues( tCount ) = tDOF->value();
                                ++tCount;
                            }
                        }

                        // send fixed IDs and values
                        comm_barrier() ;
                        send( tIDs );
                        send( tValues );
                        send( tAllIDs );
                        comm_barrier() ;

                        // receive confirmation from master
                        receive( tIDs );
                        receive( tValues );

                        // reset counter
                        tCount = 0 ;

                        // loop over all IDs
                        for( id_t tID : tIDs )
                        {
                            // get dof
                            Dof * tDof = this->dof( tID );
                            tDof->fix( tValues( tCount++ ) );
                        }
                    }
                }
            }



//------------------------------------------------------------------------------

            /**
             * Reorders DOFs to optimize matrix bandwidth for efficient solving.
             *
             * Purpose:
             *   - Applies graph reordering algorithms to minimize matrix bandwidth
             *   - Reduces fill-in during factorization and improves cache locality
             *   - Handles free and fixed DOFs independently
             *
             * Algorithm:
             *   Free DOFs:
             *     - Symmetric Reverse Cuthill-McKee (symrcm), bandwidth reduction only.
             *       Solver-internal nested dissection (STRUMPACK/ParMETIS) is applied
             *       by the solver wrapper, not here.
             *   Fixed DOFs:
             *     - symrcm (simpler, adequate for constraint equations)
             *
             * Workflow - Master Rank (rank 0):
             *   1. Initialize indices in SEPARATE mode:
             *      - Free DOFs: index = my_index = 0, 1, ..., n-1
             *      - Fixed DOFs: index = my_index = 0, 1, ..., m-1
             *   2. Build Jacobian (free-free) sparsity pattern with aLinkToSelf=false
             *   3. Apply reordering algorithm to free DOFs (updates index and sorts graph)
             *   4. Build Imposition (fixed-fixed) sparsity pattern with aLinkToSelf=false
             *   5. Apply symrcm to fixed DOFs (updates index and sorts graph)
             *   6. For serial: Update my_index to match reordered positions and return
             *   7. For parallel:
             *      - Broadcast global DOF counts to all ranks
             *      - Collect DOF IDs from worker ranks
             *      - Look up reordered index for each worker's DOF
             *      - Distribute index tables back to workers
             *
             * Workflow - Worker Ranks (rank > 0):
             *   1. Collect all local DOF IDs (free first, then fixed)
             *   2. Receive global DOF counts from master
             *   3. Send DOF IDs to master
             *   4. Receive reordered indices from master
             *   5. Apply indices to graphs:
             *      - Free DOFs: tDofIndices[0..n_local-1]
             *      - Fixed DOFs: tDofIndices[n_local..n_local+m_local-1] (NO counter reset!)
             *   6. Sort graphs by reordered index
             *
             * Index Fields After Reordering:
             *   - index: Global matrix row/column index in SEPARATE mode (0..n-1 or 0..m-1)
             *   - my_index: Local position in graph (0..n_local-1 or 0..m_local-1)
             *
             * Preconditions:
             *   - aFreeDofs and aFixedDofs populated from split_dof_container()
             *   - aGraphData contains global DOF-to-DOF connectivity
             *
             * Postconditions:
             *   - Both graphs sorted by optimized index ordering
             *   - DOF index field contains final matrix indices (SEPARATE mode)
             *   - DOF my_index field contains local graph positions (SEPARATE mode)
             *   - mDofIndexTables populated for parallel assembly (master only)
             *
             * Note on aLinkToSelf:
             *   - Set to false for reordering: METIS/SCOTCH require no self-loops
             *   - Matrix allocation later uses aLinkToSelf=true to include diagonals
             */
            void
            DofData::reorder_dofs(
                const Vector< id_t > & aGraphData,
                               Graph & aFreeDofs,
                               Graph & aFixedDofs )
            {
                if ( mCommRank == 0 )
                {
                    mNumberOfFreeDofs = 0 ;
                    mNumberOfFixedDofs = 0 ;

                    for ( graph::Vertex * tVertex : aFreeDofs )
                    {
                        Dof * tDof = reinterpret_cast< Dof * >( tVertex );
                        tDof->set_index( mNumberOfFreeDofs );
                        tDof->set_my_index( mNumberOfFreeDofs++ );
                    }
                    for ( graph::Vertex * tVertex : aFixedDofs )
                    {
                        Dof * tDof = reinterpret_cast< Dof * >( tVertex );
                        tDof->set_index( mNumberOfFixedDofs );
                        tDof->set_my_index( mNumberOfFixedDofs++ );
                    }

                    mMyNumberOfFreeDofs = mNumberOfFreeDofs ;
                    mMyNumberOfFixedDofs = mNumberOfFixedDofs ;

                    mParent->solver_data()->populate_graph(
                        aGraphData,
                        Jacobian,
                        aFreeDofs,
                        aFixedDofs,
                        false );

                    // Reorder free DOFs using symrcm (bandwidth reduction).
                    // Solver-internal nested dissection (e.g. STRUMPACK ParMETIS)
                    // is handled by the solver wrapper, not here.
                    graph::symrcm( aFreeDofs );

                    // Reorder fixed DOFs (always use symrcm)
                    if ( mNumberOfFixedDofs > 0 )
                    {
                        mParent->solver_data()->populate_graph(
                        aGraphData,
                        Imposition,
                        aFreeDofs,
                        aFixedDofs,
                        false );

                        graph::symrcm( aFixedDofs );
                    }

                    // Serial execution: Set my_index and return (no parallel communication needed)
                    if ( mCommSize == 1 )
                    {
                        index_t tCount = 0 ;
                        for ( graph::Vertex * tVertex : aFreeDofs )
                        {
                            reinterpret_cast< Dof * >( tVertex )->set_my_index( tCount++ );
                        }
                        tCount = 0 ;
                        for ( graph::Vertex * tVertex : aFixedDofs )
                        {
                            reinterpret_cast< Dof * >( tVertex )->set_my_index( tCount++ );
                        }
                        return;
                    }

                    // Parallel execution: Distribute reordered indices to worker ranks
                    comm_barrier() ;
                    broadcast( mNumberOfFreeDofs );
                    broadcast( mNumberOfFixedDofs );

                    // Collect DOF IDs from all worker ranks
                    Cell< Vector< id_t > > tDofIdTables( mCommSize, {} );
                    collect( tDofIdTables );
                    mDofIndexTables.set_size( mCommSize, {} );

                    // For each worker rank, look up reordered indices for their DOFs
                    for ( proc_t p = 1; p < mCommSize; ++p )
                    {
                        Vector< id_t    > & tDofIDs = tDofIdTables( p );
                        Vector< index_t > & tDofIndices = mDofIndexTables( p );

                        tDofIndices.set_size( tDofIDs.length() );
                        index_t tCount = 0 ;
                        for ( id_t tID : tDofIDs )
                        {
                            Dof * tDof = this->dof( tID );
                            if ( tDof->is_fixed() )
                            {
                                tDofIndices( tCount++ ) = tDof->index() + mNumberOfFreeDofs ;
                            }
                            else
                            {
                                tDofIndices( tCount++ ) = tDof->index() ;
                            }
                        }
                    }

                    // Send index tables back to worker ranks
                    comm_barrier() ;
                    distribute( mDofIndexTables );
                    for ( proc_t p = 1; p < mCommSize; ++p )
                    {
                        sort( mDofIndexTables( p ) );
                    }

                }
                else // Worker ranks (rank > 0)
                {
                    // Solver-internal reordering is handled by the solver wrapper.
                    // Worker ranks no longer participate in DofData-level nested dissection.

                    // Collect local DOF IDs to send to master
                    Vector< id_t > tDofIDs( aFreeDofs.size() + aFixedDofs.size() );
                    index_t tCount = 0 ;
                    for ( graph::Vertex * tVertex : aFreeDofs )
                    {
                        tDofIDs( tCount++ ) = tVertex->id() ;
                    }
                    for ( graph::Vertex * tVertex : aFixedDofs )
                    {
                        tDofIDs( tCount++ ) = tVertex->id() ;
                    }

                    // Send DOF IDs to master, receive reordered indices
                    comm_barrier() ;
                    broadcast( mNumberOfFreeDofs );
                    broadcast( mNumberOfFixedDofs );
                    send( tDofIDs );
                    Vector< index_t > tDofIndices ;
                    comm_barrier() ;
                    receive( tDofIndices );

                    // Apply reordered indices to graphs
                    // CRITICAL: tCount is NOT reset between free and fixed DOFs!
                    // tDofIndices contains concatenated indices: [free_0...free_n, fixed_0...fixed_m]
                    tCount = 0 ;

                    for ( graph::Vertex * tVertex : aFreeDofs )
                    {
                        tVertex->set_index( tDofIndices( tCount++ ) );
                    }
                    BELFEM_ASSERT( tCount == aFreeDofs.size() , "Number of dofs does not match" );
                    // Continue from tCount (do NOT reset to 0)
                    for ( graph::Vertex * tVertex : aFixedDofs )
                    {
                        tVertex->set_index( tDofIndices( tCount++ ) - mNumberOfFreeDofs );
                    }
                    BELFEM_ASSERT( tCount == aFreeDofs.size() + aFixedDofs.size() , "Number of dofs does not match" );

                    // Sort graphs by reordered index
                    sort( aFreeDofs, opVertexIndex );
                    sort( aFixedDofs, opVertexIndex );

                }

                // Set my_index to local graph positions for all ranks (SEPARATE mode)
                index_t tCount = 0 ;
                for ( graph::Vertex * tVertex : aFreeDofs )
                {
                    reinterpret_cast< Dof * >( tVertex )->set_my_index( tCount++ );
                }
                tCount = 0 ;
                for ( graph::Vertex * tVertex : aFixedDofs )
                {
                    reinterpret_cast< Dof * >( tVertex )->set_my_index( tCount++ );
                }
            }

//------------------------------------------------------------------------------

            uint
            DofData::num_dofs_per_element( const id_t aBlockID ) const
            {
                // get block
                ElementType tType = mMesh->block( aBlockID )->element_type() ;

                uint aN = mesh::number_of_nodes( tType );

                aN *= mParent->iwg()->number_of_dofs_per_node( aBlockID );

                aN += mParent->iwg()->number_of_dofs_per_edge( aBlockID ) *
                        mesh::number_of_edges( tType );

                aN += mParent->iwg()->number_of_dofs_per_face( aBlockID ) *
                        mesh::number_of_faces( tType );

                return aN ;
            }

//------------------------------------------------------------------------------

            uint
            DofData::num_dofs_per_facet( const id_t aSideSetID ) const
            {
                return mParent->iwg()->number_of_dofs_per_element( mParent->sideset( aSideSetID ) );
            }

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

            void
            DofData::create_dofwise_t_matrices_master()
            {
                Cell< Dof * >     tSources;
                Vector< real >    tWeights;

                Vector< real >    tCoefficients;
                Vector< index_t > tIndices;

                Matrix< real >         tNodeWeights ;

                // edge dofs whose source edge dof was not resolved yet when
                // visited ( container order ); flattened after the main loop
                Cell< Dof * > tDeferredEdgeDofs ;

                this->connect_dofs_to_mesh() ;

                mMesh->unflag_all_edges();
                mMesh->unflag_all_faces();

                for ( Dof * tDof: mHangingDOFs )
                {
                    tSources.clear();

                    // for trivial node dependency
                    if ( tDof->entity_type() == EntityType::NODE )
                    {

                        // check for trivial dependency if both dofs ase nodes
                        if (    tDof->number_of_sources() == 1
                             && tDof->source( 0 )->entity_type() == EntityType::NODE
                             && tDof->source( 0 )->mesh_basis()->number_of_dofs() == 1 )
                        {
                            BELFEM_ASSERT( ! tDof->source( 0 )->mesh_basis()->is_hanging(), "Source of dof is hanging. This should not happen" );

                            tDof->set_source( reinterpret_cast< Dof * >(
                                tDof->source( 0 )->mesh_basis()->dof( 0 ) ), tDof->mesh_basis()->weight( 0 ) );

                            continue;
                        }


                        tSources.set_size( tDof->mesh_basis()->number_of_sources(), nullptr );
                        tWeights.set_size( tDof->mesh_basis()->number_of_sources(), 0.0 );

                        uint t = tDof->type_id();

                        bool tAllSourcesFound = true;

                        for ( uint k = 0; k < tDof->mesh_basis()->number_of_sources(); ++k )
                        {
                            bool tFound = false;

                            for ( uint j = 0; j < tDof->mesh_basis()->source( k )->number_of_dofs(); ++j )
                            {
                                // grab other dof
                                Dof * tSource = reinterpret_cast< Dof * > ( tDof->mesh_basis()->source( k )->dof(
                                        j ));

                                if ( tSource->type_id() == t )
                                {
                                    tSources( k ) = tSource;
                                    tWeights( k ) = tDof->mesh_basis()->weight( k );
                                    tFound = true;
                                    break;
                                }
                            }

                            if ( ! tFound )
                            {
                                // Couldn't find matching DoF type in this source
                                // This happens at domain boundaries where source is outside formulation
                                tAllSourcesFound = false;
                                break;
                            }
                        }

                        if ( tAllSourcesFound )
                        {
                            // Normal case: all sources have matching DoF types
                            tDof->set_sources( tSources, tWeights );
                        }
                        // else: Don't set sources - DoF will be treated as non-hanging
                        // is_hanging() will return false, and remove_hanging_dofs_from_container()
                        // will move it back to regular DoFs

                        continue ;
                    } // end if node

                    if ( tDof->entity_type() == EntityType::EDGE )
                    {
                        mesh::Edge * tEdge = tDof->edge();

                        if ( tEdge->is_flagged() ) continue;

                        mesh::Basis * tSource = tEdge->source( 0 );

                        if ( tSource->entity_type() == EntityType::NODE )
                        {
                            for ( uint i = 0; i < tEdge->number_of_nodes(); ++i )
                            {
                                mesh::Node * tNode = reinterpret_cast< mesh::Node * >( tEdge->source( i ) );

                                if ( tNode->is_hanging() )
                                {
                                    for ( uint j=0; j < tNode->number_of_sources(); ++j )
                                    {
                                        BELFEM_ASSERT(  tNode->source( j )->number_of_dofs() == 1 , "Internal error: couldn't figure out how do assign dependencies of dof %lu ( edge %lu )",
                                              ( long unsigned int ) tDof->id(),
                                              ( long unsigned int ) tDof->mesh_basis()->id() );

                                        tSources.push( reinterpret_cast< Dof * >( tNode->source( j )->dof( 0 ) ) );
                                    }
                                }

                                else
                                { BELFEM_ASSERT(  tNode->number_of_dofs() == 1 ,
                                        "Internal error: couldn't figure out how do assign dependencies of dof %lu ( edge %lu [%lu->%lu] ), linked to node %lu. (Node has %u dofs, element %lu)",
                                        ( long unsigned int ) tDof->id(),
                                        ( long unsigned int ) tEdge->id(),
                                        ( long unsigned int ) tEdge->node( 0 )->id(),
                                        ( long unsigned int ) tEdge->node( 1 )->id(),
                                        ( long unsigned int ) tNode->id() ,
                                    ( unsigned int ) tNode->number_of_dofs(),
                                    ( long unsigned int ) tNode->element(0)->id() );

                                    tSources.push( reinterpret_cast< Dof * >( tNode->dof( 0 ) ) );
                                }
                            }


                            unique( tSources );

                            uint n = tSources.size();

                            // backup dof indices
                            tIndices.set_size( n );
                            for ( uint k=0; k<tSources.size(); ++k )
                            {
                                tIndices( k ) = tSources( k )->index();
                                tSources( k )->set_index( k );
                            }

                            tNodeWeights.set_size( n, tEdge->number_of_nodes(), 0.0 );

                            for ( uint i=0; i<tEdge->number_of_nodes(); ++i )
                            {
                                mesh::Node * tNode = reinterpret_cast< mesh::Node * >( tEdge->source( i ) );

                                if ( tNode->is_hanging() )
                                {
                                    for ( uint k=0; k<tNode->number_of_sources(); ++k )
                                    {
                                        tNodeWeights( tNode->source( k )->dof( 0 )->index(), i ) = tNode->weight( k );
                                    }
                                }
                                else
                                {
                                    tNodeWeights( tNode->dof( 0 )->index(), i ) = 1.0 ;
                                }
                            }

                            if ( tEdge->number_of_nodes() == 2 )
                            {
                                tCoefficients.set_size( 2, 0.0 );
                                tCoefficients( 0 ) = 1.0 ;
                                tCoefficients( 1 ) = -1.0 ;

                                tWeights = tNodeWeights * tCoefficients;

                                tDof->set_sources( tSources, tWeights );
                            }
                            else if ( tEdge->number_of_nodes() == 3 )
                            {
                                tCoefficients.set_size( 3, 0.0 );

                                // first dof
                                tCoefficients( 0 ) =  1.0 ;
                                tCoefficients( 1 ) =  1./3. ;
                                tCoefficients( 2 ) = -4./3. ;

                                tWeights = tNodeWeights * tCoefficients;

                                reinterpret_cast< Dof * >( tEdge->dof( 0 ) )->set_sources( tSources, tWeights );

                                // second dof
                                tCoefficients( 0 ) =  1./3. ;
                                tCoefficients( 1 ) =  1.0 ;
                                tCoefficients( 2 ) =  -4./3. ;

                                tWeights = tNodeWeights * tCoefficients;

                                reinterpret_cast< Dof * >( tEdge->dof( 1 ) )->set_sources( tSources, tWeights );
                            }
                            else
                            {
                                BELFEM_ERROR( false, "Invalid edge");
                            }

                            // restore dof indices
                            for ( uint k=0; k<n; ++k )
                            {
                                tSources( k )->set_index( tIndices( k ) );
                            }

                        }
                        else if ( tSource->entity_type() == EntityType::EDGE )
                        {
                            if ( tEdge->number_of_sources() == 1 ) // 1:1 edge coupling
                            {
                                mesh::Edge * tOther = reinterpret_cast< mesh::Edge * >( tSource );
                                real tWeight = tEdge->weight( 0 );

                                int n = tEdge->number_of_dofs();
                                int m = tOther->number_of_dofs();


                                switch ( n )
                                {
                                    case 1 :
                                    {
                                        switch ( m )
                                        {
                                            case 1 :
                                            {
                                                // LINE2 hanging on LINE2: 1:1 coupling
                                                Dof * tOtherDof = reinterpret_cast< Dof * >( tOther->dof( 0 ) );

                                                if ( tOtherDof->is_hanging() )
                                                {
                                                    // cascade: source edge is itself hanging
                                                    uint nSrc = tOtherDof->number_of_sources();
                                                    tSources.set_size( nSrc, nullptr );
                                                    tWeights.set_size( nSrc );
                                                    for ( uint k = 0; k < nSrc; ++k )
                                                    {
                                                        tSources( k ) = tOtherDof->source( k );
                                                        tWeights( k ) = tWeight * tOtherDof->weight( k );
                                                    }
                                                    reinterpret_cast< Dof * >( tEdge->dof( 0 ) )->set_sources( tSources, tWeights );

                                                }
                                                else if ( tOther->is_hanging() )
                                                {
                                                    // the source edge hangs at mesh level, but its dof
                                                    // was not resolved yet ( loop order ); flatten later
                                                    tDeferredEdgeDofs.push( reinterpret_cast< Dof * >( tEdge->dof( 0 ) ) );
                                                }
                                                else
                                                {
                                                    reinterpret_cast< Dof * >( tEdge->dof( 0 ) )->set_source( tOtherDof, tWeight );

                                                }
                                                break ;
                                            }
                                            case 2 :
                                            {
                                                // LINE2 hanging on LINE3: average the two source DOFs
                                                Dof * tOtherDof0 = reinterpret_cast< Dof * >( tOther->dof( 0 ) );
                                                Dof * tOtherDof1 = reinterpret_cast< Dof * >( tOther->dof( 1 ) );

                                                if ( ! tOtherDof0->is_hanging() && ! tOtherDof1->is_hanging() )
                                                {
                                                    // no cascade
                                                    tSources.set_size( 2, nullptr );
                                                    tWeights.set_size( 2 );
                                                    if ( tWeight == 1.0 )
                                                    {
                                                        tSources( 0 ) = tOtherDof0;
                                                        tSources( 1 ) = tOtherDof1;
                                                        tWeights( 0 ) = 0.5 ;
                                                        tWeights( 1 ) = 0.5 ;
                                                    }
                                                    else if ( tWeight == -1.0 )
                                                    {
                                                        tSources( 0 ) = tOtherDof1;
                                                        tSources( 1 ) = tOtherDof0;
                                                        tWeights( 0 ) = -0.5 ;
                                                        tWeights( 1 ) = -0.5 ;
                                                    }
                                                    else
                                                    {
                                                        BELFEM_ERROR( false, "invalid edge weight: %g", tWeight );
                                                    }
                                                    reinterpret_cast< Dof * >( tEdge->dof( 0 ) )->set_sources( tSources, tWeights );
                                                }
                                                else
                                                {
                                                    // cascade: expand hanging source DOFs and merge
                                                    BELFEM_ERROR( tWeight == 1.0 || tWeight == -1.0,
                                                        "invalid edge weight: %g", tWeight );

                                                    tSources.clear();
                                                    if ( tOtherDof0->is_hanging() )
                                                    {
                                                        for ( uint k = 0; k < tOtherDof0->number_of_sources(); ++k )
                                                            tSources.push( tOtherDof0->source( k ) );
                                                    }
                                                    else
                                                    {
                                                        tSources.push( tOtherDof0 );
                                                    }
                                                    if ( tOtherDof1->is_hanging() )
                                                    {
                                                        for ( uint k = 0; k < tOtherDof1->number_of_sources(); ++k )
                                                            tSources.push( tOtherDof1->source( k ) );
                                                    }
                                                    else
                                                    {
                                                        tSources.push( tOtherDof1 );
                                                    }
                                                    unique( tSources );
                                                    uint nSrc = tSources.size();

                                                    // backup and set temporary indices
                                                    tIndices.set_size( nSrc );
                                                    for ( uint k = 0; k < nSrc; ++k )
                                                    {
                                                        tIndices( k ) = tSources( k )->index();
                                                        tSources( k )->set_index( k );
                                                    }

                                                    tWeights.set_size( nSrc, 0.0 );

                                                    Dof * tFirst  = ( tWeight == 1.0 ) ? tOtherDof0 : tOtherDof1;
                                                    Dof * tSecond = ( tWeight == 1.0 ) ? tOtherDof1 : tOtherDof0;
                                                    real  tSign   = ( tWeight == 1.0 ) ? 0.5 : -0.5;

                                                    // accumulate from first source DOF
                                                    if ( tFirst->is_hanging() )
                                                    {
                                                        for ( uint k = 0; k < tFirst->number_of_sources(); ++k )
                                                            tWeights( tFirst->source( k )->index() ) += tSign * tFirst->weight( k );
                                                    }
                                                    else
                                                    {
                                                        tWeights( tFirst->index() ) += tSign;
                                                    }

                                                    // accumulate from second source DOF
                                                    if ( tSecond->is_hanging() )
                                                    {
                                                        for ( uint k = 0; k < tSecond->number_of_sources(); ++k )
                                                            tWeights( tSecond->source( k )->index() ) += tSign * tSecond->weight( k );
                                                    }
                                                    else
                                                    {
                                                        tWeights( tSecond->index() ) += tSign;
                                                    }

                                                    reinterpret_cast< Dof * >( tEdge->dof( 0 ) )->set_sources( tSources, tWeights );

                                                    // restore indices
                                                    for ( uint k = 0; k < nSrc; ++k )
                                                    {
                                                        tSources( k )->set_index( tIndices( k ) );
                                                    }
                                                }
                                                break ;
                                            }
                                            default :
                                            {
                                                BELFEM_ERROR( false, "invalid number of dofs on source edge : %u", ( unsigned int ) m );
                                            }

                                        }
                                        break ;
                                    }
                                    case 2 :
                                    {
                                        switch ( m )
                                        {
                                            case 1 :
                                            {
                                                // LINE3 hanging on LINE2: both DOFs map to same source
                                                Dof * tOtherDof = reinterpret_cast< Dof * >( tOther->dof( 0 ) );

                                                if ( tOtherDof->is_hanging() )
                                                {
                                                    // cascade: source edge is itself hanging
                                                    uint nSrc = tOtherDof->number_of_sources();
                                                    tSources.set_size( nSrc, nullptr );
                                                    tWeights.set_size( nSrc );
                                                    for ( uint k = 0; k < nSrc; ++k )
                                                    {
                                                        tSources( k ) = tOtherDof->source( k );
                                                        tWeights( k ) = tWeight * tOtherDof->weight( k );
                                                    }
                                                    reinterpret_cast< Dof * >( tEdge->dof( 0 ) )->set_sources( tSources, tWeights );
                                                    reinterpret_cast< Dof * >( tEdge->dof( 1 ) )->set_sources( tSources, tWeights );
                                                }
                                                else
                                                {
                                                    reinterpret_cast< Dof * >( tEdge->dof( 0 ) )->set_source( tOtherDof, tWeight );
                                                    reinterpret_cast< Dof * >( tEdge->dof( 1 ) )->set_source( tOtherDof, tWeight );
                                                }
                                                break;
                                            }
                                            case 2 :
                                            {
                                                // LINE3 hanging on LINE3: DOF-to-DOF mapping with orientation
                                                uint d0 = ( tWeight == 1.0 ) ? 0 : 1;
                                                uint d1 = ( tWeight == 1.0 ) ? 1 : 0;

                                                BELFEM_ERROR( tWeight == 1.0 || tWeight == -1.0,
                                                    "invalid edge weight: %g", tWeight );

                                                for ( uint d = 0; d < 2; ++d )
                                                {
                                                    uint ds = ( d == 0 ) ? d0 : d1;
                                                    Dof * tOtherDof = reinterpret_cast< Dof * >( tOther->dof( ds ) );

                                                    if ( tOtherDof->is_hanging() )
                                                    {
                                                        // cascade
                                                        uint nSrc = tOtherDof->number_of_sources();
                                                        tSources.set_size( nSrc, nullptr );
                                                        tWeights.set_size( nSrc );
                                                        for ( uint k = 0; k < nSrc; ++k )
                                                        {
                                                            tSources( k ) = tOtherDof->source( k );
                                                            tWeights( k ) = tWeight * tOtherDof->weight( k );
                                                        }
                                                        reinterpret_cast< Dof * >( tEdge->dof( d ) )->set_sources( tSources, tWeights );
                                                    }
                                                    else
                                                    {
                                                        reinterpret_cast< Dof * >( tEdge->dof( d ) )->set_source( tOtherDof, tWeight );
                                                    }
                                                }
                                                break;
                                            }
                                            default :
                                            {
                                                BELFEM_ERROR( false, "invalid number of dofs on source edge : %u", ( unsigned int ) m );
                                            }
                                        }
                                        break ;
                                    }
                                    default:
                                    {
                                        BELFEM_ERROR( false, "invalid number of dofs on target edge : %u", ( unsigned int ) n );
                                    }
                                }
                            }
                            else
                            {
                                BELFEM_ASSERT( tEdge->number_of_dofs()== 1, "Only linear edges are supported for multiple edge sources" );

                                Map< Dof * , real > tMap ;

                                for ( uint s=0; s<tEdge->number_of_sources(); ++s )
                                {
                                    real tWeight = tEdge->weight( s );
                                    mesh::Edge * tOtherEdge = reinterpret_cast< mesh::Edge * >( tEdge->source( s ) );

                                    BELFEM_ASSERT( tOtherEdge->number_of_dofs()== 1, "Only linear edges are supported for multiple edge sources" );
                                    Dof * tOtherDof = reinterpret_cast< Dof * >( tOtherEdge->dof( 0 ) );

                                    if ( tOtherDof->is_hanging() )
                                    {
                                        // cascade: source edge is itself hanging — expand into its sources
                                        uint nSrc = tOtherDof->number_of_sources();
                                        for ( uint k = 0; k < nSrc; ++k )
                                        {
                                            Dof * tOtherSource =  tOtherDof->source( k );
                                            real tCascadeWeight = tWeight * tOtherDof->weight( k );

                                            auto tIterator = tMap.find( tOtherSource );
                                            if ( tIterator == tMap.end() )
                                            {
                                                tMap[ tOtherSource ] = tCascadeWeight;
                                            }
                                            else
                                            {
                                                // accumulate
                                                tMap[ tOtherSource ] = tIterator->second + tCascadeWeight;
                                            }
                                        }
                                    }
                                    else
                                    {
                                        auto tIterator = tMap.find( tOtherDof );
                                        if ( tIterator == tMap.end() )
                                        {
                                            tMap[ tOtherDof ] = tWeight ;
                                        }
                                        else
                                        {
                                            // overwriting map entry with accumulated value
                                            tMap[ tOtherDof ]= tIterator->second + tWeight ;
                                        }
                                    }
                                }

                                uint nSrc = tMap.size();
                                tSources.set_size( nSrc, nullptr );
                                tWeights.set_size( nSrc );
                                uint tCount = 0 ;
                                for ( auto tIterator = tMap.begin(); tIterator != tMap.end(); ++tIterator )
                                {
                                    tSources( tCount ) = tIterator->first;
                                    tWeights( tCount++ ) = tIterator->second;
                                }
                                tDof->set_sources( tSources, tWeights );
                            }
                        }
                        tEdge->flag();
                        continue ;
                    } // end if edge

                    if ( tDof->entity_type() == EntityType::FACE )
                    {
                        mesh::Face * tFace = tDof->face();

                        if ( tFace->is_flagged() ) continue ;

                        for ( uint i = 0; i < tFace->number_of_nodes(); ++i )
                        {
                            mesh::Node * tNode = tFace->node( i );

                            if ( tNode->is_hanging() )
                            {
                                for ( uint j=0; j < tNode->number_of_sources(); ++j )
                                {
                                    BELFEM_ASSERT(  tNode->source( j )->number_of_dofs() == 1 , "Internal error: couldn't figure out how do assign dependencies of dof %lu ( edge %lu )",
                                          ( long unsigned int ) tDof->id(),
                                          ( long unsigned int ) tDof->mesh_basis()->id() );

                                    tSources.push( reinterpret_cast< Dof * >( tNode->source( j )->dof( 0 ) ) );
                                }
                            }
                            else
                            {
                                BELFEM_ASSERT(  tNode->number_of_dofs() == 1 , "Internal error: couldn't figure out how do assign dependencies of dof %lu ( edge %lu )",
                                    ( long unsigned int ) tDof->id(),
                                    ( long unsigned int ) tDof->mesh_basis()->id() );

                                tSources.push( reinterpret_cast< Dof * >( tFace->node( i )->dof( 0 ) ) );
                            }
                        }

                        unique( tSources );

                        uint n = tSources.size();

                        // backup dof indices
                        tIndices.set_size( n );
                        for ( uint k=0; k<tSources.size(); ++k )
                        {
                            tIndices( k ) = tSources( k )->index();
                            tSources( k )->set_index( k );
                        }

                        uint m = tFace->number_of_nodes();

                        tNodeWeights.set_size( n, m, 0.0 );

                        for ( uint i=0; i<m; ++i )
                        {
                            mesh::Node * tNode = tFace->node( i );

                            if ( tNode->is_hanging() )
                            {
                                for ( uint k=0; k<tNode->number_of_sources(); ++k )
                                {
                                    tNodeWeights( tNode->source( k )->dof( 0 )->index(), i ) = tNode->weight( k );
                                }
                            }
                            else
                            {
                                tNodeWeights( tNode->dof( 0 )->index(), i ) = 1.0 ;
                            }
                        }

                        tCoefficients.set_size( m );
                        uint tCount = 0 ;

                        for ( uint l=0; l<2; ++l )
                        {
                            for ( uint i=0; i<m; ++i )
                            {
                                tCoefficients( i ) = tFace->weight( tCount++ );
                            }

                            tWeights = tNodeWeights * tCoefficients;

                            reinterpret_cast< Dof * >( tFace->dof( l ) )->set_sources( tSources, tWeights );
                        }

                        // restore dof indices
                        for ( uint k=0; k<tSources.size(); ++k )
                        {
                            tSources( k )->set_index( tIndices( k ) );
                        }

                        tFace->flag();

                        continue ;
                    }



                    BELFEM_ERROR( false,"Don't know how to create T-Matrix");
                } // end loop over dofs

                // flatten deferred edge-on-edge chains, e.g. a periodic slave
                // rim edge hanging on an interface-hanging master edge
                for ( Dof * tDof : tDeferredEdgeDofs )
                {
                    mesh::Edge * tEdge  = tDof->edge();
                    mesh::Edge * tOther = reinterpret_cast< mesh::Edge * >( tEdge->source( 0 ) );
                    Dof * tOtherDof = reinterpret_cast< Dof * >( tOther->dof( 0 ) );
                    real tWeight = tEdge->weight( 0 );

                    BELFEM_ERROR( tOtherDof->is_hanging(),
                        "could not flatten hanging chain of edge %lu onto edge %lu",
                        ( long unsigned int ) tEdge->id(),
                        ( long unsigned int ) tOther->id() );

                    uint nSrc = tOtherDof->number_of_sources();
                    tSources.set_size( nSrc, nullptr );
                    tWeights.set_size( nSrc );
                    for ( uint k = 0; k < nSrc; ++k )
                    {
                        BELFEM_ASSERT( ! tOtherDof->source( k )->is_hanging(),
                            "hanging chain of edge %lu is deeper than one level",
                            ( long unsigned int ) tEdge->id() );

                        tSources( k ) = tOtherDof->source( k );
                        tWeights( k ) = tWeight * tOtherDof->weight( k );
                    }
                    tDof->set_sources( tSources, tWeights );

                }
            }

//------------------------------------------------------------------------------

            void
            DofData::remove_hanging_dofs_from_container()
            {
                index_t tCount = 0 ;
                for ( Dof * tDof : mDOFs )
                {
                    if( ! tDof->is_hanging() )
                    {
                        ++tCount ;
                    }
                }

                Cell< Dof * > tDofs ;
                tDofs.vector_data() = std::move( mDOFs.vector_data() );
                mDOFs.set_size( tCount, nullptr );
                tCount = 0 ;
                for ( Dof * tDof : tDofs )
                {
                    if( ! tDof->is_hanging() )
                    {
                        mDOFs( tCount++ ) = tDof ;
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            DofData::extract_abstract_dofs_from_mesh()
            {
                // save abstract dofs before we disconnect the mesh
                if ( mCommRank == 0 )
                {
                    Cell< mesh::Node * > & tAbstractNodes = mParent->iwg()->abstract_nodes() ;

                    // count dofs
                    index_t tCount = 0 ;
                    for ( mesh::Node * tNode : tAbstractNodes )
                    {
                        tCount += tNode->number_of_dofs();
                    }

                    // get the container from the dof data
                    Cell< Dof * > & tDofs = this->abstract_dofs() ;

                    tDofs.set_size( tCount, nullptr );
                    tCount = 0 ;
                    for ( mesh::Node * tNode : tAbstractNodes )
                    {
                        for ( uint d=0; d<tNode->number_of_dofs(); ++d )
                        {
                            tDofs( tCount++ ) = reinterpret_cast< Dof * >( tNode->dof( d ) );
                        }
                    }
                }
            }

//------------------------------------------------------------------------------

        } /* end namespace dofmgr */
    } /* end namespace fem */
} /* end namespace belfem */
