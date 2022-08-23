//
// Created by christian on 7/9/21.
//

#include "cl_FEM_DofMgr_DofData.hpp"
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
                mMyRank( comm_rank() )
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
                mDofIDs.clear() ;
                mDofIndices.clear() ;

                // delete the dofs
                for ( Dof * tDof : mDOFs )
                {
                    delete tDof ;
                }
                mDOFs.clear() ;

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

                if( mMyRank == mKernel->master() )
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

                    index_t tNumberOfProcs = mKernel->number_of_procs() ;

                    Cell< Vector< id_t > > tAllDofIDs( tNumberOfProcs, Vector< id_t >() );
                    Cell< Vector< id_t > > tAllEntityIDs( tNumberOfProcs, Vector< id_t >() );
                    Cell< Vector< index_t > > tAllDofTypes( tNumberOfProcs, Vector< index_t >() );

                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // Send Node-Based DOFs
                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    for( index_t p=0; p<tNumberOfProcs; ++p )
                    {
                        // get the proc ID
                        proc_t tProc = mKernel->comm_table( p );

                        // only do something if this is not the master proc
                        if ( tProc != mMyRank )
                        {
                            // collect dof entities
                            this->count_dofs_for_proc( tProc,
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
                    send( mKernel->comm_table(), tAllDofIDs );
                    send( mKernel->comm_table(), tAllEntityIDs );
                    send( mKernel->comm_table(), tAllDofTypes );

                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // Send Edge-Based DOFs
                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    for( index_t p=0; p<tNumberOfProcs; ++p )
                    {
                        // get the proc ID
                        proc_t tProc = mKernel->comm_table( p );

                        // only do something if this is not the master proc
                        if ( tProc != mMyRank )
                        {

                            // collect dof entities
                            this->count_dofs_for_proc( tProc,
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
                    send( mKernel->comm_table(), tAllDofIDs );
                    send( mKernel->comm_table(), tAllEntityIDs );
                    send( mKernel->comm_table(), tAllDofTypes );

                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // Send Face-Based DOFs
                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    for( index_t p=0; p<tNumberOfProcs; ++p )
                    {
                        // get the proc ID
                        proc_t tProc = mKernel->comm_table( p );

                        // only do something if this is not the master proc
                        if ( tProc != mMyRank )
                        {

                            // collect dof entities
                            this->count_dofs_for_proc( tProc,
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
                    send( mKernel->comm_table(), tAllDofIDs );
                    send( mKernel->comm_table(), tAllEntityIDs );
                    send( mKernel->comm_table(), tAllDofTypes );

                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // Send Cell-Based DOFs
                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    for( index_t p=0; p<tNumberOfProcs; ++p )
                    {
                        // get the proc ID
                        proc_t tProc = mKernel->comm_table( p );

                        // only do something if this is not the master proc
                        if ( tProc != mMyRank )
                        {

                            // collect dof entities
                            this->count_dofs_for_proc( tProc,
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
                    send( mKernel->comm_table(), tAllDofIDs );
                    send( mKernel->comm_table(), tAllEntityIDs );
                    send( mKernel->comm_table(), tAllDofTypes );

                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // Send LambdaDOFs
                    // - - - - - - - - - - - - - - - - - - - - - - - - - -
                    for( index_t p=0; p<tNumberOfProcs; ++p )
                    {
                        // get the proc ID
                        proc_t tProc = mKernel->comm_table( p );

                        // only do something if this is not the master proc
                        if ( tProc != mMyRank )
                        {
                            // collect dof entities
                            this->count_dofs_for_proc( tProc,
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
                    send( mKernel->comm_table(), tAllDofIDs );
                    send( mKernel->comm_table(), tAllEntityIDs );
                    send( mKernel->comm_table(), tAllDofTypes );


                }
                else
                {
                    // wait for other procs
                    comm_barrier();

                    // receive node info from master
                    receive( mKernel->master(), tNodeDofIDs );
                    receive( mKernel->master(), tNodeDofEntityIDs );
                    receive( mKernel->master(), tNodeDofTypes );

                    // wait for other procs
                    comm_barrier();

                    // receive edge info from master
                    receive( mKernel->master(), tEdgeDofIDs );
                    receive( mKernel->master(), tEdgeDofEntityIDs );
                    receive( mKernel->master(), tEdgeDofTypes );

                    // wait for other procs
                    comm_barrier();

                    // receive face info from master
                    receive( mKernel->master(), tFaceDofIDs );
                    receive( mKernel->master(), tFaceDofEntityIDs );
                    receive( mKernel->master(), tFaceDofTypes );

                    // wait for other procs
                    comm_barrier();

                    // receive cell info from master
                    receive( mKernel->master(), tCellDofIDs );
                    receive( mKernel->master(), tCellDofEntityIDs );
                    receive( mKernel->master(), tCellDofTypes );

                    // wait for other procs
                    comm_barrier();

                    receive( mKernel->master(), tLambdaDofIDs );
                    receive( mKernel->master(), tLambdaDofEntityIDs );
                    receive( mKernel->master(), tLambdaDofTypes );
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
                    index_t tIndexOnField =  tEdge->index() * tMult + i++;

                    // create a new dof
                    Dof * tDof = new Dof( tEdgeDofIDs( k ),
                                          tEdgeDofTypes( k ),
                                          tEdge ,
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
                    index_t tIndexOnField =  tFace->index() * tMult + i++;

                    // create a new dof
                    Dof * tDof = new Dof( tFaceDofIDs( k ),
                                          tFaceDofTypes( k ),
                                          tFace,
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
                    index_t tIndexOnField =  tCell->index() * tMult + i++;

                    // create a new dof
                    Dof * tDof = new Dof( tCellDofIDs( k ),
                                          tCellDofTypes( k ),
                                          tCell,
                                          tIndexOnField );

                    // set the dof index. This will be overwritten during Jacobi initialization
                    tDof->set_index( tCount );

                    // add the dof to the container
                    mDOFs( tCount++ ) = tDof ;
                }


                // create the lambda DOFs
                for( index_t k=0; k<tNumLambdaDofs; ++k )
                {
                    // create a new dof
                    Dof * tDof = new Dof( tLambdaDofIDs( k ),
                                          tLambdaDofTypes( k ),
                                          mMesh->facet( tLambdaDofEntityIDs( k ) ) );

                    // set the dof index. This will be overwritten during Jacobi initialization
                    tDof->set_index( tCount );

                    // add the dof to the container
                    mDOFs( tCount++ ) = tDof ;
                }

                // create the dof map
                this->create_dof_map();

                // tables tell the master which dofs exist on which proc
                // and in what order they are stored
                this->create_dof_table() ;

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

//---------------------------------------------------------------------------

            void
            DofData::compute_dof_offsets( IWG  * aIWG )
            {
                // wait
                comm_barrier() ;

                BELFEM_ASSERT( mParent->iwg()->is_initialized(),
                     "Equation object has not been initialized");



                if( comm_rank() == mKernel->master() )
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
                    send_same( mKernel->comm_table(), tOffsets );
                }
                else
                {
                    // receive data from master
                    Vector< index_t > tOffsets( 5 ) ;
                    receive( mKernel->master(), tOffsets );

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

//------------------------------------------------------------------


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
                            if linear interpolation is enforced, only corner nodes are flagged
    
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


                if ( mParent->enforce_linear_interpolation() )
                {
                    for( id_t tID : aIWG->selected_blocks() )
                    {

                        if ( aIWG->number_of_dofs_per_node( tID ) > 0 )
                        {
                            tHaveEntityDofs = true;
                            mMesh->block( tID )->flag_corner_nodes();
                        }
                    }

                    // also flag ghost sidesets
                    for( id_t tID : aIWG->ghost_sideset_ids() )
                    {
                        // grab sideset on mesh
                        if ( aIWG->dofs_per_node_on_sideset( tID ).length() > 0 )
                        {
                            tHaveEntityDofs = true;
                            mMesh->sideset( tID )->flag_corner_nodes();
                        }
                    }
                }
                else
                {
                    for ( id_t tID : aIWG->selected_blocks() )
                    {
                        if ( aIWG->number_of_dofs_per_node( tID ) > 0 )
                        {
                            tHaveEntityDofs = true;
                            mMesh->block( tID )->flag_nodes();
                        }
                    }

                    // also flag ghost sidesets (if they exist)
                    for( id_t tID : aIWG->ghost_sideset_ids() )
                    {
                        // grab sideset on mesh
                        if ( aIWG->dofs_per_node_on_sideset( tID ).length() > 0 )
                        {
                            tHaveEntityDofs = true;
                            mMesh->sideset( tID )->flag_all_nodes() ;
                        }
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
                        // get the number of nodes per element on this block.
                        // We must check if linearity is enforced
                        uint tNumEntities = mParent->enforce_linear_interpolation() ?
                                            mesh::number_of_corner_nodes( tBlock->element_type() ) :
                                            mesh::number_of_nodes( tBlock->element_type() );

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

                // also consider ghost dofs
                for( id_t tSideSetID: aIWG->ghost_sidesets() )
                {
                    // grab block pointer on mesh
                    mesh::SideSet * tSideSet = mMesh->sideset( tSideSetID );

                    // ask IWG about selected dofs
                    const Vector< index_t > & tSelectedDofs = aIWG->dofs_per_node_on_sideset( tSideSetID );

                    // check if any dofs are selected
                    if( tSelectedDofs.length() > 0 )
                    {
                        // get the number of edges per element on this sideset
                        uint tNumEntities = mParent->enforce_linear_interpolation() ?
                                 mesh::number_of_corner_nodes( tSideSet->element_type() )
                                : mesh::number_of_nodes( tSideSet->element_type() );

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
                        // get the number of nodes per element on this block.
                        // We must check if linearity is enforced
                        uint tNumEntities = mParams->enforce_linear() ?
                                            mesh::number_of_corner_nodes( tBlock->element_type() ) :
                                            mesh::number_of_nodes( tBlock->element_type() );


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

                /**
                 * for edge dofs such as h-phi, the ghost sidesets don't see nodes
                 * for others such as temperature, they do !
                 */
                if( ! aIWG->has_edge_dofs() )
                {
                    for( id_t tSideSetID: aIWG->ghost_sidesets() )
                    {
                        // grab block on mesh
                        mesh::SideSet * tSideSet = mMesh->sideset( tSideSetID );

                        // grab the DOFs for this entity
                        const Vector< index_t > & tSelectedDofs = aIWG->dofs_per_node( tSideSetID );

                        if ( tSelectedDofs.length() > 0 )
                        {
                            // get the number of nodes per element on this block.
                            // We must check if linearity is enforced
                            uint tNumEntities = mParams->enforce_linear() ?
                                                mesh::number_of_corner_nodes( tSideSet->element_type() ) :
                                                mesh::number_of_nodes( tSideSet->element_type() );


                            // loop over all elements on this sideset
                            for ( mesh::Facet * tFacet: tSideSet->facets() )
                            {
                                mesh::Element * tElement = tFacet->element() ;

                                // loop over all entities on this element
                                for ( uint k = 0; k < tNumEntities; ++k )
                                {
                                    // loop over all selected dofs
                                    for ( index_t d: tSelectedDofs )
                                    {
                                        // compute the unique dof ID
                                        luint tID = tElement->node( k )->id() * BELFEM_MAX_DOFTYPES + d;

                                        // grab the corresponding bitset
                                        Bitset< BELFEM_MAX_NUMPROCS > & tBitset = aProcFlags( tDofMap( tID ));

                                        // set the proc bitset
                                        tBitset.set( tElement->owner());
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

                // also flag ghost sidesets (if they exist)
                for( id_t tID : aIWG->ghost_sideset_ids() )
                {

                    if ( aIWG->dofs_per_edge_on_sideset( tID ).length() > 0 )
                    {
                        tHaveEntityDofs = true;
                        mMesh->sideset( tID )->flag_edges() ;
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

                // also consider ghost dofs
                for( id_t tSideSetID: aIWG->ghost_sidesets() )
                {

                    // grab block pointer on mesh
                    mesh::SideSet * tSideSet = mMesh->sideset( tSideSetID );

                    // ask IWG about selected dofs
                    const Vector< index_t > & tSelectedDofs = aIWG->dofs_per_edge_on_sideset( tSideSetID );

                    // check if any dofs are selected
                    if( tSelectedDofs.length() > 0 )
                    {
                        // get the number of edges per element on this sideset
                        uint tNumEntities = mesh::number_of_edges( tSideSet->element_type() );

                        // now we loop over all elements on the block and the selected number oedges
                        for( mesh::Facet * tFacet : tSideSet->facets() )
                        {
                            // loop over all mesh entities on this element
                            for( uint k=0; k<tNumEntities; ++k )
                            {
                                // get the temporary index of this entity
                                index_t tIndex = tEntityMap( tFacet->edge( k )->id() );

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

                // also consider ghost dofs
                for( id_t tSideSetID: aIWG->ghost_sidesets() )
                {
                    // grab sideset on mesh
                    mesh::SideSet * tSideSet = mMesh->sideset( tSideSetID );

                    // grab the DOFs for this entity
                    const Vector< index_t > & tSelectedDofs = aIWG->dofs_per_edge_on_sideset( tSideSetID );

                    if ( tSelectedDofs.length() > 0 )
                    {
                        // get the number of edges per element
                        uint tNumEntities = mesh::number_of_edges( tSideSet->element_type() );

                        // loop over all elements on this block
                        for( mesh::Facet * tFacet : tSideSet->facets() )
                        {
                            mesh::Element * tElement = tFacet->element() ;

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

                // also flag ghost sidesets (if they exist)
                for( id_t tID : aIWG->ghost_sideset_ids() )
                {
                    if ( aIWG->dofs_per_face_on_sideset( tID ).length() > 0 )
                    {
                        tHaveEntityDofs = true;
                        mMesh->sideset( tID )->flag_faces() ;
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

                // also consider ghost dofs
                if( mMesh->number_of_dimensions() == 3 )
                {
                    for ( id_t tSideSetID: aIWG->ghost_sidesets() )
                    {
                        // grab block pointer on mesh
                        mesh::SideSet * tSideSet = mMesh->sideset( tSideSetID );

                        // ask IWG about selected dofs
                        const Vector< index_t > & tSelectedDofs = aIWG->dofs_per_face_on_sideset( tSideSetID );

                        // check if any dofs are selected
                        if ( tSelectedDofs.length() > 0 )
                        {
                            // now we loop over all elements on the block and the selected number oedges
                            for ( mesh::Facet * tFacet: tSideSet->facets() )
                            {
                                // get the temporary index of this entity
                                index_t tIndex = tEntityMap( tFacet->face()->id() );

                                // grab the corresponding bitset
                                Bitset< BELFEM_MAX_DOFTYPES > & tBitset = tDofFlags( tIndex );

                                // flip the bits that represent the selected dofs
                                for ( index_t s: tSelectedDofs )
                                {
                                    for ( index_t i = 0; i < tMultiplicity; ++i )
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

                // also consider ghost sidesets
                for( id_t tSideSetID : aIWG->ghost_sideset_ids() )
                {
                    // grab block on mesh
                    mesh::SideSet * tSideSet = mMesh->sideset( tSideSetID );

                    // grab the DOFs for this entity
                    const Vector< index_t > & tSelectedDofs = aIWG->dofs_per_face_on_sideset( tSideSetID );

                    if ( tSelectedDofs.length() > 0 )
                    {
                        // get the number of faces per element
                        uint tNumEntities = mesh::number_of_faces( tSideSet->element_type() );

                        // loop over all elements on this block
                        for( mesh::Facet * tFacet : tSideSet->facets() )
                        {
                            mesh::Element * tElement = tFacet->element() ;

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
                for( mesh::Facet * tFacet : mMesh->connectors() )
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
                for( uint s=0; s<2; ++s )
                {
                    Cell< mesh::Facet * > & tFacets = ( s == 0 ) ? mMesh->facets() : mMesh->connectors();

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
                    const index_t                              aProc,
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
            DofData::create_dof_table()
            {
                // initialize counter
                index_t tCount = 0 ;

                if ( mMyRank != mKernel->master() )
                {
                    // any other proc collects its DOF IDs and sends them to the master
                    index_t tNumberOfDOFs = mDOFs.size();

                    // collect IDs of DOFs
                    Vector< id_t > tIDs( tNumberOfDOFs );

                    for ( Dof * tDof : mDOFs )
                    {
                        tIDs( tCount++ ) = tDof->id() ;
                    }

                    // send IDs to master
                    send( mKernel->master(), tIDs );
                }
                else
                {
                    // get number of procs
                    uint tNumberOfProcs = mKernel->comm_table().length() ;

                    // initialize container
                    mDofIDs.set_size( tNumberOfProcs, {} );

                    // collect IDs from all other procs
                    receive( mKernel->comm_table(), mDofIDs );

                    // get ref to my IDs
                    Vector< id_t > & tMyIDs = mDofIDs( 0 );
                    tMyIDs.set_size( mDOFs.size() );

                    for ( Dof * tDof : mDOFs )
                    {
                        tMyIDs( tCount++ ) = tDof->id() ;
                    }

                    // now we can create the dof table that contains the dof indices
                    // per proc
                    mDofIndices.set_size( tNumberOfProcs, Vector< index_t >() );

                    // loop over all procs
                    for ( uint p = 0; p < tNumberOfProcs; ++p )
                    {
                        // get ID container
                        Vector< id_t > & tIDs = mDofIDs( p );

                        // get index list from dof table
                        Vector< index_t > & tIndices = mDofIndices( p );

                        // get number of DOFs for thos proc
                        index_t tNumberOfDOFs = tIDs.length();

                        // allocate memory for indices
                        tIndices.set_size( tNumberOfDOFs );

                        // populate Indices using map
                        for ( index_t k = 0; k < tNumberOfDOFs; ++k )
                        {
                            tIndices( k ) = mDofMap( tIDs( k ) )->index() ;
                        }
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            DofData::init_dof_values()
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
                                    tValues( tDof->dof_index_on_field() ) = tDof->value();
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

            void
            DofData::init_dirichlet_bcs()
            {

                uint tNumProcs = mKernel->number_of_procs() ;

                if( tNumProcs > 1 )
                {
                    if ( mMyRank == mKernel->master() )
                    {
                        // cell with dof indices of fixed dofs
                        Cell< Vector< id_t > > tAllIDs( tNumProcs, {} );

                        // cell with dof values
                        Cell< Vector< real > > tAllValues( tNumProcs, {} );

                        receive( mKernel->comm_table(), tAllIDs );
                        receive( mKernel->comm_table(), tAllValues );

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
                            Vector< id_t > & tDofIDs = mDofIDs( p );

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
                        send( mKernel->comm_table(), tAllIDs ) ;
                        send( mKernel->comm_table(), tAllValues );
                    }
                    else
                    {
                        // count fixed dofs
                        index_t tCount = 0;

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
                        send( mKernel->master(), tIDs );
                        send( mKernel->master(), tValues );
                        comm_barrier() ;

                        // receive confirmation from master
                        receive( mKernel->master(), tIDs );
                        receive( mKernel->master(), tValues );

                        // reset counter
                        tCount = 0 ;

                        // this->free_all_dofs() ;

                        // loop over all IDs
                        for( id_t tID : tIDs )
                        {
                            // get dof
                            Dof * tDof = this->dof( tID );
                            tDof->fix( tValues( tCount ) );
                            ++tCount;
                        }
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            DofData::compute_dof_indices()
            {
                if (  mMyRank != mKernel->master() )
                {
                    // get indices
                    Vector< index_t > tIndices ;
                    receive( mKernel->master(), tIndices );

                    BELFEM_ASSERT( tIndices.length() == mDOFs.size(),
                                  "Length of indices does not match. ( is %lu, expect %lu )",
                                  ( long unsigned int ) tIndices.length(),
                                  ( long unsigned int ) mDOFs.size() );

                    // receive global counters
                    Vector< index_t > tCounters( 2 );
                    receive( mKernel->master(), tCounters );
                    mNumberOfFreeDofs = tCounters( 0 );
                    mNumberOfFixedDofs = tCounters( 1 );

                    // reset counter
                    index_t tCount = 0 ;

                    // write indices into dofs
                    for( Dof * tDof : mDOFs )
                    {
                        tDof->set_index( tIndices( tCount++ ) );
                    }
                }
                else
                {
                    // step 1: count fixed and free dofs
                    mNumberOfFreeDofs = 0;
                    mNumberOfFixedDofs = 0;

                    // loop over all DOFs, note that on master,
                    // myindex and index are the same
                    for ( Dof * tDof : mDOFs )
                    {

                        if ( tDof->is_fixed() )
                        {
                            tDof->set_index( mNumberOfFixedDofs++ );
                        }
                        else
                        {
                            tDof->set_index( mNumberOfFreeDofs++ );
                        }
                    }

                    // get size of communication table
                    uint tNumProcs = mKernel->comm_table().length() ;

                    // this container holds the global indices for all procs
                    Cell< Vector< index_t > > tAllIndices( tNumProcs, {} );

                    // populate indices
                    for( uint p=1; p<tNumProcs; ++p )
                    {
                        // get id vector
                        const Vector< id_t > & tIDs = mDofIDs( p );

                        // get index vector
                        Vector< index_t >    & tIndices = tAllIndices( p );

                        // allocate vector
                        tIndices.set_size( tIDs.length() );

                        // initialize counter
                        index_t tCount = 0 ;

                        // loop over all IDs
                        for( id_t tID : tIDs )
                        {
                            // get index
                            tIndices( tCount++ ) = this->dof( tID )->index() ;
                        }
                    }

                    // send indices to procs
                    send( mKernel->comm_table(), tAllIndices );

                    // send global counters to other procs
                    Vector< index_t > tCounters( 2 );
                    tCounters( 0 ) = mNumberOfFreeDofs ;
                    tCounters( 1 ) = mNumberOfFixedDofs ;
                    send_same( mKernel->comm_table(), tCounters );
                }
            }

//------------------------------------------------------------------------------

            uint
            DofData::num_dofs_per_element( const id_t aBlockID ) const
            {
                // get block
                ElementType tType = mMesh->block( aBlockID )->element_type() ;

                uint aN = mParent->enforce_linear_interpolation() ?
                        mesh::number_of_corner_nodes( tType ) :
                        mesh::number_of_nodes( tType );

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
        } /* end namespace dofmgr */
    } /* end namespace fem */
} /* end namespace belfem */