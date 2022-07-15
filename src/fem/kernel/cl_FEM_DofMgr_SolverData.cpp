//
// Created by christian on 7/15/21.
//


#include "commtools.hpp"
#include "cl_Logger.hpp"
#include "cl_FEM_DofMgr_SolverData.hpp"
#include "cl_FEM_DofMgr_DofData.hpp"
#include "cl_FEM_DofMgr_BlockData.hpp"
#include "cl_FEM_DofMgr_SideSetData.hpp"
#include "cl_Mesh.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_DofManager.hpp"
#include "fn_max.hpp"
#include "fn_unique.hpp"
#include "fn_norm.hpp"
#include "cl_Timer.hpp"
#include "meshtools.hpp"

namespace belfem
{
    namespace fem
    {
        namespace dofmgr
        {
//------------------------------------------------------------------------------

            SolverData::SolverData( DofManager * aParent, DofData * aDofData ,
                                    BlockData * aBlockData,
                                    SideSetData * aSideSetData ) :
                    mParent( aParent ),
                    mKernel( aParent->parent() ),
                    mDofData( aDofData ),
                    mBlockData( aBlockData ),
                    mSideSetData( aSideSetData ),
                    mMyRank( aParent->rank() ),
                    mDOFs( aDofData->dofs() ),
                    mNumberOfFreeDofs( aDofData->number_of_free_dofs() ),
                    mNumberOfFixedDofs( aDofData->number_of_fixed_dofs() )
            {

            }

//------------------------------------------------------------------------------

            SolverData::~SolverData()
            {
                this->reset() ;

                // delete solver if it exists
                if( mSolver != nullptr )
                {
                    delete mSolver ;
                }
            }

//------------------------------------------------------------------------------

            void SolverData::reset()
            {
                mMyNumberOfFreeDofs = gNoIndex ;
                mMyNumberOfFixedDofs = gNoIndex ;

                if( mJacobian != nullptr )
                {
                    delete mJacobian ;
                    mJacobian = nullptr ;
                }
                if( mDirichletMatrix != nullptr )
                {
                    delete mDirichletMatrix ;
                    mDirichletMatrix = nullptr ;
                }

                mJacobianTable.clear() ;
                mDirichletTable.clear() ;
            }

//------------------------------------------------------------------------------

            void
            SolverData::allocate_matrices()
            {
                // restore factory settings
                this->reset() ;

                // local element-to-dof adjacency
                Vector< id_t > tElementWiseData ;
                this->compute_element_dof_connectivity( tElementWiseData );

                // local dof-to-element adjacency
                Vector< id_t > tDofWiseData ;
                this->compute_dof_element_connectivity( tDofWiseData );

                // create dof wise data
                Vector< id_t > tData ;

                if( mKernel->number_of_procs() > 1 )
                {
                    Cell< Vector< id_t > > tConnectivities( mKernel->number_of_procs(),
                                                            Vector< id_t >());
                    Vector< id_t > & tConnectivity = mMyRank == mKernel->master() ? tConnectivities( 0 ) : tData;
                    this->compute_dof_dof_connectivity( tDofWiseData, tElementWiseData, tConnectivity );

                    comm_barrier();

                    if ( mMyRank == mKernel->master() )
                    {
                        receive( mKernel->comm_table(), tConnectivities );
                        this->unite_dofs( tConnectivities, tData );
                    }
                    else
                    {
                        send( mKernel->master(), tData );
                    }
                }
                else
                {
                    this->compute_dof_dof_connectivity( tDofWiseData, tElementWiseData, tData );
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // write local node indices
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                mMyNumberOfFixedDofs = 0 ;
                mMyNumberOfFreeDofs = 0 ;

                for( Dof * tDof : mDOFs )
                {
                    if( tDof->is_fixed() )
                    {
                        tDof->set_my_index( mMyNumberOfFixedDofs++ );
                    }
                    else
                    {
                        tDof->set_my_index( mMyNumberOfFreeDofs++ );
                    }
                }

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // create the graph and initialize the matrices
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                Cell< graph::Vertex * > tGraph ;
                if( mMyNumberOfFixedDofs > 0 )
                {
                    this->populate_graph( tData, true, tGraph );
                    mDirichletMatrix = new SpMatrix( tGraph, SpMatrixType::CSR,
                                                     mNumberOfFreeDofs, mNumberOfFixedDofs );

                    tGraph.clear() ;
                }

                this->populate_graph( tData, false, tGraph );

                BELFEM_ASSERT( mSolver != nullptr, "no solver created" );

                mJacobian =  new SpMatrix( tGraph,
                                           mSolver->type() == SolverType::PETSC ?
                                           SpMatrixType::CSR : SpMatrixType::CSC,
                                           mNumberOfFreeDofs, mNumberOfFreeDofs );

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // allocate RHS
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // allocate right hand side
                if( mParent->iwg()->num_rhs_cols() <= 1 )
                {
                    mRhsVector.set_size( mMyNumberOfFreeDofs, 0.0 );
                }
                else
                {
                    mRhsMatrix.set_size( mMyNumberOfFreeDofs,
                                         mParent->iwg()->num_rhs_cols(), 0.0 );
                }

                comm_barrier() ;

                if( mMyRank == mKernel->master() )
                {
                    // get dof counters from clients
                    receive( mKernel->comm_table(), mNumberOfFreeDofsPerProc, mMyNumberOfFreeDofs );
                    receive( mKernel->comm_table(), mNumberOfFixedDofsPerProc, mMyNumberOfFixedDofs );
                }
                else
                {
                    // send counters to master
                    send( mKernel->master(), mMyNumberOfFreeDofs );
                    send( mKernel->master(), mMyNumberOfFixedDofs );
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::create_assembly_tables()
            {
                // get master proc
                proc_t tMaster = mKernel->master();

                if ( mMyRank == tMaster )
                {
                    const Vector< proc_t > & tComm = mKernel->comm_table();

                    Cell< Vector< int > > tJacobianRows;
                    Cell< Vector< int > > tJacobianCols;
                    Cell< Vector< int > > tDirichletRows;
                    Cell< Vector< int > > tDirichletCols;

                    // get data
                    receive( tComm, tJacobianRows );
                    receive( tComm, tJacobianCols );

                    uint tNumberOfProcs = tComm.length();

                    // allocate memory
                    mJacobianTable.set_size( tNumberOfProcs, {} );

                    // create jacobian indices
                    for ( uint p = 1; p < tNumberOfProcs; ++p )
                    {
                        // get index
                        Vector< index_t > & tIndex = mJacobianTable( p );

                        // get rows
                        Vector< int > & tRows = tJacobianRows( p );

                        // get  cols
                        Vector< int > & tCols = tJacobianCols( p );

                        // get number of nonzeros in this matrix
                        index_t tNNZ = tRows.length();

                        tIndex.set_size( tNNZ );

                        // loop over all entries
                        for ( index_t k = 0; k < tNNZ; ++k )
                        {
                            // compute index
                            tIndex( k ) = mJacobian->index( tRows( k ), tCols( k ) );
                        }
                    }

                    receive( tComm, tDirichletRows );
                    receive( tComm, tDirichletCols );

                    if ( mDirichletMatrix != NULL )
                    {
                        mDirichletTable.set_size( tNumberOfProcs, {} );

                        // create jacobian indices
                        for ( uint p = 1; p < tNumberOfProcs; ++p )
                        {
                            // get index
                            Vector< index_t > & tIndex = mDirichletTable( p );

                            // get rows
                            Vector< int > & tRows = tDirichletRows( p );

                            // get  cols
                            Vector< int > & tCols = tDirichletCols( p );

                            // get number of nonzeros in this matrix
                            index_t tNNZ = tRows.length();

                            if( tNNZ > 0 )
                            {
                                tIndex.set_size( tNNZ );

                                // loop over all entries
                                for ( index_t k = 0; k < tNNZ; ++k )
                                {
                                    // compute index
                                    tIndex( k ) = mDirichletMatrix->index( tRows( k ), tCols( k ) );
                                }
                            }
                        }
                    }


                }
                else
                {
                    // prepare data
                    mJacobian->create_coo_indices();

                    // send rows
                    send( tMaster, mJacobian->number_of_nonzeros(), mJacobian->rows() );

                    // send columns
                    send( tMaster, mJacobian->number_of_nonzeros(), mJacobian->cols() );

                    if ( mDirichletMatrix == NULL )
                    {

                        index_t tZero = 0;

                        // send zero to master ( for rows )
                        send( tMaster, tZero );
                        // send zero to master ( for cols )
                        send( tMaster, tZero );
                    }
                    else
                    {
                        // prepare data
                        mDirichletMatrix->create_coo_indices();

                        send( tMaster, mDirichletMatrix->number_of_nonzeros(), mDirichletMatrix->rows() );

                        send( tMaster, mDirichletMatrix->number_of_nonzeros(), mDirichletMatrix->cols() );
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::populate_graph( const Vector< id_t > & aData,
                                        const bool      aFixedFlag,
                                        Cell< graph::Vertex * > & aGraph )
            {
                aGraph.set_size( mMyNumberOfFreeDofs, nullptr );

                index_t tPivot = 0 ;

                // total number of dofs
                index_t tNumberOfDofs = aData( tPivot++ );

                index_t tNumDofs ;

                Cell< Dof * > tDofs( 1024, nullptr );
                index_t tCount ;

                // loop over all dofs
                for( index_t k=0; k<tNumberOfDofs; ++k )
                {
                    // grab dof
                    Dof * tA = mDofData->dof( aData( tPivot++ ) );
                    tA->reset_vertex_container() ;

                    tNumDofs = aData( tPivot++ );

                    BELFEM_ASSERT( tNumDofs < 1024, "Too many dofs" );

                    // check if dof is fidex
                    if( tA->is_fixed() )
                    {
                        tPivot += tNumDofs ;
                    }
                    else
                    {

                        tCount = 0 ;

                        for( index_t i=0; i<tNumDofs; ++i )
                        {
                            // get other dof
                            Dof * tB = mDofData->dof( aData( tPivot++ ) );

                            // check flag
                            if( tB->is_fixed() == aFixedFlag )
                            {
                                tDofs( tCount++ ) = tB ;
                            }
                        }

                        // allocate memory
                        tA->init_vertex_container( tCount );

                        for( index_t i=0; i<tCount; ++i )
                        {
                            tA->insert_vertex( tDofs( i )  );
                        }

                        aGraph( tA->my_index() ) = tA ;
                    }
                }



#if !defined( NDEBUG ) || defined( DEBUG )
                // check that graph is ok
                for( index_t k=0; k<aGraph.size(); ++k )
                {
                    BELFEM_ERROR( aGraph( k ) != nullptr, "entry %lu in graph is empty",
                                 ( long unsigned int ) k );
                }
#endif

            }

//------------------------------------------------------------------------------

            void
            SolverData::collect_fields( Cell< mesh::Field * > & aFields )
            {
                // get IWG
                IWG * tIWG = mParent->iwg() ;
                BELFEM_ERROR( tIWG != nullptr, "no equation was set" );

                // grab field data
                const Cell< string > & tFieldLabels = tIWG->all_fields() ;

                // how many fields exist on the mesh
                uint tNumFields = tFieldLabels.size() ;

                // count fields ( need to account for multiplicities )
                uint tCount = 0 ;
                for ( uint f = 0; f < tNumFields; ++f )
                {
                    switch ( mParent->mesh()->field( tFieldLabels( f ))->entity_type())
                    {
                        case ( EntityType::EDGE ) :
                        {
                            tCount += tIWG->edge_multiplicity();
                            break;
                        }
                        case ( EntityType::FACE ) :
                        {
                            tCount += tIWG->face_multiplicity();
                            break;
                        }
                        case ( EntityType::FACET ) :
                        {
                            tCount += tIWG->lambda_multiplicity();
                            break;
                        }
                        default :
                        {
                            tCount += 1;
                            break;
                        }
                    }
                }

                // allocate container
                aFields.set_size( tCount, nullptr );

                // reset counter
                tCount = 0 ;

                for ( uint f = 0; f < tNumFields; ++f )
                {
                    // grab field
                    mesh::Field * tField = mParent->mesh()->field( tFieldLabels( f ));

                    uint tMultiplicity;

                    switch ( tField->entity_type())
                    {
                        case ( EntityType::EDGE ) :
                        {
                            tMultiplicity = tIWG->edge_multiplicity();
                            break;
                        }
                        case ( EntityType::FACE ) :
                        {
                            tMultiplicity = tIWG->face_multiplicity();
                            break;
                        }
                        case ( EntityType::FACET ) :
                        {
                            tMultiplicity = tIWG->lambda_multiplicity();
                            break;
                        }
                        default :
                        {
                            tMultiplicity = 1;
                            break;
                        }
                    }


                    for( uint i=0; i<tMultiplicity; ++i )
                    {
                        aFields( tCount++ ) = tField ;
                    }

                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::compute_element_dof_connectivity( Vector< id_t > & aData )
            {
                // contains
                // number of elements
                // element id
                // number of dofs
                // dof index

                // set temporary dof index
                index_t tCount = 0 ;
                for( Dof * tDof : mDOFs )
                {
                    tDof->set_my_index( tCount++ );
                }

                // counter for number of elements
                index_t tElemCount = 0 ;

                // determine memory size
                tCount = 1 ;

                for( Block * tBlock : mBlockData->blocks() )
                {
                    tElemCount += tBlock->number_of_elements();

                    tCount += tBlock->number_of_elements() * (
                            mParent->iwg()->number_of_dofs_per_element( tBlock ) + 2 );
                }

                for( id_t tID : mParent->iwg()->selected_sidesets() )
                {
                    SideSet * tSideSet = mSideSetData->sideset( tID );

                    tElemCount += tSideSet->number_of_elements();
                    tCount += tSideSet->number_of_elements() * (
                            mParent->iwg()->number_of_dofs_per_element( tSideSet ) + 2 );
                }

                aData.set_size( tCount );

                // reset the counter
                tCount = 0 ;

                aData( tCount++ ) = tElemCount ;

                for( Block * tBlock : mBlockData->blocks() )
                {
                    id_t tNumDofsPerElement
                            = mParent->iwg()->number_of_dofs_per_element( tBlock );

                    for ( Element * tElement: tBlock->elements())
                    {
                        aData( tCount++ ) = tElement->id();
                        aData( tCount++ ) = tNumDofsPerElement;

                        BELFEM_ASSERT( tElement->number_of_dofs() == tNumDofsPerElement,
                                      "number of dofs do not match for element %lu on block %lu  is %lu but expect %lu",
                                      ( long unsigned int ) tElement->id(),
                                      ( long unsigned int ) tBlock->id(),
                                      ( long unsigned int ) tElement->number_of_dofs(),
                                      ( long unsigned int ) tNumDofsPerElement );

                        for ( uint d=0; d<tNumDofsPerElement; ++d )
                        {
                            aData( tCount++ ) = tElement->dof( d )->my_index() ;
                        }
                    }
                }

                for( id_t tID : mParent->iwg()->selected_sidesets() )
                {
                    SideSet * tSideSet = mSideSetData->sideset( tID );

                    uint tNumDofsPerElement
                            = mParent->iwg()->number_of_dofs_per_element( tSideSet );

                    for ( Element * tElement: tSideSet->elements() )
                    {

                        aData( tCount++ ) = tElement->id();
                        aData( tCount++ ) = tNumDofsPerElement;

                        /*if(  tElement->id() == 25 )
                        {
                            std::cout << "sideset " << tSideSet->id() << " " << ( uint ) tSideSet->domain_type() << std::endl ;
                            mParent->iwg()->print_dofs( tElement );
                            exit( 0 );
                        }*/

                        BELFEM_ASSERT( tElement->number_of_dofs() == tNumDofsPerElement,
                                      "number of dofs do not match for element %lu on sideset %lu  is %u but expect %u",
                                      ( long unsigned int ) tElement->id(),
                                      ( long unsigned int ) tSideSet->id(),
                                      ( unsigned int ) tElement->number_of_dofs(),
                                      ( unsigned int ) tNumDofsPerElement );

                        for ( uint d=0; d<tNumDofsPerElement; ++d )
                        {
                            aData( tCount++ ) = tElement->dof( d )->my_index();
                        }
                    }
                }

                BELFEM_ASSERT( tCount == aData.length(), "something went wrong");

            }

//------------------------------------------------------------------------------

            void
            SolverData::compute_dof_element_connectivity( Vector< id_t > & aData )
            {
                index_t tNumberOfDofs = mDOFs.size() ;

                // count elements per dof
                Vector< id_t > tWork( tNumberOfDofs, 0 );

                index_t tCount = 0 ;

                // set temporary index for dofs
                for( Dof * tDof : mDOFs )
                {
                    tDof->set_my_index( tCount++ );
                }

                // contains
                // number of dofs
                // dof index
                // number of elements
                // element ids

                tCount = mDOFs.size() + 1 ;

                for( Block * tBlock : mBlockData->blocks() )
                {
                    for( Element * tElement : tBlock->elements() )
                    {
                        tCount += tElement->number_of_dofs() ;

                        for( uint d=0; d<tElement->number_of_dofs(); ++d )
                        {
                            ++tWork( tElement->dof( d )->my_index() ) ;
                        }
                    }
                }

                for( id_t tID : mParent->iwg()->selected_sidesets() )
                {
                    SideSet * tSideSet = mSideSetData->sideset( tID );

                    for( Element * tElement : tSideSet->elements() )
                    {
                        tCount += tElement->number_of_dofs() ;

                        for( uint d=0; d<tElement->number_of_dofs(); ++d )
                        {
                            ++tWork( tElement->dof( d )->my_index() ) ;
                        }
                    }
                }

                aData.set_size( tCount );
                tCount = 0;
                aData( tCount++ ) = tNumberOfDofs ;

                id_t tSwap ;

                for( index_t k=0; k<tNumberOfDofs; ++k )
                {
                    // tWork(k) contains number of elements
                    aData( tCount++ ) = tWork( k );

                    // save in temporary index
                    tSwap = tWork( k );

                    // overwrite work with memory offset
                    tWork( k ) = tCount ;

                    // increment memory counter
                    tCount += tSwap ;
                }

                // populate vector
                tCount = 0 ;
                aData( tCount++ ) = mDOFs.size();

                for( Block * tBlock : mBlockData->blocks() )
                {
                    for( Element * tElement : tBlock->elements() )
                    {
                        for( uint d=0; d<tElement->number_of_dofs(); ++d )
                        {
                            aData( tWork( tElement->dof( d )->my_index() )++ ) = tElement->id() ;
                        }
                    }
                }

                for( id_t tID : mParent->iwg()->selected_sidesets() )
                {
                    SideSet * tSideSet = mSideSetData->sideset( tID );

                    for( Element * tElement : tSideSet->elements() )
                    {
                        tCount += tElement->number_of_dofs() ;

                        for( uint d=0; d<tElement->number_of_dofs(); ++d )
                        {
                            aData( tWork( tElement->dof( d )->my_index() )++ )
                                    = tElement->id() ;
                        }
                    }
                }

            }

//------------------------------------------------------------------------------

            void
            SolverData::compute_dof_dof_connectivity(
                    const Vector< id_t > & aDofWiseData,
                    const Vector< id_t > & aElementWiseData,
                    Vector< id_t > & aConnectivity )
            {

                // a memory pointer
                index_t tPivot = 0 ;

                // create a temporary map that keeps track of the element offsets
                Map< id_t, index_t > tOffsets ;

                // get the number of elements
                index_t tNumElems = aElementWiseData( tPivot++ );
                index_t tNumDofs ;

                id_t tElementID ;

                // loop over all elements
                for( index_t e=0; e<tNumElems; ++e )
                {
                    // get the element ID
                    tElementID = aElementWiseData( tPivot++ );

                    // remember the memory position
                    tOffsets[ tElementID ] = tPivot ;

                    // get the number of DOFs for this element
                    tNumDofs = aElementWiseData( tPivot++ );

                    // jump pivot to next element ID
                    tPivot += tNumDofs ;
                }

                // reset the pivot
                tPivot = 0 ;
                index_t tOffset = 0 ;

                // number of DOFs on this proc
                index_t tNumberOfDofs = aDofWiseData( tPivot++ );

                // count number of dofs per dof
                Vector< index_t > tCount( tNumberOfDofs, 0 );

                for( index_t k=0; k<tNumberOfDofs; ++k )
                {
                    // get the number of elements for this dof
                    tNumElems = aDofWiseData( tPivot++ );

                    // loop over all elements of this dof
                    for( index_t e=0; e<tNumElems; ++e )
                    {
                        // get the element ID
                        tElementID = aDofWiseData( tPivot++ );

                        // get the offset in the other vector
                        tOffset = tOffsets[ tElementID ];

                        // get the number of DOFs for this element
                        tCount( k ) += aElementWiseData( tOffset++ );
                    }
                }

                // allocate memory
                // allocate list with IDs
                Cell< Vector< id_t > > tAllIDs( tNumberOfDofs,
                                                Vector< id_t >() );
                for( index_t k=0; k<tNumberOfDofs; ++k )
                {
                    tAllIDs( k ).set_size( tCount( k ) );
                }

                tPivot = 1 ;

                index_t j ;

                index_t tMemCount = 1 ;

                for( index_t k=0; k<tNumberOfDofs; ++k )
                {
                    // get the number of elements for this dof
                    tNumElems = aDofWiseData( tPivot++ );

                    // get ID list
                    Vector< id_t > & tIDs = tAllIDs( k );

                    if( tIDs.length() > 0 )
                    {
                        j = 0;

                        // loop over all elements of this dof
                        for ( index_t e = 0; e < tNumElems; ++e )
                        {
                            // get the element ID
                            tElementID = aDofWiseData( tPivot++ );

                            // get the offset in the other vector
                            tOffset = tOffsets[ tElementID ];

                            // get the number of DOFs for this element
                            tNumDofs = aElementWiseData( tOffset++ );

                            for ( index_t i = 0; i < tNumDofs; ++i )
                            {
                                tIDs( j++ ) = aElementWiseData( tOffset++ );
                            }
                        }
                        BELFEM_ASSERT( j == tCount( k ),"memory error" );

                        unique( tIDs );
                    }
                    tMemCount += 2 + tIDs.length();
                }

                aConnectivity.set_size( tMemCount );
                tPivot = 0 ;
                aConnectivity( tPivot++ ) = tNumberOfDofs ;

                for( index_t k=0; k<tNumberOfDofs; ++k )
                {

                    Vector< id_t > & tIDs = tAllIDs( k );

                    aConnectivity( tPivot++ ) = mDOFs( k )->id() ;

                    tNumDofs = tIDs.length() ;

                    aConnectivity( tPivot++ ) = tIDs.length() ;

                    for( index_t i=0; i<tNumDofs; ++i )
                    {
                        aConnectivity( tPivot++ ) = mDOFs( tIDs( i ) )->id() ;
                    }
                }

                BELFEM_ASSERT( tPivot == aConnectivity.length(), "memory error");
            }


//------------------------------------------------------------------------------

            void
            SolverData::unite_dofs(
                    const Cell< Vector<id_t > > & aConnectivities,
                    Vector< id_t > & aConnectivity )
            {
                index_t tNumberOfDofs = mDOFs.size() ;

                // maximum number of dofs, may include multiple entries
                Vector< index_t > tCount( mDOFs.size(), 0 ) ;

                const uint tNumProcs = aConnectivities.size() ;

                index_t tPivot  ;
                index_t tNumDofs ;

                // loop over all procs and count memory needs
                for( uint p=0; p<tNumProcs; ++p )
                {
                    // get vector with proc-wise data
                    const Vector< id_t > & tData = aConnectivities( p );

                    // reset pivot
                    tPivot = 0 ;

                    // get number of dofs
                    index_t tNumberOfDofsPerProc = tData( tPivot++ );

                    for( index_t k=0; k<tNumberOfDofsPerProc; ++k )
                    {
                        // get dof
                        Dof * tDof = mDofData->dof( tData( tPivot++ ));

                        // get number of dofs per dof
                        tNumDofs = tData( tPivot++ );

                        // add to counter
                        tCount( tDof->my_index() ) += tNumDofs ;

                        // jump
                        tPivot += tNumDofs;
                    }
                }

                // allocate list with IDs
                Cell< Vector< id_t > > tAllIDs( tNumberOfDofs,
                                                Vector< id_t >() );
                for( index_t k=0; k<tNumberOfDofs; ++k )
                {
                    tAllIDs( k ).set_size( tCount( k ) );
                }
                tCount.fill( 0 );

                // loop over all procs and combine dofs
                for( uint p=0; p<tNumProcs; ++p )
                {
                    // get vector with proc-wise data
                    const Vector< id_t > & tData = aConnectivities( p );

                    // reset pivot
                    tPivot = 0 ;

                    // get number of dofs on this proc
                    index_t tNumberOfDofsPerProc = tData( tPivot++ );

                    // populate cell
                    for( index_t k=0; k<tNumberOfDofsPerProc; ++k )
                    {
                        // get dof index
                        index_t tIndex = mDofData->dof( tData( tPivot++ ))->my_index() ;

                        // get ID list
                        Vector< id_t > & tIDs = tAllIDs( tIndex );

                        // get number of dofs per dof
                        tNumDofs = tData( tPivot++ );

                        // add dofs to ID list
                        for( index_t i=0; i<tNumDofs; ++i )
                        {
                            tIDs( tCount( tIndex )++ ) = tData( tPivot++ );
                        }

                    }
                }

                // unify data
                // memory counter
                tPivot = 2 * tNumberOfDofs + 1 ;

                for( index_t k=0; k<tNumberOfDofs; ++k )
                {
                    unique( tAllIDs( k ) );
                    tPivot += tAllIDs( k ).length() ;
                }

                aConnectivity.set_size( tPivot );
                tPivot = 0 ;

                aConnectivity( tPivot++ ) = tNumberOfDofs ;

                for( index_t k=0; k<tNumberOfDofs; ++k )
                {
                    Vector< id_t > & tIDs = tAllIDs( k );

                    aConnectivity( tPivot++ ) = mDOFs( k )->id();
                    tNumDofs = tIDs.length() ;

                    aConnectivity( tPivot++ ) = tNumDofs ;

                    for( index_t i=0; i<tNumDofs; ++i )
                    {
                        aConnectivity( tPivot++ ) = tIDs( i );
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::reset_matrices()
            {
                BELFEM_ASSERT( mJacobian != nullptr,
                              "Jacobian Matrix was not initialized");

                mJacobian->fill( 0.0 );

                if ( mDirichletMatrix != nullptr )
                {
                    mDirichletMatrix->fill( 0.0 );
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::reset_rhs_vector()
            {
                mRhsVector.fill( 0.0 );
            }

//------------------------------------------------------------------------------

            void
            SolverData::reset_rhs_matrix()
            {
                mRhsMatrix.fill( 0.0 );
            }

//------------------------------------------------------------------------------

            void
            SolverData::assemble_jacobian(
                    Element * aElement,
                    const Matrix< real > & aJacobian )
            {
                SpMatrix & J       =  *mJacobian;
                SpMatrix & D       =  *mDirichletMatrix;

                // get dimension of element Jacobian
                uint tN = aElement->number_of_dofs() ;

                // add element jacobian to system matrices
                for ( uint i = 0; i < tN; ++i )
                {
                    Dof * tRow = aElement->dof( i );
                    if ( !tRow->is_fixed() )
                    {
                        for ( uint j = 0; j < tN; ++j )
                        {
                            Dof * tCol = aElement->dof( j );

                            if ( tCol->is_fixed() )
                            {
                                D( tRow->index(), tCol->index() ) -= aJacobian( i, j );
                            }
                            else
                            {
                                J( tRow->index(), tCol->index() ) += aJacobian( i, j );
                            }
                        }
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::assemble_jacobian_and_rhs( Element * aElement,
                                                   const Matrix< real > & aJacobian,
                                                   const Vector< real > & aResidual )
            {
                SpMatrix & J       =  *mJacobian;
                SpMatrix & D       =  *mDirichletMatrix;

                // get dimension of element Jacobian
                uint tN = aElement->number_of_dofs() ;

                // add element jacobian to system matrices
                for ( uint i = 0; i < tN; ++i )
                {
                    Dof * tRow = aElement->dof( i );
                    if ( !tRow->is_fixed() )
                    {
                        for ( uint j = 0; j < tN; ++j )
                        {
                            Dof * tCol = aElement->dof( j );

                            if ( !tCol->is_fixed() )
                            {
                                J( tRow->index(), tCol->index() ) += aJacobian( i, j );
                            }
                            else
                            {
                                D( tRow->index(), tCol->index() ) -= aJacobian( i, j );
                            }
                        }
                    }
                }



                // add residual to vector
                for ( uint i = 0; i < tN; ++i )
                {
                    Dof * tRow = aElement->dof( i );
                    if ( !tRow->is_fixed() )
                    {
                        mRhsVector( tRow->my_index() ) += aResidual( i );
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::asseble_rhs( Element * aElement,
                                     const Vector< real > & aRHS )
            {
                // get dimension of element Jacobian
                uint tN = aElement->number_of_dofs() ;

                // add residual to vector
                for ( uint i = 0; i < tN; ++i )
                {
                    Dof * tRow = aElement->dof( i );
                    if ( !tRow->is_fixed() )
                    {
                        mRhsVector( tRow->my_index() ) += aRHS( i );
                    }
                }
            }

            //------------------------------------------------------------------------------

            void
            SolverData::assemble_volume_loads( Element * aElement,
                                               const Vector< real > & aRHS )
            {

                // get dimension of element Jacobian
                uint tN = aElement->number_of_dofs() ;

                // add residual to vector
                for ( uint i = 0; i < tN; ++i )
                {
                    Dof * tRow = aElement->dof( i );
                    if ( !tRow->is_fixed() )
                    {
                        mVolumeLoads( tRow->my_index() ) += aRHS( i );
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::asseble_surface_loads( Element * aElement,
                                               const Vector< real > & aRHS )
            {

                // get dimension of element Jacobian
                uint tN = aElement->number_of_dofs() ;

                // add residual to vector
                for ( uint i = 0; i < tN; ++i )
                {
                    Dof * tRow = aElement->dof( i );
                    if ( !tRow->is_fixed() )
                    {
                        mConvection( tRow->my_index() ) += aRHS( i );
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::asseble_rhs(
                    Element * aElement,
                    const Matrix< real > & aRHS )
            {
                // get dimension of element Jacobian
                uint tN = aElement->number_of_dofs() ;
                uint tM = aRHS.n_cols() ;

                // add residual to vector
                for ( uint i = 0; i < tN; ++i )
                {
                    Dof * tRow = aElement->dof( i );
                    if ( !tRow->is_fixed() )
                    {
                        for( uint j=0; j<tM; ++j )
                        {
                            mRhsMatrix( tRow->my_index(), j ) += aRHS( i, j );
                        }
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::collect_jacobian()
            {
                if ( mMyRank == mKernel->master() )
                {

                    // get number of procs
                    uint tNumberOfProcs = mKernel->comm_table().length();

                    Cell< Vector< real > > tJacobian;

                    receive( mKernel->comm_table(), tJacobian );

                    // assemble jacobian
                    for ( uint p = 1; p < tNumberOfProcs; ++p )
                    {
                        // get data
                        Vector< real > & tData = tJacobian( p );

                        // get indices
                        Vector< index_t > & tIndices = mJacobianTable( p );

                        // get number of nonzeros
                        index_t tNNZ = tData.length();

                        BELFEM_ASSERT( tIndices.length() == tData.length(),
                                      "Jacobian Matrix from proc %u has wrong number of nonzeros ( is %lu, expect %lu )",
                                      ( unsigned int ) mKernel->comm_table( p ),
                                      ( unsigned int ) tData.length(),
                                      ( unsigned int ) tIndices.length() );
                        // loop over all entries
                        for ( index_t k = 0; k < tNNZ; ++k )
                        {
                            // add entries
                            mJacobian->data( tIndices( k ) ) += tData( k );
                        }

                    }

                    // tidy up memory a bit
                    tJacobian.clear();

                    Cell< Vector< real > > tDirichlet( tNumberOfProcs, {} );

                    // assemble Dirichlet matrix
                    if ( mNumberOfFixedDofs > 0 )
                    {
                        // collect data from clients
                        for ( uint p = 1; p < tNumberOfProcs; ++p )
                        {
                            if ( mNumberOfFixedDofsPerProc( p ) > 0 )
                            {
                                // get data
                                Vector< real > & tData = tDirichlet( p );

                                receive( mKernel->comm_table( p ), tData );
                            }
                        }
                        // assemble matrix
                        for ( uint p = 1; p < tNumberOfProcs; ++p )
                        {
                            // get data
                            Vector< real > & tData = tDirichlet( p );

                            // get indices
                            Vector< index_t > & tIndices = mDirichletTable( p );

                            // get number of nonzeros
                            index_t tNNZ = tData.length() ;

                            BELFEM_ASSERT( tIndices.length() == tData.length(),
                                          "Dirichlet Matrix from proc %u has wrong number of nonzeros ( is %lu, expect %lu )",
                                          ( unsigned int ) mKernel->comm_table( p ),
                                          ( unsigned int ) tData.length(),
                                          ( unsigned int ) tIndices.length() );

                            // loop over all entries
                            for ( index_t k = 0; k < tNNZ ; ++k )
                            {
                                // add entries
                                mDirichletMatrix->data( tIndices( k ) ) += tData( k );
                            }
                        }
                    }
                }
                else
                {
                    // send my data to master
                    send( mKernel->master(),
                          mJacobian->number_of_nonzeros(),
                          mJacobian->data());

                    if ( mMyNumberOfFixedDofs > 0 )
                    {
                        send( mKernel->master(),
                              mDirichletMatrix->number_of_nonzeros(),
                              mDirichletMatrix->data() );
                    }

                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::collect_rhs_vector()
            {
                this->collect_vector( mRhsVector );
            }

//------------------------------------------------------------------------------

            void
            SolverData::collect_vector( Vector< real > & aVector )
            {
                if ( mMyRank == mKernel->master() )
                {
                    const Vector< proc_t > & tComm = mKernel->comm_table() ;

                    Cell< Vector< real > > tAllVectors;
                    receive( tComm, tAllVectors );

                    // assemble system
                    for ( uint p = 1; p < tComm.length(); ++p )
                    {
                        // get dof table for this proc
                        const Vector< index_t > & tDOFs = mDofData->dof_indices( p );

                        Vector< real > & tVector = tAllVectors( p );

                        // sometimes, the vector may be of zero length
                        // eg, if a proc does not have a wetted surface
                        if ( tVector.length() > 0 )
                        {
                            // get number of dofs fixme: check this
                            index_t tNumDOFs = tDOFs.length();

                            index_t tCount = 0 ;

                            // loop over all dofs
                            for ( index_t i = 0; i < tNumDOFs; ++i )
                            {
                                // get dof
                                Dof * tDOF = mDOFs( tDOFs( i ) );

                                if ( !tDOF->is_fixed() )
                                {
                                    aVector( tDOF->index() ) += tVector( tCount++ );
                                }
                            }
                        }
                    }
                }
                else
                {
                    // send vector to master
                    send( mKernel->master(), aVector );
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::collect_rhs_matrix()
            {
                // number of cols in rhs matrix
                uint tNumCols = mParent->iwg()->num_rhs_cols();

                if ( mMyRank == mKernel->master() )
                {
                    Cell< Matrix< real > > tAllRHS;
                    receive( mKernel->comm_table(), tAllRHS );

                    // assemble system
                    for ( uint p = 1; p < mKernel->comm_table().length(); ++p )
                    {
                        // get dof table for this proc
                        const Vector< index_t > & tDOFs = mDofData->dof_indices( p );

                        // get the current data
                        Matrix< real > & tRHS = tAllRHS( p );

                        // get number of dofs
                        index_t tNumDOFs = tDOFs.length();

                        // loop over all dofs
                        for ( index_t j = 0; j < tNumCols; ++j )
                        {
                            for ( index_t i = 0; i < tNumDOFs; ++i )
                            {
                                // get dof
                                Dof * tDOF = mDOFs( tDOFs( i ));

                                if ( !tDOF->is_fixed() )
                                {
                                    mRhsMatrix( tDOF->index(), j ) += tRHS( i, j );
                                }
                            }
                        }
                    }
                }
                else
                {
                    // send vector to master
                    send( mKernel->master(), mRhsMatrix );
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::update_field_values()
            {
                if ( mMyRank == mKernel->master() )
                {
                    // allocate vector
                    if ( mFieldValues.length() != mNumberOfFreeDofs )
                    {
                        mFieldValues.set_size( mNumberOfFreeDofs, 0.0 );
                    }

                    // collect values for free dofs
                    for ( Dof * tDof: mDOFs )
                    {
                        if ( !tDof->is_fixed() )
                        {
                            mFieldValues( tDof->index() ) = tDof->value();
                        }
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::set_solver( const SolverType aSolver )
            {
                // delete solver if it exists already
                if( mSolver != nullptr )
                {
                    delete mSolver ;
                }

                // create a new solver
                mSolver = new Solver( aSolver );

                mSolver->set_symmetry_mode( mParent->iwg()->symmetry_mode() );
            }

//------------------------------------------------------------------------------

            void
            SolverData::solve()
            {
                // get pointer to the equation
                IWG * tIWG = mParent->iwg() ;
                BELFEM_ERROR( tIWG != nullptr, "no equation was set" );

                Cell< mesh::Field * > tFields;
                this->collect_fields( tFields );

                if ( mKernel->is_master() )
                {

                    Timer tTimer;

                    // right hand side
                    Vector< real > tFixedValues( mNumberOfFixedDofs );

                    index_t tCount = 0;

                    // loop over all dofs
                    for ( Dof * tDof: mDOFs )
                    {
                        if ( tDof->is_fixed() )
                        {
                            tFixedValues( tCount++ ) = tDof->value();
                        }
                    }

                    // add loads over boundary
                    if( mConvection.length() > 0 )
                    {
                        mRhsVector += mConvection ;
                    }

                    // add volume loads
                    if( mVolumeLoads.length() > 0 )
                    {
                        mRhsVector += mVolumeLoads ;
                    }

                    if ( mNumberOfFixedDofs != 0 )
                    {
                        BELFEM_ASSERT( tIWG->num_rhs_cols() == 1,
                                      "Can only impose values of RHS is a vector, not a matrix!");

                        mDirichletMatrix->multiply( tFixedValues, mRhsVector, 1.0, 1.0 );
                    }

                    if(  tIWG->num_rhs_cols() == 1 ) // right hand side is vector
                    {
                        // compute the norm of the rhs vector
                        mRhsNorm = norm( mRhsVector );

                        switch( tIWG->mode() )
                        {
                            case( IwgMode::Direct ) :
                            {
                                // wait for other procs
                                comm_barrier() ;

                                // solve the system
                                mSolver->solve( *mJacobian, mLhsVector, mRhsVector ) ;

                                // write values into field
                                for ( Dof * tDof: mDOFs )
                                {
                                    if ( tDof->is_fixed() )
                                    {
                                        tFields( tDof->type_id() )->value(
                                                tDof->dof_index_on_field() ) = tDof->value();
                                    }
                                    else
                                    {
                                        tFields( tDof->type_id() )->value(
                                                tDof->dof_index_on_field() )
                                                = mLhsVector( tDof->index() );
                                    }
                                }

                                break ;
                            }
                            case( IwgMode::Iterative ) :
                            {

                                BELFEM_ASSERT( mFieldValues.length() == mRhsVector.length(),
                                              "Length of Field values and RHS vector do not match ( %lu vs. %lu, free dofs: %lu )",
                                              ( long unsigned int ) mFieldValues.length(),
                                              ( long unsigned int ) mRhsVector.length(),
                                              ( long unsigned int ) mNumberOfFreeDofs );


                                switch( tIWG->algorithm() )
                                {
                                    case( SolverAlgorithm::NewtonRaphson ) :
                                    {
                                        // compute the residual as r = A * x - b and write it into RHS vector
                                        mJacobian->multiply( mFieldValues, mRhsVector, 1.0, -1.0 );

                                        // wait for other procs
                                        comm_barrier() ;

                                        // solve the system
                                        mSolver->solve( *mJacobian, mLhsVector, mRhsVector ) ;

                                        for ( Dof * tDof: mDOFs )
                                        {
                                            // update DOF values
                                            if ( !tDof->is_fixed() )
                                            {
                                                tDof->value() -= tIWG->omega() * mLhsVector( tDof->index() );
                                            }

                                            // update value in field
                                            tFields( tDof->type_id() )->value(
                                                    tDof->dof_index_on_field() ) = tDof->value();
                                        }

                                        break ;
                                    }
                                    case( SolverAlgorithm::Picard ) :
                                    {
                                        // wait for other procs
                                        comm_barrier() ;

                                        // solve the system
                                        mSolver->solve( *mJacobian, mLhsVector, mRhsVector ) ;

                                        real tA = 1. - tIWG->omega() ;
                                        real tB = tIWG->omega() ;

                                        for ( Dof * tDof: mDOFs )
                                        {
                                            // update DOF values
                                            if ( !tDof->is_fixed() )
                                            {
                                                tDof->value() *= tA ;
                                                tDof->value() += tB * mLhsVector( tDof->index() );
                                            }

                                            // update value in field
                                            tFields( tDof->type_id() )->value(
                                                    tDof->dof_index_on_field() ) = tDof->value() ;
                                        }

                                        // compute the residual as r = A * x - b and write it into RHS vector
                                        mJacobian->multiply( mFieldValues, mRhsVector, 1.0, -1.0 );


                                        break ;
                                    }
                                    default :
                                    {
                                        BELFEM_ERROR( false, "Unsupported Solve Algorithm");
                                    }
                                }

                                break ;
                            }
                            default :
                            {
                                BELFEM_ERROR( false, "Undefinded IWG Mode");
                            }
                        } // end IWG mode
                    }
                    else // right hand side is matrix
                    {
                        BELFEM_ASSERT( tIWG->mode() == IwgMode::Direct, "IWG must be direct when rhs is matrix" );

                        // perform a sanity check
                        const Vector< id_t> & tBlocks = tIWG->selected_blocks();
                        uint tNumDofsPerNode = tIWG->number_of_dofs_per_node( tBlocks( 0 ) );

                        for( uint b=1; b<tBlocks.length(); ++b )
                        {
                            BELFEM_ERROR( tNumDofsPerNode == tIWG->number_of_dofs_per_node( tBlocks( b ) ),
                                         "All selected blocks must have the same number of dofs per node" );
                        }
                        for( id_t b: tBlocks )
                        {
                            BELFEM_ERROR( tIWG->number_of_dofs_per_edge( b ) == 0,
                                         "Edge DOFs not supported if using an RHS matrix" );
                        }
                        // note: this only works for node fields
                        // wait for other procs
                        comm_barrier() ;
                        mSolver->solve( *mJacobian, mLhsMatrix, mRhsMatrix );

                        uint k = 0;
                        for( index_t j=0; j< tIWG->num_rhs_cols(); ++j )
                        {
                            for ( index_t i = 0; i <  tIWG->number_of_dofs_per_node(); ++i )
                            {
                                // get field
                                mesh::Field * tField = mParent->mesh()->field( tIWG->all_fields()( k++ ) );

                                // sanity check
                                BELFEM_ERROR( tField->entity_type() == EntityType::NODE,
                                             "RHS matrices require all fields to be nodal but field %s is not a node field.",
                                             tField->label().c_str() );

                                // get field data
                                Vector< real > & tData = tField->data();

                                // loop over all nodes
                                for( mesh::Node *tNode : mParent->mesh()->nodes() )
                                {
                                    // get dof
                                    Dof * tDOF = mDofData->dof( mDofData->node_dof_id( tNode->id(), i ) );

                                    // write data of dof into field
                                    tData( tDOF->dof_index_on_field() ) = mLhsMatrix( tDOF->index(), j );
                                }
                            }
                        }
                    }

                    message( 4, "    ... time for solving system of equations    : %u ms\n",
                             ( unsigned int ) tTimer.stop());
                }
                else if ( mSolver->type() == SolverType::MUMPS ||
                          mSolver->type() == SolverType::PETSC )
                {
                    if ( tIWG->num_rhs_cols() == 1 )
                    {
                        // wait for other procs
                        comm_barrier() ;
                        mSolver->solve( *mJacobian, mLhsVector, mRhsVector ) ;
                    }
                    else
                    {
                        // wait for other procs
                        comm_barrier() ;
                        mSolver->solve( *mJacobian, mLhsMatrix, mRhsMatrix ) ;
                    }
                }
                else
                {
                    // wait for other procs
                    comm_barrier() ;
                }
            }

//------------------------------------------------------------------------------

            real
            SolverData::residual( const uint aIteration )
            {
                if( mKernel->is_master() )
                {
                    // compute value
                    // note that this vector now contains the error r=A*x-b
                    // while the value of mRhsNorm was computed before with the real rhs vectror
                    real aResidual = norm( mRhsVector ) / mRhsNorm ;

                    // distribute data
                    Vector< real > tResidual( mKernel->comm_table().length(), aResidual );

                    // catch error
                    if( ( aResidual == 0  && aIteration == 0 ) || aResidual > 1E12 || std::isnan( aResidual ) )
                    {

                        // mParent->save_system("error.hdf5");

                        BELFEM_ERROR( false,
                                     "ITERATION SCHEME FAILED: R2 = %8.2g \nYOU CAN TRY THE FOLLIOWING THINGS: \n    * CHECK YOUR BOUNDARY CONDITIONS\n    * MAKE A BETTER CONDITIONED MESH\n    * USE A DIFFERENT PRECONDITIONER OR A DIRECT SOLVER\n    * DECREASE THE TIMESTEP\n    * DECREASE THE RELAXATION PARAMETER\n", aResidual );
                    }

                    send( mKernel->comm_table(), tResidual );

                    return aResidual;
                }
                else
                {
                    // get value from master
                    real aResidual = BELFEM_REAL_MAX ;

                    receive( mKernel->master(), aResidual );

                    return  aResidual;
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::save_system( const string & aPath )
            {
#ifdef BELFEM_HDF5

                if( mParent->parent()->is_master() )
                {

                    HDF5 tFile( aPath, FileMode::NEW );

                    this->save_system( tFile ) ;

                    tFile.close() ;
                }
#endif
            }

//------------------------------------------------------------------------------

            void
            SolverData::load_system( const string & aPath )
            {
#ifdef BELFEM_HDF5

                if( mParent->parent()->is_master() )
                {

                    HDF5 tFile( aPath, FileMode::OPEN_RDONLY );

                    this->load_system( tFile ) ;

                    tFile.close() ;
                }
#endif
            }

//------------------------------------------------------------------------------

#ifdef BELFEM_HDF5
            void
            SolverData::load_system( HDF5 & aFile )
            {
                herr_t tError = 0 ;

                hid_t tGroup = aFile.select_group("Matrix");
                mJacobian->load( tGroup, tError );
                aFile.close_active_group() ;
                if( mRhsVector.length() > 0 )
                {
                    aFile.load_data( "LHS", mLhsVector );
                    aFile.load_data( "RHS", mRhsVector );
                }
                else
                {
                    aFile.load_data( "LHS", mLhsMatrix );
                    aFile.load_data( "RHS", mRhsMatrix );
                }
                if( mConvection.length() > 0 )
                {
                    aFile.load_data( "SurfaceLoads", mConvection );
                }
                if( mVolumeLoads.length() > 0 )
                {
                    aFile.load_data( "VolumeLoads", mVolumeLoads );
                }
            }

            void
            SolverData::save_system( HDF5 & aFile )
            {
                herr_t tError = 0 ;

                hid_t tGroup = aFile.create_group("Matrix");
                mJacobian->save( tGroup, tError );
                aFile.close_active_group() ;
                if( mRhsVector.length() > 0 )
                {
                    aFile.save_data( "LHS", mLhsVector );
                    aFile.save_data( "RHS", mRhsVector );
                }
                else
                {
                    aFile.save_data( "LHS", mLhsMatrix );
                    aFile.save_data( "RHS", mRhsMatrix );
                }
                if( mConvection.length() > 0 )
                {
                    aFile.save_data( "SurfaceLoads", mConvection );
                }
                if( mVolumeLoads.length() > 0 )
                {
                    aFile.save_data( "VolumeLoads", mVolumeLoads );
                }

            }

#endif

//------------------------------------------------------------------------------
        } /* end namespace dofmgr */
    } /* end namespace fem */
} /* end namespace belfem */