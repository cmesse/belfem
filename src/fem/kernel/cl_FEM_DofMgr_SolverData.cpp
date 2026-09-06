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
#include <iterator>
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
#include "fn_FEM_anderson_mixing.hpp"
#include "cl_Timer.hpp"

#include "fn_entity_type.hpp"
#include "fn_matrix_type.hpp"
//#include "op_Graph_Vertex_Index.hpp"
#include "petsctools.hpp"

namespace belfem
{
    namespace fem
    {
        namespace dofmgr
        {
            string
            to_string( const MatrixType aType )
            {
                switch ( aType )
                {
                    case System : return "System" ;
                    case Jacobian : return "Jacobian" ;
                    case Dirichlet : return "Dirichlet" ;
                    case Enforcement : return "Enforcement" ;
                    case Imposition : return "Imposition" ;
                    case FullMass : return "FullMass" ;
                    case FullStiffness : return "FullStiffness" ;
                    default : return "unknown" ;
                }
            }

//------------------------------------------------------------------------------

            SolverData::SolverData( DofManager * aParent, DofData * aDofData ,
                                    BlockData * aBlockData,
                                    SideSetData * aSideSetData ) :
                    mParent( aParent ),
                    mKernel( aParent->parent() ),
                    mDofData( aDofData ),
                    mBlockData( aBlockData ),
                    mSideSetData( aSideSetData ),
                    mCommRank( comm_rank() ),
                    mCommSize( comm_size() ),
                    mDOFs( aDofData->dofs() ),
                    mNumberOfFreeDofs( aDofData->number_of_free_dofs() ),
                    mNumberOfFixedDofs( aDofData->number_of_fixed_dofs() ),
                    //mNumberOfHangingDofs( aDofData->number_of_hanging_dofs() ),
                    mMyNumberOfFreeDofs( aDofData->my_number_of_free_dofs() ),
                    mMyNumberOfFixedDofs( aDofData->my_number_of_fixed_dofs() ),
                    mMyNumberOfHangingDofs( aDofData->my_number_of_hanging_dofs() )
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

                // delete the Anderson history registers if they exist
                if ( mAndersonX != nullptr )
                {
                    delete mAndersonX ;
                    delete mAndersonR ;
                }
            }

//------------------------------------------------------------------------------

            void SolverData::reset()
            {

                if( mSystemMatrix != nullptr )
                {
                    delete mSystemMatrix ;
                    mSystemMatrix = nullptr ;
                }
                if( mEnforcementMatrix != nullptr )
                {
                    delete mEnforcementMatrix ;
                    mEnforcementMatrix = nullptr ;
                }
                if( mDirichletMatrix != nullptr )
                {
                    delete mDirichletMatrix ;
                    mDirichletMatrix = nullptr ;
                }
                if ( mImpositionMatrix != nullptr )
                {
                    delete mImpositionMatrix ;
                    mImpositionMatrix = nullptr ;
                }
                if ( mFullMassMatrix != nullptr )
                {
                    delete mFullMassMatrix ;
                    mFullMassMatrix = nullptr ;
                }
                if ( mFullStiffnessMatrix != nullptr )
                {
                    delete mFullStiffnessMatrix ;
                    mFullStiffnessMatrix = nullptr ;
                }
                if ( mJacobianMatrix != nullptr )
                {
                    delete mJacobianMatrix ;
                    mJacobianMatrix = nullptr ;
                }

                mSystemTable.clear() ;
                mJacobianTable.clear() ;
                mEnforcementTable.clear() ;
                mDirichletTable.clear() ;
                mImpositionTable.clear() ;

                mFullMatrixTable.clear() ;

                this->reset_convection() ;
            }

//------------------------------------------------------------------------------

            void
            SolverData::reset_convection()
            {
                mConvection.fill( 0.0 );
            }

//------------------------------------------------------------------------------

            void
            SolverData::extract_graph_from_mesh( Vector< id_t > & aGraphData )
            {
                // restore factory settings
                this->reset() ;

                // local element-to-dof adjacency
                Vector< id_t > tElementWiseData ;
                this->compute_element_dof_connectivity( tElementWiseData );

                // local dof-to-element adjacency
                Vector< id_t > tDofWiseData ;
                this->compute_dof_element_connectivity( tDofWiseData );

                if( mKernel->number_of_procs() > 1 )
                {
                    Cell< Vector< id_t > > tConnectivities( mKernel->number_of_procs(),
                                                            Vector< id_t >());
                    Vector< id_t > & tConnectivity = mKernel->is_master() ? tConnectivities( 0 ) : aGraphData;
                    this->compute_dof_dof_connectivity( tDofWiseData, tElementWiseData, tConnectivity );

                    comm_barrier();

                    if ( mKernel->is_master() )
                    {
                        collect( tConnectivities );

                        this->unite_dofs( tConnectivities, aGraphData );
                    }
                    else
                    {
                        send( aGraphData );
                    }
                }
                else
                {
                    this->compute_dof_dof_connectivity( tDofWiseData, tElementWiseData, aGraphData );
                }

            }

//------------------------------------------------------------------------------

            void
            SolverData::use_jedi_force( const bool aFlag  )
            {
                mUseJediForce = aFlag ;
            }

//------------------------------------------------------------------------------

            void
            SolverData::use_full_force( const bool aFlag )
            {
                mUseFullForce = aFlag ;
            }

//------------------------------------------------------------------------------

            /**
             * Allocates sparse matrices for the finite element system.
             *
             * Purpose:
             *   - Creates sparse matrix structures with optimized sparsity patterns
             *   - Supports both JEDI Force matrices (J, D, E, I) and Full matrices (M, K)
             *   - Allocates RHS vectors/matrices for solution
             *
             * Matrix Types:
             *   JEDI Force:
             *     - Jacobian (J): n×n matrix (free-free interactions)
             *     - Dirichlet (D): n×m matrix (free-fixed interactions, negative sign in assembly)
             *     - Enforcement (E): m×n matrix (fixed-free, computes reaction forces) [if mUseJediForce]
             *     - Imposition (I): m×m matrix (fixed-fixed self-coupling) [if mUseJediForce]
             *
             *   Full Force (if mUseFullForce):
             *     - FullMass (M): (n+m)×(n+m) matrix (all DOFs, eigenvalue analysis)
             *     - FullStiffness (K): (n+m)×(n+m) matrix (all DOFs, eigenvalue analysis)
             *
             * Workflow:
             *   1. If mUseFullForce: allocate FullMass / FullStiffness in the combined index space
             *   2. Allocate the Dirichlet matrix (always)
             *   3. If mUseJediForce: allocate Enforcement and Imposition
             *   4. Allocate the System matrix from the System graph (must be last), then
             *      the Jacobian as its child (shared pointers and indices)
             *   5. Allocate RHS vector or matrix
             *
             * Index Space Handling for Full Matrices:
             *   - Free DOFs normally use indices 0..n-1 (SEPARATE mode)
             *   - Fixed DOFs normally use indices 0..m-1 (SEPARATE mode)
             *   - Full matrices require combined indexing: free=0..n-1, fixed=n..n+m-1
             *   - Temporarily shift fixed DOF indices, merge graphs, create matrices, restore
             *
             * Graph Requirements:
             *   - aFreeDofs: Contains free DOF vertices with reordered indices 0..n-1
             *   - aFixedDofs: Contains fixed DOF vertices with reordered indices 0..m-1
             *   - For matrix rows: Graph passed to SpMatrix must match row DOF type
             *     • Jacobian (J) → aFreeDofs (free rows × free cols)
             *     • Dirichlet (D) → aFreeDofs (free rows × fixed cols)
             *     • Enforcement (E) → aFixedDofs (fixed rows × free cols)
             *     • Imposition (I) → aFixedDofs (fixed rows × fixed cols)
             *
             * Preconditions:
             *   - DOFs must be reordered (via reorder_dofs())
             *   - aGraphData contains complete DOF-to-DOF connectivity
             *   - Solver must be initialized
             *
             * Postconditions:
             *   - All required sparse matrices allocated with correct sparsity patterns
             *   - RHS vectors/matrices allocated
             */
            void
            SolverData::allocate_matrices(
                const Vector< id_t > & aGraphData ,
                Graph & aFreeDofs,
                Graph & aFixedDofs )
            {
                // full-force matrices are opt-in
                if ( mUseFullForce )
                {
                    this->populate_graph( aGraphData, MatrixType::FullMass, aFreeDofs, aFixedDofs );

                    // Temporarily shift fixed DOF indices to combined space
                    for ( graph::Vertex * tVertex : aFixedDofs )
                    {
                        Dof * tDof = reinterpret_cast< Dof * >( tVertex );
                        tDof->set_index( tDof->index() + mNumberOfFreeDofs );
                        tDof->set_my_index( tDof->my_index() + mMyNumberOfFreeDofs );
                    }

                    // Merge graphs into combined graph for full matrices
                    Graph tGraph( mMyNumberOfFreeDofs + mMyNumberOfFixedDofs, nullptr );

                    index_t tCount = 0 ;
                    for ( graph::Vertex * tVertex : aFreeDofs )
                    {
                        tGraph( tCount++ ) = tVertex;
                    }
                    for ( graph::Vertex * tVertex : aFixedDofs )
                    {
                        tGraph( tCount++ ) = tVertex;
                    }

                    index_t n = mNumberOfFreeDofs + mNumberOfFixedDofs ;

                    // Create full mass and stiffness matrices
                    mFullMassMatrix = new SpMatrix( tGraph, SpMatrixType::CSR,
                                                     n, n, false );

                    // matrix shares pointers and indices with full mass matrix
                    mFullStiffnessMatrix = new SpMatrix( mFullMassMatrix );

                    tGraph.clear();

                    // Shift fixed DOF indices back to SEPARATE mode
                    tCount = 0 ;
                    for ( graph::Vertex * tVertex : aFixedDofs )
                    {
                        Dof * tDof = reinterpret_cast< Dof * >( tVertex );
                        tDof->set_index( tDof->index() - mNumberOfFreeDofs );
                        tDof->set_my_index( tDof->my_index() - mMyNumberOfFreeDofs );
                    }
#ifdef DEBUG
                    // Verify all DOF indices are in valid ranges after shift-back
                    for ( graph::Vertex * tVertex : aFreeDofs )
                    {
                        Dof * tDof = reinterpret_cast< Dof * >( tVertex );
                        BELFEM_ASSERT( reinterpret_cast< Dof * >( tVertex )->my_index() < mMyNumberOfFreeDofs,
                            "Free DOF %lu has invalid my_index %lu (should be < %lu) after FullForce allocation on rank %u",
                            (long unsigned int)tDof->id(),
                            (long unsigned int)tDof->my_index(),
                            (long unsigned int)mMyNumberOfFreeDofs,
                            (unsigned int)mCommRank );
                    }

                    for ( graph::Vertex * tVertex : aFixedDofs )
                    {
                        Dof * tDof = reinterpret_cast< Dof * >( tVertex );
                        BELFEM_ASSERT( tDof->my_index() < mMyNumberOfFixedDofs,
                            "Fixed DOF %lu has invalid my_index %lu (should be < %lu) after shift-back on rank %u",
                            (long unsigned int)tDof->id(),
                            (long unsigned int)tDof->my_index(),
                            (long unsigned int)mMyNumberOfFixedDofs,
                            (unsigned int)mCommRank );
                    }

#endif
                }

                this->populate_graph( aGraphData, MatrixType::Dirichlet, aFreeDofs, aFixedDofs , false);

                mDirichletMatrix = new SpMatrix(
                    aFreeDofs,
                    SpMatrixType::CSR,
                    mNumberOfFreeDofs,
                    mNumberOfFixedDofs,
                    false);

                if ( mUseJediForce )
                {
                    this->populate_graph( aGraphData,
                        Enforcement,
                        aFreeDofs,
                        aFixedDofs );

                    mEnforcementMatrix = new SpMatrix(
                        aFixedDofs, SpMatrixType::CSR,
                          mNumberOfFixedDofs,
                          mNumberOfFreeDofs,
                          false );

                    this->populate_graph(
                        aGraphData,
                        MatrixType::Imposition,
                        aFreeDofs,
                        aFixedDofs );

                    mImpositionMatrix = new SpMatrix(
                        aFixedDofs,
                        SpMatrixType::CSR,
                       mNumberOfFixedDofs,
                       mNumberOfFixedDofs,
                       false );
                }

                // J-Matrix must be created last because the dof-dof connectivities remain
                this->populate_graph(
                    aGraphData,
                    System,
                    aFreeDofs,
                    aFixedDofs );

                BELFEM_ASSERT( mSolver != nullptr, "no solver created" );

                mSystemMatrix =  new SpMatrix(
                    aFreeDofs,
                    matrix_type( mSolver->type() ),
                    mNumberOfFreeDofs,
                    mNumberOfFreeDofs,
                    false );

                // jacobian matrix inherits pointers and indices from
                // system matrix
                mJacobianMatrix =  new SpMatrix( mSystemMatrix );

                // - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // allocate RHS
                // - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // allocate right hand side
                if( mParent->iwg()->num_rhs_cols() <= 1 )
                {
                    mRhsVector.set_size( mMyNumberOfFreeDofs, 0.0 );
                    if( mParent->iwg()->has_convection() )
                    {
                        mConvection.set_size( mMyNumberOfFreeDofs, 0.0 );
                    }
                }
                else
                {
                    mRhsMatrix.set_size( mMyNumberOfFreeDofs,
                                         mParent->iwg()->num_rhs_cols(), 0.0 );
                }

                // ( a per-proc counter exchange lived here until 2026-08-15;
                //   the collected tables were never read — removed on both
                //   the master and the worker side )
            }

 //------------------------------------------------------------------------------

            void
            SolverData::compute_memory()
            {
                const uint tNumSlots = 12 ;
                Vector< size_t > tMyMemory( tNumSlots, 0 );

                // slot 0: mesh
                tMyMemory( 0 ) = mParent->mesh()->memory();

                // slot 1: dof graph
                size_t tCount = 0 ;
                for ( Dof * tDof : mDOFs )
                {
                    tCount += tDof->memory();
                }
                tMyMemory( 1 ) = tCount ;

                // slot 2: dense vectors and matrices
                tCount  = mLhsVector.length()   * sizeof( real );
                tCount += mLhsMatrix.n_rows()   * mLhsMatrix.n_cols() * sizeof( real );
                tCount += mRhsVector.length()   * sizeof( real );
                tCount += mRhsMatrix.n_rows()   * mRhsMatrix.n_cols() * sizeof( real );
                tCount += mFieldValues.length() * sizeof( real );
                tCount += mConvection.length()  * sizeof( real );
                tCount += mVolumeLoads.length() * sizeof( real );
                tMyMemory( 2 ) = tCount ;

                // slot 3: assembly tables
                tCount = 0 ;
                for ( uint i = 0; i < mSystemTable.size(); ++i )
                    tCount += mSystemTable( i ).length() * sizeof( index_t );
                for ( uint i = 0; i < mJacobianTable.size(); ++i )
                    tCount += mJacobianTable( i ).length() * sizeof( index_t );
                for ( uint i = 0; i < mEnforcementTable.size(); ++i )
                    tCount += mEnforcementTable( i ).length() * sizeof( index_t );
                for ( uint i = 0; i < mDirichletTable.size(); ++i )
                    tCount += mDirichletTable( i ).length() * sizeof( index_t );
                for ( uint i = 0; i < mImpositionTable.size(); ++i )
                    tCount += mImpositionTable( i ).length() * sizeof( index_t );
                for ( uint i = 0; i < mFullMatrixTable.size(); ++i )
                    tCount += mFullMatrixTable( i ).length() * sizeof( index_t );
                tMyMemory( 3 ) = tCount ;

                // slot 4: reset buffers
                tCount  = mRhsVector0.length()         * sizeof( real );
                tCount += mRhsMatrix0.n_rows()         * mRhsMatrix0.n_cols() * sizeof( real );
                tCount += mJacobianValues0.length()    * sizeof( real );
                tCount += mDirichletValues0.length()   * sizeof( real );
                tCount += mEnforcementValues0.length() * sizeof( real );
                tCount += mImpositionValues0.length()  * sizeof( real );
                tMyMemory( 4 ) = tCount ;

                // slots 5-11: sparse matrices
                if ( mSystemMatrix != nullptr )
                {
                    tMyMemory( 5 ) = mSystemMatrix->memory();
                }
                if ( mJacobianMatrix != nullptr )
                {
                    tMyMemory( 6 ) = mJacobianMatrix->memory();
                }
                if ( mEnforcementMatrix != nullptr )
                {
                    tMyMemory( 7 ) = mEnforcementMatrix->memory();
                }
                if ( mDirichletMatrix != nullptr )
                {
                    tMyMemory( 8 ) = mDirichletMatrix->memory();
                }
                if ( mImpositionMatrix != nullptr )
                {
                    tMyMemory( 9 ) = mImpositionMatrix->memory();
                }
                if ( mFullMassMatrix != nullptr )
                {
                    tMyMemory( 10 ) = mFullMassMatrix->memory();
                }
                if ( mFullStiffnessMatrix != nullptr )
                {
                    tMyMemory( 11 ) = mFullStiffnessMatrix->memory();
                }
                comm_barrier();

                if ( mCommRank == 0 )
                {
                    Cell< Vector< size_t > > tAllMemory( mCommSize, {} );
                    collect( tAllMemory );

                    for ( proc_t p = 1; p < mCommSize; ++p )
                    {
                        tMyMemory += tAllMemory( p );
                    }

                    Cell< string > tType = { "Mesh           ",
                                             "Dof Graph      ",
                                             "Vectors        ",
                                             "Assembly Tables",
                                             "Reset Buffers  ",
                                             "System         ",
                                             "Jacobian       ",
                                             "Enforcement    ",
                                             "Dirichlet      ",
                                             "Imposition     ",
                                             "FullMass       ",
                                             "FullStiffness  "};


                    // name the dof fields so the report says WHICH system it
                    // belongs to — two kernels print this block back to back,
                    // and unlabeled counts from different runs have been
                    // compared as if they were the same quantity (2026-08-15)
                    const Cell< string > & tDofFields = mParent->iwg()->dof_fields() ;
                    string tFieldList ;
                    for ( uint f = 0; f < tDofFields.size() && f < 4; ++f )
                    {
                        if ( f > 0 ) tFieldList += ", " ;
                        tFieldList += tDofFields( f );
                    }
                    if ( tDofFields.size() > 4 ) tFieldList += ", ..." ;

                    message( InfoLevel::Default, "\n    Number of Degrees of Freedom ( %s ): ",
                        tFieldList.c_str() );
                    message( InfoLevel::Default, "            * unknowns   (free)    : %lu",
                        ( long unsigned int ) mNumberOfFreeDofs );
                    message( InfoLevel::Default, "            * prescribed (fixed)   : %lu",
                        ( long unsigned int ) mNumberOfFixedDofs );
                    message( InfoLevel::Default, "            * condensed  (hanging) : %lu",
                        ( long unsigned int ) mDofData->number_of_hanging_dofs() );

                    for ( uint k = 0; k < tNumSlots; ++k )
                    {
                        if ( k == 0 )
                        {
                            message( InfoLevel::Default, "\n    General Memory Needs: " );
                        }
                        if ( k == 5 )
                        {
                            message( InfoLevel::Default, "\n    Matrix Memory Needs: " );
                        }
                        if ( tMyMemory( k ) == 0 ) continue;

                        string tUnit ;
                        long unsigned int tMem ;
                        if ( tMyMemory( k ) < 1048576 )
                        {
                            tMem = tMyMemory( k ) / 1024 ;
                            tUnit = "KiB" ;
                        }
                        else if ( tMyMemory( k ) < 1073741824 )
                        {
                            tMem = tMyMemory( k ) / 1048576 ;
                            tUnit = "MiB" ;
                        }
                        else
                        {
                            tMem = tMyMemory( k ) / 1073741824 ;
                            tUnit = "GiB" ;
                        }

                        message( InfoLevel::Default, "            * %s : %lu %s ", tType( k ).c_str(), tMem, tUnit.c_str() );
                    }
                    message( InfoLevel::Default, "\n" );
                }
                else
                {
                    send( tMyMemory );
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::create_assembly_tables()
            {
                const uint tNumMatrices = 6 ;
                if ( mCommRank == 0 )
                {

                    comm_barrier() ;
                    Cell< Cell< Vector< int_t > > > tAllIndices( 2*tNumMatrices, {} ) ;

                    for ( uint k=0; k<2*tNumMatrices; ++k )
                    {
                        comm_barrier() ;
                        Cell< Vector< int_t > > & tTable = tAllIndices( k );
                        tTable.set_size( mCommSize, {} );

                        collect( tTable );
                    }

                    index_t tCount = 0 ;

                    for ( uint m=0; m<tNumMatrices; ++m )
                    {
                        // get the matrix
                        SpMatrix * tMatrix = this->matrix( static_cast< MatrixType >( m ) );

                        // Always consume the rows/cols indices, even if matrix doesn't exist
                        Cell< Vector< int_t > > & tAllRows = tAllIndices( tCount++ );
                        Cell< Vector< int_t > > & tAllCols = tAllIndices( tCount++ );

                        if ( tMatrix == nullptr ) continue ;
                        tMatrix->set_indexing_base( SpMatrixIndexingBase::Cpp );

                        Cell< Vector< index_t > > & tTables = this->tables( static_cast< MatrixType >( m ) );
                        tTables.set_size( mCommSize, {} );

                        for ( proc_t p = 1; p < mCommSize; ++p )
                        {
                            // get the local table
                            Vector< index_t > & tTable = tTables( p );

                            // get rows
                            Vector< int_t > & tRows = tAllRows( p );

                            // get cols
                            Vector< int_t > & tCols = tAllCols( p );

                            // get number of nonzeros in this matrix
                            index_t tNNZ = tRows.length();

                            tTable.set_size( tNNZ );

                            // loop over all entries
                            for ( index_t k = 0; k < tNNZ; ++k )
                            {
                                // compute index
                                tTable( k ) = tMatrix->index( tRows( k ), tCols( k ) );
                            }
                        }
                    }
                }
                else
                {
                    for ( uint m=0; m<tNumMatrices; ++m )
                    {
                        SpMatrix * tMatrix = this->matrix( static_cast< MatrixType >( m ) );
                        if ( tMatrix != nullptr )
                        {
                            tMatrix->set_indexing_base( SpMatrixIndexingBase::Cpp );
                            tMatrix->create_coo_indices();
                        }
                    }
                    comm_barrier() ;
                    for ( uint m=0; m<tNumMatrices; ++m )
                    {
                        SpMatrix * tMatrix = this->matrix( static_cast< MatrixType >( m ) );


                        if ( tMatrix == nullptr )
                        {
                            int_t tZero = 0;

                            // send zero to master ( for rows )
                            comm_barrier() ;
                            send( tZero );

                            // send zero to master ( for cols )
                            comm_barrier() ;
                            send( tZero );
                        }
                        else
                        {
                            comm_barrier() ;
                            send( tMatrix->rows(), tMatrix->number_of_nonzeros(), 0 );
                            comm_barrier() ;
                            send( tMatrix->cols(), tMatrix->number_of_nonzeros(), 0 );
                        }
                    }
                    for ( uint m=0; m<tNumMatrices; ++m )
                    {
                        SpMatrix * tMatrix = this->matrix( static_cast< MatrixType >( m ) );
                        if ( tMatrix != nullptr )
                        {
                            tMatrix->free_coo_indices();
                        }
                    }
                }
            }

//------------------------------------------------------------------------------

            /**
             * Builds vertex connectivity for sparse matrix sparsity patterns.
             *
             * Purpose:
             *   - Creates graph connectivity by linking DOF vertices based on element topology
             *   - Filters connections by DOF type (free/fixed) to match specific matrix structure
             *   - Supports both JEDI Force matrices (J, D, E, I) and Full matrices (M, K)
             *   - Optionally excludes self-connections for graph reordering algorithms
             *
             * Matrix Type Mapping (rows × columns):
             *   - Jacobian (J):     free × free   (standard stiffness, used for solving)
             *   - Dirichlet (D):    free × fixed  (BC enforcement, negative sign in assembly)
             *   - Enforcement (E):  fixed × free  (computes reaction forces at constraints)
             *   - Imposition (I):   fixed × fixed (self-coupling of prescribed DOFs)
             *   - FullMass (M):     all × all     (eigenvalue analysis, combined index space)
             *   - FullStiffness (K):all × all     (eigenvalue analysis, combined index space)
             *
             * Algorithm:
             *   1. Determine row/column DOF types from matrix type
             *   2. For each DOF in connectivity data:
             *      a. Check if DOF should be a row vertex (matches tRowIsFixed)
             *      b. Loop over all connected DOFs from elements
             *      c. Filter connections by column type (matches tColIsFixed)
             *      d. Mark connections in bitset using my_index + offset
             *      e. [Optional] Remove self-connection if aLinkToSelf=false
             *      f. Convert bitset to index list
             *      g. Insert vertex pointers into DOF's connectivity list
             *
             * Index Space Handling:
             *   - Free DOFs: my_index = 0..n-1, offset = 0
             *   - Fixed DOFs: my_index = 0..m-1, offset = tMyNumberOfFreeDofs
             *   - Bitset size: (n+m) to accommodate both types with offset
             *   - For full matrices: uses combined space, fixed DOFs already offset
             *
             * Bitset Mechanism:
             *   - Efficiently tracks unique connections (prevents duplicates)
             *   - Position i in bitset represents DOF with my_index i (or i-offset for fixed)
             *   - bitset.where() extracts set positions into index list
             *
             * Graph Vertex Container:
             *   - Each Dof inherits from graph::Vertex
             *   - Vertex has container of connected vertices (adjacency list)
             *   - SpMatrix constructor reads this connectivity to build sparsity pattern
             *
             * Parameters:
             *   @param aData         DOF-to-DOF connectivity from extract_graph_from_mesh()
             *                        Format: [num_dofs, dof_id, num_connected, dof_id, ...]
             *   @param aMatrixType   Which matrix to populate (J, D, E, I, FullMass, FullStiffness)
             *   @param aFreeDofs     Graph vertices for free DOFs (indices 0..n-1)
             *   @param aFixedDofs    Graph vertices for fixed DOFs (indices 0..m-1)
             *   @param aLinkToSelf   Include diagonal connections (default true)
             *                        - true: Include self-connections (for matrix allocation)
             *                        - false: Exclude self-connections (for METIS/SCOTCH reordering)
             *
             * Preconditions:
             *   - aData contains complete DOF connectivity from extract_graph_from_mesh()
             *   - aFreeDofs and aFixedDofs contain reordered DOF vertices
             *   - Each DOF's my_index is set to its position in respective graph (0..n-1 or 0..m-1)
             *
             * Postconditions:
             *   - Each row DOF has vertex container populated with column DOF vertices
             *   - Connectivity matches requested matrix type structure
             *   - Ready for SpMatrix construction
             *
             * Usage Example:
             *   // For Jacobian (free × free):
             *   populate_graph(graphData, Jacobian, freeDofs, fixedDofs, true);
             *   // Result: Each free DOF points to all connected free DOFs
             *
             *   // For Dirichlet (free × fixed):
             *   populate_graph(graphData, Dirichlet, freeDofs, fixedDofs, false);
             *   // Result: Each free DOF points to all connected fixed DOFs, no self-links
             */
            void
            SolverData::populate_graph( const Vector< id_t > & aData,
                                        const MatrixType aMatrixType,
                                        Graph                & aFreeDofs,
                                        Graph                & aFixedDofs,
                                        const bool aLinkToSelf )
            {

                bool tRowIsFixed = ( aMatrixType == Enforcement || aMatrixType == Imposition );
                bool tColIsFixed = ( aMatrixType == Dirichlet || aMatrixType == Imposition );

                bool tUseFullMatrix =  aMatrixType == FullMass || aMatrixType == FullStiffness ;

                index_t tPivot = 0 ;

                index_t tMyNumberOfFreeDofs = aFreeDofs.size() ;

                // total number of dofs
                index_t tNumberOfDofs = aData( tPivot++ );

                index_t tNumDofsPerDof ;

                DynamicBitset tBitset( aFreeDofs.size() + aFixedDofs.size() );
                Cell< index_t > tIndices ;

                for( index_t k=0; k<tNumberOfDofs; ++k )
                {
                    // grab dof
                    Dof * tA = mDofData->dof( aData( tPivot++ ) );
                    tA->reset_vertex_container();
                    tNumDofsPerDof = aData( tPivot++ );

                    BELFEM_ASSERT( ! tA->is_hanging(), "dofs should not be hanging in graph creation" );

                    // check if dof is taken into account
                    if ( tA->is_fixed() == tRowIsFixed || tUseFullMatrix )
                    {
                        tBitset.reset();

                        for( index_t i=0; i<tNumDofsPerDof; ++i )
                        {
                            // get other dof
                            Dof * tB = mDofData->dof( aData( tPivot++ ) );

                            // check flag
                            index_t tOff = tB->is_fixed() ? tMyNumberOfFreeDofs : 0 ;

                            if ( tB->is_fixed() == tColIsFixed || tUseFullMatrix )
                            {
                                tBitset.set( tB->my_index() + tOff );
                            }
                        }

                        // for metis or scotch, a dof must not link to itself!
                        index_t tOff = tA->is_fixed() ? tMyNumberOfFreeDofs : 0 ;
                        if ( ! aLinkToSelf ) tBitset.reset( tA->my_index() + tOff );

                        tBitset.where( tIndices );

                        tA->reset_vertex_container() ;
                        tA->init_vertex_container( tIndices.size() );

                        for ( index_t i : tIndices )
                        {
                            if ( i < tMyNumberOfFreeDofs )
                            {
                                tA->insert_vertex(aFreeDofs(i));
                            }
                            else
                            {
                                tA->insert_vertex(aFixedDofs(i-tMyNumberOfFreeDofs));
                            }
                        }

                        BELFEM_ASSERT( tA->number_of_vertices() == tIndices.size(),
                            "Expected %lu vertices but got %u (MatrixType=%d)",
                            (long unsigned int)tIndices.size(), tA->number_of_vertices(), (int)aMatrixType );
                    }
                    else
                    {
                        tPivot += tNumDofsPerDof ;
                    }
                }
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

                // counter for number of elements
                index_t tElemCount = 0 ;

                // determine memory size
                index_t tCount = 1 ;

                for( Block * tBlock : mBlockData->blocks() )
                {
                    tElemCount += tBlock->number_of_elements();

                    for( Element * tElement : tBlock->elements() )
                    {
                        tCount += tElement->number_of_dofs() + 2 ;
                    }
                }

                for( id_t tID : mParent->iwg()->selected_sidesets() )
                {
                    SideSet * tSideSet = mSideSetData->sideset( tID );

                    if( ! tSideSet->is_active() )
                    {
                        continue ;
                    }

                    tElemCount += tSideSet->number_of_elements();
                    for( Element * tElement : tSideSet->elements() )
                    {
                        tCount += tElement->number_of_dofs() + 2 ;
                    }
                }

                aData.set_size( tCount );

                // reset the counter
                tCount = 0 ;

                aData( tCount++ ) = tElemCount ;

                for( Block * tBlock : mBlockData->blocks() )
                {
#ifdef DEBUG
                    uint tNumDofsPerElement
                            = mParent->iwg()->number_of_dofs_per_element( tBlock );
#endif
                    for ( Element * tElement: tBlock->elements())
                    {
                        aData( tCount++ ) = tElement->id();
                        aData( tCount++ ) = tElement->number_of_dofs() ;
#ifdef DEBUG
                        BELFEM_ASSERT( tElement->number_of_local_dofs() == tNumDofsPerElement ,
                                      "number of local dofs do not match for element %lu on block %lu  is %lu but expect %lu",
                                      ( long unsigned int ) tElement->id(),
                                      ( long unsigned int ) tBlock->id(),
                                      ( long unsigned int ) tElement->number_of_local_dofs(),
                                      ( long unsigned int ) tNumDofsPerElement );
#endif
                        for ( uint d=0; d<tElement->number_of_dofs(); ++d )
                        {
                            aData( tCount++ ) = tElement->dof( d )->id() ;
                        }
                    }
                }

                for( id_t tID : mParent->iwg()->selected_sidesets() )
                {
                    SideSet * tSideSet = mSideSetData->sideset( tID );

                    if( ! tSideSet->is_active() )
                    {
                        continue ;
                    }
#ifdef DEBUG
                    uint tNumDofsPerElement
                            = mParent->iwg()->number_of_dofs_per_element( tSideSet );
#endif
                    for ( Element * tElement: tSideSet->elements() )
                    {

                        aData( tCount++ ) = tElement->id();
                        aData( tCount++ ) = tElement->number_of_dofs();
#ifdef DEBUG
                        BELFEM_ASSERT( tElement->number_of_local_dofs() == tNumDofsPerElement,
                                      "number of dofs do not match for element %lu ( master: %lu ) on sideset %lu  is %u but expect %u",
                                      ( long unsigned int ) tElement->id(),
                                      ( long unsigned int ) tElement->master()->id(),
                                      ( long unsigned int ) tSideSet->id(),
                                      ( unsigned int ) tElement->number_of_dofs(),
                                      ( unsigned int ) tNumDofsPerElement );
#endif
                        for ( uint d=0; d<tElement->number_of_dofs(); ++d )
                        {
                            aData( tCount++ ) = tElement->dof( d )->id();
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
                            /*// debug output
                            if ( tElement->dof( d )->is_hanging() )
                            {
                                std::cout << "#ELEMENT " << tElement->id() << " " << to_string( tElement->element()->type() ) << std::endl;

                                for ( uint i=0; i<tElement->number_of_dofs(); ++i )
                                {
                                    std::cout << "  #GDOF : " << tElement->dof( i )->id() << " " << tElement->dof( i )->is_hanging() << std::endl;
                                }
                                std::cout << std::endl;
                                for ( uint i=0; i<tElement->number_of_local_dofs(); ++i )
                                {
                                    std::cout << "  #LDOF : " << tElement->local_dof( i )->id() << " " << tElement->local_dof( i )->is_hanging() << std::endl;
                                }
                            }*/
                            BELFEM_ASSERT( ! tElement->dof( d )->is_hanging(), "Element %lu on block %lu has a hanging dof that is denoted as real dof: %lu (%s %lu)" ,
                                ( long unsigned int ) tElement->id(),
                                ( long unsigned int ) tElement->element()->block_id(),
                                ( long unsigned int ) tElement->dof( d )->id(),
                                to_string( tElement->dof( d )->entity_type() ).c_str(),
                                ( long unsigned int ) tElement->dof( d )->mesh_basis()->id() );

                            ++tWork( tElement->dof( d )->my_index() ) ;
                        }
                    }
                }


                for( id_t tID : mParent->iwg()->selected_sidesets() )
                {
                    SideSet * tSideSet = mSideSetData->sideset( tID );

                    if( ! tSideSet->is_active() )
                    {
                        continue ;
                    }

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

                    if( ! tSideSet->is_active() )
                    {
                        continue ;
                    }

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
                        tOffset = tOffsets( tElementID );

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
                    Vector< index_t > & tIDs = tAllIDs( k );

                    if( tIDs.length() > 0 )
                    {
                        j = 0;

                        // loop over all elements of this dof
                        for ( index_t e = 0; e < tNumElems; ++e )
                        {
                            // get the element ID
                            tElementID = aDofWiseData( tPivot++ );

                            // get the offset in the other vector
                            tOffset = tOffsets( tElementID );

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
                        aConnectivity( tPivot++ ) = tIDs( i );
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
                index_t tNumDofs = 0 ;

                // set temporary dof index
                for ( Dof * tDof : mDOFs )
                {
                    tDof->set_my_index( tNumDofs++ );
                }


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
                        Dof * tDof = mDofData->dof( tData( tPivot++ ) );

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
            SolverData::reset_matrices( const bool aFullForce )
            {
                BELFEM_ASSERT( mSystemMatrix != nullptr,
                              "Jacobian Matrix was not initialized");

                if ( aFullForce )
                {
                    if ( mFullMassMatrix != nullptr )
                    {
                        mFullMassMatrix->fill( 0.0 );
                    }
                    if ( mFullStiffnessMatrix != nullptr )
                    {
                        mFullStiffnessMatrix->fill( 0.0 );
                    }
                    return;
                }

                if( mUseResetValues && mParent->parent()->is_master() )
                {
                    std::copy( mSystemValues0.data(),
                               mSystemValues0.data()+mSystemMatrix->number_of_nonzeros(),
                               mSystemMatrix->data() );

                    std::copy( mJacobianValues0.data(),
                               mJacobianValues0.data()+mJacobianMatrix->number_of_nonzeros(),
                               mJacobianMatrix->data() );

                    if ( mEnforcementMatrix != nullptr )
                    {
                        std::copy( mEnforcementValues0.data(),
                                   mEnforcementValues0.data()+mEnforcementMatrix->number_of_nonzeros(),
                                   mEnforcementMatrix->data() );
                    }
                    if( mDirichletMatrix != nullptr )
                    {
                        std::copy( mDirichletValues0.data(),
                                   mDirichletValues0.data()+mDirichletMatrix->number_of_nonzeros(),
                                   mDirichletMatrix->data() );
                    }
                    if ( mImpositionMatrix != nullptr )
                    {
                        std::copy( mImpositionValues0.data(),
                                   mImpositionValues0.data()+mImpositionMatrix->number_of_nonzeros(),
                                   mImpositionMatrix->data() );
                    }
                }
                else
                {
                    mSystemMatrix->fill( 0.0 );
                    mJacobianMatrix->fill( 0.0 );

                    if ( mEnforcementMatrix != nullptr )
                    {
                        mEnforcementMatrix->fill( 0.0 );
                    }
                    if ( mDirichletMatrix != nullptr )
                    {
                        mDirichletMatrix->fill( 0.0 );
                    }

                    if ( mImpositionMatrix != nullptr )
                    {
                        mImpositionMatrix->fill( 0.0 );
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::reset_rhs_vector()
            {
                if( mUseResetValues && mParent->parent()->is_master() &&
                    mRhsVector0.length() > 0 )
                {
                    mRhsVector = mRhsVector0 ;
                }
                else
                {
                    mRhsVector.fill( 0.0 );
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::reset_rhs_matrix()
            {
                if( mUseResetValues && mParent->parent()->is_master() &&
                    mRhsMatrix0.n_cols() * mRhsMatrix0.n_rows() > 0 )
                {
                    mRhsMatrix = mRhsMatrix0 ;
                }
                else
                {
                    mRhsMatrix.fill( 0.0 );
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::assemble_jacobian(
                    Element * aElement,
                    const Matrix< real > & aJacobian )
            {
                SpMatrix & A       =  *mSystemMatrix;
                SpMatrix & E       =  *mEnforcementMatrix;
                SpMatrix & D       =  *mDirichletMatrix;
                SpMatrix & I       =  *mImpositionMatrix;
                SpMatrix & J       =  *mJacobianMatrix;

                // get dimension of element Jacobian
                uint tN = aElement->number_of_dofs() ;

                // cache dof indices
                Cell< index_t > & idx = mWorkDofIndices ;
                idx.clear();
                for ( uint i = 0; i < tN; ++i )
                {
                    idx.push( aElement->dof( i )->index() );
                }


#if !defined( NDEBUG ) || defined( DEBUG )
                for ( uint i = 0; i < tN; ++i )
                {
                    for ( uint j = 0; j < tN; ++j )
                    {
                        BELFEM_ASSERT( ! std::isnan( aJacobian( i, j ) ),
                            "Detected NaN for Element %lu", ( long unsigned int ) aElement->id() );
                    }
                }
#endif

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
                                D( idx( i ), idx( j ) ) -= aJacobian( i, j );
                            }
                            else
                            {
                                // A and J are built with the same sparsity.
                                // hence we only need to compute the position once
                                // and use the raw pointer for accessing
                                // A( idx( i ), idx( j ) ) += aJacobian( i, j );
                                // J( idx( i ), idx( j ) ) += aJacobian( i, j );
                                index_t tPos = A.index( idx( i ), idx( j  ) );

                                BELFEM_ASSERT( tPos < ( index_t ) A.number_of_nonzeros(),
                                    "entry ( %lu, %lu ) not in sparsity pattern",
                                    ( long unsigned int ) idx( i ),
                                    ( long unsigned int ) idx( j ) );

                                A.data()[ tPos ] += aJacobian( i, j );
                                J.data()[ tPos ] += aJacobian( i, j );

                            }
                        }
                    }
                    else if ( mUseJediForce )
                    {
                        for ( uint j = 0; j < tN; ++j )
                        {
                            Dof * tCol = aElement->dof( j );

                            if ( tCol->is_fixed() )
                            {
                                I( idx( i ), idx( j ) ) += aJacobian( i, j );
                            }
                            else
                            {
                                E( idx( i ), idx( j ) ) += aJacobian( i, j );
                            }
                        }
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::assemble_newton(
                    Element * aElement,
                    const Matrix< real > & adJdx)
            {

                //Add the derivative term to the jacobian
                SpMatrix & J       =  *mJacobianMatrix;

                // get dimension of element Jacobian
                uint tN = aElement->number_of_dofs() ;

                // cache dof indices
                Cell< index_t > & idx = mWorkDofIndices ;
                idx.clear();
                for ( uint i = 0; i < tN; ++i )
                {
                    idx.push( aElement->dof( i )->index() );
                }


                if ( adJdx.n_cols() > 0 )
                {
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
                                    J( idx( i ), idx( j ) ) += adJdx( i, j ) ;
                                }
                            }
                        }
                    }
                }

            }

//------------------------------------------------------------------------------

            void
            SolverData::assemble_full_matrices( Element * aElement,
                   const Matrix< real > & aMass,
                   const Matrix< real > & aStiffness )
            {

                real * M = mFullMassMatrix->data();
                real * K = mFullStiffnessMatrix->data();

                // get dimension of element Jacobian
                uint tN = aElement->number_of_dofs() ;


                // cache dof indices
                Cell< index_t > & idx = mWorkDofIndices ;
                idx.clear();
                for ( uint i = 0; i < tN; ++i )
                {
                    idx.push( aElement->dof( i )->index() );
                }

                index_t k ;
#ifdef DEBUG
                index_t nnz = mFullMassMatrix->number_of_nonzeros() ;
#endif
                for ( uint i = 0; i < tN; ++i )
                {
                    for ( uint j = 0; j < tN; ++j )
                    {
                        // we compute the index once and use it for both matrices
                        k = mFullMassMatrix->index( idx( i ), idx( j ) );
#ifdef DEBUG
                        BELFEM_ASSERT( k < nnz, "memory error" );
#endif
                        M[ k ] += aMass( i, j );
                        K[ k ] += aStiffness( i, j );
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::assemble_rhs( Element * aElement,
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
            SolverData::assemble_rhs(
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
            SolverData::collect_matrices( const bool aFullForce )
            {
                // where the historical stray message struck -- catch a leftover
                // BEFORE the innocent size collect below consumes it, and
                // name this boundary instead of an MPI_ERR_TRUNCATE
                comm_drain_check( "SolverData::collect_matrices" );

                comm_barrier() ;

                uint a = aFullForce ? 5 : 0 ;
                uint b = aFullForce ? 7 : 5 ;

                if ( mCommRank == 0 )
                {
                    Cell< Vector< real > > tAllData( mCommSize, {} );

                    for ( uint m=a; m<b; ++m )
                    {
                        MatrixType tType = static_cast< MatrixType >( m );
                        SpMatrix * tMatrix = this->matrix( tType );

                        tAllData.clear() ;
                        tAllData.set_size( mCommSize, {} );
                        comm_barrier() ;
                        collect( tAllData );
                        if ( tMatrix == nullptr ) continue ;
                        const Cell< Vector< index_t > > & tTables = this->tables( tType );

                        for ( proc_t p=1; p<mCommSize; ++p )
                        {
                            Vector< real > & tData = tAllData( p );

                            const Vector< index_t > & tTable = tTables( p );

                            index_t tNNZ = tData.length() ;
                            BELFEM_ASSERT( tTable.length() == tData.length(),
                                   "%s Matrix from proc %u has wrong number of nonzeros ( is %lu, expect %lu )",
                                   to_string( tType ).c_str(),
                                   ( unsigned int ) p,
                                   ( unsigned int ) tData.length(),
                                   ( unsigned int ) tTable.length() );

                            for ( index_t k=0; k<tNNZ; ++k )
                            {
                                tMatrix->data( tTable( k ) ) += tData( k );
                            }
                        }

                    }
                }
                else for ( uint m=a; m<b; ++m )
                {
                    SpMatrix * tMatrix = this->matrix( static_cast< MatrixType >( m ) );
                    comm_barrier() ;
                    if ( tMatrix == nullptr )
                    {
                        int_t tZero = 0 ;
                        send( tZero );
                    }
                    else
                    {
                        send( tMatrix->data(), tMatrix->number_of_nonzeros(), 0 );
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
                comm_barrier() ;

                if ( mCommRank == 0 )
                {
                    Cell< Vector< real > > tAllVectors;
                    collect( tAllVectors );

                    // assemble system
                    for ( proc_t p = 1; p < mCommSize; ++p )
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
                    send( aVector );
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::collect_rhs_matrix()
            {
                // number of cols in rhs matrix
                uint tNumCols = mParent->iwg()->num_rhs_cols();

                if ( mCommRank == 0 )
                {
                    Cell< Matrix< real > > tAllRHS;
                    collect( tAllRHS );

                    // assemble system
                    for ( proc_t p = 1; p < mCommSize; ++p )
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
                    send( mRhsMatrix );
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::update_field_values()
            {
                if ( mCommRank == 0 )
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
            SolverData::set_solver( const SolverParameters & aParams )
            {
                // delete solver if it exists already
                if( mSolver != nullptr )
                {
                    delete mSolver ;
                }

                // create a new solver
                mSolver = new Solver( aParams );

                mSolver->set_symmetry_mode( mParent->iwg()->symmetry_mode() );


            }

//------------------------------------------------------------------------------

            void
            SolverData::solve()
            {
                // get pointer to the equation
                IWG * tIWG = mParent->iwg() ;
                BELFEM_ERROR( tIWG != nullptr, "no equation was set" );

                // the iterative vector-RHS path runs through the two-phase
                // split ( certified exit ); this fused entry is the
                // composition, so unrestructured callers keep their contract
                if ( tIWG->mode() == IwgMode::Iterative
                     && tIWG->num_rhs_cols() == 1 )
                {
                    this->compute_residual();
                    this->solve_from_residual();
                    return ;
                }

                // clear a recorded soft failure from the previous solve
                // ( rank-uniform: the flag is local, the recording is
                // replicated by the wrapper )
                BELFEM_ASSERT( mSolver != nullptr, "no solver was set" );
                mSolver->wrapper()->clear_failure() ;

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
                                mSolver->solve( *mSystemMatrix, mLhsVector, mRhsVector ) ;
                                if (std::isnan(mRhsNorm))
                                {
                                    std::cout << "# NAN values caught, reseting the time step" << std::endl ;
                                    return  ;
                                }

                                // a soft solver failure must not write the
                                // garbage LHS into dofs/fields ( the contract
                                // is armed on the wrapper, so the Direct path
                                // must gate too )
                                if ( mSolver->wrapper()->failed() )
                                {
                                    return ;
                                }

                                // write values into field
                                for ( Dof * tDof: mDOFs )
                                {
                                    if ( ! tDof->is_fixed() )
                                    {
                                        tDof->value() = mLhsVector( tDof->index() );
                                    }

                                    tFields( tDof->type_id() )->value(
                                            tDof->dof_index_on_field() ) = tDof->value();
                                }

                                break ;
                            }
                            case( IwgMode::Iterative ) :
                            {
                                // unreachable: intercepted at the top of
                                // solve() and routed through the two-phase
                                // compute_residual / solve_from_residual pair
                                BELFEM_ERROR( false,
                                    "iterative vector-RHS solve must run through the two-phase path" );
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
                        mSolver->solve( *mSystemMatrix, mLhsMatrix, mRhsMatrix );

                        // soft-failure gate, cf. the vector paths above
                        if ( mSolver->wrapper()->failed() )
                        {
                            return ;
                        }

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

                                    // write data into dof
                                    tDOF->value() = mLhsMatrix( tDOF->index(), j );

                                    // write data of dof into field
                                    tData( tDOF->dof_index_on_field() ) = tDOF->value() ;
                                }

                            }
                        }
                    }

                    this->compute_hanging_dofs( tFields );

                    message( InfoLevel::Verbose, "    ... time for solving system of equations    : %u ms\n",
                             ( unsigned int ) tTimer.stop());
                }
                else if ( mSolver->type() == SolverType::MUMPS ||
                          mSolver->type() == SolverType::STRUMPACK ||
                          mSolver->type() == SolverType::PETSc )
                {

                    // note: it doesn't matter what matrix we pass here because
                    // the solver wrapper redistributes the matrix from the main proc
                    if ( tIWG->num_rhs_cols() == 1 )
                    {
                        // wait for other procs
                        comm_barrier() ;

                        mSolver->solve( *mSystemMatrix, mLhsVector, mRhsVector ) ;
                    }
                    else
                    {
                        // wait for other procs
                        comm_barrier() ;
                        mSolver->solve( *mSystemMatrix, mLhsMatrix, mRhsMatrix ) ;
                    }

                }
                else
                {
                    // wait for other procs
                    comm_barrier() ;
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::compute_residual()
            {
                IWG * tIWG = mParent->iwg() ;
                BELFEM_ERROR( tIWG != nullptr, "no equation was set" );

                BELFEM_ERROR( tIWG->mode() == IwgMode::Iterative
                              && tIWG->num_rhs_cols() == 1,
                    "compute_residual() serves only the iterative vector-RHS path - use solve()" );

                // one-shot per assembly: the load terms below are ADDED into
                // the RHS, so a second pass would double-add them
                BELFEM_ERROR( ! mResidualReady,
                    "compute_residual() called twice on one assembly - reassemble first" );

                this->collect_fields( mCollectedFields );

                if ( mKernel->is_master() )
                {
                    // scrape the fixed-dof values into retained scratch
                    // ( no per-iteration allocation on the solve path )
                    if ( mFixedValuesScratch.length() != mNumberOfFixedDofs )
                    {
                        mFixedValuesScratch.set_size( mNumberOfFixedDofs );
                    }

                    index_t tCount = 0;
                    for ( Dof * tDof: mDOFs )
                    {
                        if ( tDof->is_fixed() )
                        {
                            mFixedValuesScratch( tCount++ ) = tDof->value();
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
                        mDirichletMatrix->multiply( mFixedValuesScratch, mRhsVector, 1.0, 1.0 );
                    }

                    // ||b|| of the load-adjusted right hand side: residual()
                    // keeps reporting ||A x - b|| / ||b||
                    mRhsNorm = norm( mRhsVector );

                    BELFEM_ASSERT( mFieldValues.length() == mRhsVector.length(),
                                  "Length of Field values and RHS vector do not match ( %lu vs. %lu, free dofs: %lu )",
                                  ( long unsigned int ) mFieldValues.length(),
                                  ( long unsigned int ) mRhsVector.length(),
                                  ( long unsigned int ) mNumberOfFreeDofs );

                    switch( tIWG->algorithm() )
                    {
                        case( SolverAlgorithm::NewtonRaphson ) :
                        {
                            // preserve b DURABLY: the post-update recompute in
                            // solve_from_residual() restores b from here.
                            // mRhsBackup keeps its Picard meaning ( = r ) and
                            // must never carry b
                            if ( mRhsOriginal.length() != mRhsVector.length() )
                            {
                                mRhsOriginal.set_size( mRhsVector.length() );
                            }
                            mRhsOriginal = mRhsVector ;

                            // r = A * x - b of the COMMITTED state
                            mSystemMatrix->multiply( mFieldValues, mRhsVector, 1.0, -1.0 );

                            // pre-update residual of the ENTRY state under the
                            // CURRENT assembly ( honest backtracking reference )
                            mPreUpdateResidual = mRhsNorm > BELFEM_EPS ?
                                norm( mRhsVector ) / mRhsNorm : BELFEM_QUIET_NAN ;

                            break ;
                        }
                        case( SolverAlgorithm::Picard ) :
                        {
                            // no line search on Picard: the reference is unused
                            // and must not go stale ( Newton sets it fresh )
                            mPreUpdateResidual = BELFEM_QUIET_NAN ;

                            // increment form: r = A * x - b ( incremental-
                            // iterative scheme, Bathe 2016 §8.4.1 ). The solve
                            // in phase 2 then runs on A * delta = r
                            mSystemMatrix->multiply( mFieldValues, mRhsVector, 1.0, -1.0 );

                            break ;
                        }
                        default :
                        {
                            BELFEM_ERROR( false, "Unsupported Solve Algorithm");
                        }
                    }
                }

                // rank-uniform: every rank flips the handshake, so the
                // one-shot guard and the phase-2 precondition behave
                // identically across the communicator
                mResidualReady = true ;
            }

//------------------------------------------------------------------------------

            void
            SolverData::solve_from_residual()
            {
                IWG * tIWG = mParent->iwg() ;
                BELFEM_ERROR( tIWG != nullptr, "no equation was set" );
                BELFEM_ASSERT( mSolver != nullptr, "no solver was set" );

                // phase-1 handshake: mRhsVector must hold r, not b
                BELFEM_ERROR( mResidualReady,
                    "solve_from_residual() without a prior compute_residual()" );
                mResidualReady = false ;

                // clear a recorded soft failure from the previous solve — a
                // residual-only head exit must never read as solve_failed
                mSolver->wrapper()->clear_failure() ;

                if ( mKernel->is_master() )
                {
                    Timer tTimer;

                    switch( tIWG->algorithm() )
                    {
                        case( SolverAlgorithm::NewtonRaphson ) :
                        {
                            // wait for other procs
                            comm_barrier() ;

                            // solve the system
                            mSolver->solve( *mJacobianMatrix, mLhsVector, mRhsVector ) ;
                            if (std::isnan(mRhsNorm))
                            {
                                std::cout << "# NAN values caught, reseting the time step" << std::endl ;
                                return  ;
                            }

                            // a soft solver failure must not write
                            // the garbage LHS into dofs/fields; the
                            // controller queries solve_failed()
                            if ( mSolver->wrapper()->failed() )
                            {
                                return ;
                            }

                            // the update refreshes mFieldValues alongside the dofs:
                            // the residual recompute below multiplies with it, and
                            // a stale ( assembly-time ) vector would silently report
                            // the PRE-update residual despite the recompute
                            for ( Dof * tDof: mDOFs )
                            {
                                // update DOF values
                                if ( ! tDof->is_fixed() )
                                {
                                    tDof->value() -= tIWG->omega() * mLhsVector( tDof->index() );
                                    mFieldValues( tDof->index() ) = tDof->value() ;
                                }

                                // update value in field
                                mCollectedFields( tDof->type_id() )->value(
                                        tDof->dof_index_on_field() ) = tDof->value();
                            }

                            // recompute the residual at the UPDATED iterate, with the
                            // same lagged-operator semantics as the Picard branch;
                            // without this, residual() reports the pre-update residual
                            // and every controller decision ( promotion, stagnation
                            // guard, relaxation ) runs one iteration behind. b comes
                            // back from the DURABLE backup — never from mRhsBackup
                            mRhsVector = mRhsOriginal ;
                            mSystemMatrix->multiply( mFieldValues, mRhsVector, 1.0, -1.0 );

                            break ;
                        }
                        case( SolverAlgorithm::Picard ) :
                        {
                            // preserve r: residual() reads it after the solve
                            mRhsBackup = mRhsVector ;

                            // wait for other procs
                            comm_barrier() ;
                            // solve the system
                            mSolver->solve( *mSystemMatrix, mLhsVector, mRhsVector ) ;
                            if (std::isnan(mRhsNorm))
                            {
                                std::cout << "# NAN values caught, reseting the time step" << std::endl ;
                                return  ;
                            }

                            // a soft solver failure must not write
                            // the garbage LHS into dofs/fields
                            if ( mSolver->wrapper()->failed() )
                            {
                                return ;
                            }

                            if ( mAndersonDepth > 0 )
                            {
                                // Anderson-mixed update ( opt-in )
                                this->anderson_update( tIWG, mCollectedFields );
                            }
                            else
                            {

                            real tOmega = tIWG->omega() ;

                            for ( Dof * tDof: mDOFs )
                            {
                                // update DOF values: x - omega * delta
                                if ( !tDof->is_fixed() )
                                {
                                    tDof->value() -= tOmega * mLhsVector( tDof->index() );
                                }

                                // update value in field
                                mCollectedFields( tDof->type_id() )->value(
                                        tDof->dof_index_on_field() ) = tDof->value() ;
                            }

                            }

                            // restore the pre-update residual r = A * x - b for residual():
                            // identical content to the recompute the absolute form did here.
                            // mFieldValues is deliberately NOT refreshed: Picard reports
                            // the PRE-update residual ( always-accept semantics, Messe 2023 §4 )
                            mRhsVector = mRhsBackup ;

                            break ;
                        }
                        default :
                        {
                            BELFEM_ERROR( false, "Unsupported Solve Algorithm");
                        }
                    }

                    // hanging dofs follow their sources after every accepted update
                    this->compute_hanging_dofs( mCollectedFields );

                    message( InfoLevel::Verbose, "    ... time for solving system of equations    : %u ms\n",
                             ( unsigned int ) tTimer.stop());
                }
                else if ( mSolver->type() == SolverType::MUMPS ||
                          mSolver->type() == SolverType::STRUMPACK ||
                          mSolver->type() == SolverType::PETSc )
                {
                    // wait for other procs
                    comm_barrier() ;

                    // note: it doesn't matter what matrix we pass here because
                    // the solver wrapper redistributes the matrix from the main proc
                    mSolver->solve( *mSystemMatrix, mLhsVector, mRhsVector ) ;
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
                real aResidual = BELFEM_REAL_MAX ;

                if( mKernel->is_master() )
                {
                    // compute value
                    // note that this vector now contains the error r=A*x-b
                    // while the value of mRhsNorm was computed before with the real rhs vector

                    real tRhsNorm = norm( mRhsVector );

                    mAbsoluteResidual = tRhsNorm ;

                    aResidual = tRhsNorm/ mRhsNorm ;

                    // catch case if rhs is zero
                    if ( mRhsNorm < BELFEM_EPS || tRhsNorm < BELFEM_EPS )
                    {
                        if ( std::abs( aResidual ) < 1. + BELFEM_EPS )
                        {
                            aResidual = BELFEM_EPS;
                        }
                    }


                    //this->write_residuals_to_mesh() ;

                    // catch error
                    /*if( ( aResidual == 0  && aIteration == 0 ) || aResidual > 1E12 || std::isnan( aResidual ) )
                    {

                        mParent->save_system("error.hdf5");

                        BELFEM_ERROR( false,
                                     "ITERATION SCHEME FAILED: R2 = %8.2g \nYOU CAN TRY THE FOLLOWING THINGS: \n    * CHECK YOUR BOUNDARY CONDITIONS\n    * MAKE A BETTER CONDITIONED MESH\n    * USE A DIFFERENT PRECONDITIONER OR A DIRECT SOLVER\n    * DECREASE THE TIMESTEP\n    * DECREASE THE RELAXATION PARAMETER\n", aResidual );
                    }*/
                }

                comm_barrier();

                broadcast( aResidual );
                broadcast( mAbsoluteResidual );
                broadcast( mFixedPointResidual );
                broadcast( mPreUpdateResidual );
                return  aResidual;
            }

//------------------------------------------------------------------------------

            void
            SolverData::remember_initialization_values( const bool aSaveRHS )
            {
                if( aSaveRHS )
                {
                    if ( mRhsVector.length() > 0 )
                    {
                        mRhsVector0 = mRhsVector;
                    }
                    if ( mRhsMatrix.n_cols() * mRhsMatrix.n_rows() > 0 )
                    {
                        mRhsMatrix0 = mRhsMatrix;
                    }
                }

                mSystemValues0.set_size( mSystemMatrix->number_of_nonzeros() );
                std::copy( mSystemMatrix->data(),
                           mSystemMatrix->data() + mSystemMatrix->number_of_nonzeros(),
                           mSystemValues0.data() );

                mJacobianValues0.set_size( mJacobianMatrix->number_of_nonzeros() );
                std::copy( mJacobianMatrix->data(),
                           mJacobianMatrix->data() + mJacobianMatrix->number_of_nonzeros(),
                           mJacobianValues0.data() );
                
                if ( mEnforcementMatrix != nullptr )
                {
                    mEnforcementValues0.set_size( mEnforcementMatrix->number_of_nonzeros() );
                    std::copy( mEnforcementMatrix->data(),
                               mEnforcementMatrix->data() + mEnforcementMatrix->number_of_nonzeros(),
                               mEnforcementValues0.data() );
                }

                if( mDirichletMatrix != nullptr )
                {
                    mDirichletValues0.set_size( mDirichletMatrix->number_of_nonzeros() );
                    std::copy( mDirichletMatrix->data(),
                               mDirichletMatrix->data() + mDirichletMatrix->number_of_nonzeros(),
                               mDirichletValues0.data() );
                }

                if ( mImpositionMatrix != nullptr )
                {
                    mImpositionValues0.set_size( mImpositionMatrix->number_of_nonzeros() );
                    std::copy( mImpositionMatrix->data(),
                               mImpositionMatrix->data() + mImpositionMatrix->number_of_nonzeros(),
                               mImpositionValues0.data() );
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
                mSystemMatrix->load( tGroup, tError );
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
                if( mFieldValues.length() > 0 )
                {
                    aFile.load_data( "FieldValues", mFieldValues );
                    for ( Dof * tDof : mDOFs )
                    {
                        if( ! tDof->is_fixed() )
                        {
                            tDof->value() = mFieldValues( tDof->index() );
                        }
                    }
                }
            }

            void
            SolverData::save_system( HDF5 & aFile )
            {
                herr_t tError = 0 ;

                hid_t tGroup = aFile.create_group("Matrix");
                mSystemMatrix->save( tGroup, tError );
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
                if( mFieldValues.length() > 0 )
                {
                    aFile.save_data( "FieldValues", mFieldValues );
                }
            }

#endif
//------------------------------------------------------------------------------

            void
            SolverData::compute_hanging_dofs( Cell< mesh::Field * > & aFields )
            {
                if( mParent->is_master() )
                {
                    for ( Dof * tDof: mDofData->hanging_dofs() )
                    {
                        real tValue = 0.0;
                        for ( uint k = 0; k < tDof->number_of_sources(); ++k )
                        {
                            tValue += tDof->weight( k ) * tDof->source( k )->value();
                        }
                        tDof->value() = tValue;

                        aFields( tDof->type_id() )->value( tDof->dof_index_on_field()) = tValue;
                    }
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::set_anderson_depth( const uint aDepth )
            {
                mAndersonDepth = aDepth ;
                mAndersonStaged = false ;

                // the registers exist only while mixing is active: a
                // ShiftRegister must not be reserved with zero capacity
                if ( mAndersonX != nullptr )
                {
                    delete mAndersonX ;
                    delete mAndersonR ;
                    mAndersonX = nullptr ;
                    mAndersonR = nullptr ;
                }
                if ( aDepth > 0 )
                {
                    mAndersonX = new ShiftRegister< Vector< real > >( aDepth );
                    mAndersonR = new ShiftRegister< Vector< real > >( aDepth );
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::anderson_commit()
            {
                // controller code runs on every rank; only the master holds
                // data. Do NOT branch controller logic on anything computed
                // here — it would diverge across ranks.
                if ( mAndersonDepth == 0 || ! mKernel->is_master() )
                {
                    return ;
                }
                if ( mAndersonStaged )
                {
                    mAndersonX->push( mAndersonXStage );
                    mAndersonR->push( mAndersonRStage );
                    mAndersonStaged = false ;
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::anderson_discard()
            {
                mAndersonStaged = false ;
            }

//------------------------------------------------------------------------------

            void
            SolverData::anderson_clear()
            {
                mAndersonStaged = false ;
                if ( mAndersonX != nullptr )
                {
                    mAndersonX->clear() ;
                    mAndersonR->clear() ;
                }
            }

//------------------------------------------------------------------------------

            void
            SolverData::anderson_update( IWG * aIWG, Cell< mesh::Field * > & aFields )
            {
                // master only; mFieldValues holds the assembly-time iterate
                // x_k ( updated in compute_jacobian_and_rhs ), mLhsVector
                // holds the increment delta = x_k - G( x_k ) of the
                // increment-form solve. The fixed-point residual
                // r = G - x = -delta lives in the compact free-dof space; it
                // is distinct from the linear residual A x - b that
                // residual() reports.
                const index_t tN = mFieldValues.length() ;

                if ( mAndersonXNew.length() != tN )
                {
                    mAndersonXStage.set_size( tN );
                    mAndersonRStage.set_size( tN );
                    mAndersonXNew.set_size( tN );
                    mAndersonDeltaR.set_size( tN, mAndersonDepth );
                    mAndersonRhs.set_size( tN );
                    mAndersonGamma.set_size( mAndersonDepth );
                    mAndersonColNorm.set_size( mAndersonDepth );

                    // a changed system size invalidates any old history
                    this->anderson_clear() ;
                }

                // stage the pair BEFORE the update overwrites the dofs
                for ( index_t k = 0; k < tN; ++k )
                {
                    mAndersonXStage( k ) = mFieldValues( k );
                    mAndersonRStage( k ) = -mLhsVector( k );
                }

                // fixed-point residual ||G(x)-x|| / ||x||: the honest
                // per-iterate quality of the Picard map under mixing. The
                // reported epsilon stays the pre-update force residual
                // ( see the write-back below ); rank-uniform after the
                // broadcast in residual()
                real tXNorm = norm( mAndersonXStage );
                real tRNorm = norm( mAndersonRStage );
                mFixedPointResidual = tXNorm > BELFEM_EPS ?
                    tRNorm / tXNorm : tRNorm ;

                // beta = the live relaxation the controller line search
                // manages ( O2 decision: backtracking must damp this step )
                const uint tUsed = anderson_mixing_step(
                    mAndersonXStage,
                    mAndersonRStage,
                    *mAndersonX,
                    *mAndersonR,
                    aIWG->omega(),
                    mAndersonDeltaR,
                    mAndersonRhs,
                    mAndersonWork,
                    mAndersonGamma,
                    mAndersonColNorm,
                    mAndersonXNew );

                // O3: a pair whose difference columns just failed the solve
                // must not re-enter the history ( it would poison the next
                // attempt ); the bootstrap pair of an empty window always
                // stages, otherwise the window could never fill
                mAndersonStaged = ( tUsed > 0 ) || ( mAndersonX->size() == 0 );

                // write the mixed iterate, same masking as the legacy loop.
                // mFieldValues is deliberately NOT refreshed: it must keep the
                // assembly-time iterate, so that residual() reports the
                // pre-update force residual || A(x_k) x_k - b(x_k) || / ||b||
                // -- the published convergence criterion ( Messe et al. 2023
                // §4, Eq. 10-11 ), identical to the depth-0 path. Since the
                // increment-form solve, that residual is computed
                // BEFORE the solve and restored from mRhsBackup afterwards;
                // do NOT re-insert a post-update multiply here -- with
                // mRhsVector holding r it would compute A x - r = b and pin
                // epsilon at 1. Per-iterate mixing quality is exposed as
                // fixed_point_residual() instead.
                for ( Dof * tDof: mDOFs )
                {
                    if ( ! tDof->is_fixed() )
                    {
                        tDof->value() = mAndersonXNew( tDof->index() );
                    }
                    aFields( tDof->type_id() )->value(
                            tDof->dof_index_on_field() ) = tDof->value() ;
                }
            }

//------------------------------------------------------------------------------
        } /* end namespace dofmgr */
    } /* end namespace fem */
} /* end namespace belfem */
