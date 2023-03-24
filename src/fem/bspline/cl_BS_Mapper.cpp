//
// Created by Christian Messe on 29.04.20.
//
#include "commtools.hpp"
#include "cl_BS_Mapper.hpp"
#include "fn_linspace.hpp"
#include "cl_Element_Factory.hpp"
#include "meshtools.hpp"
#include "fn_intpoints.hpp"
#include "cl_IF_InterpolationFunctionFactory.hpp"
#include "fn_trans.hpp"
#include "fn_inv.hpp"
#include "cl_Solver.hpp"
#include "filetools.hpp"
#include "cl_HDF5.hpp"

namespace belfem
{
    namespace bspline
    {
//------------------------------------------------------------------------------

        Mapper::Mapper( const uint aNumberOfDimensions,
                        const uint aOrder,
                        const Vector <index_t> & aNumberOfElementsPerDimension,
                        const Vector< real > & aMinPoint,
                        const Vector< real > & aMaxPoint ) :
                mNumberOfDimensions( aNumberOfDimensions ),
                mOrder( aOrder ),
                mNumberOfElementsPerDimension( aNumberOfElementsPerDimension ),
                mOffset( aMinPoint )
        {
            BELFEM_ERROR( comm_size() == 1, "The Mapper can not be run in parallel mode" );

            this->create_t_matrix();
            this->create_axis( aMinPoint, aMaxPoint );
            this->create_mesh() ;
            this->create_integration_points() ;
            this->create_integration_grid();
            this->create_basis();
            this->create_bspline_elements();
            this->create_basis_adjency();
            this->create_interpolation_function();
            this->compute_element_matrices();
            this->create_jacobian();
        }

//------------------------------------------------------------------------------

        Mapper::~Mapper()
        {
            for( Vector< real > * tField : mFields )
            {
                delete tField ;
            }

            delete mJacobian ;

            delete mInterpolationFunction;

            // delete elements
            for( Element * tElement : mElements )
            {
                delete tElement ;
            }

            // delete basis
            for( Basis * tBasis : mBasis )
            {
                delete tBasis ;
            }

            delete mMesh ;
            delete mTMatrix;
        }

//------------------------------------------------------------------------------

        const Matrix< real > &
        Mapper::integration_grid() const
        {
            return mIntegrationGrid ;
        }

//------------------------------------------------------------------------------

        Vector< real > &
        Mapper::create_field( const string & aLabel )
        {
            // create the new field
            Vector< real > * aField = new Vector< real >( mGridSize );

            // safe field into memory
            mFieldLables.push( aLabel );
            mFields.push( aField );

            return *aField ;
        }
//------------------------------------------------------------------------------
        Vector< real > &
        Mapper::field( const index_t aIndex )
        {
            return *mFields( aIndex );
        }

//------------------------------------------------------------------------------

        void
        Mapper::compute_node_values()
        {
            uint tNumberOfFields = mFields.size() ;

            for( uint k=0; k<tNumberOfFields; ++k )
            {
                this->compute_dofs( k );
            }
        }

//------------------------------------------------------------------------------

        void
        Mapper::create_t_matrix()
        {
            mTMatrix = new TMatrix(mNumberOfDimensions, mOrder );

            mNumberOfNodesPerElement = mesh::number_of_nodes( mTMatrix->lagrange_type() );
            mNumberOfBasisPerElement = mesh::number_of_nodes( mTMatrix->bspline_type() );

        }

//------------------------------------------------------------------------------
        void
        Mapper::create_axis( const Vector< real >   & aMinPoint,
                             const Vector< real >   & aMaxPoint )
        {
            mNumberOfNodesPerDimension.set_size( mNumberOfDimensions );
            mElementLength.set_size( mNumberOfDimensions );

            mNumberOfNodes = 1 ;
            mNumberOfElements = 1 ;
            for( uint k=0; k<mNumberOfDimensions; ++k )
            {
                mNumberOfNodesPerDimension( k ) = mOrder * mNumberOfElementsPerDimension( k ) + 1 ;
                mNumberOfNodes *= mNumberOfNodesPerDimension( k ) ;
                mNumberOfElements *= mNumberOfElementsPerDimension( k );
                mElementLength( k ) = ( aMaxPoint( k ) - aMinPoint( k ) ) / mNumberOfElementsPerDimension( k );
            }

            linspace(
                    aMinPoint( 0 ),
                    aMaxPoint( 0 ),
                    mNumberOfNodesPerDimension( 0 ),
                    mAxisX );

            if( mNumberOfDimensions > 1 )
            {
                linspace(
                        aMinPoint( 1 ),
                        aMaxPoint( 1 ),
                        mNumberOfNodesPerDimension( 1 ),
                        mAxisY);

                if( mNumberOfDimensions > 2 )
                {
                    linspace(
                            aMinPoint( 2 ),
                            aMaxPoint( 2 ),
                            mNumberOfNodesPerDimension( 2 ),
                            mAxisZ );
                }
            }
        }
//------------------------------------------------------------------------------

        void
        Mapper::create_mesh()
        {
            mMesh = new Mesh( mNumberOfDimensions, 0 );

            this->create_nodes();
            this->create_elements();
            mMesh->finalize();
            this->create_global_variables();
        }

//------------------------------------------------------------------------------

        void
        Mapper::create_nodes()
        {
            // create the nodes
            Cell< mesh::Node * > & tNodes = mMesh->nodes() ;

            tNodes.set_size( mNumberOfNodes, nullptr );

            index_t tCount = 0;

            if( mNumberOfDimensions == 1 )
            {
                for( index_t i=0; i<mNumberOfNodesPerDimension( 0 ); ++i )
                {
                    tNodes( tCount ) = new mesh::Node( tCount+1, mAxisX( i ) );
                    ++tCount ;
                }
            }
            else if ( mNumberOfDimensions == 2 )
            {
                for( index_t j=0; j<mNumberOfNodesPerDimension( 1 ); ++j )
                {
                    for ( index_t i = 0; i < mNumberOfNodesPerDimension( 0 ); ++i )
                    {
                        tNodes( tCount ) = new mesh::Node( tCount + 1,
                                mAxisX( i ),
                                mAxisY( j ) );
                        ++tCount ;
                    }
                }
            }
            else if ( mNumberOfDimensions == 3 )
            {
                for( index_t k=0; k<mNumberOfNodesPerDimension( 2 ); ++k )
                {
                    for ( index_t j = 0; j < mNumberOfNodesPerDimension( 1 ); ++j )
                    {
                        for ( index_t i = 0; i < mNumberOfNodesPerDimension( 0 ); ++i )
                        {
                            tNodes( tCount ) = new mesh::Node( tCount + 1,
                                                                 mAxisX( i ),
                                                                 mAxisY( j ),
                                                                 mAxisZ( k ) );
                            ++tCount ;
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Mapper::create_elements()
        {
            // container of nodes on the mesh
            Cell< mesh::Node * > & tNodes = mMesh->nodes() ;

            // create a block on the mesh
            Cell< mesh::Block *> & tBlocks = mMesh->blocks();
            tBlocks.set_size( 1, nullptr );
            mesh::Block * tBlock = new mesh::Block( 1, mNumberOfElements ) ;
            tBlocks( 0 ) = tBlock;

            // get the element type
            const ElementType tType = mTMatrix->lagrange_type() ;

            // get the node index
            const Matrix< index_t > & tIndex = mTMatrix->node_index() ;

            mesh::ElementFactory tFactory;

            if( mNumberOfDimensions == 1 )
            {
                index_t tOff = 0;

                id_t tCount = 1;
                for ( index_t I = 0; I < mNumberOfElementsPerDimension( 0 ); ++I )
                {
                    // create a new element
                    mesh::Element * tElement = tFactory.create_element(
                            tType,
                            tCount++ );

                    // link element with nodes
                    for( uint n=0; n<mNumberOfNodesPerElement; ++n )
                    {
                        tElement->insert_node( tNodes( tOff + tIndex( 0, n ) ), n );
                    }

                    // add element to block
                    tBlock->insert_element( tElement );

                    // increment anchor offset
                    tOff += mOrder ;


                }
            }
            else if ( mNumberOfDimensions == 2 )
            {
                index_t tOff = 0;
                index_t tOff0 = 0;

                id_t tCount = 1;
                for( index_t j=0; j<mNumberOfElementsPerDimension( 1 ); ++j )
                {
                    tOff = tOff0 ;
                    for ( index_t i = 0; i < mNumberOfElementsPerDimension( 0 ); ++i )
                    {
                        // create a new element
                        mesh::Element * tElement = tFactory.create_element(
                                tType,
                                tCount++ );
                        // link element with nodes
                        for( uint n=0; n<mNumberOfNodesPerElement; ++n )
                        {
                            tElement->insert_node( tNodes(
                                    tOff + tIndex( 1, n ) * mNumberOfNodesPerDimension( 0 ) + tIndex( 0, n ) ), n );
                        }

                        // add element to block
                        tBlock->insert_element( tElement );

                        // increment anchor offset
                        tOff += mOrder ;
                    }
                    tOff0 += mOrder*mNumberOfNodesPerDimension( 0 );
                }
            }
            else if ( mNumberOfDimensions == 3 )
            {
                index_t tOff = 0;
                index_t tOff1 = 0;
                id_t tCount = 1;
                for( index_t k=0; k<mNumberOfElementsPerDimension( 2 ); ++k )
                {
                    index_t tOff0 = tOff1;
                    for ( index_t j = 0; j < mNumberOfElementsPerDimension( 1 ); ++j )
                    {
                        tOff = tOff0;
                        for ( index_t i = 0; i < mNumberOfElementsPerDimension( 0 ); ++i )
                        {
                            // create a new element
                            mesh::Element * tElement = tFactory.create_element(
                                    tType,
                                    tCount++ );

                            // link element with nodes
                            for ( uint n = 0; n < mNumberOfNodesPerElement; ++n )
                            {
                                tElement->insert_node( tNodes(
                                        tOff + mNumberOfNodesPerDimension( 0 ) *
                                                ( tIndex( 2, n ) * mNumberOfNodesPerDimension( 1 ) + tIndex( 1, n ) )
                                              + tIndex( 0, n )), n );
                            }

                            // add element to block
                            tBlock->insert_element( tElement );

                            // increment anchor offset
                            tOff += mOrder;
                        }
                        tOff0 += mOrder*mNumberOfNodesPerDimension( 0 );
                    }
                    tOff1 += mOrder*mNumberOfNodesPerDimension( 0 ) * mNumberOfNodesPerDimension( 1 );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Mapper::create_integration_points()
        {
            // we can afford the classic scheme ( product combination )
            // since most computation speed is lost during the ray cast
            if( mNumberOfDimensions == 1 )
            {
                intpoints( IntegrationScheme::GAUSSCLASSIC,
                           GeometryType::LINE,
                           2*mOrder + 1,
                           mIntegrationWeights,
                           mIntegrationPoints );
            }
            else if ( mNumberOfDimensions == 2 )
            {
                intpoints( IntegrationScheme::GAUSSCLASSIC,
                           GeometryType::QUAD,
                           2*mOrder + 1,
                           mIntegrationWeights,
                           mIntegrationPoints );

            }
            else if ( mNumberOfDimensions == 3 )
            {
                intpoints( IntegrationScheme::GAUSSCLASSIC,
                           GeometryType::HEX,
                           2*mOrder + 1,
                           mIntegrationWeights,
                           mIntegrationPoints );
            }

        }

//------------------------------------------------------------------------------

        void
        Mapper::create_global_variables()
        {

            id_t tCount = 1;
            // global variables to be understood by the interpolator object
            mMesh->create_global_variable(
                    "numElemsX",
                    ( real ) mNumberOfElementsPerDimension( 0 ),
                    tCount++ );

            mMesh->create_global_variable(
                    "Xmin",
                    mAxisX( 0 ),
                    tCount++ );

            mMesh->create_global_variable(
                    "Xmax",
                    mAxisX( mAxisX.length() - 1 ),
                    tCount++ );

            mMesh->create_global_variable(
                    "DeltaX",
                    mElementLength( 0 ),
                    tCount++ );

            if( mNumberOfDimensions > 1 )
            {
                mMesh->create_global_variable(
                        "numElemsY",
                        ( real ) mNumberOfElementsPerDimension( 1 ),
                        tCount++ );

                mMesh->create_global_variable(
                        "Ymin",
                        mAxisY( 0 ),
                        tCount++ );

                mMesh->create_global_variable(
                        "Ymax",
                        mAxisY( mAxisY.length() - 1 ),
                        tCount++ );

                mMesh->create_global_variable(
                        "DeltaY",
                        mElementLength( 1 ),
                        tCount++ );

                if( mNumberOfDimensions > 2 )
                {
                    mMesh->create_global_variable(
                            "numElemsZ",
                            ( real ) mNumberOfElementsPerDimension( 2 ),
                            tCount++ );

                    mMesh->create_global_variable(
                            "Zmin",
                            mAxisZ( 0 ),
                            tCount++ );

                    mMesh->create_global_variable(
                            "Zmax",
                            mAxisZ( mAxisZ.length() - 1 ),
                            tCount++ );

                    mMesh->create_global_variable(
                            "DeltaZ",
                            mElementLength( 2 ),
                            tCount++ );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Mapper::create_integration_grid()
        {
            uint tNumPointsPerElement = mIntegrationWeights.length() ;

            mGridSize = tNumPointsPerElement * mNumberOfElements ;
            mIntegrationGrid.set_size( mNumberOfDimensions,
                                       mGridSize );

            Cell < mesh::Element * > & tElements = mMesh->block( 1 )->elements();

            Vector< real > tCenter( mNumberOfDimensions );

            index_t tCount = 0;

            // loop over all elements
            for( mesh::Element * tElement : tElements )
            {
                // compute center of element
                for( uint i=0; i<mNumberOfDimensions; ++i )
                {
                    tCenter( i ) = tElement->node( 0 )->x( i ) + 0.5 * mElementLength( i );
                }

                for( uint k=0; k<tNumPointsPerElement; ++k )
                {
                    for( uint i=0; i<mNumberOfDimensions; ++i )
                    {
                        mIntegrationGrid( i, tCount )
                            = tCenter( i ) + 0.5 * mElementLength( i ) * mIntegrationPoints( i, k );
                    }

                    // increment counter
                    ++tCount ;
                }
            }

        }

//------------------------------------------------------------------------------

        void
        Mapper::create_basis()
        {
            // determine number of basis in this mesh
            mNumberOfBasis = 1;
            mNumberOfBasisPerDimension.set_size( mNumberOfDimensions );
            for ( uint i = 0; i < mNumberOfDimensions; ++i )
            {
                mNumberOfBasisPerDimension( i ) = mNumberOfElementsPerDimension( i ) + mOrder;
                mNumberOfBasis *= mNumberOfBasisPerDimension( i );
            }

            // allocate memory
            mBasis.set_size( mNumberOfBasis, nullptr );

            // vector with dofs
            mDOFs.set_size( mNumberOfBasis );

            // vector with RHS
            mRHS.set_size( mNumberOfBasis );

            for ( index_t k = 0; k<mNumberOfBasis; ++k )
            {
                // create dof
                mBasis( k ) = new Basis( k+1 );

                // set index and ID of this dof
                mBasis( k )->set_index( k );
            }
        }

//------------------------------------------------------------------------------

        void
        Mapper::create_bspline_elements()
        {
            // elements on mesh
            Cell < mesh::Element * > & tElements = mMesh->block( 1 )->elements();

            mElements.set_size( tElements.size(), nullptr );

            // create elements and link them to the elements on the mesh
            for( index_t k=0; k<mNumberOfElements; ++k )
            {
                mElements( k ) = new Element( tElements( k ) );
            }

            // get basis indices
            const Matrix< index_t > & tIndex = mTMatrix->basis_index();

            index_t tCount = 0;

            if( mNumberOfDimensions == 1 )
            {
                index_t tOff = 0;

                for( index_t i=0; i<mNumberOfElementsPerDimension( 0 ); ++i )
                {
                    // create a new element
                    Element * tElement = new Element( tElements( tCount ) );

                    // link element with basis
                    for( uint n=0; n<mNumberOfBasisPerElement; ++n )
                    {
                        tElement->set_basis( mBasis( tOff + tIndex( 0, n ) ), n );
                    }

                    tOff += 1;

                    // add element to container
                    mElements( tCount++ ) = tElement;
                }
            }
            else if ( mNumberOfDimensions == 2 )
            {
                index_t tOff0 = 0;

                for( index_t j=0; j<mNumberOfElementsPerDimension( 1 ); ++j )
                {
                    index_t tOff = tOff0;

                    for( index_t i=0; i<mNumberOfElementsPerDimension( 0 ); ++i )
                    {
                        // create a new element
                        Element * tElement = new Element( tElements( tCount ) );

                        // link element with basis
                        for( uint n=0; n<mNumberOfBasisPerElement; ++n )
                        {
                            tElement->set_basis( mBasis(
                                    tOff + tIndex( 1, n ) * mNumberOfBasisPerDimension( 0 ) + tIndex( 0, n ) ), n );
                        }

                        tOff += 1;

                        // add element to container
                        mElements( tCount++ ) = tElement;
                    }

                    tOff0 += mNumberOfBasisPerDimension( 0 );
                }
            }
            else if ( mNumberOfDimensions == 3 )
            {
                index_t tOff = 0;
                index_t tOff1 = 0;
                for( index_t k=0; k<mNumberOfElementsPerDimension( 2 ); ++k )
                {
                    index_t tOff0 = tOff1;
                    for ( index_t j = 0; j < mNumberOfElementsPerDimension( 1 ); ++j )
                    {
                        tOff = tOff0;
                        for ( index_t i = 0; i < mNumberOfElementsPerDimension( 0 ); ++i )
                        {
                            // create a new element
                            Element * tElement = new Element( tElements( tCount ) );

                            // link element with nodes
                            for ( uint n = 0; n < mNumberOfBasisPerElement; ++n )
                            {
                                tElement->set_basis( mBasis(
                                        tOff + mNumberOfBasisPerDimension( 0 ) *
                                               ( tIndex( 2, n ) * mNumberOfBasisPerDimension( 1 ) + tIndex( 1, n ) )
                                        + tIndex( 0, n )), n );
                            }

                            // increment anchor offset
                            tOff += 1;
                        }
                        tOff0 += mNumberOfBasisPerDimension( 0 );
                    }
                    tOff1 += mNumberOfBasisPerDimension( 0 ) * mNumberOfBasisPerDimension( 1 );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Mapper::create_basis_adjency()
        {

            // step 1: count elements per basis
            for( Element * tElement : mElements )
            {
                // loop over all basis of this element
                for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                {
                    tElement->basis( k )->increment_element_counter();
                }
            }

            // step 2: initialize element containers
            for( Basis * tBasis : mBasis )
            {
                tBasis->init_element_container();
            }

            // step 3: insert elements into basis
            for( Element * tElement : mElements )
            {
                // loop over all basis of this element
                for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                {
                    tElement->basis( k )->insert_element( tElement );
                }
            }

            // step 4: create graph
            for( Basis * tBasis : mBasis )
            {
                tBasis->link_basis();
            }

        }

//------------------------------------------------------------------------------

        void
        Mapper::create_interpolation_function()
        {
            fem::InterpolationFunctionFactory tFactory ;

            mInterpolationFunction = tFactory.create_lagrange_function(
                    mTMatrix->lagrange_type() );
        }

//------------------------------------------------------------------------------

        void
        Mapper::create_jacobian()
        {
            mJacobian = new SpMatrix( reinterpret_cast< Cell< graph::Vertex * > & >( mBasis ) );

            // convert pointer to ref
            SpMatrix & tM = *mJacobian ;

            // index vector
            Vector< index_t > tIndex( mNumberOfBasisPerElement );

            // assemble jacobian
            for( Element * tElement : mElements )
            {
                // populate element index
                for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                {
                    tIndex( k ) = tElement->basis( k )->index() ;
                }

                // add mass matrix to element
                for( uint j=0; j<mNumberOfBasisPerElement; ++j )
                {
                    for( uint i=0; i<mNumberOfBasisPerElement; ++i )
                    {
                        tM( tIndex( i ), tIndex( j ) ) += mMassMatrix( i, j );
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Mapper::compute_element_matrices()
        {
            uint tNumberOfIntegrationPoints = mIntegrationWeights.length();

            // allocate memory
            Matrix< real > tM( mNumberOfNodesPerElement, mNumberOfNodesPerElement, 0.0 );

            // for RHS
            Matrix< real > tR( mNumberOfNodesPerElement, tNumberOfIntegrationPoints );

            // determinant of geometry jacobian
            real tDetJ = 1;
            for( uint i=0; i<mNumberOfDimensions; ++i )
            {
                tDetJ *= 0.5 * mElementLength( i ) ;
            }

            // the N Matrix
            Matrix< real > tN( 1, mNumberOfNodesPerElement );

            // vector with point data
            Vector< real > tXi( mNumberOfDimensions );

            // loop over all integration points
            for( uint k=0; k<tNumberOfIntegrationPoints; ++k )
            {
                // populate point
                for( uint i=0; i<mNumberOfDimensions; ++i )
                {
                    tXi( i ) = mIntegrationPoints( i, k );
                }

                // compute function
                mInterpolationFunction->N( tXi, tN );

                // add value to mass matrix
                tM += mIntegrationWeights( k ) * trans( tN ) * tN * tDetJ ;

                for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                {
                    tR( i, k ) = mIntegrationWeights( k ) * tN( 0, i ) * tDetJ ;
                }
            }

            mMassMatrix = trans( mTMatrix->Lagrange() ) * tM * mTMatrix->Lagrange();
            mRhsMatrix  = trans( mTMatrix->Lagrange() ) * tR ;
        }

//------------------------------------------------------------------------------

        void
        Mapper::compute_dofs( const index_t aFieldIndex )
        {
            // reset vectors
            mDOFs.fill( 0.0 );
            mRHS.fill( 0.0 );

            // get field
            Vector< real > & tIntegrationField = *mFields( aFieldIndex );

            // init counter
            index_t tCount = 0;

            // number of integration points per element
            uint tNumberOfIntegrationPoints = mIntegrationWeights.length() ;

            // matrix with values
            Vector< real > tValues( tNumberOfIntegrationPoints );

            // RHS vector
            Vector< real > tRHS( mNumberOfBasisPerElement );

            // loop over all elements
            for ( Element * tElement : mElements )
            {
                // populate values vector
                for( uint i=0; i<tNumberOfIntegrationPoints; ++i )
                {
                   tValues( i ) =  tIntegrationField( tCount++ );
                }

                tRHS = mRhsMatrix * tValues ;

                // add values to RHS
                for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                {
                    mRHS( tElement->basis( k )->index() ) += tRHS( k );
                }
            }

            // create a solver
            Solver tSolver( SolverType::UMFPACK );

            // solve system
            tSolver.solve( *mJacobian, mDOFs, mRHS ) ;

            // delete the solver
            tSolver.free();

            // create the nodal field
            Vector< real > & tField = mMesh->create_field( mFieldLables( aFieldIndex ) );

            // transformation matrix
           const Matrix< real > & tT = mTMatrix->Lagrange() ;

            // DOFs and node values
            Vector< real > tBasis( mNumberOfBasisPerElement );
            Vector< real > tNodes( mNumberOfNodesPerElement );

            // unflag all nodes on the mesh
            mMesh->unflag_all_nodes();

            // now, we can compute the node values
            for( Element * tElement : mElements )
            {
                // populate basis
                for( uint k=0; k<mNumberOfBasisPerElement; ++k )
                {
                    tBasis( k ) = mDOFs( tElement->basis( k )->index() );
                }

                // compute nodal values
                tNodes = tT * tBasis ;

                // write values onto mesh
                for( uint i=0; i<mNumberOfNodesPerElement; ++i )
                {
                    if( ! tElement->element()->node( i )->is_flagged() )
                    {
                        tField( tElement->element()->node( i )->index() ) = tNodes( i );
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Mapper::write_field_to_database(
                const string & aLabel,
                const string & aDatabase,
                const uint     aDimensions )
        {

            // check if file exists
            FileMode tMode = file_exists( aDatabase ) ? FileMode::OPEN_RDWR : FileMode::NEW ;

            // open database
            HDF5 tFile( aDatabase, tMode );

            // create a new group
            tFile.create_group( aLabel );

            // write the dimension
            uint tDimension = aDimensions == 0 ? mNumberOfDimensions : aDimensions ;
            tFile.save_data( "dimension", tDimension );

            // write the number of points
            Vector<  uint > tNumPoints( tDimension );
            for( uint k=0; k<tDimension; ++k )
            {
                tNumPoints( k ) = mNumberOfNodesPerDimension( k );
            }
            tFile.save_data( "numpoints", tNumPoints );

            // write the offset
            Vector< double > tOffset( tDimension );
            for( uint k=0; k<tDimension; ++k )
            {
                tOffset( k ) = mOffset( k );
            }
            tFile.save_data( "offset", tOffset );

            // write the interpolation order
            tFile.save_data("order", mOrder );

            // write the stepsize
            Vector< double > tStep( tDimension );
            tStep( 0 ) = mAxisX( mOrder ) - mAxisX( 0 );

            index_t tMemory = mNumberOfNodesPerDimension( 0 );

            if( tDimension > 1 )
            {
                tStep( 1 ) = mAxisY( mOrder ) - mAxisY( 0 );
                tMemory *= mNumberOfNodesPerDimension( 1 );
                if( tDimension > 2 )
                {
                    tStep( 2 ) = mAxisZ( mOrder ) - mAxisZ( 0 );
                    tMemory *= mNumberOfNodesPerDimension( 2 );
                }
            }
            tFile.save_data("step", tStep );

            // write the dataset
            Vector< real > & tData = mMesh->field_data( aLabel );
            Vector< double > tValues( tMemory );
            for( index_t k=0; k<tMemory; ++k )
            {
                tValues( k ) = tData( k );
            }
            tFile.save_data("values", tValues );

            // close the group
            tFile.close_active_group() ;

            // close the file
            tFile.close() ;
        }

//------------------------------------------------------------------------------
    }
}