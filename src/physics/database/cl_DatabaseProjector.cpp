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

#include "commtools.hpp"
#include "cl_Logger.hpp"
#include "fn_intpoints.hpp"
#include "cl_DatabaseProjector.hpp"

#include "cl_IF_InterpolationFunctionFactory.hpp"
#include "cl_TensorMeshFactory.hpp"
#include "fn_trans.hpp"
#include "cl_Solver.hpp"
#include "fn_matrix_type.hpp"
#include "cl_Timer.hpp"
#include "fn_unique.hpp"

namespace belfem
{
    namespace database
    {
        Projector::Projector( Mesh * aMesh ) :
        mCommRank( comm_rank() ),
        mCommSize( comm_size() ),
        mConfig( aMesh->tensorconf() ), // <-- will trow error if not tensormesh
        mMesh( aMesh )
        {

            // note: we don't want STRUMPACK here because the
            //       matrices are too small
#ifdef BELFEM_MUMPS
                SolverType tSolverType = SolverType::MUMPS;
#elif BELFEM_PARDISO
                SolverType tSolverType = SolverType::PARDISO;
#elif BELFEM_SUPERLU
                SolverType tSolverType = SolverType::SUPERLU;
#elif BELFEM_SUITESPARSE
                SolverType tSolverType = SolverType::UMFPACK ;
#else
            BELFEM_ERROR( false, "Need either one of MUMPS, PARDISO, SUPERLU, or SUITESPARSE" );
#endif

            mSolver = new Solver( tSolverType );

            if ( mCommRank == 0 )
            {
                this->allocate_system_matrix( aMesh );
                this->compute_element_matrices();
            }
            else
            {
                // the workers join a distributed solve with an empty matrix
                // ( the master supplies the system ); the object must exist,
                // because project() binds a reference to it before the solve
                mA = new SpMatrix();
            }
        }

        Projector::~Projector()
        {
            if ( mSolver != nullptr )
            {
                delete mSolver ;
            }
            if ( mA != nullptr )
            {
                delete mA ;
            }
        }

        void
        Projector::allocate_system_matrix( Mesh * aMesh )
        {
            if ( mCommRank == 0 )
            {
                TensorMeshFactory tFactory;
                if ( aMesh->number_of_control_points() == 0 )
                {
                    tFactory.create_bsplines( aMesh );
                    mT = tFactory.t_matrix();
                }
                else
                {
                    mT = tFactory.t_matrix( aMesh->tensorconf() );
                }

                // create the graph
                Cell< mesh::ControlPoint * > & tControlPoints = aMesh->control_points();
                Cell< mesh::Element * > & tElements = aMesh->elements();

                // make sure than indices are sane
                index_t tNumPoints = 0 ;
                index_t tNumPointsPerElem = mesh::number_of_nodes( mConfig->element_type() );

                for ( mesh::ControlPoint * tControlPoint : tControlPoints )
                {
                    tControlPoint->set_index( tNumPoints++ );
                }

                // create graph
                Graph tGraph( tNumPoints, nullptr );
                for ( uint k=0; k<tNumPoints; ++k )
                {
                    tGraph( k ) = new graph::Vertex();
                    tGraph( k )->set_index( k );
                }
                Vector< uint > tNumElemsPerPoint( tNumPoints, 0 );
                for ( mesh::Element * tElement : tElements )
                {
                    for ( uint k=0; k<tElement->number_of_control_points(); ++k )
                    {
                        ++tNumElemsPerPoint( tElement->control_point( k )->index() ) ;
                    }
                }
                Cell< Cell< mesh::Element * > > tElementsPerPoint( tNumPoints, {}) ;
                for ( index_t k=0; k<tNumPoints; ++k )
                {
                    tElementsPerPoint( k ).reserve( tNumElemsPerPoint( k ) );
                }
                tNumElemsPerPoint.fill( 0 );
                for ( mesh::Element * tElement : tElements )
                {
                    for ( uint k=0; k<tElement->number_of_control_points(); ++k )
                    {

                        tElementsPerPoint( tElement->control_point( k )->index() ).push( tElement );
                    }
                }
                for ( index_t k=0; k<tNumPoints; ++k )
                {
                    Cell< index_t > tIndices( tNumPointsPerElem * tElementsPerPoint( k ).size(), 0 );

                    index_t tCount = 0 ;
                    for ( mesh::Element * tElement : tElementsPerPoint( k ) )
                    {
                        for ( uint j=0; j<tNumPointsPerElem; ++j )
                        {
                            tIndices( tCount++ ) = tElement->control_point( j )->index();
                        }
                    }
                    unique( tIndices );
                    graph::Vertex * tVertex = tGraph( k ) ;
                    tCount = tIndices.size() ;

                    tVertex->init_vertex_container( tCount );
                    for ( uint i=0; i<tCount; ++i )
                    {
                        tVertex->insert_vertex( tGraph( tIndices( i ) ) );
                    }
                }

                // create the matrix
                mA = new SpMatrix( tGraph, matrix_type( mSolver->type() ) );

                // tidy up graph
                for ( graph::Vertex * tVertex : tGraph )
                {
                    delete tVertex ;
                }
                tGraph.clear();
            }
        }

        void
        Projector::compute_element_matrices()
        {
            fem::InterpolationFunctionFactory tFactory;

            fem::InterpolationFunction * tFun = tFactory.create_lagrange_function( mConfig->element_type() );

            Vector< real > w;
            Matrix< real > xi;

            uint tOrder = mConfig->order() * mConfig->order() ;

            intpoints(
                IntegrationScheme::GAUSSCLASSIC,
                mesh::geometry_type( mConfig->element_type() ),
                tOrder,
                w,
                xi );

            uint tNumNodes = mesh::number_of_nodes( mConfig->element_type() );

            mMel.set_size( tNumNodes, tNumNodes, 0.0 );

            Matrix< real > N( 1, tNumNodes) ;

            uint n = w.length();

            for ( uint k=0; k<n; ++k )
            {
                tFun->N( xi.col( k ), N  );
                mMel += w(k) * trans( N ) * N ;
            }

            for ( uint d=0; d<mConfig->dimension(); ++d )
            {
                mMel *= 0.5 * mConfig->element_step( d ) ;
            }

            mBel = trans( mT ) * mMel ;
            mAel = mBel * mT ;

            delete tFun;
        }

        void
        Projector::project( const string & aField, Vector< real > & aResult )
        {
            Vector< real > X ;
            Vector< real > Y ;
            if ( mCommRank == 0 )
            {
                // get the field data
                Vector< real > F = mMesh->field( aField )->data() ;

                // allocate the rhs vector
                X.set_size( mMesh->number_of_control_points(), 0.0 );
                Y.set_size( mMesh->number_of_control_points(), 0.0 );
                // assemble the matrix
                mA->fill( 0.0 );

                Cell< mesh::Element * > & tElements = mMesh->elements();

                uint n = mesh::number_of_nodes( mConfig->element_type() );
                Vector< real > Fel( n );
                Vector< index_t > idx( n );
                Vector< real > Yel( n );

                SpMatrix & A = *mA;

                for ( mesh::Element * tElement : tElements )
                {
                    // populate indices
                    for ( uint i=0; i<n; ++i )
                    {
                        idx( i ) = tElement->control_point( i )->index();
                    }

                    // assemble matrix
                    for ( uint i=0; i<n; ++i )
                    {
                        for ( uint j=0; j<n; ++j )
                        {
                            A( idx( i ), idx( j ) ) += mAel( i, j );
                        }
                    }

                    // populate rhs
                    for ( uint i=0; i<n; ++i )
                    {
                        Fel( i ) = F( tElement->node( i )->index() );
                    }

                    Yel = mBel * Fel ;

                    // assemble rhs
                    for ( uint i=0; i<n; ++i )
                    {
                        Y( idx( i ) ) += Yel( i );
                    }
                }

                comm_barrier();
                Timer tTimer;
                message( InfoLevel::Default, "\n    creating lookup table for %s ...", mLabel.size() == 0 ? aField.c_str() : mLabel.c_str() );
                mSolver->solve( A, X, Y );
                message( InfoLevel::Default, "    ... elapsed time %4.1f s.\n",
                         ( float ) tTimer.stop() * 0.001 );

                // project back
                mMesh->unflag_all_nodes();
                Vector< real > Xel( n );

                aResult.set_size( mMesh->number_of_nodes(), 0.0 );

                for ( mesh::Element * tElement : tElements )
                {
                    // populate indices
                    for ( uint i=0; i<n; ++i )
                    {
                        Xel( i ) = X( tElement->control_point( i )->index() );
                    }

                    Fel = mT * Xel ;

                    for ( uint i=0; i<n; ++i )
                    {
                        if ( ! tElement->node( i )->is_flagged() )
                        {
                            aResult( tElement->node( i )->index() ) = Fel( i );
                            tElement->node( i )->flag();
                        }
                    }
                }

                /*Vector< real > & tZ = mMesh->create_field( aField + "_proj" );
                tZ = aResult ;
                string tFile = aField + ".exo";
                mMesh->save( tFile );*/

                comm_barrier();
                share( aResult );
            }
            else
            {
                comm_barrier();

                if ( mSolver->wrapper()->uses_mpi() )
                {
                    SpMatrix & A = *mA;
                    mSolver->solve( A, X, Y );
                }

                comm_barrier();
                receive( aResult );
            }
        }

    }
}
