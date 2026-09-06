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

#include <iostream>

#include "typedefs.hpp"
#include "cl_Communicator.hpp"
#include "cl_Logger.hpp"
#include "banner.hpp"
#include "cl_Mesh_GmshReader.hpp"
#include "cl_Mesh_ExodusWriter.hpp"
#include "cl_Mesh.hpp"
#include "cl_Mesh_OrientationChecker.hpp"
#include "cl_Visualizer.hpp"
#include "cl_Element_Factory.hpp"

using namespace belfem;
using namespace mesh;

Communicator gComm;
Logger       gLog( 5 );

int main( int    argc,
          char * argv[] )
{
    // create communicator
    gComm.init( argc, argv );


    Matrix< index_t > tTetTopo = { { 0, 4, 6, 7, 10, 31, 15, 16, 23, 29 },
                                   { 4, 1, 5, 8, 11, 12, 33, 23, 19, 26 },
                                   { 6, 5, 2, 9, 32, 13, 14, 30, 27, 21 },
                                   { 7, 8, 9, 3, 24, 27, 29, 17, 20, 22 } };

    Matrix< index_t > tPyraTopo  = { { 5, 8, 7, 6, 4, 25, 24, 28, 32, 33, 23, 22, 31, 34 },
                                     { 5, 6, 7, 8, 9, 32, 28, 24, 25, 26, 30, 29, 27, 34 } };


    // tet geometries:
    real tA = 1.0 ;
    real tB = 0.5 * std::sqrt( 3. ) * tA ;
    real tH = tA * std::sqrt( 6 ) / 3.0 ;

    ElementFactory tFactory ;

    Mesh * tMesh = new Mesh( 3, 0, true );

    Cell< Node * > & tNodes = tMesh->nodes() ;
    tNodes.set_size( 4, nullptr );

    id_t tNodeID = 0 ;
    id_t tElementID = 0 ;

    tNodes( 0 ) = new Node( ++tNodeID, -0.5*tA, -tB/3.0, -tH/4.0 );
    tNodes( 1 ) = new Node( ++tNodeID,  0.5*tA, -tB/3.0, -tH/4.0 );
    tNodes( 2 ) = new Node( ++tNodeID,  0, 2.0*tB/3.0, -tH/4.0 );
    tNodes( 3 ) = new Node( ++tNodeID,  0 , 0, 0.75 * tH );

    Block * tBlock = new Block( 1, 1 );
    Element * tElement = tFactory.create_element( ElementType::TET4, ++tElementID );
    for( uint k=0; k<4; ++k )
    {
        tElement->insert_node( tNodes( k ), k );
    }
    tBlock->elements()( 0 ) = tElement ;
    tMesh->blocks().push( tBlock );

    tMesh->finalize() ;
    tMesh->create_edges();
    tMesh->create_faces();
    tMesh->unfinalize() ;

    // create nodes
    tElement->flag_edges() ;
    index_t tCount = 0 ;
    for( Edge * tEdge : tMesh->edges() )
    {
        if( tEdge->is_flagged() )
        {
            tEdge->set_index( tCount++ );
        }
    }

    Cell< Node * > tOldNodes ;
    tOldNodes.vector_data() = std::move( tNodes.vector_data() );
    tNodes.set_size( tOldNodes.size() + tCount, nullptr );
    index_t tNumNodes = tOldNodes.size() ;
    tCount = tNumNodes ;

    std::memcpy( tNodes.data(), tOldNodes.data(), tOldNodes.size() * sizeof( Node * ) );
    for( Edge * tEdge : tMesh->edges() )
    {
        if( tEdge->is_flagged() )
        {
            real tX = 0.5 * ( tEdge->node( 0 )->x() + tEdge->node( 1 )->x() ) ;
            real tY = 0.5 * ( tEdge->node( 0 )->y() + tEdge->node( 1 )->y() ) ;
            real tZ = 0.5 * ( tEdge->node( 0 )->z() + tEdge->node( 1 )->z() ) ;

            Node * tNode = new Node( ++tNodeID, tX, tY, tZ );
            tNodes( tCount++ ) = tNode ;
        }
    }

    Cell< Node * > tMyNodes( 10, nullptr );
    tCount = 0 ;
    for( uint k=0; k<4; ++k )
    {
        tMyNodes( tCount++ ) = tElement->node( k );
    }
    for( uint k=0; k<6; ++k )
    {
        tMyNodes( tCount++ ) = tNodes( tNumNodes + tElement->edge( k )->index() );
    }

    Block * tBlock1 = new Block( 2, 4 );


    for( uint e=0; e<4; ++e )
    {
        Element * tNewElement = tFactory.create_element(  ElementType::TET4, ++tElementID );
        for( uint k=0; k<4; ++k )
        {
            tNewElement->insert_node( tMyNodes( tTetTopo( e, k )), k );
        }
        tBlock1->elements()( e )  = tNewElement ;
    }

    tMesh->blocks().push( tBlock1 );

    Block * tBlock2 = new Block( 3, 2 );
    for( uint e=0; e<2; ++e )
    {
        Element * tNewElement = tFactory.create_element(  ElementType::PYRA5, ++tElementID );
        for( uint k=0; k<5; ++k )
        {
            tNewElement->insert_node( tMyNodes( tPyraTopo( e, k )), k );
        }
        tBlock2->elements()( e )  = tNewElement ;
    }
    tMesh->blocks().push( tBlock2 );
    tMesh->finalize() ;

    tMesh->save( "tet.exo");

    delete tMesh ;

    return  gComm.finalize();
}
