//
// Created by Christian Messe on 05.04.21.
//

#include "cl_TensorMeshFactory.hpp"
#include "assert.hpp"
#include "fn_linspace.hpp"
#include "cl_Element_Factory.hpp"
#include "commtools.hpp"

namespace belfem
{
    Mesh *
    TensorMeshFactory::create_tensor_mesh(
            const Vector< uint > & aNumElems,
            const Vector< real > & aMinPoint,
            const Vector< real > & aMaxPoint,
            const uint aOrder )
    {
        BELFEM_ERROR( aOrder <= 2, "Only 1st and 2nd order tensor meshes are implemented" );

        uint tDimension = aNumElems.length() ;

        BELFEM_ERROR( aMinPoint.length() == tDimension,
                     "dimension of min point does not match");
        BELFEM_ERROR( aMaxPoint.length() == tDimension,
                     "dimension of max point does not match");
        BELFEM_ERROR( tDimension == 2 || tDimension == 3,
                     "dimension mist be 2 or 3");


        Mesh * aMesh = new Mesh( tDimension, 0 );

        if( comm_rank() == 0 )
        {
            if( aOrder == 1 )
            {
                if( tDimension == 2 )
                {
                    this->create_2d_mesh_linear( aNumElems, aMinPoint, aMaxPoint, aMesh );
                }
                else if ( tDimension == 3 )
                {
                    this->create_3d_mesh_linear( aNumElems, aMinPoint, aMaxPoint, aMesh );
                }
            }
            else if ( aOrder == 2 )
            {
                if( tDimension == 2 )
                {
                    this->create_2d_mesh_quadratic( aNumElems, aMinPoint, aMaxPoint, aMesh );
                }
                else if ( tDimension == 3 )
                {
                    this->create_3d_mesh_linear( aNumElems, aMinPoint, aMaxPoint, aMesh );
                }
            }

            aMesh->finalize() ;
        }

        return aMesh ;
    }

//------------------------------------------------------------------------------

    void
    TensorMeshFactory::create_2d_mesh_linear(
            const Vector< uint > & aNumElems,
            const Vector< real > & aMinPoint,
            const Vector< real > & aMaxPoint,
                  Mesh           * aMesh )
    {
        index_t tNumElemsI = aNumElems( 0 );
        index_t tNumElemsJ = aNumElems( 1 );

        uint tNumNodesI = tNumElemsI + 1 ;
        uint tNumNodesJ = tNumElemsJ + 1 ;

        // create node coordinates
        Vector< real > tX( tNumNodesI );
        Vector< real > tY( tNumNodesJ );

        linspace( aMinPoint( 0 ), aMaxPoint( 0 ), tNumNodesI, tX );
        linspace( aMinPoint( 1 ), aMaxPoint( 1 ), tNumNodesJ, tY );

        // get node container
        Cell< mesh::Node * > & tNodes = aMesh->nodes() ;

        // allocate node container
        tNodes.set_size( tNumNodesI * tNumNodesJ, nullptr );

        // initialize counter
        id_t tCount = 0 ;

        // create nodes
        for( index_t j=0; j<tNumNodesJ; ++j )
        {
            for ( index_t i=0; i<tNumNodesI; ++i )
            {
                tNodes( tCount ) = new mesh::Node( tCount+1, tX( i ), tY( j ) );
                ++tCount ;
            }
        }

        // get block container
        Cell< mesh::Block * > & tBlocks = aMesh->blocks() ;
        tBlocks.set_size( 1, nullptr );

        // create a block
        mesh::Block * tBlock = new mesh::Block( 1,
                                            tNumElemsI * tNumElemsJ );

        // link block
        tBlocks( 0 ) = tBlock ;

        // create element factory
        mesh::ElementFactory tFactory ;

        // reset counter
        tCount = 0 ;
        index_t tOff = 0 ;
        for( index_t j=0; j<tNumElemsJ; ++j )
        {
            for ( index_t i=0; i<tNumElemsI; ++i )
            {
                // create element
                mesh::Element * tElement = tFactory.create_element( ElementType::QUAD4, tCount + 1 );

                // link element with nodes
                tElement->insert_node( tNodes( tOff + i ), 0 );
                tElement->insert_node( tNodes( tOff + i + 1 ), 1 );
                tElement->insert_node( tNodes(tOff + i + 1 + tNumNodesI ), 2 );
                tElement->insert_node( tNodes( tOff + i + tNumNodesI ), 3 );

                // add element to container
                tBlock->insert_element( tElement );
            }

            // increment offset
            tOff += tNumNodesI ;
        }
    }

//------------------------------------------------------------------------------

    void
    TensorMeshFactory::create_2d_mesh_quadratic(
            const Vector< uint > & aNumElems,
            const Vector< real > & aMinPoint,
            const Vector< real > & aMaxPoint,
            Mesh           * aMesh )
    {
        index_t tNumElemsI = aNumElems( 0 );
        index_t tNumElemsJ = aNumElems( 1 );

        uint tNumNodesI = 2 * tNumElemsI + 1 ;
        uint tNumNodesJ = 2 * tNumElemsJ + 1 ;

        // create node coordinates
        Vector< real > tX( tNumNodesI );
        Vector< real > tY( tNumNodesJ );

        linspace( aMinPoint( 0 ), aMaxPoint( 0 ), tNumNodesI, tX );
        linspace( aMinPoint( 1 ), aMaxPoint( 1 ), tNumNodesJ, tY );

        // get node container
        Cell< mesh::Node * > & tNodes = aMesh->nodes() ;

        // allocate node container
        tNodes.set_size( tNumNodesI * tNumNodesJ, nullptr );

        // initialize counter
        id_t tCount = 0 ;

        // create nodes
        for( index_t j=0; j<tNumNodesJ; ++j )
        {
            for ( index_t i=0; i<tNumNodesI; ++i )
            {
                tNodes( tCount ) = new mesh::Node( tCount+1, tX( i ), tY( j ) );
                ++tCount ;
            }
        }

        // get block container
        Cell< mesh::Block * > & tBlocks = aMesh->blocks() ;
        tBlocks.set_size( 1, nullptr );

        // create a block
        mesh::Block * tBlock = new mesh::Block( 1,
                                                tNumElemsI * tNumElemsJ );

        // link block
        tBlocks( 0 ) = tBlock ;

        // create element factory
        mesh::ElementFactory tFactory ;

        // reset counter
        tCount = 0 ;
        index_t tOff = 0 ;
        for( index_t j=0; j<tNumElemsJ; ++j )
        {
            for ( index_t i=0; i<tNumElemsI; ++i )
            {
                // create element
                mesh::Element * tElement = tFactory.create_element( ElementType::QUAD9, tCount + 1 );

                // link element with nodes
                tElement->insert_node( tNodes( tOff + 2*i ), 0 );
                tElement->insert_node( tNodes( tOff + 2*i + 2 ), 1 );
                tElement->insert_node( tNodes( tOff + 2*i + 2 * tNumNodesI + 2 ), 2 );
                tElement->insert_node( tNodes( tOff + 2*i + 2 * tNumNodesI ), 3 );
                tElement->insert_node( tNodes( tOff + 2*i + 1 ), 4 );
                tElement->insert_node( tNodes( tOff + 2*i + tNumNodesI + 2 ), 5 );
                tElement->insert_node( tNodes( tOff + 2*i + 2 * tNumNodesI + 1 ), 6 );
                tElement->insert_node( tNodes( tOff + 2*i + tNumNodesI ), 7 );
                tElement->insert_node( tNodes( tOff + 2*i + tNumNodesI+1 ), 8 );

                // add element to container
                tBlock->insert_element( tElement );
            }

            // increment offset
            tOff += 2 * tNumNodesI ;
        }
    }


//------------------------------------------------------------------------------

    void
    TensorMeshFactory::create_3d_mesh_linear(
            const Vector< uint > & aNumElems,
            const Vector< real > & aMinPoint,
            const Vector< real > & aMaxPoint,
                  Mesh           * aMesh )
    {
        index_t tNumElemsI = aNumElems( 0 );
        index_t tNumElemsJ = aNumElems( 1 );
        index_t tNumElemsK = aNumElems( 2 );

        uint tNumNodesI = tNumElemsI + 1 ;
        uint tNumNodesJ = tNumElemsJ + 1 ;
        uint tNumNodesK = tNumElemsK + 1 ;

        // create node coordinates
        Vector< real > tX( tNumNodesI );
        Vector< real > tY( tNumNodesJ );
        Vector< real > tZ( tNumNodesK );

        linspace( aMinPoint( 0 ), aMaxPoint( 0 ), tNumNodesI, tX );
        linspace( aMinPoint( 1 ), aMaxPoint( 1 ), tNumNodesJ, tY );
        linspace( aMinPoint( 2 ), aMaxPoint( 2 ), tNumNodesK, tZ );

        // get node container
        Cell< mesh::Node * > & tNodes = aMesh->nodes() ;

        // allocate node container
        tNodes.set_size( tNumNodesI * tNumNodesJ * tNumNodesK, nullptr );

        // initialize counter
        id_t tCount = 0 ;

        // create nodes
        for( index_t k=0; k<tNumNodesK; ++k )
        {
            for( index_t j=0; j<tNumNodesJ; ++j )
            {
                for ( index_t i = 0; i < tNumNodesI; ++i )
                {
                    tNodes( tCount ) = new mesh::Node( tCount + 1, tX( i ), tY( j ), tZ( k ) );
                    ++tCount;
                }
            }
        }

        // get block container
        Cell< mesh::Block * > & tBlocks = aMesh->blocks() ;
        tBlocks.set_size( 1, nullptr );

        // create a block
        mesh::Block * tBlock = new mesh::Block( 1,
                                            tNumElemsI * tNumElemsJ * tNumElemsK );

        // link block
        tBlocks( 0 ) = tBlock ;

        // create element factory
        mesh::ElementFactory tFactory ;

        // reset counter
        tCount = 0 ;
        index_t tOff ;
        index_t tOff1;
        index_t tOff2 = 0 ;
        index_t tNumNodesIJ = tNumNodesI * tNumNodesJ ;

        for( index_t k=0; k<tNumElemsK; ++k )
        {

            tOff1 = 0 ;
            for( index_t j=0; j<tNumElemsJ; ++j )
            {

                for ( index_t i=0; i<tNumElemsI; ++i )
                {
                    // create element
                    mesh::Element * tElement = tFactory.create_element( ElementType::HEX8,
                                                                                 tCount + 1 );

                    // compute offset
                    tOff = tOff2 + tOff1 ;

                    // link element with nodes
                    tElement->insert_node( tNodes( tOff + i ), 0 );
                    tElement->insert_node( tNodes( tOff + i + 1 ), 1 );
                    tElement->insert_node( tNodes(tOff + i + 1 + tNumNodesI) , 2 );
                    tElement->insert_node( tNodes(tOff + i + tNumNodesI ), 3 );
                    tElement->insert_node( tNodes( tOff + i + tNumNodesIJ ), 4 );
                    tElement->insert_node( tNodes( tOff + i + 1 + tNumNodesIJ ), 5 );
                    tElement->insert_node( tNodes(tOff + i + 1 + tNumNodesI + tNumNodesIJ ), 6 );
                    tElement->insert_node( tNodes(tOff + i + tNumNodesI + tNumNodesIJ ), 7 );

                    // add element to container
                    tBlock->insert_element( tElement );
                }

                // increment offset
                tOff1 += tNumNodesI ;
            }
            tOff2 += tNumNodesIJ ;
        }
    }

//------------------------------------------------------------------------------

    void
    TensorMeshFactory::create_3d_mesh_quadratic(
            const Vector< uint > & aNumElems,
            const Vector< real > & aMinPoint,
            const Vector< real > & aMaxPoint,
            Mesh           * aMesh )
    {
        index_t tNumElemsI = aNumElems( 0 );
        index_t tNumElemsJ = aNumElems( 1 );
        index_t tNumElemsK = aNumElems( 2 );

        uint tNumNodesI = 2*tNumElemsI + 1 ;
        uint tNumNodesJ = 2*tNumElemsJ + 1 ;
        uint tNumNodesK = 2*tNumElemsK + 1 ;

        // create node coordinates
        Vector< real > tX( tNumNodesI );
        Vector< real > tY( tNumNodesJ );
        Vector< real > tZ( tNumNodesK );

        linspace( aMinPoint( 0 ), aMaxPoint( 0 ), tNumNodesI, tX );
        linspace( aMinPoint( 1 ), aMaxPoint( 1 ), tNumNodesJ, tY );
        linspace( aMinPoint( 2 ), aMaxPoint( 2 ), tNumNodesJ, tZ );

        // get node container
        Cell< mesh::Node * > & tNodes = aMesh->nodes() ;

        // allocate node container
        tNodes.set_size( tNumNodesI * tNumNodesJ * tNumNodesK, nullptr );

        // initialize counter
        id_t tCount = 0 ;

        // create nodes
        for( index_t k=0; k<tNumNodesK; ++k )
        {
            for( index_t j=0; j<tNumNodesJ; ++j )
            {
                for ( index_t i = 0; i < tNumNodesI; ++i )
                {
                    tNodes( tCount ) = new mesh::Node( tCount + 1, tX( i ), tY( j ), tZ( k ) );
                    ++tCount;
                }
            }
        }

        // get block container
        Cell< mesh::Block * > & tBlocks = aMesh->blocks() ;
        tBlocks.set_size( 1, nullptr );

        // create a block
        mesh::Block * tBlock = new mesh::Block( 1,
                                                tNumElemsI * tNumElemsJ * tNumElemsK );

        // link block
        tBlocks( 0 ) = tBlock ;

        // create element factory
        mesh::ElementFactory tFactory ;

        // reset counter
        tCount = 0 ;
        index_t tOff ;
        index_t tOff1;
        index_t tOff2 = 0 ;
        index_t tNumNodesIJ = tNumNodesI * tNumNodesJ ;

        for( index_t k=0; k<tNumElemsK; ++k )
        {

            tOff1 = 0 ;
            for( index_t j=0; j<tNumElemsJ; ++j )
            {

                for ( index_t i=0; i<tNumElemsI; ++i )
                {
                    // create element
                    mesh::Element * tElement = tFactory.create_element( ElementType::HEX27,
                                                                                 tCount + 1 );

                    // compute offset
                    tOff = tOff2 + tOff1 ;

                    // link element with nodes
                    tElement->insert_node( tNodes( tOff + 2*i ), 0 );
                    tElement->insert_node( tNodes( tOff + 2*i + 2 ), 1 );
                    tElement->insert_node( tNodes( tOff + 2*i + 2 * tNumNodesI + 2 ), 2 );
                    tElement->insert_node( tNodes( tOff + 2*i + 2 * tNumNodesI ), 3 );
                    tElement->insert_node( tNodes( tOff + 2*i + 2 + 2*tNumNodesIJ ), 4 );
                    tElement->insert_node( tNodes( tOff + 2*i + 2 * tNumNodesI + 2 + 2*tNumNodesIJ ), 5 );
                    tElement->insert_node( tNodes( tOff + 2*i + 2 * tNumNodesI + 2*tNumNodesIJ ), 6 );
                    tElement->insert_node( tNodes( tOff + 2*i + 2*tNumNodesIJ ), 7 );
                    tElement->insert_node( tNodes( tOff + 2*i + 1 ), 8 );
                    tElement->insert_node( tNodes( tOff + 2*i + tNumNodesI + 2 ), 9 );
                    tElement->insert_node( tNodes( tOff + 2*i + 2 * tNumNodesI + 1 ), 10 );
                    tElement->insert_node( tNodes( tOff + 2*i + tNumNodesI ), 11 );
                    tElement->insert_node( tNodes( tOff + 2*i + 2 + tNumNodesIJ ), 12 );
                    tElement->insert_node( tNodes( tOff + 2*i + 2 * tNumNodesI + 2 + tNumNodesIJ ), 13 );
                    tElement->insert_node( tNodes( tOff + 2*i + 2 * tNumNodesI + tNumNodesIJ ), 14 );
                    tElement->insert_node( tNodes( tOff + 2*i + tNumNodesIJ ), 15 );
                    tElement->insert_node( tNodes( tOff + 2*i + tNumNodesI + 2 + 2*tNumNodesIJ ), 16 );
                    tElement->insert_node( tNodes( tOff + 2*i + 2 * tNumNodesI + 1 + 2*tNumNodesIJ ), 17 );
                    tElement->insert_node( tNodes( tOff + 2*i + tNumNodesI + 2*tNumNodesIJ ), 18);
                    tElement->insert_node( tNodes( tOff + 2*i + 1 + 2*tNumNodesIJ ), 19 );
                    tElement->insert_node( tNodes( tOff + 2*i + tNumNodesI+1 + tNumNodesIJ ), 20 );
                    tElement->insert_node( tNodes( tOff + 2*i + tNumNodesI+1 ), 21 );
                    tElement->insert_node( tNodes( tOff + 2*i + tNumNodesI+1 + 2*tNumNodesIJ ), 22 );
                    tElement->insert_node( tNodes( tOff + 2*i + tNumNodesI + tNumNodesIJ ), 23);
                    tElement->insert_node( tNodes( tOff + 2*i + tNumNodesI + 2 + tNumNodesIJ ), 24 );
                    tElement->insert_node( tNodes( tOff + 2*i + 1 + tNumNodesIJ ), 25 );
                    tElement->insert_node( tNodes( tOff + 2*i + 2 * tNumNodesI + 1 + tNumNodesIJ ), 26 );

                    // add element to container
                    tBlock->insert_element( tElement );
                }

                // increment offset
                tOff1 += 2*tNumNodesI ;
            }
            tOff2 += 2*tNumNodesIJ ;
        }
    }

//------------------------------------------------------------------------------
}