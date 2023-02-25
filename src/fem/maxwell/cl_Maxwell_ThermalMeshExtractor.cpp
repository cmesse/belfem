//
// Created by christian on 2/7/23.
//
#include "cl_Maxwell_ThermalMeshExtractor.hpp"
#include "commtools.hpp"
#include "assert.hpp"
#include "cl_Pipette.hpp"
#include "fn_unique.hpp"
#include "cl_Element_Factory.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        ThermalMeshExtractor::ThermalMeshExtractor(
                Mesh * aInputMesh ) :
                mInputMesh( * aInputMesh ),
                mGhostBlockIDs( aInputMesh->ghost_block_ids() )
        {

        }

//-----------------------------------------------------------------------------

        Mesh *
        ThermalMeshExtractor::create_mesh()
        {
            if( comm_rank() == 0 )
            {
                if ( mInputMesh.number_of_dimensions() == 2 )
                {
                    return this->create_2d_mesh();
                }
                else
                {
                    return this->create_3d_mesh();
                }
            }
            else
            {
                return new Mesh( mInputMesh.number_of_dimensions(), mInputMesh.master() );
            }
        }

//-----------------------------------------------------------------------------

        Mesh *
        ThermalMeshExtractor::create_2d_mesh()
        {
            // create a new mesh object
            Mesh * aMesh = new Mesh( mInputMesh.number_of_dimensions(),
                                     mInputMesh.master() ) ;

            // create the nodes
            this->create_nodes( aMesh );

            // create resistor objects
            this->create_resistors( aMesh );

            // finalize the mesh
            aMesh->finalize() ;

            // compute the volumes ( not needed in 2D )
            //this->compute_element_volumes( aMesh );

            // copy element ownerships from input mesh
            // onto node ownerships from output mesh
            this->set_ownerships( &mInputMesh, aMesh );

            return aMesh ;
        }

//-----------------------------------------------------------------------------

        Mesh *
        ThermalMeshExtractor::create_3d_mesh()
        {
            BELFEM_ERROR( false, "3D feature not imlemented yet!" );

            return nullptr ;
        }

//-----------------------------------------------------------------------------

        void
        ThermalMeshExtractor::create_nodes( Mesh * aMesh )
        {
            // count nodes on new mesh ( also elements on old mesh )
            index_t tCount = 0 ;

            // loop over all elements on original mesh
            for( id_t b: mGhostBlockIDs )
            {
                tCount += mInputMesh.block( b )->number_of_elements() ;
            }

            // allocate node container
            Cell< mesh::Node * > & tNodes = aMesh->nodes() ;
            tNodes.set_size( tCount, nullptr );

            // reset counter
            tCount = 0 ;

            //uint tNumLayers = mGhostBlockIDs.length() ;
            //uint tCount = 0 ;

            // loop over all elements on original mesh
            for( id_t b: mGhostBlockIDs )
            {
                // grab all elements on original block
                Cell< mesh::Element * > & tElements = mInputMesh.block( b )->elements() ;

                // check how many nodes exist per element
                uint tNumNodesPerElement = mesh::number_of_nodes( mInputMesh.block( b )->element_type() );

                // loop over all elements on this block
                for( mesh::Element * tElement : tElements )
                {
                    // node coordinates
                    real tX = 0 ;
                    real tY = 0 ;
                    real tZ = 0 ;

                    for( uint k = 0 ; k<tNumNodesPerElement; ++k )
                    {
                        tX += tElement->node( k )->x() ;
                        tY += tElement->node( k )->y() ;
                        tZ += tElement->node( k )->z() ;
                    }

                    tX /= tNumNodesPerElement ;
                    tY /= tNumNodesPerElement ;
                    tZ /= tNumNodesPerElement ;

                    // add index to map
                    mNodeIndexMap[ tElement->id() ] = tCount ;

                    // create the new node and add node to container
                    tNodes( tCount++ ) = new mesh::Node( tElement->id(), tX, tY, tZ );
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        ThermalMeshExtractor::compute_element_volumes( Mesh * aMesh )
        {
            // create volume field
            Vector< real > & tVolumes = aMesh->create_field( "Volumes" );

            // create pipetter
            mesh::Pipette tPipette ;
            // loop over all elements on original mesh
            for( id_t b: mGhostBlockIDs )
            {
                // grab all elements on original block
                Cell< mesh::Element * > & tElements = mInputMesh.block( b )->elements() ;

                // set the element tupe
                tPipette.set_element_type( mInputMesh.block( b )->element_type() );

                // loop over all elements on this block
                for( mesh::Element * tElement : tElements )
                {
                    std::cout << tElement->id() << " " << tPipette.measure( tElement ) << std::endl ;

                    tVolumes( aMesh->node( tElement->id() )->index() )
                        = tPipette.measure( tElement );
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        ThermalMeshExtractor::create_resistors( Mesh * aMesh )
        {
            mInputMesh.unflag_all_elements() ;
            mInputMesh.unflag_all_nodes() ;

            uint tNumDim = mInputMesh.number_of_dimensions() ;

            // - - - - - - - - - - - - - - - - - - - - - - -
            // step 1: count how many connectors exist
            // - - - - - - - - - - - - - - - - - - - - - - -

            index_t tCount = 0 ;

            // loop over all elements on original mesh
            for( id_t b: mGhostBlockIDs )
            {
                // grab all elements on original block
                Cell< mesh::Element * > & tElements = mInputMesh.block( b )->elements() ;

                // loop over all elements on this block
                for( mesh::Element * tElement : tElements )
                {
                    // flag all nodes of this element
                    tElement->flag_nodes() ;

                    for( uint e = 0 ; e<tElement->number_of_elements() ; ++e )
                    {
                        // grab element
                        mesh::Element * tOther = tElement->element( e );

                        // init node counter
                        uint tNodeCount = 0 ;

                        // loop over all corner nodes of other element
                        for( uint k=0 ; k<tOther->number_of_corner_nodes(); ++k )
                        {
                            if( tOther->node( k )->is_flagged() )
                            {
                                // increment node counter
                                ++tNodeCount ;
                            }
                        }

                        if( tNodeCount >= tNumDim )
                        {
                            ++tCount ;
                        }
                    }

                    // unflag all nodes of this element
                    tElement->unflag_nodes() ;
                }
            }

            // - - - - - - - - - - - - - - - - - - - - - - -
            // step 2: populate resistor IDs
            // - - - - - - - - - - - - - - - - - - - - - - -

            Vector< index_t > tIDs( tCount );

            // reset counter
            tCount = 0 ;

            // get number of nodes
            index_t tNumNodes = aMesh->number_of_nodes() ;

            // loop over all elements on original mesh
            for( id_t b: mGhostBlockIDs )
            {
                // grab all elements on original block
                Cell< mesh::Element * > & tElements = mInputMesh.block( b )->elements() ;

                // loop over all elements on this block
                for( mesh::Element * tElement : tElements )
                {
                    // flag all nodes of this element
                    tElement->flag_nodes() ;

                    // index of first node on new mesh
                    index_t tIndexA = mNodeIndexMap( tElement->id() );

                    for( uint e = 0 ; e<tElement->number_of_elements() ; ++e )
                    {
                        // grab element
                        mesh::Element * tOther = tElement->element( e );

                        index_t tIndexB = mNodeIndexMap( tOther->id() );

                        // init node counter
                        uint tNodeCount = 0 ;

                        // loop over all corner nodes of other element
                        for( uint k=0 ; k<tOther->number_of_corner_nodes(); ++k )
                        {
                            if( tOther->node( k )->is_flagged() )
                            {
                                // increment node counter
                                ++tNodeCount ;
                            }
                        }

                        if( tNodeCount >= tNumDim )
                        {
                            // compute an ID
                            if( tIndexA > tIndexB )
                            {
                                tIDs( tCount++ ) = tNumNodes * tIndexA + tIndexB ;
                            }
                            else
                            {
                                tIDs( tCount++ ) = tNumNodes * tIndexB + tIndexA ;
                            }
                        }
                    }

                    // unflag all nodes of this element
                    tElement->unflag_nodes() ;
                }
            }

            // make IDs unique
            unique( tIDs );

            // create a block on output mesh
            mesh::Block * tBlock = new mesh::Block( 1 , tIDs.length() );

            // grab the nodes from the mesh
            Cell< mesh::Node * > & tNodes = aMesh->nodes() ;

            // grab element container from block
            Cell< mesh::Element * > & tElements = tBlock->elements() ;

            // create an element factory
            mesh::ElementFactory tFactory ;

            // reset the counter
            tCount = 0 ;

            // create container on mesh
            for( index_t tID : tIDs )
            {
                index_t tB = tID % tNumNodes ;
                index_t tA = ( tID - tB ) / tNumNodes ;

                // create the element
                mesh::Element * tElement = tFactory.create_lagrange_element( ElementType::LINE2, tCount+1 );

                // connect element with nodes
                tElement->insert_node( tNodes( tA ), 0 );
                tElement->insert_node( tNodes( tB ), 1 );

                // add element to container
                tElements( tCount++ ) = tElement ;
            }

            // add the block to the mesh
            aMesh->blocks().set_size( 1, nullptr );
            aMesh->blocks()( 0 ) = tBlock ;
        }

//-----------------------------------------------------------------------------

        void
        ThermalMeshExtractor::set_ownerships(  Mesh * aInputMesh,  Mesh * aOutputMesh )
        {
            if( comm_rank() == 0 )
            {
                Cell< mesh::Node * > & tNodes = aOutputMesh->nodes();

                // get node owners from input mesh
                for ( mesh::Node * tNode: tNodes )
                {
                    tNode->set_owner( aInputMesh->element( tNode->id())->owner());
                }

                // generate element owners
                Cell< mesh::Element * > & tElements = aOutputMesh->elements();
                proc_t tNumProcs = comm_size();
                for ( mesh::Element * tElement: tElements )
                {
                    proc_t tOwner = tNumProcs;
                    for ( uint k = 0; k < tElement->number_of_nodes(); ++k )
                    {
                        if ( tElement->node( k )->owner() < tOwner )
                        {
                            tOwner = tElement->node( k )->owner();
                        }
                        tElement->set_owner( tOwner );
                    }
                }
            }
            comm_barrier();
        }

//-----------------------------------------------------------------------------
    }
}
