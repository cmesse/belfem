//
// Created by christian on 10/12/21.
//
#include "commtools.hpp"
#include "fn_FEM_compute_element_current.hpp"
#include "cl_IWG_Maxwell.hpp"
#include "cl_FEM_Kernel.hpp"

namespace belfem
{
    namespace fem
    {
//----------------------------------------------------------------------------

        void
        compute_element_current( DofManager * aField )
        {
            // get the mesh
            Mesh * tMesh = aField->mesh() ;

            // get the equation object
            IWG_Maxwell * tFormulation = reinterpret_cast< IWG_Maxwell * >( aField->iwg() );

            index_t tCount = 0 ;

            // count elements that belong to this block
            for( id_t tID : tFormulation->current_blocks() )
            {
                tCount += aField->block( tID )->number_of_elements() ;
            }

            // Vectors for the element IDs
            Vector< id_t > tMyIDs( tCount );

            Vector< real > tMyJx;
            Vector< real > tMyJy ;
            Vector< real > tMyJz ;

            if( tMesh->number_of_dimensions() == 2 )
            {
                // allocate memory
                tMyJz.set_size( tCount );

                // reset counter
                tCount = 0 ;

                // loop over all blocks
                for( id_t tID : tFormulation->current_blocks() )
                {
                    // get block
                    Block * tBlock = aField->block( tID );
                    tFormulation->link_to_group( tBlock );

                    // loop over all elements of this block
                    for( Element * tElement : tBlock->elements() )
                    {
                        // compute current and write result into field on mesh
                        tFormulation->compute_element_current_2d( tElement, tMyJz( tCount ) );

                        // write element ID into field
                        tMyIDs( tCount++ ) = tElement->id() ;
                    }
                }

            }
            else if ( tMesh->number_of_dimensions() == 3 )
            {
                // allocate memory
                tMyJx.set_size( tCount );
                tMyJy.set_size( tCount );
                tMyJz.set_size( tCount );

                // reset counter
                tCount = 0 ;

                // loop over all blocks
                for( id_t tID : tFormulation->current_blocks() )
                {
                    // get block
                    Block * tBlock = aField->block( tID );
                    tFormulation->link_to_group( tBlock );

                    // loop over all elements of this block
                    for( Element * tElement : tBlock->elements() )
                    {
                        // compute current and write result into field on mesh
                        tFormulation->compute_element_current_3d( tElement,
                                                                  tMyJx( tCount ), tMyJy( tCount ), tMyJz( tCount ) );

                        // write element ID into field
                        tMyIDs( tCount++ ) = tElement->id() ;
                    }
                }
            }

            // wait for all procs
            comm_barrier();

            if( aField->is_master() )
            {
                // write my own data into the mesh
                Cell< Vector< id_t > > tAllIDs ;


                receive( aField->parent()->comm_table(), tAllIDs );

                if( tMesh->number_of_dimensions() == 2 )
                {
                    Cell< Vector< real > > tAllJz ;
                    receive( aField->parent()->comm_table(), tAllJz );

                    // get fields from mesh
                    Vector< real > & tJz = tMesh->field_exists("ElementCurrent") ?
                                           tMesh->field_data( "ElementCurrent") :
                                           tMesh->create_field("ElementCurrent" , EntityType::ELEMENT );

                    // loop over all procs
                    for( proc_t p : aField->parent()->comm_table() )
                    {
                        Vector< id_t > & tIDs = p == aField->rank() ? tMyIDs : tAllIDs( p );
                        Vector< real > & tPJz = p == aField->rank() ? tMyJz : tAllJz( p );

                        // reset counter
                        tCount = 0 ;

                        // write data onto mesh
                        for( id_t tID : tIDs )
                        {
                            tJz( tMesh->element( tID )->index() ) = tPJz( tCount++ );
                        }
                    }
                }
                else if( tMesh->number_of_dimensions() == 3 )
                {
                    Cell< Vector< real > > tAllJx ;
                    Cell< Vector< real > > tAllJy ;
                    Cell< Vector< real > > tAllJz ;

                    receive( aField->parent()->comm_table(), tAllJx );
                    receive( aField->parent()->comm_table(), tAllJy );
                    receive( aField->parent()->comm_table(), tAllJz );

                    // get fields from mesh
                    Vector< real > & tJx = tMesh->field_exists("ElementCurrentx") ?
                                           tMesh->field_data( "ElementCurrentx") :
                                           tMesh->create_field("ElementCurrentx" , EntityType::ELEMENT );

                    // get fields from mesh
                    Vector< real > & tJy = tMesh->field_exists("ElementCurrenty") ?
                                           tMesh->field_data( "ElementCurrenty") :
                                           tMesh->create_field("ElementCurreny" , EntityType::ELEMENT );

                    // get fields from mesh
                    Vector< real > & tJz = tMesh->field_exists("ElementCurrentz") ?
                                           tMesh->field_data( "ElementCurrentz") :
                                           tMesh->create_field("ElementCurrentz" , EntityType::ELEMENT );

                    // loop over all procs
                    for( proc_t p : aField->parent()->comm_table() )
                    {
                        Vector< id_t > & tIDs = p == aField->rank() ? tMyIDs : tAllIDs( p );
                        Vector< real > & tPJx = p == aField->rank() ? tMyJx : tAllJx( p );
                        Vector< real > & tPJy = p == aField->rank() ? tMyJy : tAllJy( p );
                        Vector< real > & tPJz = p == aField->rank() ? tMyJz : tAllJz( p );

                        // reset counter
                        tCount = 0 ;
                        index_t tIndex;

                        // write data onto mesh
                        for( id_t tID : tIDs )
                        {
                            tIndex = tMesh->element( tID )->index();
                            tJx( tIndex  ) = tPJx( tCount );
                            tJy( tIndex  ) = tPJy( tCount );
                            tJz( tIndex  ) = tPJz( tCount );
                            ++tCount ;
                        }
                    }
                }
            }
            else
            {
                // send data to master
                send( aField->parent()->master(), tMyIDs );

                // send currents
                if( tMesh->number_of_dimensions() == 2 )
                {
                    send( aField->parent()->master(), tMyJz );
                }
                else if( tMesh->number_of_dimensions() == 3 )
                {
                    send( aField->parent()->master(), tMyJx );
                    send( aField->parent()->master(), tMyJy );
                    send( aField->parent()->master(), tMyJz );
                }

            }

            // wait for all procs
            comm_barrier();
        }

//----------------------------------------------------------------------------
    }
}