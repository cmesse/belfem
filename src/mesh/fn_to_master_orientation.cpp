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

#include "fn_to_master_orientation.hpp"
#include "assert.hpp"
namespace belfem
{
    namespace mesh
    {
        void
        to_master_orientation( Facet * aFacet,
            Cell< Node * > & aSlaveOrientation,
            Cell< Node * > & aMasterOrientation )
        {
            BELFEM_ASSERT( aSlaveOrientation.size() == aFacet->number_of_nodes(), "Container not populated" );

            BELFEM_ASSERT( aFacet->slave() != nullptr , "Facet %lu does not have a slave element",
                ( long unsigned int ) aFacet->id() );

            aMasterOrientation.set_size( aSlaveOrientation.size(), nullptr );

            switch ( aFacet->element()->type() )
            {
                case( ElementType::LINE2 ) :
                {
                    aMasterOrientation( 0 ) =  aSlaveOrientation( 1 );
                    aMasterOrientation( 1 ) =  aSlaveOrientation( 0 );
                    break ;
                }
                case( ElementType::LINE3 ) :
                {
                    aMasterOrientation( 0 ) =  aSlaveOrientation( 1 );
                    aMasterOrientation( 1 ) =  aSlaveOrientation( 0 );
                    aMasterOrientation( 2 ) =  aSlaveOrientation( 2 );
                    break ;
                }
                case( ElementType::TRI3 ) :
                {
                    switch ( aFacet->orientation_on_slave() )
                    {
                        case( 1 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 1 );
                            break ;
                        }
                        case( 2 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 2 );
                            break ;
                        }
                        case( 3 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 0 );
                            break ;
                        }

                        default:
                        {
                            BELFEM_ERROR( false ,
                                "Invalid orientation of facet %lu on slave element %lu ( master is %lu )",
                            ( long unsigned int ) aFacet->id(),
                            ( long unsigned int ) aFacet->slave()->id(),
                            ( long unsigned int ) aFacet->master()->id() );
                        }
                    }
                    break ;
                }
                case( ElementType::TRI6 ) :
                {
                    switch ( aFacet->orientation_on_slave() )
                    {
                        case( 1 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 5 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 4 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 3 );
                            break ;
                        }
                        case( 2 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 5 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 4 );
                            break ;
                        }
                        case( 3 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 4 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 5 );
                            break ;
                        }
                        default:
                        {
                            BELFEM_ERROR( false ,
                                "Invalid orientation of facet %lu on slave element %lu ( master is %lu )",
                            ( long unsigned int ) aFacet->id(),
                            ( long unsigned int ) aFacet->slave()->id(),
                            ( long unsigned int ) aFacet->master()->id() );
                        }
                    }
                    break ;
                }
                case( ElementType::QUAD4 ) :
                {
                    switch ( aFacet->orientation_on_slave() )
                    {
                        case( 1 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 1 );
                            break ;
                        }
                        case( 2 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 2 );
                            break ;
                        }
                        case( 3 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 3 );
                            break ;
                        }
                        case( 4 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 0 );
                            break ;
                        }
                        default:
                        {
                            BELFEM_ERROR( false ,
                                "Invalid orientation of facet %lu on slave element %lu ( master is %lu )",
                            ( long unsigned int ) aFacet->id(),
                            ( long unsigned int ) aFacet->slave()->id(),
                            ( long unsigned int ) aFacet->master()->id() );
                        }
                    }
                    break ;
                }
                case( ElementType::QUAD8 ) :
                {
                    switch ( aFacet->orientation_on_slave() )
                    {
                        case( 1 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 7 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 6 );
                            aMasterOrientation( 6 ) =  aSlaveOrientation( 5 );
                            aMasterOrientation( 7 ) =  aSlaveOrientation( 4 );
                            break ;
                        }
                        case( 2 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 4 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 7 );
                            aMasterOrientation( 6 ) =  aSlaveOrientation( 6 );
                            aMasterOrientation( 7 ) =  aSlaveOrientation( 5 );
                            break ;
                        }
                        case( 3 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 5 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 4 );
                            aMasterOrientation( 6 ) =  aSlaveOrientation( 7 );
                            aMasterOrientation( 7 ) =  aSlaveOrientation( 6 );
                            break ;
                        }
                        case( 4 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 6 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 5 );
                            aMasterOrientation( 6 ) =  aSlaveOrientation( 4 );
                            aMasterOrientation( 7 ) =  aSlaveOrientation( 7 );
                            break ;
                        }
                        default:
                        {
                            BELFEM_ERROR( false ,
                                "Invalid orientation of facet %lu on slave element %lu ( master is %lu )",
                            ( long unsigned int ) aFacet->id(),
                            ( long unsigned int ) aFacet->slave()->id(),
                            ( long unsigned int ) aFacet->master()->id() );
                        }
                    }
                    break ;
                }
                case( ElementType::QUAD9 ) :
                {
                    switch ( aFacet->orientation_on_slave() )
                    {
                        case( 1 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 7 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 6 );
                            aMasterOrientation( 6 ) =  aSlaveOrientation( 5 );
                            aMasterOrientation( 7 ) =  aSlaveOrientation( 4 );
                            aMasterOrientation( 8 ) =  aSlaveOrientation( 8 );
                            break ;
                        }
                        case( 2 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 4 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 7 );
                            aMasterOrientation( 6 ) =  aSlaveOrientation( 6 );
                            aMasterOrientation( 7 ) =  aSlaveOrientation( 5 );
                            aMasterOrientation( 8 ) =  aSlaveOrientation( 8 );
                            break ;
                        }
                        case( 3 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 5 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 4 );
                            aMasterOrientation( 6 ) =  aSlaveOrientation( 7 );
                            aMasterOrientation( 7 ) =  aSlaveOrientation( 6 );
                            aMasterOrientation( 8 ) =  aSlaveOrientation( 8 );
                            break ;
                        }
                        case( 4 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 6 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 5 );
                            aMasterOrientation( 6 ) =  aSlaveOrientation( 4 );
                            aMasterOrientation( 7 ) =  aSlaveOrientation( 7 );
                            aMasterOrientation( 8 ) =  aSlaveOrientation( 8 );
                            break ;
                        }
                        default:
                        {
                            BELFEM_ERROR( false ,
                                "Invalid orientation of facet %lu on slave element %lu ( master is %lu )",
                            ( long unsigned int ) aFacet->id(),
                            ( long unsigned int ) aFacet->slave()->id(),
                            ( long unsigned int ) aFacet->master()->id() );
                        }
                    }
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false ,
                               "Invalid Element type for facet %lu",
                           ( long unsigned int ) aFacet->id() );
                }
            }
        }

        void
        to_master_orientation( Face * aFace,
            Cell< Node * > & aSlaveOrientation,
            Cell< Node * > & aMasterOrientation )
        {
            BELFEM_ASSERT( aFace->slave() != nullptr , "Face %lu does not have a slave element",
                ( long unsigned int ) aFace->id() );

            BELFEM_ASSERT( aSlaveOrientation.size() == aFace->number_of_nodes(), "Container not populated" );

            aMasterOrientation.set_size( aSlaveOrientation.size(), nullptr );

            switch ( aSlaveOrientation.size() )
            {
                case( 3 ) : // tri3
                {
                    switch ( aFace->orientation_on_slave() )
                    {
                        case( 1 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 1 );
                            break ;
                        }
                        case( 2 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 2 );
                            break ;
                        }
                        case( 3 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 0 );
                            break ;
                        }

                        default:
                        {
                            BELFEM_ERROR( false ,
                                "Invalid orientation of face %lu on slave element",
                            ( long unsigned int ) aFace->id() );
                        }
                    }
                    break ;
                }
                case( 6 ) : // tri6
                {
                    switch ( aFace->orientation_on_slave() )
                    {
                        case( 1 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 5 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 4 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 3 );
                            break ;
                        }
                        case( 2 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 5 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 4 );
                            break ;
                        }
                        case( 3 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 4 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 5 );
                            break ;
                        }
                        default:
                        {
                            BELFEM_ERROR( false ,
                                "Invalid orientation of face %lu on slave element",
                            ( long unsigned int ) aFace->id() );
                        }
                    }
                    break ;
                }
                case( 4 ) : // quad4
                {
                    switch ( aFace->orientation_on_slave() )
                    {
                        case( 1 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 1 );
                            break ;
                        }
                        case( 2 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 2 );
                            break ;
                        }
                        case( 3 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 3 );
                            break ;
                        }
                        case( 4 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 0 );
                            break ;
                        }
                        default:
                        {
                            BELFEM_ERROR( false ,
                                "Invalid orientation of face %lu on slave element",
                            ( long unsigned int ) aFace->id() );
                        }
                    }
                    break ;
                }
                case( 8 ) : // quad8
                {
                    switch ( aFace->orientation_on_slave() )
                    {
                        case( 1 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 7 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 6 );
                            aMasterOrientation( 6 ) =  aSlaveOrientation( 5 );
                            aMasterOrientation( 7 ) =  aSlaveOrientation( 4 );
                            break ;
                        }
                        case( 2 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 4 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 7 );
                            aMasterOrientation( 6 ) =  aSlaveOrientation( 6 );
                            aMasterOrientation( 7 ) =  aSlaveOrientation( 5 );
                            break ;
                        }
                        case( 3 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 5 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 4 );
                            aMasterOrientation( 6 ) =  aSlaveOrientation( 7 );
                            aMasterOrientation( 7 ) =  aSlaveOrientation( 6 );
                            break ;
                        }
                        case( 4 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 6 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 5 );
                            aMasterOrientation( 6 ) =  aSlaveOrientation( 4 );
                            aMasterOrientation( 7 ) =  aSlaveOrientation( 7 );
                            break ;
                        }
                        default:
                        {
                            BELFEM_ERROR( false ,
                                "Invalid orientation of face %lu on slave element",
                            ( long unsigned int ) aFace->id() );
                        }
                    }
                    break ;
                }
                case( 9 ) : // quad9
                {
                    switch ( aFace->orientation_on_slave() )
                    {
                        case( 1 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 7 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 6 );
                            aMasterOrientation( 6 ) =  aSlaveOrientation( 5 );
                            aMasterOrientation( 7 ) =  aSlaveOrientation( 4 );
                            aMasterOrientation( 8 ) =  aSlaveOrientation( 8 );
                            break ;
                        }
                        case( 2 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 4 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 7 );
                            aMasterOrientation( 6 ) =  aSlaveOrientation( 6 );
                            aMasterOrientation( 7 ) =  aSlaveOrientation( 5 );
                            aMasterOrientation( 8 ) =  aSlaveOrientation( 8 );
                            break ;
                        }
                        case( 3 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 5 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 4 );
                            aMasterOrientation( 6 ) =  aSlaveOrientation( 7 );
                            aMasterOrientation( 7 ) =  aSlaveOrientation( 6 );
                            aMasterOrientation( 8 ) =  aSlaveOrientation( 8 );
                            break ;
                        }
                        case( 4 ) :
                        {
                            aMasterOrientation( 0 ) =  aSlaveOrientation( 3 );
                            aMasterOrientation( 1 ) =  aSlaveOrientation( 2 );
                            aMasterOrientation( 2 ) =  aSlaveOrientation( 1 );
                            aMasterOrientation( 3 ) =  aSlaveOrientation( 0 );
                            aMasterOrientation( 4 ) =  aSlaveOrientation( 6 );
                            aMasterOrientation( 5 ) =  aSlaveOrientation( 5 );
                            aMasterOrientation( 6 ) =  aSlaveOrientation( 4 );
                            aMasterOrientation( 7 ) =  aSlaveOrientation( 7 );
                            aMasterOrientation( 8 ) =  aSlaveOrientation( 8 );
                            break ;
                        }
                        default:
                        {
                            BELFEM_ERROR( false ,
                                "Invalid orientation of face %lu on slave element",
                            ( long unsigned int ) aFace->id() );
                        }
                    }
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false ,
                               "Invalid number of nodes for face %lu : %u",
                           ( long unsigned int ) aFace->id(),  ( unsigned int ) aFace->number_of_nodes() );
                }
            }
        }

        void
        to_master_orientation( Facet * aFacet,
               Cell< Edge * > & aSlaveOrientation,
               Cell< Edge * > & aMasterOrientation )
        {
            BELFEM_ASSERT( aFacet->slave() != nullptr , "Facet %lu does not have a slave element",
        ( long unsigned int ) aFacet->id() );

            aMasterOrientation.set_size( aSlaveOrientation.size() , nullptr );
             switch ( aSlaveOrientation.size() )
            {
                case 1 : // line ( 2D facet ) : only one edge, so the master
                         // ordering is the identity. The relative direction is
                         // not expressed here -- the caller reads it off the
                         // node order, which the Node overload has already put
                         // into master orientation.
                {
                    aMasterOrientation( 0 ) = aSlaveOrientation( 0 );
                    break ;
                }
                case 3 : // triangle
                {
                    switch( aFacet->orientation_on_slave() )
                    {
                        case 1 :
                        {
                            aMasterOrientation( 0 ) = aSlaveOrientation( 2 );
                            aMasterOrientation( 1 ) = aSlaveOrientation( 1 );
                            aMasterOrientation( 2 ) = aSlaveOrientation( 0 );
                            break ;
                        }
                        case 2 :
                        {
                            aMasterOrientation( 0 ) = aSlaveOrientation( 0 );
                            aMasterOrientation( 1 ) = aSlaveOrientation( 2 );
                            aMasterOrientation( 2 ) = aSlaveOrientation( 1 );
                            break ;
                        }
                        case 3 :
                        {
                            aMasterOrientation( 0 ) = aSlaveOrientation( 1 );
                            aMasterOrientation( 1 ) = aSlaveOrientation( 0 );
                            aMasterOrientation( 2 ) = aSlaveOrientation( 2 );
                            break ;
                        }
                        default:
                        {
                            BELFEM_ERROR( false ,
                                "Invalid orientation of facet %lu on slave element",
                            ( long unsigned int ) aFacet->id() );
                        }
                    }
                    break ;
                }
                case 4 :
                {
                    // quadrangle
                    switch( aFacet->orientation_on_slave() )
                    {
                        case 1 :
                        {
                            aMasterOrientation( 0 ) = aSlaveOrientation( 3 );
                            aMasterOrientation( 1 ) = aSlaveOrientation( 2 );
                            aMasterOrientation( 2 ) = aSlaveOrientation( 1 );
                            aMasterOrientation( 3 ) = aSlaveOrientation( 0 );
                            break ;
                        }
                        case 2 :
                        {
                            aMasterOrientation( 0 ) = aSlaveOrientation( 0 );
                            aMasterOrientation( 1 ) = aSlaveOrientation( 3 );
                            aMasterOrientation( 2 ) = aSlaveOrientation( 2 );
                            aMasterOrientation( 3 ) = aSlaveOrientation( 1 );
                            break ;
                        }
                        case 3 :
                        {
                            aMasterOrientation( 0 ) = aSlaveOrientation( 1 );
                            aMasterOrientation( 1 ) = aSlaveOrientation( 0 );
                            aMasterOrientation( 2 ) = aSlaveOrientation( 3 );
                            aMasterOrientation( 3 ) = aSlaveOrientation( 2 );
                            break ;
                        }
                        case 4 :
                        {
                            aMasterOrientation( 0 ) = aSlaveOrientation( 2 );
                            aMasterOrientation( 1 ) = aSlaveOrientation( 1 );
                            aMasterOrientation( 2 ) = aSlaveOrientation( 0 );
                            aMasterOrientation( 3 ) = aSlaveOrientation( 3 );
                            break ;
                        }
                        default:
                        {
                            BELFEM_ERROR( false ,
                                "Invalid orientation of facet %lu on slave element",
                            ( long unsigned int ) aFacet->id() );
                        }
                    }
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false ,
                        "Invalid number of edges for facet %lu : %u",
                    ( long unsigned int ) aFacet->id(),
                    ( unsigned int ) aSlaveOrientation.size() );
                }
            }

        }

        void
        to_master_orientation( Face * aFace ,
           Cell< Edge * > & aSlaveOrientation,
           Cell< Edge * > & aMasterOrientation )
        {
            BELFEM_ASSERT( aFace->slave() != nullptr , "Facet %lu does not have a slave element",
            ( long unsigned int ) aFace->id() );

            aMasterOrientation.set_size( aSlaveOrientation.size() , nullptr );
            switch ( aSlaveOrientation.size() )
            {
                case 3 : // triangle
                {
                    switch( aFace->orientation_on_slave() )
                    {
                        case 1 :
                        {
                            aMasterOrientation( 0 ) = aSlaveOrientation( 2 );
                            aMasterOrientation( 1 ) = aSlaveOrientation( 1 );
                            aMasterOrientation( 2 ) = aSlaveOrientation( 0 );
                            break ;
                        }
                        case 2 :
                        {
                            aMasterOrientation( 0 ) = aSlaveOrientation( 0 );
                            aMasterOrientation( 1 ) = aSlaveOrientation( 2 );
                            aMasterOrientation( 2 ) = aSlaveOrientation( 1 );
                            break ;
                        }
                        case 3 :
                        {
                            aMasterOrientation( 0 ) = aSlaveOrientation( 1 );
                            aMasterOrientation( 1 ) = aSlaveOrientation( 0 );
                            aMasterOrientation( 2 ) = aSlaveOrientation( 2 );
                            break ;
                        }
                        default:
                        {
                            BELFEM_ERROR( false ,
                                "Invalid orientation of face %lu on slave element",
                            ( long unsigned int ) aFace->id() );
                        }
                    }
                    break ;
                }
                case 4 :
                {
                    // quadrangle
                    switch( aFace->orientation_on_slave() )
                    {
                        case 1 :
                        {
                            aMasterOrientation( 0 ) = aSlaveOrientation( 3 );
                            aMasterOrientation( 1 ) = aSlaveOrientation( 2 );
                            aMasterOrientation( 2 ) = aSlaveOrientation( 1 );
                            aMasterOrientation( 3 ) = aSlaveOrientation( 0 );
                            break ;
                        }
                        case 2 :
                        {
                            aMasterOrientation( 0 ) = aSlaveOrientation( 0 );
                            aMasterOrientation( 1 ) = aSlaveOrientation( 3 );
                            aMasterOrientation( 2 ) = aSlaveOrientation( 2 );
                            aMasterOrientation( 3 ) = aSlaveOrientation( 1 );
                            break ;
                        }
                        case 3 :
                        {
                            aMasterOrientation( 0 ) = aSlaveOrientation( 1 );
                            aMasterOrientation( 1 ) = aSlaveOrientation( 0 );
                            aMasterOrientation( 2 ) = aSlaveOrientation( 3 );
                            aMasterOrientation( 3 ) = aSlaveOrientation( 2 );
                            break ;
                        }
                        case 4 :
                        {
                            aMasterOrientation( 0 ) = aSlaveOrientation( 2 );
                            aMasterOrientation( 1 ) = aSlaveOrientation( 1 );
                            aMasterOrientation( 2 ) = aSlaveOrientation( 0 );
                            aMasterOrientation( 3 ) = aSlaveOrientation( 3 );
                            break ;
                        }
                        default:
                        {
                            BELFEM_ERROR( false ,
                                "Invalid orientation of face %lu on slave element",
                            ( long unsigned int ) aFace->id() );
                        }
                    }
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false ,
                        "Invalid number of edges for face %lu",
                    ( long unsigned int ) aFace->id() );
                }
            }
        }
    }
}
