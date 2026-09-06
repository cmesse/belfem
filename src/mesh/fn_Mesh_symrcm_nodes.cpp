/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#include "fn_Mesh_symrcm_nodes.hpp"
#include "fn_Graph_symrcm.hpp"
#include "op_Node_Index.hpp"
#include "cl_DynamicBitset.hpp"
namespace belfem
{
    namespace mesh
    {
        void
        symrcm( Cell< Node * > & aNodes )
        {
            Graph tGraph( aNodes.size() , nullptr );

            index_t tCount = 0 ;
            for ( Node * tNode : aNodes )
            {
                tNode->set_index( tCount );
                tGraph( tCount++ ) = tNode ;
                tNode->unflag() ;
            }

            DynamicBitset * tBitset = new DynamicBitset( aNodes.size() );
            Cell< index_t > tIndices ;
            for ( Node * tNode : aNodes )
            {
                Node * tOrg = tNode->original();

                tBitset->reset();
                for ( uint k=0; k<tOrg->number_of_nodes(); ++k )
                {
                    tBitset->set( tOrg->node(k)->index() );
                }
                for ( uint d=0; d<tOrg->number_of_duplicates(); ++d )
                {
                    Node * tDup = tOrg->duplicate( d );
                    for ( uint k=0; k<tDup->number_of_nodes(); ++k )
                    {
                        tBitset->set( tDup->node(k)->index() );
                    }
                }
                if ( tOrg->is_periodic() )
                {
                    Node * tPeriodic = tOrg->periodic();

                    for ( uint k=0; k<tPeriodic->number_of_nodes(); ++k )
                    {
                        tBitset->set( tPeriodic->node(k)->index() );
                    }

                    for ( uint d=0; d<tPeriodic->number_of_duplicates(); ++d )
                    {
                        Node * tDup = tPeriodic->duplicate( d );
                        for ( uint k=0; k<tDup->number_of_nodes(); ++k )
                        {
                            tBitset->set( tDup->node(k)->index() );
                        }
                    }
                }
                tBitset->reset( tNode->index() );
                tBitset->where( tIndices );
                tNode->init_vertex_container( tIndices.size() );
                for ( index_t k : tIndices )
                {
                    tNode->insert_vertex( tGraph( k ) );
                }
            }
            delete tBitset;
            graph::symrcm( tGraph );
            for ( Node * tNode : aNodes )
            {
                tNode->reset_vertex_container();
            }
            sort( aNodes, opNodeIndex );
        }
    }
}
