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

#ifndef BELFEM_CL_EDGEFACTORY_HPP
#define BELFEM_CL_EDGEFACTORY_HPP

#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Cell.hpp"
#include "cl_Element.hpp"
#include "cl_Edge.hpp"
#include "cl_Node.hpp"
#include "cl_Graph_Vertex.hpp"
#include "cl_Mesh.hpp"
#include "cl_Map.hpp"
namespace belfem
{
    namespace mesh
    {
//-----------------------------------------------------------------------

        class EdgeFactory
        {
            const proc_t mCommRank ;

            // ref to mesh we work on
            Mesh & mMesh ;

            // number of nodes on mesh
            const key_t mNumberOfNodes ;

            // container for relevant elements
            Cell< Element * > mElements ;

            // map for edges
            Map< key_t, Edge * > mMap ;

            // max order of elements
            uint mElementOrder = 0;

//-----------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------

            EdgeFactory( Mesh & aMesh );

//-----------------------------------------------------------------------

            EdgeFactory( Mesh * aMesh );

//-----------------------------------------------------------------------

            ~EdgeFactory() = default ;

//----------------------------------------------------------------------

            void
            create_edges( const Vector< id_t > aNedelecBlocks = Vector< id_t >(),
                          const Vector< id_t > aNedelecSideSets = Vector< id_t >(),
                          const bool aCreateEdgesOnAllSideSets = true);

//----------------------------------------------------------------------

            // print result for debugging
            void
            print();

//----------------------------------------------------------------------
        private:
//----------------------------------------------------------------------

            // a routine for the lazy that collects all block IDs
            void
            get_all_block_ids( Vector< id_t > & aBlockIDs );

//----------------------------------------------------------------------

            void
            collect_elements(
                    const Vector< id_t > & aBlockIDs,
                    const Vector< id_t > & aSideSetIDs,
                    const bool aCreateEdgesOnAllSideSets );

//----------------------------------------------------------------------

            // create edge IDs
            void
            create_edge_keys( Vector< key_t > & aIDs );

//----------------------------------------------------------------------

            key_t
            edge_key( Element * aElement, const uint aEdgeIndex,
                     Cell< Node * > & aNodes );

//----------------------------------------------------------------------

            void
            create_edges_on_master(  const Vector< key_t > & aKeys );

//----------------------------------------------------------------------

            void
            link_elements_to_edges();

//----------------------------------------------------------------------

            void
            link_curves_to_edges();

//----------------------------------------------------------------------

            void
            compute_edge_ownerships();

//----------------------------------------------------------------------

            /**
             * complete the node list with midside nodes
             */
            void
            grab_nodes( const key_t aKey, Cell< Node* > & aNodes );

//-----------------------------------------------------------------------

            void
            set_edge_ids();

//-----------------------------------------------------------------------
        };

//-----------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_EDGEFACTORY_HPP
