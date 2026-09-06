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

#ifndef BELFEM_CL_MESH_PARTITIONER_HPP
#define BELFEM_CL_MESH_PARTITIONER_HPP



#include "graph_typedefs.hpp"
#include "cl_Mesh.hpp"


namespace belfem
{
    namespace mesh
    {


        class Partitioner
        {
            Mesh * mMesh;

            const metis_t mNumberOfPartitions;

            Vector< metis_t > mPartition;

            const bool mForceContiguousPartitions = true ;

            const bool mResetVertexContainers = true ;
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Partitioner( Mesh * aMesh,
                    const uint aNumberOfPartitions,
                    bool aSetProcOwnerships = true,
                    bool aForceContiguousPartitions = true,
                    bool aResetVertexContainers = true );

//------------------------------------------------------------------------------

            ~Partitioner();

//------------------------------------------------------------------------------

            /**
             * returns a vector containing the partitioning data
             */
             const Vector < metis_t > &
             partition();

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            create_graph( Graph & aGraph );

//------------------------------------------------------------------------------

            void
            run_metis( Graph & aGraph );

//------------------------------------------------------------------------------

            void
            set_element_owners();

//------------------------------------------------------------------------------

            void
            fix_facet_related_ownerships( );

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_MESH_PARTITIONER_HPP
