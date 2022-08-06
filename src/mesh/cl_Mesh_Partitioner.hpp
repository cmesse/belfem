//
// Created by Christian Messe on 27.10.19.
//

#ifndef BELFEM_CL_MESH_PARTITIONER_HPP
#define BELFEM_CL_MESH_PARTITIONER_HPP


#ifdef BELFEM_SUITESPARSE
// compiler settings
#ifdef BELFEM_GCC
#ifndef __GNUC__
#define __GNUC__
#endif
#elif BELFEM_INTEL
#ifndef __ICC
#define __ICC
#endif
#endif

#include <metis.h>

// bugfix to avoid error with BLAZE
#ifdef BELFEM_BLAZE
#undef abs
#undef iabs
#endif

#endif

#include "cl_Mesh.hpp"


namespace belfem
{
    namespace mesh
    {
#ifdef BELFEM_SUITESPARSE
        typedef idx_t metis_t;
#else
        typedef int metis_t;
#endif

        class Partitioner
        {
            Mesh * mMesh;

            metis_t   mNumberOfPartitions;
            metis_t   mNumberOfElements = 0;
            metis_t   mAdjacencySize = 0;
            metis_t * mElementPointers;
            metis_t * mAdjacency;

            Vector< metis_t > mPartition;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Partitioner( Mesh * aMesh,
                    const uint aNumberOfPartitions,
                    bool aSetProcOwnerships = true );

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
            create_graph();

//------------------------------------------------------------------------------

            void
            run_metis();

//------------------------------------------------------------------------------

            string
            metis_status( const int aStatus );

//------------------------------------------------------------------------------

            void
            set_element_owners();

//------------------------------------------------------------------------------

            void
            fix_facet_related_ownerships();

//------------------------------------------------------------------------------

            void
            fix_ghost_sideset_related_ownerships() ;

//------------------------------------------------------------------------------

            void
            fix_cut_related_ownerships() ;

//------------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_MESH_PARTITIONER_HPP
