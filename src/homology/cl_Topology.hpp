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

#ifndef CL_TOPOLOGY_HPP
#define CL_TOPOLOGY_HPP
#include "cl_DynamicBitset.hpp"
#include "cl_Map.hpp"
#include "cl_Mesh.hpp"
#include "en_DomainType.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------
        /**
         * @class Topology
         * @brief Analyzes and classifies mesh topology for finite element simulations
         * 
         * The Topology class provides functionality to analyze mesh structures and
         * categorize different regions based on their domain types. It identifies
         * blocks, boundaries, and interfaces within the mesh, and is particularly
         * useful for electromagnetic simulations where regions need to be classified
         * as phi (electromagnetic) or non-phi domains.
         * 
         * Key features:
         * - Categorizes mesh regions by domain type (blocks, boundaries, interfaces)
         * - Identifies electromagnetic domains (phi vs non-phi blocks)
         * - Detects sideset types based on neighboring block relationships
         * - Manages topological mappings between mesh IDs and domain classifications
         * 
         * This class is essential for preprocessing meshes to determine which regions
         * require different physics solvers or boundary conditions.
         */
        class Topology
        {
            const proc_t mCommRank ;
            Mesh * mMesh ;

            Vector< id_t > mPhiInterfaceIDs ;
            Vector< id_t > mPhiBoundaryIDs ;
            Vector< id_t > mPhiPeriodicIDs ;

            Map< DomainType, Vector< id_t > * > mTypeMap ;

            Map< id_t, DomainType > mBlockTypes ;

            // needed for the cut factory
            Vector< id_t > mPhiBlockIDs;
            Vector< id_t > mNonPhiBlockIDs;

            // enrichment entities to skip when scanning a reloaded mesh
            // ( filled by run_on_enriched_mesh, empty on the fresh path )
            Map< id_t, index_t > mEnrichedBlocks ;
            Map< id_t, index_t > mEnrichedSideSets ;

        public:

            Topology( Mesh * aMesh );

            ~Topology();

            void
            run();

            /**
             * variant of run() for a mesh that was restored from a .bfm file:
             * the fresh path builds the maps BEFORE the enrichment factories
             * run, so the thin-shell layer/buffer blocks and the tape/ghost
             * sidesets must not enter the maps, and the sideset types
             * (already final on the reloaded mesh) must not be re-detected
             */
            void
            run_on_enriched_mesh();

            const Vector< id_t > &
            groups( const DomainType aType ) const;

            const Vector< id_t > &
            phi_block_ids() const ;

            const Vector< id_t >  &
            non_phi_block_ids() const ;

            const Vector< id_t > &
            phi_interface_ids() const ;

            const Vector< id_t > &
            phi_boundary_ids() const ;

            const Vector< id_t > &
            phi_periodic_ids() const ;

            DomainType
            sideset_type( const DomainType aMaster, const DomainType aSlave );

            DomainType
            sideset_type( const DomainType aMaster );

            const Map< id_t, DomainType > &
            block_map() const ;

            void
            synchronize_maps();

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            collect_enrichment_ids();

            void
            collect_block_and_sideset_types();

            void
            detect_sideset_types();

            void
            select_blocks();

            void
            select_sidesets();

            void
            update_block_map();

        };

//------------------------------------------------------------------------------

        inline const Vector< id_t > &
        Topology::groups( const DomainType aType ) const
        {
            return *mTypeMap( aType );
        }

        inline const Map< id_t, DomainType > &
        Topology::block_map() const
        {
            return mBlockTypes ;
        }

        inline const Vector< id_t > &
        Topology::phi_block_ids() const
        {
            return mPhiBlockIDs ;
        }

        inline const Vector< id_t > &
        Topology::phi_boundary_ids() const
        {
            return mPhiBoundaryIDs ;
        }

        inline const Vector< id_t > &
        Topology::phi_periodic_ids() const
        {
            return mPhiPeriodicIDs ;
        }

        inline const Vector< id_t > &
        Topology::phi_interface_ids() const
        {
            return mPhiInterfaceIDs ;
        }

        inline const Vector< id_t > &
        Topology::non_phi_block_ids() const
        {
            return mNonPhiBlockIDs ;
        }


//------------------------------------------------------------------------------
    }
}
#endif //CL_TOPOLOGY_HPP
