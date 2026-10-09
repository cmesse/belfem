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

#ifndef CL_INTERFACEPROCESSOR_HPP
#define CL_INTERFACEPROCESSOR_HPP

#include "typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "cl_DynamicBitset.hpp"
#include "cl_Topology.hpp"

namespace belfem
{
    namespace mesh
    {
        // how an interface set ties its duplicates to their originals
        enum class InterfaceTreatment
        {
            TieWeight1,   // default (ferro-air, conductor-air, conductor-ferro): duplicate hangs on its original with weight 1 (viz split only)
            Decouple      // coil-touching: duplicate carries no sources, excluded from the solve
        };

        class InterfaceSet
        {
             Mesh * mMesh ;

             DynamicBitset * mBitset ;

             Cell< Element * > mElements ;
             Map< id_t, Node * > mOriginals ;

             Map< id_t, Node * > mDuplicates ;
             uint mNumberOfSideSets = 0 ;

             InterfaceTreatment mTreatment = InterfaceTreatment::TieWeight1 ;

             // domain type of the master block of the first sideset's first facet
             DomainType mMasterType = DomainType::UNDEFINED ;

//-----------------------------------------------------------------------------
             Cell< SideSet * > mSideSets ;
//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            InterfaceSet( Mesh * aMesh, const string aBitsetHex );

            ~InterfaceSet();

            void
            increment_sideset_counter();

            void
            allocate_sideset_container();

            void
            add_sideset( SideSet * aSideSet );

            InterfaceTreatment
            treatment() const ;

            DomainType
            master_type() const ;

            // a tie-weight set with conductor masters shares its copies
            bool
            shares_duplicates() const ;

            void
            collect_nodes_and_elements();

            /**
             * takes the copy of each collected node from aSharedDuplicates,
             * keyed by the node's own id, and allocates and registers it there
             * on a miss. The map does not own the copies
             */
            void
            duplicate_nodes( id_t & aMaxNodeID, Map< id_t, Node * > & aSharedDuplicates );

            void
            relink_elements();

            index_t
            number_of_nodes() const ;

            /**
             * a set with private copies appends them to aDuplicates; a
             * tie-weight set registers the original of each copy in
             * aOriginalMap
             */
            void
            add_duplicates( index_t & aCount, Cell< Node * > & aDuplicates, Map< id_t, Node * > & aOriginalMap );

            Map< id_t, Node * > &
            duplicate_map();

            Map< id_t, Node * > &
            original_map();

        };

        class InterfaceProcessor
        {
            Mesh * mMesh ;
            Topology * mTopology ;

            Cell< Node * > & mAbstractNodes ;

            const uint mNumOriginalSidesets ;
            id_t mMaxNodeID ;

            DynamicBitset * mIsAirBlock = nullptr ;
            DynamicBitset * mIsFerroBlock = nullptr ;

            Map< id_t, index_t > mBlockIndices ;
            Cell< InterfaceSet * > mSets ;

            // link sidesetr IDs to set
            Map< id_t, InterfaceSet * > mSetsMap ;

            // one conductor-side copy per node, shared by every tie-weight
            // set with conductor masters and keyed by the id of the node it
            // duplicates; the mesh owns the copies
            Map< id_t, Node * > mSharedDuplicates ;

//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            InterfaceProcessor(
                Mesh * aMesh,
                Topology * aTopology,
                Cell< Node * > & aAbstractNodes,
                const uint aNumOriginalSideSets,
                const id_t aMaxNodeID );

            ~InterfaceProcessor();

//-----------------------------------------------------------------------------
        private:
//-----------------------------------------------------------------------------

            void
            create_interface_sets();

            void
            connect_sidesets_to_interface_sets();

            void
            populate_interface_sets();

            void
            duplicate_nodes();

            void
            pair_cross_set_periodic_duplicates();

            void
            relink_elements();

            void
            add_duplicate_nodes_to_mesh();

            void
            unify_duplicates( Cell< Node * > & aDuplicates );

            bool
            pair_within( InterfaceSet * aSet, InterfaceSet * aOther ) const ;

        };

        inline InterfaceTreatment
        InterfaceSet::treatment() const
        {
            return mTreatment ;
        }

        inline DomainType
        InterfaceSet::master_type() const
        {
            return mMasterType ;
        }

        inline bool
        InterfaceSet::shares_duplicates() const
        {
            return mTreatment == InterfaceTreatment::TieWeight1 && mMasterType == DomainType::Conductor ;
        }

        inline Map< id_t, Node * > &
        InterfaceSet::duplicate_map()
        {
            return mDuplicates ;
        }

        inline Map< id_t, Node * > &
        InterfaceSet::original_map()
        {
            return mOriginals ;
        }

    }
}
#endif //CL_INTERFACEPROCESSOR_HPP
