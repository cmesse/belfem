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

             // one master block per sideset, handed over by the processor
             // from a duplicated facet of that sideset ( never read from
             // the first facet here ); their elements are relinked
             Cell< id_t > mMasterBlocks ;

             // domain type of the master blocks; one type per sideset
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

            void
            add_master_block( const id_t aBlockID );

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

            DynamicBitset mIsAirBlock ;
            DynamicBitset mIsFerroBlock ;

            Map< id_t, index_t > mBlockIndices ;
            Cell< InterfaceSet * > mSets ;

            // link sideset IDs to set
            Map< id_t, InterfaceSet * > mSetsMap ;

            // the master block of each admitted sideset, from its first
            // duplicated facet
            Map< id_t, id_t > mMasterBlockOfSideSet ;

            // one conductor-side copy per node, shared by every tie-weight
            // set with conductor masters and keyed by the id of the node it
            // duplicates; the mesh owns the copies
            //
            // Copy numbering: copies are allocated set by set, in the sets'
            // order ( the sorted keys hex( master blocks )|hex( slave
            // blocks ) of the block-index bitsets, a fixed function of the
            // mesh's block order ), and within a set in ascending id of the
            // original; a shared copy is allocated by the first set that
            // reaches the original. They are appended to the mesh in
            // ascending copy id, which is the allocation order. Both follow
            // from the block order, the set order and the original ids
            // alone, never from the iteration order of a Map
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

            // whether the interface at aFacet is one this processor
            // duplicates: a slave exists, the two sides are not both air,
            // and one side is air or ferro
            bool
            duplicates_at( Facet * aFacet ) const ;

            // whether a coil block lies on either side of aFacet
            bool
            touches_coil( Facet * aFacet ) const ;

            // why aFacet is no interface facet: no slave, both sides air,
            // or neither side air or ferro
            const char *
            not_interface_reason( Facet * aFacet ) const ;

            // every sideset with a facet this processor duplicates at must
            // have masters of one domain type on all its facets that have a
            // slave, its facets must agree on touching a coil, every one of
            // its facets must be duplicated at, and if its copies are private
            // ( non-conductor masters ) they must lie on one master block;
            // checked in that order, each check over the whole sideset, so
            // the class of a refusal does not depend on the facet order;
            // runs before any sideset is admitted and before any interface
            // copy is allocated
            void
            refuse_mixed_sidesets() const ;

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

            // aSet holds a copy of aNode but no admitted set holds a copy of
            // its periodic partner aPartner; aOtherType and aOtherTreatment
            // are the excluded sets that hold one, or nullptr. Always refuses
            void
            refuse_unpaired_periodic_copy( InterfaceSet * aSet,
                                           Node * aNode,
                                           Node * aPartner,
                                           InterfaceSet * aOtherType,
                                           InterfaceSet * aOtherTreatment ) const ;

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
