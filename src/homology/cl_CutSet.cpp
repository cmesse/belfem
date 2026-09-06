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

#include "cl_CutSet.hpp"

#include "cl_CutData.hpp"

namespace belfem
{
    namespace mesh
    {
        CutSet::CutSet(
                Mesh * aMesh,
                Cell< Node * >          & aNodeOriginals,
                const string & aHexString,
                const index_t  aNumberOfCuts ):
             mMesh( aMesh ),
             mNodeOriginals( aNodeOriginals )
        {
            mBitset = new DynamicBitset( aNumberOfCuts );
            mBitset->set_from_hex( aHexString );
            mBitset->lock();
            mNodeBitset = new DynamicBitset( aNodeOriginals.size() );

        }

        CutSet::~CutSet()
        {
            delete mBitset;
        }


        void
        CutSet::create_duplicates( id_t & aMaxNodeID, Cell< Node * > & aAbstractNodes )
        {
            // create containers
            uint tNumAbstractDofs = mBitset->count();

            Cell< Node * > tSources( tNumAbstractDofs + 1, nullptr );

            // all weights are 1 at any Lagrange order: the jump [phi] = I is
            // constant over the cut, so each duplicate ( corner or midside )
            // ties 1:1 to its original plus the abstract current dofs
            Vector< real > tWeights( tNumAbstractDofs + 1, 1.0 );

            index_t tCount = 0 ;
            for ( uint c=0; c<aAbstractNodes.size(); ++c )
            {
                if ( mBitset->test( c ) )
                {
                    tSources( tCount++ ) = aAbstractNodes( c );
                }
            }

            if ( mMesh->has_periodicity() )
            {
                for( Node * tOrg : mNodeOriginals )
                {
                    tOrg->unflag( 6 );
                }
            }

            // tier 1, test non-periodic nodes
            for( Node * tOrg : mNodeOriginals )
            {
                if ( tOrg->is_periodic() ) continue;
                if ( mNodeBitset->test( tOrg->index() ) )
                {
                    Node * tDup = new Node( ++aMaxNodeID, tOrg->x(), tOrg->y(), tOrg->z() );

                    tSources( tCount ) = tOrg ;
                    tDup->set_sources( tSources, tWeights );

                    mNodeDuplicates[ tOrg->id() ] = tDup ;
                }

            }

            if ( mMesh->has_periodicity() )
            {
                for( Node * tOrg : mNodeOriginals )
                {
                    if ( ! tOrg->is_periodic() )
                    {
                        tOrg->flag( 6 );
                    }
                }

                Periodicity * tPeriodicity = mMesh->periodicity();

                for( Node * tOrg : mNodeOriginals )
                {
                    // either periodic or processed
                    if ( tOrg->is_flagged( 6 ) ) continue ;

                    Node * tOrgA = tOrg->id() < tOrg->periodic()->id() ? tOrg : tOrg->periodic() ;
                    Node * tOrgB = tOrg->id() < tOrg->periodic()->id() ? tOrg->periodic() : tOrg ;

                    // Step 6 three-way branch ( replaces the too-strict symmetry
                    // assert / the debug-only 4c skip; correct in BOTH builds ).
                    // This gets every pair PAST relink ( all dups exist ). A prior
                    // "Step 6c" was proposed to explicitly register FRAGMENTATION's
                    // extra single-partner dups via a verdict table, then REJECTED
                    // as unsafe ( see todo/closed/periodic_thin_cut_continuity_fix.md,
                    // "Why Pairing Was Rejected" ) — registering could clobber a
                    // node's single mPeriodic slot and corrupt slave T-matrices.
                    // Those unbacked seam dups don't need backup registration
                    // because the periodic rebuild keys facet/edge corners by
                    // original identity instead ( cl_Mesh_PeriodicityFactory.cpp ),
                    // a different mechanism than this branch.
                    //   both bits   -> paired duplicates ( periodic continuity )
                    //   exactly one -> one-sided duplicate, this CutSet's pattern
                    //                  sources, no tie/backup. Correct for a true
                    //                  one-sided JUMP ( jump rides through the
                    //                  original p_A<->p_B tie ).
                    //   neither     -> nothing to duplicate in this CutSet
                    const bool tBitA = mNodeBitset->test( tOrgA->index() );
                    const bool tBitB = mNodeBitset->test( tOrgB->index() );

                    if ( tBitA && tBitB )
                    {
                        BELFEM_ASSERT( tOrgA != tOrgB, "Node %lu is its own periodic partner",
                                       ( long unsigned int ) tOrgA->id() );

                        Node * tDupA = new Node( ++aMaxNodeID, tOrgA->x(), tOrgA->y(), tOrgA->z() );
                        Node * tDupB = new Node( ++aMaxNodeID, tOrgB->x(), tOrgB->y(), tOrgB->z() );
                        tSources( tCount ) = tOrgA ;
                        tDupA->set_sources( tSources, tWeights );
                        tSources( tCount ) = tOrgB ;
                        tDupB->set_sources( tSources, tWeights );

                        tDupA->set_periodic( tDupB );
                        tDupB->set_periodic( tDupA );

                        tPeriodicity->add_node_pair_to_backup( tDupA, tDupB );

                        mNodeDuplicates[ tOrgA->id() ] = tDupA ;
                        mNodeDuplicates[ tOrgB->id() ] = tDupB ;
                    }
                    else if ( tBitA != tBitB )
                    {
                        // one-sided: duplicate the member only, this CutSet's
                        // pattern signature, no periodic tie, no backup pair.
                        Node * tMember = tBitA ? tOrgA : tOrgB ;

                        Node * tDup = new Node( ++aMaxNodeID, tMember->x(), tMember->y(), tMember->z() );
                        tSources( tCount ) = tMember ;
                        tDup->set_sources( tSources, tWeights );

                        mNodeDuplicates[ tMember->id() ] = tDup ;
                    }

                    // gate both partners so the periodic pair is processed exactly once
                    tOrgA->flag( 6 );
                    tOrgB->flag( 6 );
                }

                for( Node * tOrg : mNodeOriginals )
                {
                    tOrg->unflag( 6 );
                }
            }
        }

    }
}
