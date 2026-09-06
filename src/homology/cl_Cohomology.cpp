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

#include <iostream>
#include <string>

#include "cl_Cohomology.hpp"
#include "cl_DynamicBitset.hpp"
#include "cl_Logger.hpp"

#include "fn_Graph_spfa.hpp"
#include "fn_inv.hpp"
#include "fn_trans.hpp"

namespace belfem
{
    namespace mesh
    {
//-----------------------------------------------------------------------------

        Cohomology::Cohomology( SimplicialComplex * aSimplicialComplex, Mesh * aMesh ) :
                mSimplicialComplex( aSimplicialComplex ),
                mMesh( aMesh )
        {
            mGenerators.set_size(4,Cell< Cochain * >());
            mOrders.set_size(4,Cell< int >());

            mD = mSimplicialComplex->createMatrixFromCoboundaryMap();
            this->cohomologyGroupOfChainComplex();
            this->generatorsOfCohomology();
            this->clean();
        }

//------------------------------------------------------------------------------

        Cohomology::Cohomology( SimplicialComplex * aSimplicialComplex, Mesh * aMesh, BeltedTree * aBTree ) :
                mSimplicialComplex( aSimplicialComplex ),
                mMesh( aMesh )
        {
            mGenerators.set_size(4,Cell< Cochain * >());
            mOrders.set_size(4,Cell< int >());

            aBTree->compute_cohomology();

            //Copy the belted tree cohomology cochain in the cohomology
            for (Cochain * tCochain : aBTree->get_cohomology())
            {
                Cochain * tGenerator = new Cochain(1,mMesh, false, false);
                tGenerator->addCochainToCochain(tCochain,1);
                mGenerators(1).push(tGenerator);
            }
            this->clean();
        }

//------------------------------------------------------------------------------

        Cohomology::~Cohomology()
        {
            for(uint i = 0; i < mGenerators.size(); i++)
            {
                for (uint j = 0; j < mGenerators(i).size(); j++)
                {
                    delete mGenerators(i)(j);
                }
            }
            mGenerators.clear();
            mOrders.clear();
            mV.clear();
            mW.clear();
            mU.clear();
            mB.clear();
            mD.clear();
            if ( mProgress != nullptr ) delete mProgress ;
        }

//-----------------------------------------------------------------------------

        void
        Cohomology::cohomologyGroupOfChainComplex()
        {
            mV.clear();
            mW.clear();

            // loop over all dimensions
            for(uint k = 0; k < 4; k++)
            {
                // populate the matrices W and V (kernel and image of the coboundary)
                auto [w, v] = kernelImage(mD(k));
                mW[k] = w;
                mV[k+1] = v;
            }
            mV[0] = Matrix< int >(mW[0].n_rows(),1,0);

            //Compute the quotient groups
            this->quotientGroup();
        }

//-----------------------------------------------------------------------------

        void
        Cohomology::quotientGroup()
        {
            uint n;

            //loop over all dimensions
            for(uint k = 0; k < 4; k++)
            {
                n = mV[k].n_cols();
                mB[k].set_size(mW[k].n_cols(),n,0);

                //Solve W\V to get the A matrix
                for (uint i = 0; i < n; i++)
                {
                    Matrix< int > tM( mV(k).n_rows(), 1 );
                    tM.set_col( 0, mV(k).col(i) );
                    Matrix < int > tv = SolveInt(mW[k],tM);

                    for (uint j = 0; j < mW[k].n_cols(); j++)
                    {
                        mB[k](j,i) = tv(j,0);
                    }
                }

                // Compute the Smith normal form of W\V
                auto [tQ, tQ_, tR, tR_, s, t] = smithForm(mB[k]);
                ms[k] = s;
                mt[k] = t;

                //Get the U matrix, with the cohomology generators as the last columns.
                mU[k] = mW[k];
                mU(k)*=tQ;
            }
        }

//-----------------------------------------------------------------------------

        void
        Cohomology::generatorsOfCohomology()
        {
            uint tCount;
            uint tCount2;

            //loop over all dimensions
            for (uint k = 0; k < 4; k++)
            {
                tCount = 0;
                if (mSimplicialComplex->get_kcochainMap(k).size() == 0)
                {
                    continue ;
                }
                // loop over the last columns of U
                for (uint j = ms[k]+1; j < mU[k].n_cols()+1; j++)
                {

                    //Get the order (0 means infinity... Maybe should change that)
                    if(j > mt[k])
                    {
                        mOrders(k).push(0);
                    }
                    else
                    {
                        mOrders(k).push(mB[k](j-1,j-1));
                    }

                    // Create the generator from column data
                    Cochain * tCochain = new Cochain(k,mMesh, false, false);
                    mGenerators(k).push(tCochain);
                    tCount2 = 0;
                    for(const auto& [tKey, tCochain2] : mSimplicialComplex->get_kcochainMap(k))
                    {
                        mGenerators(k)(tCount)->addCochainToCochain(tCochain2,mU[k](tCount2,j-1));
                        tCount2++;
                    }
                    tCount++;
                }
            }
        }

//-----------------------------------------------------------------------------

        void Cohomology::clean()
        {
            this->clean_spfa();
        }

//-----------------------------------------------------------------------------

        namespace cohomology
        {
            /**
             * add aWeight times the coboundary of a representative node to the
             * generator, restricted to the complex edge domain, under the
             * periodic master + slave quotient: boundary node indices fold
             * slave -> master via is_flagged(), and only slave entities are
             * flagged ( inherited from the retired greedy cleaner )
             */
            static void
            fire_node_coboundary(
                    Cochain             * aGenerator,
                    Node                * aNode,
                    const int             aWeight,
                    const DynamicBitset & aComplexEdges )
            {
                for ( uint j = 0 ; j < aNode->number_of_edges(); ++j )
                {
                    Edge * tEdge = aNode->edge( j );
                    if ( ! aComplexEdges.test( tEdge->index() ) ) continue ;

                    aGenerator->addSimplexToCochain( tEdge->index(),
                            tEdge->node( 0 )->index() == aNode->index() ? -aWeight : aWeight );

                    Node * tN0 = tEdge->node( 0 );
                    Node * tN1 = tEdge->node( 1 );
                    aGenerator->add_simplex_to_boundary(
                            tN0->is_flagged() ? tN0->periodic()->index() : tN0->index(), -aWeight );
                    aGenerator->add_simplex_to_boundary(
                            tN1->is_flagged() ? tN1->periodic()->index() : tN1->index(),  aWeight );
                }

                // the periodic counterpart is the same node in the quotient:
                // its adjacency must fire too; a self-paired node has no
                // second adjacency
                if ( aNode->is_periodic() && aNode->periodic() != aNode )
                {
                    Node * tPartner = aNode->periodic() ;

                    for ( uint j = 0 ; j < tPartner->number_of_edges(); ++j )
                    {
                        Edge * tEdge = tPartner->edge( j );
                        if ( ! aComplexEdges.test( tEdge->index() ) ) continue ;
                        if ( tEdge->is_flagged() ) continue ;

                        aGenerator->addSimplexToCochain( tEdge->index(),
                                tEdge->node( 0 )->index() == tPartner->index() ? -aWeight : aWeight );

                        Node * tN0 = tEdge->node( 0 );
                        Node * tN1 = tEdge->node( 1 );
                        aGenerator->add_simplex_to_boundary(
                                tN0->is_flagged() ? tN0->periodic()->index() : tN0->index(), -aWeight );
                        aGenerator->add_simplex_to_boundary(
                                tN1->is_flagged() ? tN1->periodic()->index() : tN1->index(),  aWeight );
                    }
                }
            }
        }

//-----------------------------------------------------------------------------

        namespace cohomology
        {
            /**
             * classic greedy rectification sweeps on one generator: fire the
             * outward node of every non-unit edge once per sweep, until unit
             * or stalled. Fires only where |c| >= 2, so a sparse generator
             * stays sparse — this is the preferred rectifier once
             * feasibility is proven. Returns false on stall ( exact count
             * plateau, inherited from the retired greedy cleaner ) or on the
             * sweep cap ( insurance against oscillation the plateau detector
             * cannot see ).
             */
            static bool
            rectify_greedy_sweeps(
                    Cochain             * aGenerator,
                    Cell< Edge * >      & aEdges,
                    const DynamicBitset & aComplexEdges )
            {
                // freeze both fuses on the PRE-greedy state: cascading can
                // inflate the live map, so live sizes are not a baseline
                // ( Grok P2: plateau detection alone is incomplete — a
                // non-constant fire-count cycle never trips it )
                const index_t tSupport0 = aGenerator->getSimplicesMap().size() ;
                index_t tNumNonUnit0 = 0 ;
                for ( const auto & [ tIndex, tCoeff ] : aGenerator->getSimplicesMap() )
                {
                    if ( tCoeff > 1 || tCoeff < -1 )
                    {
                        ++tNumNonUnit0 ;
                    }
                }

                if ( tNumNonUnit0 == 0 )
                {
                    return true ; // already unit
                }

                const index_t tSweepCap    = 64 > 8 * tNumNonUnit0 ? 64 : 8 * tNumNonUnit0 ;
                const index_t tSupportFuse = 4 * tSupport0 ;

                // seed with the map size, exactly like the retired greedy
                // cleaner: the plateau test keeps its battle-tested semantics
                index_t tNumNonUnit = tSupport0 ;
                index_t tNumNonUnitOld ;
                index_t tSweep = 0 ;

                while ( tNumNonUnit > 0 )
                {
                    tNumNonUnitOld = tNumNonUnit ;
                    tNumNonUnit = 0 ;

                    for ( auto & [ tIndex, tCoeff ] : aGenerator->getSimplicesMap() )
                    {
                        if ( tCoeff > 1 || tCoeff < -1 )
                        {
                            ++tNumNonUnit ;

                            // fire the node for which the edge points outward,
                            // folding a slave onto its master
                            Node * tNode = aEdges( tIndex )->node( tCoeff > 0 ? 0 : 1 );
                            if ( tNode->is_flagged() )
                            {
                                tNode = tNode->periodic() ;
                            }

                            fire_node_coboundary( aGenerator, tNode, 1, aComplexEdges );
                        }
                    }

                    // multi-trigger bail-out: exact plateau, sweep cap, or
                    // support-growth fuse ( catches cascade floods early )
                    if ( tNumNonUnit == tNumNonUnitOld
                         || ++tSweep > tSweepCap
                         || aGenerator->getSimplicesMap().size() > tSupportFuse )
                    {
                        return false ;
                    }
                }

                return true ;
            }
        }

//-----------------------------------------------------------------------------

        // Rectifies non-unit coefficients ( |c| >= 2 ) to unit ones, keeping
        // each generator in its equivalence class, or fails with a
        // negative-cycle certificate if no unit representative exists.
        // Strategy ( certify - sparse - dense fallback ):
        //   1. SPFA difference-constraint solve proves feasibility or yields
        //      an edge-ID certificate of the obstruction.
        //   2. Feasible generators are rectified by the classic greedy sweeps,
        //      which fire only at non-unit edges and preserve the original
        //      sparse support ( a pure feasibility theta write-back is unit
        //      but DENSE: off-tree edges land at +-1 generically ).
        //   3. If greedy stalls on a proven-feasible generator ( never
        //      observed ), fall back to the dense theta write-back after
        //      re-solving on the mutated ( same-class ) coefficients.
        // Theory and plan: src/homology/doc/thin_cut_nonunit_rectification.md,
        // todo/thin_cut_nonunit_rectification_implementation.md
        void
        Cohomology::clean_spfa()
        {
            BELFEM_ERROR( mSimplicialComplex != nullptr,
                    "clean_spfa() requires a simplicial complex" );

            mMesh->update_node_indices(); // <-- might not be necessary but is cheap

            // flag the entities on the slave periodic side; rep() below
            // folds flagged slaves onto masters
            mMesh->unflag_everything() ;

            if ( mGenerators( 1 ).size() == 0 )
            {
                return ; // flags already cleared on every exit path
            }
            const DynamicBitset & tComplexEdges = mSimplicialComplex->original_edges() ;
            if ( mMesh->has_periodicity() )
            {
                for ( Node * tNode : mMesh->periodicity()->slave_nodes() )
                {
                    tNode->flag() ;
                }
                for ( Edge * tEdge : mMesh->periodicity()->slave_edges() )
                {
                    tEdge->flag() ;
                }
            }

            Cell< Edge* > & tEdges = mMesh->edges() ;
            Cell< Node* > & tNodes = mMesh->nodes() ;
            const index_t tNumNodes = mMesh->number_of_nodes() ;

            // the cohomology edge domain: all edges of the complex, including
            // those with zero coefficient ( their unit constraints still
            // couple the graph )
            Cell< index_t > tDomainEdges ;
            tComplexEdges.where( tDomainEdges );
            const index_t tNumDomainEdges = tDomainEdges.size() ;

            // two constraint arcs per domain edge e = ( i, j ) with
            // u = rep( i ), v = rep( j ), c = c( e ) :
            //     theta( v ) - theta( u ) <= 1 + c    ( arc 2k     : u -> v )
            //     theta( u ) - theta( v ) <= 1 - c    ( arc 2k + 1 : v -> u )
            // which is | c + theta( u ) - theta( v ) | <= 1
            Cell< index_t > tArcTails(   2 * tNumDomainEdges, 0 );
            Cell< index_t > tArcHeads(   2 * tNumDomainEdges, 0 );
            Cell< int64_t > tArcWeights( 2 * tNumDomainEdges, 0 );

            Cell< int64_t > tTheta ;
            Cell< index_t > tCycle ;

            uint tGeneratorCount = 0 ;

            uint tNumGenerators = mGenerators( 1 ).size();

            InfoLevel tLevel = static_cast< InfoLevel >( gLog.info_level() );

            if ( ( tLevel == InfoLevel::Default || tLevel == InfoLevel::Detailed ) && tNumGenerators > 4 )
            {
                if ( mProgress != nullptr ) delete mProgress ;
                mProgress = new Progressbar( 2 * tNumGenerators );
            }

            for ( Cochain * tGenerator : mGenerators( 1 ) )
            {
                index_t tNumNonUnit = 0 ;

                // constraint arcs from the CURRENT coefficients ( reused by
                // the dense fallback after greedy has mutated them )
                auto tBuildArcs = [&]()
                {
                    tNumNonUnit = 0 ;
                    for ( index_t k = 0; k < tNumDomainEdges; ++k )
                    {
                        Edge * tEdge = tEdges( tDomainEdges( k ) );

                        const int64_t tC = tGenerator->getCoefficient( tEdge->index() );
                        if ( tC > 1 || tC < -1 )
                        {
                            ++tNumNonUnit ;
                        }

                        Node * tN0 = tEdge->node( 0 );
                        Node * tN1 = tEdge->node( 1 );
                        const index_t tU = tN0->is_flagged() ? tN0->periodic()->index() : tN0->index() ;
                        const index_t tV = tN1->is_flagged() ? tN1->periodic()->index() : tN1->index() ;

                        tArcTails(   2 * k     ) = tU ;
                        tArcHeads(   2 * k     ) = tV ;
                        tArcWeights( 2 * k     ) = ( int64_t ) 1 + tC ;
                        tArcTails(   2 * k + 1 ) = tV ;
                        tArcHeads(   2 * k + 1 ) = tU ;
                        tArcWeights( 2 * k + 1 ) = ( int64_t ) 1 - tC ;
                    }
                };
                tBuildArcs();

                const bool tFeasible = graph::spfa_difference_constraints(
                        tNumNodes, tArcTails, tArcHeads, tArcWeights, tTheta, tCycle );

                if ( ! tFeasible )
                {
                    // report the obstruction loop as mesh edge IDs so the
                    // throat can be located and refined in the mesher
                    const uint tMaxPrint = 32 ;
                    std::string tString ;


                    Cell< Edge * > tTroubleMakers ;

                    for ( index_t a = 0; a < tCycle.size() && a < tMaxPrint; ++a )
                    {
                        Edge * tEdge = tEdges( tDomainEdges( tCycle( a ) / 2 ) ) ;
                        tTroubleMakers.push( tEdge ) ;
                        tString += "(" ;
                        tString += std::to_string( tEdge->node( 0 )->id() );
                        tString += ")->(" ;
                        tString += std::to_string( tEdge->node( 1 )->id() );
                        tString += ") " ;
                    }


                    // find closest element
                    Cell< Element * > tCandidates ;
                    for ( Edge * tEdge : tTroubleMakers )
                    {
                        for ( uint k=0; k<2; ++k )
                        {
                            Node * tNode = tEdge->node( k );
                            for ( uint e=0; e<tNode->number_of_elements(); ++e )
                            {
                                tCandidates.push( tNode->element( e ) );
                            }
                        }
                    }

                    unique( tCandidates );

                    // note: we will trow an error below
                    // hence it doesn't matter if we corrupt the mesh flags
                    for ( Element * tElement : tCandidates )
                    {
                        tElement->unflag_edges();
                        tElement->unflag();
                    }
                    for ( Edge * tEdge : tTroubleMakers )
                    {
                        tEdge->flag();
                    }
                    for ( Element * tElement : tCandidates )
                    {
                        uint tCount =  0 ;
                        tElement->set_level( 0 );
                        if ( ! tElement->has_edges() ) continue;

                        for ( uint e=0; e<tElement->number_of_edges(); ++e )
                        {
                           if ( tElement->edge( e )->is_flagged() )
                           {
                               ++tCount ;
                           }
                        }
                        tElement->set_level( tCount );
                    }

                    if ( tCycle.size() > tMaxPrint )
                    {
                        tString += "..." ;
                    }

                    sort( tCandidates.begin(), tCandidates.end(), []( Element * a, Element * b ) { return a->level() < b->level(); } );

                    tString += " near element " + std::to_string( tCandidates.last()->id() );

                    mMesh->save( "error.exo");

                    BELFEM_ERROR( false,
                            "BELFEM found cohomology generator %u, but it cannot be used:\n"
                            "it is a thick cut, while the solver supports only thin cuts: a surface of element\n"
                            "faces with a single potential jump. Thick cuts as basis functions are not implemented yet.\n"
                            "No representative of this generator has only -1, 0, and 1 as coefficients on this mesh,\n"
                            "so no thin cut exists here. The loop of %lu edges between nodes %s winds around the\n"
                            "conductor more times than it has edges (see error.exo).\n"
                            "Refine the mesh along this loop and rerun.",
                            tGeneratorCount,
                            ( long unsigned int ) tCycle.size(),
                            tString.c_str() );
                }

                // sparse rectification first: greedy sweeps are safe now that
                // feasibility is proven, and they fire only where |c| >= 2,
                // preserving the original sparse support ( a feasibility theta
                // write-back is unit but DENSE: off-tree edges land at +-1 )
                const bool tSparse = cohomology::rectify_greedy_sweeps(
                        tGenerator, tEdges, tComplexEdges );

                if ( ! tSparse )
                {
                    // proven-feasible stall ( never observed ): fall back to
                    // the dense theta write-back. Greedy mutated the
                    // coefficients by exact coboundaries — class and
                    // feasibility intact — so theta must be re-solved on the
                    // current state before firing
                    gLog.message( InfoLevel::Verbose,
                            "    WARNING: greedy rectification stalled on feasible generator %u - using dense fallback",
                            tGeneratorCount );

                    tBuildArcs();
                    const bool tRefeasible = graph::spfa_difference_constraints(
                            tNumNodes, tArcTails, tArcHeads, tArcWeights, tTheta, tCycle );

                    BELFEM_ERROR( tRefeasible,
                            "clean_spfa: generator %u became infeasible during greedy sweeps ( internal defect )",
                            tGeneratorCount );

                    // rectify c <- c + d theta : fire every representative
                    // node with weight -theta
                    for ( index_t n = 0; n < tNumNodes; ++n )
                    {
                        Node * tNode = tNodes( n );

                        // slaves fold onto their masters; a self-paired node is
                        // its own representative and must still fire
                        if ( tNode->is_flagged() && tNode->periodic() != tNode )
                        {
                            continue ;
                        }

                        const int64_t tW = -tTheta( tNode->index() );
                        if ( tW == 0 )
                        {
                            continue ;
                        }

                        // cochain coefficients are int: guard the firing weight
                        BELFEM_ERROR( tW <= 0x3FFFFFFF && tW >= -( int64_t ) 0x3FFFFFFF,
                                "clean_spfa: potential %ld on node %lu exceeds the coefficient range",
                                ( long int ) -tW,
                                ( long unsigned int ) tNode->index() );

                        cohomology::fire_node_coboundary(
                                tGenerator, tNode, ( int ) tW, tComplexEdges );
                    }
                }

                // safety gate ( release-active ): the rectified generator must
                // be unit everywhere; closedness is preserved by construction,
                // since only node coboundaries were added ( dd = 0 )
                for ( const auto & [ tIndex, tCoeff ] : tGenerator->getSimplicesMap() )
                {
                    BELFEM_ERROR( tCoeff <= 1 && tCoeff >= -1,
                            "clean_spfa: generator %u keeps non-unit coefficient %d on edge index %lu after rectification",
                            tGeneratorCount,
                            tCoeff,
                            ( long unsigned int ) tIndex );
                }

                gLog.message( InfoLevel::Verbose,
                        "    cohomology generator %u : rectified %lu non-unit coefficients ( %s, support %lu )",
                        tGeneratorCount,
                        ( long unsigned int ) tNumNonUnit,
                        tSparse ? "sparse greedy" : "dense fallback",
                        ( long unsigned int ) tGenerator->getSimplicesMap().size() );

                if ( mProgress != nullptr )
                {
                    mProgress->step();
                }
                ++tGeneratorCount ;
            }

            // pocket census + Tier-A removal on the rectified generators
            // ( runs after ALL generators are rectified, so the
            // other-generator overlap check sees final coefficients )
            this->remove_cut_pockets( true );

            mMesh->unflag_everything() ;

            if ( mProgress != nullptr )
            {
                mProgress->finish() ;
                delete mProgress ;
                mProgress = nullptr ;
            }
        }

//-----------------------------------------------------------------------------

        // Pocket census + removal ( rules: todo/cut_pocket_removal_rules.md ).
        // A pocket is a connected component K of the zero-coefficient quotient
        // graph whose entire boundary is cancelled by one exact coboundary
        // flip s * delta( 1_K ). Flips preserve the cohomology class
        // unconditionally; only Tier-A pure pockets ( no structure contact,
        // no periodic-cap contact, no other-generator overlap, small side )
        // are fired. Everything is logged for diagnostics.
        void
        Cohomology::remove_cut_pockets( const bool aFireTierA )
        {
            const DynamicBitset & tComplexEdges = mSimplicialComplex->original_edges() ;

            Cell< Edge* > & tEdges = mMesh->edges() ;
            Cell< Node* > & tNodes = mMesh->nodes() ;
            const index_t tNumNodes = mMesh->number_of_nodes() ;

            Cell< index_t > tDomainEdges ;
            tComplexEdges.where( tDomainEdges );
            const index_t tNumDomainEdges = tDomainEdges.size() ;

            if ( tNumDomainEdges == 0 )
            {
                return ;
            }

            // conservative Tier-A structure predicate ( rules O1 ): a node on
            // ANY sideset ( interfaces, shells, outer boundaries, caps, cuts )
            // counts as structure contact; duplicates mark their originals too
            DynamicBitset tStructureNodes( tNumNodes );
            for ( SideSet * tSideSet : mMesh->sidesets() )
            {
                for ( Facet * tFacet : tSideSet->facets() )
                {
                    for ( uint k = 0; k < tFacet->number_of_nodes(); ++k )
                    {
                        tStructureNodes.set( tFacet->node( k )->index() );
                        tStructureNodes.set( tFacet->node( k )->original()->index() );
                    }
                }
            }

            // quotient fold, identical to the solver and the firing helper
            auto tRep = []( Node * aNode ) -> index_t
            {
                return aNode->is_flagged() ? aNode->periodic()->index() : aNode->index() ;
            };

            // active quotient representatives
            DynamicBitset tActive( tNumNodes );
            for ( index_t k = 0; k < tNumDomainEdges; ++k )
            {
                Edge * tEdge = tEdges( tDomainEdges( k ) );
                tActive.set( tRep( tEdge->node( 0 ) ) );
                tActive.set( tRep( tEdge->node( 1 ) ) );
            }
            Cell< index_t > tActiveList ;
            tActive.where( tActiveList );
            const index_t tNumActive = tActiveList.size() ;

            // scratch
            DynamicBitset   tVisited( tNumNodes );
            DynamicBitset   tInK( tNumNodes );
            DynamicBitset   tEdgeSeen( mMesh->number_of_edges() );
            Cell< index_t > tQueue( tNumActive, 0 );
            Cell< index_t > tComponent ;
            Cell< index_t > tCandidates ; // boundary candidate edge indices

            uint tGeneratorCount = 0 ;

            for ( Cochain * tGenerator : mGenerators( 1 ) )
            {
                index_t tNumFlipped = 0 ;

                // one scan pass over the current zero-graph components.
                // aFire: fire the first eligible Tier-A pure pocket and return
                // true ( caller rebuilds, Rule 5 ).
                auto tScan = [&]( const bool aFire ) -> bool
                {
                    tVisited.reset() ;

                    for ( index_t tRoot : tActiveList )
                    {
                        if ( tVisited.test( tRoot ) )
                        {
                            continue ;
                        }

                        // BFS over zero-coefficient edges; nonzero edges are
                        // boundary candidates ( re-checked after K is complete )
                        tComponent.clear() ;
                        tCandidates.clear() ;
                        bool tCapContact = false ;
                        bool tStructureContact = false ;

                        index_t tHead = 0, tTail = 0 ;
                        tQueue( tTail++ ) = tRoot ;
                        tVisited.set( tRoot );
                        tInK.set( tRoot );

                        while ( tHead < tTail )
                        {
                            const index_t tU = tQueue( tHead++ );
                            tComponent.push( tU );
                            Node * tUNode = tNodes( tU );

                            if ( tUNode->is_periodic() )    tCapContact = true ;
                            if ( tStructureNodes.test( tU ) ) tStructureContact = true ;

                            // master + slave-partner adjacency, as in
                            // fire_node_coboundary()
                            for ( uint tPass = 0; tPass < 2; ++tPass )
                            {
                                Node * tSide = tUNode ;
                                if ( tPass == 1 )
                                {
                                    if ( ! tUNode->is_periodic() || tUNode->periodic() == tUNode )
                                    {
                                        break ;
                                    }
                                    tSide = tUNode->periodic() ;
                                }

                                for ( uint j = 0; j < tSide->number_of_edges(); ++j )
                                {
                                    Edge * tEdge = tSide->edge( j );
                                    if ( ! tComplexEdges.test( tEdge->index() ) ) continue ;
                                    if ( tPass == 1 && tEdge->is_flagged() ) continue ;

                                    const int tC = tGenerator->getCoefficient( tEdge->index() );

                                    if ( tC == 0 )
                                    {
                                        const index_t tV0 = tRep( tEdge->node( 0 ) );
                                        const index_t tV1 = tRep( tEdge->node( 1 ) );
                                        const index_t tV  = ( tV0 == tU ) ? tV1 : tV0 ;

                                        if ( ! tVisited.test( tV ) )
                                        {
                                            tVisited.set( tV );
                                            tInK.set( tV );
                                            tQueue( tTail++ ) = tV ;
                                        }
                                    }
                                    else if ( ! tEdgeSeen.test( tEdge->index() ) )
                                    {
                                        tEdgeSeen.set( tEdge->index() );
                                        tCandidates.push( tEdge->index() );
                                    }
                                }
                            }
                        }

                        // classify the boundary: candidates with both reps in
                        // K are interior nonzero edges, not boundary
                        index_t tNumBoundary = 0 ;
                        bool tPlusOK  = true ;
                        bool tMinusOK = true ;
                        bool tOverlap = false ;

                        for ( index_t tEdgeIndex : tCandidates )
                        {
                            Edge * tEdge = tEdges( tEdgeIndex );
                            const index_t tU0 = tRep( tEdge->node( 0 ) );
                            const index_t tU1 = tRep( tEdge->node( 1 ) );
                            const bool tIn0 = tInK.test( tU0 );
                            const bool tIn1 = tInK.test( tU1 );

                            if ( tIn0 == tIn1 )
                            {
                                continue ; // interior nonzero edge
                            }

                            ++tNumBoundary ;
                            const int tC = tGenerator->getCoefficient( tEdgeIndex );

                            // flip with sign s changes this edge by -s if the
                            // K side is node( 0 ), by +s if it is node( 1 );
                            // cancellation requires c == s resp. c == -s
                            const int tRequired = tIn0 ? tC : -tC ;
                            if ( tRequired != 1  ) tPlusOK  = false ;
                            if ( tRequired != -1 ) tMinusOK = false ;

                            // footprint conservatism: the outside endpoint also
                            // counts for structure / cap contact
                            const index_t tOut = tIn0 ? tU1 : tU0 ;
                            if ( tStructureNodes.test( tOut ) )   tStructureContact = true ;
                            if ( tNodes( tOut )->is_periodic() )  tCapContact = true ;

                            // overlap with any other generator's support
                            for ( Cochain * tOther : mGenerators( 1 ) )
                            {
                                if ( tOther != tGenerator &&
                                     tOther->getCoefficient( tEdgeIndex ) != 0 )
                                {
                                    tOverlap = true ;
                                    break ;
                                }
                            }
                        }

                        const bool tPure  = ( tNumBoundary > 0 ) && ( tPlusOK || tMinusOK );
                        const bool tSmall = 2 * tComponent.size() <= tNumActive ;
                        const bool tTierA = tPure && tSmall &&
                                            ! tStructureContact && ! tCapContact && ! tOverlap ;

                        if ( aFire && tTierA )
                        {
                            const int tSign = tPlusOK ? 1 : -1 ;
                            for ( index_t tIndex : tComponent )
                            {
                                cohomology::fire_node_coboundary(
                                        tGenerator, tNodes( tIndex ), tSign, tComplexEdges );
                            }

                            ++tNumFlipped ;

                            // reset scratch and rebuild components ( Rule 5 )
                            for ( index_t tIndex : tComponent )  tInK.reset( tIndex );
                            for ( index_t tIndex : tCandidates ) tEdgeSeen.reset( tIndex );
                            return true ;
                        }

                        for ( index_t tIndex : tComponent )  tInK.reset( tIndex );
                        for ( index_t tIndex : tCandidates ) tEdgeSeen.reset( tIndex );
                    }

                    return false ;
                };

                if ( aFireTierA )
                {
                    // fire one Tier-A pocket per pass, rebuilding the
                    // components after each accepted flip, to a fixed point
                    while ( tScan( true ) ) {}
                }

                // unit gate re-run ( rules, Rule 6 )
                for ( const auto & [ tIndex, tCoeff ] : tGenerator->getSimplicesMap() )
                {
                    BELFEM_ERROR( tCoeff <= 1 && tCoeff >= -1,
                            "remove_cut_pockets: generator %u has non-unit coefficient %d on edge index %lu after pocket removal",
                            tGeneratorCount,
                            tCoeff,
                            ( long unsigned int ) tIndex );
                }

                gLog.message( InfoLevel::Verbose,
                        "    cohomology generator %u : removed %lu cut pockets",
                        tGeneratorCount,
                        ( long unsigned int ) tNumFlipped );

                if ( mProgress != nullptr )
                {
                    mProgress->step();
                }
                ++tGeneratorCount ;
            }
        }

//-----------------------------------------------------------------------------

        void
        Cohomology::check()
        {

            if (mMesh->number_of_dimensions() == 3)
            {
                for (uint i = 0 ; i < this->get_Generators()(1).size(); ++i)
                {
                    auto tMap = this->get_Generators()(1)(i)->getSimplicesMap() ;
                    for (Face * tFace : mMesh->faces())
                    {
                        int tCirc = 0;
                        for(uint j = 0; j < tFace->number_of_edges(); ++j)
                        {
                            //std::cout << (tFace->edge_direction(j)?1.0:-1.0)*tMap[tFace->edge(j)->index()] << std::endl ;
                            tCirc += (tFace->edge_direction(j) ? 1 : -1)*tMap[tFace->edge(j)->index()] ;
                        }
                        if (tCirc!= 0 && tFace->is_flagged())
                        {
                            BELFEM_ERROR(false, "Invalid cohomology computed") ;
                        }
                    }
                }
            }
            else
            {
                for (uint i = 0 ; i < this->get_Generators()(1).size(); ++i)
                {
                    auto tMap = this->get_Generators()(1)(i)->getSimplicesMap() ;
                    for (Element * tFace : mMesh->elements())
                    {
                        real tCirc = 0.0;
                        for(uint j = 0; j < tFace->number_of_edges(); ++j)
                        {
                            tCirc += (tFace->edge_direction(j)?1.0:-1.0)*tMap[tFace->edge(j)->index()] ;
                        }
                        if (tCirc!= 0 && tFace->is_flagged())
                        {
                            BELFEM_ERROR(false, "Invalid cohomology computed") ;
                        }
                    }
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        Cohomology::create_kGeneratorsField(const uint k, Mesh * aMesh, string aFieldName )
        {
            uint tCount = 0;
            Cell < Element * > & tElements = aMesh->elements();
            for(const auto& tCochain : mGenerators(k))
            {
                tCount++;
                aMesh->create_field( aFieldName+std::to_string(tCount), EntityType::ELEMENT);
                OrderedMap< index_t, int > & tSimplicesMap = tCochain->getSimplicesMap();
                for( const auto& [tID, tCoeff] :  tSimplicesMap)
                {
                    aMesh->field_data(aFieldName+std::to_string(tCount))(tElements(tID)->index()) = tCoeff;
                    if ( mMesh->edges()(tID)->is_periodic() )
                    {
                        aMesh->field_data(aFieldName+std::to_string(tCount))(tElements(mMesh->edges()(tID)->periodic()->index())->index()) = tCoeff;
                    }
                }
            }
        }

//-----------------------------------------------------------------------------

        Cell< Matrix< int > > &
        Cohomology::get_CoboundaryMatrix()
        {
            return mD;
        }

//-----------------------------------------------------------------------------

        Cell <Cell< Cochain * >> &
        Cohomology::get_Generators()
        {
            return mGenerators;
        }

//-----------------------------------------------------------------------------

        void
        Cohomology::updatekGeneratorsFromHomology(Cell< Chain * > & tkGenerators, const uint k)
        {
            const uint n = mGenerators(k).size();
            const uint m = tkGenerators.size();
            Matrix< int > tTransInv = Matrix< int >(n,n,0);

            // ---------------------------------------------------------------
            // Step 1: Compute coupling matrix between cohomology generators
            //         and suggested homology generators (n x m)
            // ---------------------------------------------------------------
            Matrix<int> tTemp = Matrix< int >(n,m,0);
            for(uint i = 0; i < n; i++)
            {
                Cochain * tCoGen = mGenerators(k)(i);
                for(uint j = 0; j < m; j++)
                {
                    tTemp(i,j) = tCoGen->operator()(tkGenerators(j));
                }
            }
            // ---------------------------------------------------------------
            // Step 2: Identify unique columns in the coupling matrix.
            //         When terminals are shared between conditions (e.g.
            //         T1->T3, T2->T4, T3->T4), multiple homology generators
            //         reference the same terminal surface, producing identical
            //         coupling columns. We deduplicate to avoid rank deficiency.
            // ---------------------------------------------------------------
            Cell<uint> tColToUnique(m, 0);        // maps each column to its unique group index
            Cell<uint> tUniqueColRepresentative;   // original column index for each unique group
            uint tNumUnique = 0;

            for(uint j = 0; j < m; j++)
            {
                bool tFound = false;
                for(uint u = 0; u < tUniqueColRepresentative.size(); u++)
                {
                    uint tOrigIdx = tUniqueColRepresentative(u);
                    bool tMatch = true;
                    for(uint i = 0; i < n; i++)
                    {
                        if(tTemp(i, j) != tTemp(i, tOrigIdx))
                        {
                            tMatch = false;
                            break;
                        }
                    }
                    if(tMatch)
                    {
                        tColToUnique(j) = u;
                        tFound = true;
                        break;
                    }
                }
                if(!tFound)
                {
                    tColToUnique(j) = tNumUnique;
                    tUniqueColRepresentative.push(j);
                    tNumUnique++;
                }
            }

            // ---------------------------------------------------------------
            // Step 3: Build terminal coupling matrix F (n x tNumUnique)
            //         from the unique columns of tTemp
            // ---------------------------------------------------------------
            Matrix<int> tF = Matrix< int >(n, tNumUnique, 0);
            for(uint u = 0; u < tNumUnique; u++)
            {
                uint tOrigIdx = tUniqueColRepresentative(u);
                for(uint i = 0; i < n; i++)
                {
                    tF(i, u) = tTemp(i, tOrigIdx);
                }
            }

            // ---------------------------------------------------------------
            // Step 4: Build condition matrix D (aNumConditions x tNumUnique)
            //
            //   In 3D, generators are ordered [in0, out0, in1, out1, ...]:
            //     condition i -> D(i, unique(2i)) += +1, D(i, unique(2i+1)) -= 1
            //
            //   In 2D, generators are ordered [cond0, cond1, ...]:
            //     condition i -> D(i, unique(i)) += +1
            //
            //   The += / -= handles the case where shared terminals
            //   map to the same unique column, correctly summing contributions.
            // ---------------------------------------------------------------
            const bool tIs3D = (mMesh->number_of_dimensions() == 3);
            const uint aNumConditions = tIs3D ? m / 2 : m;

            // Each condition consumes one cohomology generator, so asking for
            // more conditions than the mesh topology provides cannot be
            // satisfied. Caught here because the assembly loop in step 6 would
            // otherwise walk off the end of the n x n transformation matrix and
            // surface as an out-of-bounds abort from the matrix backend.
            BELFEM_ERROR( aNumConditions <= n,
                "The model defines %u current or voltage conditions, but the topology only provides %u "
                "cohomology generator%s, so at most %u condition%s can be imposed. Every galvanically "
                "connected conductor carries ONE condition: conductors joined by a bulk domain, by "
                "solder or by a shared terminal count as one, and a periodic model spends one further "
                "generator on the free axial loop. Either merge the terminals that belong to the same "
                "conductor into a single bracketed group, or separate the conductors in the topology. "
                "Note that outside a bracket every id, including each member of an a:b range, becomes "
                "its own condition.",
                ( unsigned int ) aNumConditions,
                ( unsigned int ) n,
                n == 1 ? "" : "s",
                ( unsigned int ) n,
                n == 1 ? "" : "s" );

            Matrix<int> tD = Matrix< int >(aNumConditions, tNumUnique, 0);
            for(uint i = 0; i < aNumConditions; i++)
            {
                if(tIs3D)
                {
                    uint tInGenIdx  = 2 * i;
                    uint tOutGenIdx = 2 * i + 1;
                    tD(i, tColToUnique(tInGenIdx))  += 1;
                    tD(i, tColToUnique(tOutGenIdx)) -= 1;
                }
                else
                {
                    tD(i, tColToUnique(i)) += 1;
                }
            }

            // ---------------------------------------------------------------
            // Step 5: Compute A = D * F^{dagger} via Smith Normal Form
            //         smithForm modifies tF in-place to become the diagonal Sigma
            // ---------------------------------------------------------------
            auto [tQ, tQ_, tR, tR_, s, t] = smithForm(tF);

            // t is the rank of the terminal coupling matrix. Step 6 writes one
            // row per condition and then one row per free cut, n - t of them,
            // so it stays inside the n x n transformation exactly while
            // aNumConditions <= t. Rank deficiency means the conditions are not
            // independent even though there are few enough of them, which the
            // count check above cannot see.
            BELFEM_ERROR( aNumConditions <= t,
                "The %u current or voltage conditions are not independent: they couple to the "
                "cohomology generators with rank %u, so only %u of them carry distinct "
                "information. Two conditions that drive the same conductor, an input and an "
                "output terminal that are not connected through a conductor, or terminals that "
                "share a surface all produce this.",
                ( unsigned int ) aNumConditions,
                ( unsigned int ) t,
                ( unsigned int ) t );

            Matrix<int> tFTrans = trans(tF);
            Matrix<int> tPseudoInv = Matrix<int>(
                tD.matrix_data() * tR.matrix_data() * tFTrans.matrix_data() * tQ_.matrix_data());

            // ---------------------------------------------------------------
            // Step 6: Build the full transformation matrix
            // ---------------------------------------------------------------

            // First rows: terminal-constrained generators from pseudo-inverse
            uint tCount = 0;
            for(uint i = 0; i < tPseudoInv.n_rows(); i++)
            {
                for (uint j = 0; j < n ; ++ j)
                {
                    tTransInv(tCount,j) = tPseudoInv(i,j);
                }
                tCount++;
            }

            // Add generators from the left null space of F (free cuts not linked to any terminal)
            for(uint i = t ; i < tQ_.n_cols(); ++ i)
            {
                for (uint j = 0; j < n ; ++ j)
                {
                    tTransInv(tCount,j) = tQ_(i,j);
                }
                tCount++;
            }

            // Complete the basis with the kernel of the assembled rows (also free cuts).
            // Basis choice: this takes the kernel of the PARTIALLY assembled
            // transformation; taking the kernel of the final pseudo-inverse is
            // also possible, and both choices can misbehave in rare free-cut
            // configurations.
            auto [w,v] = kernelImage(tTransInv);
            BELFEM_ERROR(w.n_cols()+tCount==n,
                         "Inconsistency in the definition of the terminals. Either some constraints are linearly dependent, "
                         "or some input/output terminals are not connected through a conductor. Please review the terminal definition");

            for(uint i = 0; i < w.n_cols(); i++)
            {
                for (uint j = 0; j < n ; ++ j)
                {
                    tTransInv(tCount,j) = w(j,i);
                }
                tCount++;
            }

            // ---------------------------------------------------------------
            // Step 7: Create new generators as linear combinations of the
            //         initial ones using the transformation matrix
            // ---------------------------------------------------------------
            Cell< Cochain * > tNewGenerators = Cell< Cochain * >();
            for(uint i = 0; i < n; i++)
            {
                Cochain * tGen = new Cochain(k, mMesh, false, false);
                for(uint j = 0; j < n; j++)
                {
                    tGen->addCochainToCochain(mGenerators(k)(j),tTransInv(i,j));
                }
                tNewGenerators.push(tGen);
            }

            for(uint i = 0; i < n; i++)
            {
                delete mGenerators(k)(i);
            }
            mGenerators(k) = tNewGenerators;

            mflagProp = true;
        }

//-----------------------------------------------------------------------------

        Matrix< real >
        Cohomology::coefficient_TMatrix(Cell< Chain * > & tkGenerators, const uint k)
        {
            Cell < int > tCoeff = Cell < int >(tkGenerators.size(),0);
            const uint n = mGenerators(k).size();
            BELFEM_ASSERT(tkGenerators.size() == n, "Error, not the same number of k-generators (%u vs %u)",
                ( unsigned int ) tkGenerators.size(), ( unsigned int ) n  );
            Matrix< real > tTransMat = Matrix< real >(mGenerators(k).size(),mGenerators(k).size(),0);

            //Populate the transformation matrix
            for(uint i = 0; i < n; i++)
            {
                Cochain * tCoGen = mGenerators(k)(i);
                for(uint j = 0; j < n; j++)
                {
                    tTransMat(j,i) = tCoGen->operator()(tkGenerators(j));
                }
            }

            Matrix< real > tTMatrix = inv(tTransMat);


            return tTMatrix;
        }

//-----------------------------------------------------------------------------
    }
}
//------------------------------------------------------------------------------
