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

#include "cl_SimplicialComplex.hpp"

#include "cl_Timer.hpp"
#include "cl_Logger.hpp"

//#define PERFORMANCE_CHECK

namespace belfem {
//------------------------------------------------------------------------------

    namespace mesh {
        SimplicialComplex::SimplicialComplex( Mesh *aMesh, bool aPeriodicity ) :
            mMesh( aMesh )
        {
            mChainsMap.set_size(4, Map<index_t, Chain *>());
            mCochainsMap.set_size(4, Map<index_t, Cochain *>());
            this->create_complex( aMesh, aPeriodicity );
        }

//------------------------------------------------------------------------------

        SimplicialComplex::~SimplicialComplex()
        {
            this->reset();
        }

//------------------------------------------------------------------------------

        void
        SimplicialComplex::reset()
        {
            // delete pointers
            for (const auto &[tKey, tChain]: mChainsMap(0))
            {
                delete tChain;
            }

            for (const auto &[tKey, tCochain]: mCochainsMap(0))
            {
                delete tCochain;
            }

            // delete map
            mChainsMap(0).clear();
            mCochainsMap(0).clear();

            // delete pointers
            for (const auto &[tKey, tChain]: mChainsMap(1))
            {
                delete tChain;
            }

            for (const auto &[tKey, tCochain]: mCochainsMap(1))
            {
                delete tCochain;
            }

            // delete map
            mChainsMap(1).clear();
            mCochainsMap(1).clear();

            // delete pointers
            for (const auto &[tKey, tChain]: mChainsMap(2))
            {
                delete tChain;
            }

            for (const auto &[tKey, tCochain]: mCochainsMap(2))
            {
                delete tCochain;
            }

            // delete map
            mChainsMap(2).clear();
            mCochainsMap(2).clear();

            // delete pointers
            for (const auto &[tKey, tChain]: mChainsMap(3))
            {
                delete tChain;
            }

            for (const auto &[tKey, tCochain]: mCochainsMap(3))
            {
                delete tCochain;
            }

            // delete map
            mChainsMap(3).clear();
            mCochainsMap(3).clear();

            delete mOriginalEdges ;
        }

//------------------------------------------------------------------------------

        void
        SimplicialComplex::create_complex( Mesh * aMesh, bool aPeriodicity )
        {
            // restore factory settings
            this->reset();

            uint tDim = aMesh->number_of_dimensions();

            //First, unflag stuff on periodic slave (if they exist)
            if (mMesh->has_periodicity() && aPeriodicity)
            {
                for (Node * tNode : mMesh->periodicity()->slave_nodes())
                {
                    tNode->unflag() ;
                }
                for (Edge * tEdge : mMesh->periodicity()->slave_edges())
                {
                    tEdge->unflag() ;
                }
                for (Face * tFace : mMesh->periodicity()->slave_faces())
                {
                    tFace->unflag() ;
                }
            }

            // loop over all nodes on mesh
            for ( Node *tNode: aMesh->nodes() )
            {
                if ( tNode->is_flagged() )
                {
                    Chain *tChain = new Chain(0, aMesh, false, false);
                    tChain->addSimplexToChain(tNode->index(), 1);

                    Cochain *tCochain = new Cochain(0, aMesh, false, false);
                    tCochain->addSimplexToCochain(tNode->index(), 1);

                    for (uint i = 0; i < tNode->number_of_edges(); i++)
                    {
                        if (tNode->edge(i)->is_flagged())
                        {
                            if (tNode->edge(i)->node(0)->index() == tNode->index())
                            {
                                tCochain->add_simplex_to_coboundary(tNode->edge(i)->index(), -1);
                                tChain->add_simplex_to_coboundary(tNode->edge(i)->index(), -1);
                            }
                            else
                            {
                                tCochain->add_simplex_to_coboundary(tNode->edge(i)->index(), 1);
                                tChain->add_simplex_to_coboundary(tNode->edge(i)->index(), 1);
                            }
                        }
                    }

                    // Add the periodic side if it exists
                    if ( tNode->is_periodic() && aPeriodicity )
                    {
                        /*if (tNode->id() == 2571)
                        {
                            std::cout << "Clone : " << tNode->periodic()->id() << std::endl;
                            std::cout << tNode->periodic()->number_of_edges() << std::endl ;
                            for (uint i = 0; i < tNode->periodic()->number_of_edges(); i++)
                            {
                                std::cout << "Edge " << i << " : " << tNode->periodic()->edge(i)->node(0)->id() << "->" <<  tNode->periodic()->edge(i)->node(1)->id()<< std::endl;
                            }
                        }*/
                        for (uint i = 0; i < tNode->periodic()->number_of_edges(); i++)
                        {
                            if (tNode->periodic()->edge(i)->is_flagged() )
                            {
                                if (tNode->periodic()->edge(i)->node(0)->index() == tNode->periodic()->index())
                                {
                                    tCochain->add_simplex_to_coboundary(tNode->periodic()->edge(i)->index(), -1);
                                    tChain->add_simplex_to_coboundary(tNode->periodic()->edge(i)->index(), -1);
                                }
                                else
                                {
                                    tCochain->add_simplex_to_coboundary(tNode->periodic()->edge(i)->index(), 1);
                                    tChain->add_simplex_to_coboundary(tNode->periodic()->edge(i)->index(), 1);
                                }
                            }
                        }
                    }


                    // add entry to map
                    mChainsMap(0)[tNode->index()] = tChain;
                    mCochainsMap(0)[tNode->index()] = tCochain;
                }
            }

            Cell< Edge * > & tEdges = aMesh->edges() ;

            // loop over all edges on mesh
            for (Edge *tEdge: tEdges )
            {
                if (tEdge->is_flagged())
                {
                    Chain *tChain = new Chain(1, aMesh, false, false);
                    tChain->addSimplexToChain(tEdge->index(), 1);

                    Cochain *tCochain = new Cochain(1, aMesh, false, false);
                    tCochain->addSimplexToCochain(tEdge->index(), 1);

                    if (tEdge->node(0)->is_flagged())
                    {
                        tChain->add_simplex_to_boundary(tEdge->node(0)->index(), -1);
                        tCochain->add_simplex_to_boundary(tEdge->node(0)->index(), -1);
                    }
                    else if (tEdge->node(0)->is_periodic() && aPeriodicity && tEdge->node(0)->periodic()->is_flagged())
                    {
                        tChain->add_simplex_to_boundary(tEdge->node(0)->periodic()->index(), -1);
                        tCochain->add_simplex_to_boundary(tEdge->node(0)->periodic()->index(), -1);
                    }
                    if (tEdge->node(1)->is_flagged())
                    {
                        tChain->add_simplex_to_boundary(tEdge->node(1)->index(), 1);
                        tCochain->add_simplex_to_boundary(tEdge->node(1)->index(), 1);
                    }
                    else if (tEdge->node(1)->is_periodic() && aPeriodicity && tEdge->node(1)->periodic()->is_flagged())
                    {
                        tChain->add_simplex_to_boundary(tEdge->node(1)->periodic()->index(), 1);
                        tCochain->add_simplex_to_boundary(tEdge->node(1)->periodic()->index(), 1);
                    }

                    if (tDim == 2)
                    {
                        for (uint i = 0; i < tEdge->number_of_elements(); ++i)
                        {
                            int tMul = 0;
                            Element *tElement = tEdge->element(i);

                            if (tElement->is_flagged())
                            {
                                if(tElement->number_of_edges()!= 1)
                                {
                                    for (uint j = 0; j < tElement->number_of_edges(); j++)
                                    {
                                        if (tElement->edge(j)->index() == tEdge->index())
                                        {
                                            tMul = (tElement->edge_direction(j) ? 1 : -1);
                                            break;
                                        }
                                    }
                                    tCochain->add_simplex_to_coboundary(tElement->index(), tMul);
                                    tChain->add_simplex_to_coboundary(tElement->index(), tMul);
                                }
                            }
                        }
                    }
                    else
                    {
                        for (uint i = 0; i < tEdge->number_of_faces(); ++i)
                        {
                            int tMul = 0;
                            Face *tFace = tEdge->face(i);

                            if (tFace->is_flagged())
                            {
                                for (uint j = 0; j < tFace->number_of_edges(); j++)
                                {
                                    if (tFace->edge(j)->index() == tEdge->index())
                                    {
                                        tMul = (tFace->edge_direction(j) ? 1 : -1);
                                        break;
                                    }
                                }
                                tCochain->add_simplex_to_coboundary(tFace->index(), tMul);
                                tChain->add_simplex_to_coboundary(tFace->index(), tMul);
                            }
                        }
                    }

                    // Add the periodic side if it exists
                    if ( tEdge->is_periodic() && aPeriodicity )
                    {
                        for (uint i = 0; i < tEdge->periodic()->number_of_faces(); ++i)
                        {
                            int tMul = 0;
                            Face *tFace = tEdge->periodic()->face(i);

                            if (tFace->is_flagged())
                            {
                                for (uint j = 0; j < tFace->number_of_edges(); j++)
                                {
                                    if (tFace->edge(j)->index() == tEdge->periodic()->index())
                                    {
                                        tMul = (tFace->edge_direction(j) ? 1 : -1);
                                        break;
                                    }
                                }
                                tCochain->add_simplex_to_coboundary(tFace->index(), tMul);
                                tChain->add_simplex_to_coboundary(tFace->index(), tMul);
                            }
                        }
                    }

                    // add entry to map
                    mChainsMap(1)[tEdge->index()] = tChain;
                    mCochainsMap(1)[tEdge->index()] = tCochain;

                }
            }

            if (tDim == 3)
            {
                for (Face *tFace: aMesh->faces())
                {
                    if (tFace->is_flagged())
                    {
                        Chain *tChain = new Chain(2, aMesh, false, false);
                        tChain->addSimplexToChain(tFace->index(), 1);
                        for (uint i = 0; i < tFace->number_of_edges(); ++i)
                        {
                            if (tFace->edge(i)->is_flagged())
                            {
                                tChain->add_simplex_to_boundary(tFace->edge(i)->index(),
                                                                (tFace->edge_direction(i) ? 1 : -1));
                            }
                            else if (tFace->edge(i)->is_periodic() && aPeriodicity &&  tFace->edge(i)->periodic()->is_flagged() )
                            {
                                tChain->add_simplex_to_boundary(tFace->edge(i)->periodic()->index(),
                                                                (tFace->edge_direction(i) ? 1 : -1));
                            }
                        }

                        Cochain *tCochain = new Cochain(2, aMesh, false, false);
                        tCochain->addSimplexToCochain(tFace->index(), 1);
                        for (uint i = 0; i < tFace->number_of_edges(); ++i)
                        {
                            if (tFace->edge(i)->is_flagged())
                            {
                                tCochain->add_simplex_to_boundary(tFace->edge(i)->index(),
                                                                  (tFace->edge_direction(i) ? 1 : -1));
                            }
                            else if (tFace->edge(i)->is_periodic() && aPeriodicity &&  tFace->edge(i)->periodic()->is_flagged() )
                            {
                                tCochain->add_simplex_to_boundary(tFace->edge(i)->periodic()->index(),
                                                                (tFace->edge_direction(i) ? 1 : -1));
                            }
                        }

                        if (tFace->master())
                        {
                            if (tFace->master()->is_flagged())
                            {
                                tCochain->add_simplex_to_coboundary(tFace->master()->index(), 1);
                                tChain->add_simplex_to_coboundary(tFace->master()->index(), 1);
                            }
                        }
                        if (tFace->slave())
                        {
                            if (tFace->slave()->is_flagged())
                            {
                                tCochain->add_simplex_to_coboundary(tFace->slave()->index(), -1);
                                tChain->add_simplex_to_coboundary(tFace->slave()->index(), -1);
                            }
                        }

                        // Add the periodic side if it exists
                        if ( tFace->is_periodic() && aPeriodicity )
                        {
                            if (tFace->periodic()->master())
                            {
                                if (tFace->periodic()->master()->is_flagged())
                                {
                                    tCochain->add_simplex_to_coboundary(tFace->periodic()->master()->index(), 1);
                                    tChain->add_simplex_to_coboundary(tFace->periodic()->master()->index(), 1);
                                }
                            }
                            if (tFace->periodic()->slave())
                            {
                                if (tFace->periodic()->slave()->is_flagged())
                                {
                                    tCochain->add_simplex_to_coboundary(tFace->periodic()->slave()->index(), -1);
                                    tChain->add_simplex_to_coboundary(tFace->periodic()->slave()->index(), -1);
                                }
                            }
                        }

                        // add entry to map
                        mChainsMap(2)[tFace->index()] = tChain;
                        mCochainsMap(2)[tFace->index()] = tCochain;

                    }
                }
            }

            // loop over all elements on mesh
            for (Element *tElement: aMesh->elements())
            {
                if (tElement->is_flagged() && tElement->dimension() == tDim)
                {
                    Chain *tChain = new Chain(tDim, aMesh, false, false);
                    tChain->addSimplexToChain(tElement->index(), 1);
                    if (tDim == 3)
                    {
                        for (uint i = 0; i < tElement->number_of_faces(); ++i)
                        {
                            if (tElement->face(i)->is_flagged())
                            {
                                tChain->add_simplex_to_boundary(tElement->face(i)->index(),
                                                                (tElement->face(i)->master()->index() == tElement->index() ? 1: -1));
                            }
                            else if ( tElement->face(i)->is_periodic() && aPeriodicity && tElement->face(i)->periodic()->is_flagged())
                            {
                                tChain->add_simplex_to_boundary(tElement->face(i)->periodic()->index(),
                                                                (tElement->face(i)->slave()->index() == tElement->index() ? 1: -1));
                            }
                        }

                    }
                    else
                    {
                        for (uint i = 0; i < tElement->number_of_edges(); ++i)
                        {
                            if (tElement->edge(i)->is_flagged())
                            {
                                tChain->add_simplex_to_boundary(tElement->edge(i)->index(),
                                                                (tElement->edge_direction(i) ? 1 : -1));
                            }
                        }
                    }

                    Cochain *tCochain = new Cochain(tDim, aMesh, false, false);
                    tCochain->addSimplexToCochain(tElement->index(), 1);
                    if (tDim == 3)
                    {
                        for (uint i = 0; i < tElement->number_of_faces(); ++i)
                        {
                            if (tElement->face(i)->is_flagged())
                            {
                                tCochain->add_simplex_to_boundary(tElement->face(i)->index(),
                                                                  (tElement->face(i)->master()->index() == tElement->index() ? 1: -1));
                            }
                            else if ( tElement->face(i)->is_periodic() && aPeriodicity && tElement->face(i)->periodic()->is_flagged())
                            {
                                tCochain->add_simplex_to_boundary(tElement->face(i)->periodic()->index(),
                                                                (tElement->face(i)->slave()->index() == tElement->index() ? 1: -1));
                            }
                        }
                    }
                    else
                    {
                        for (uint i = 0; i < tElement->number_of_edges(); ++i)
                        {
                            if (tElement->edge(i)->is_flagged())
                            {
                                tCochain->add_simplex_to_boundary(tElement->edge(i)->index(),
                                                                  (tElement->edge_direction(i) ? 1 : -1));
                            }
                        }
                    }

                    // add entry to map
                    mChainsMap(tDim)[ tElement->index() ] = tChain;
                    mCochainsMap(tDim)[ tElement->index() ] = tCochain;

                }
            }
            // Store original 1-simplex indices before any reduction
            mOriginalEdges = new DynamicBitset( aMesh->number_of_edges() ) ;
            for ( const auto & [tKey, tCochain] : mCochainsMap( 1 ) )
            {
                mOriginalEdges->set( tKey ) ;
            }
        }


//------------------------------------------------------------------------------

        // CCR algorithm as described in Computational Homology from T. Kaczynski et al.
        void
        SimplicialComplex::reduce_complexCCR()
        {
            Timer tTimer;

            gLog.message( InfoLevel::Detailed, "    reducing the complex (ccr) ... " );

        // loop over all dimensions
            for ( int p = 3; p > 0; --p )
            {
                this->pGeneralizedCombine(p);
            }

            gLog.message( InfoLevel::Detailed, "    ... reduction time %.3f s\n", tTimer.stop() * 1e-3 );

            for ( uint k = 0; k < 4; ++k )
            {
                gLog.message( InfoLevel::Detailed, "    number of %u-chains: %u",
                              ( unsigned int ) k,
                              ( unsigned int ) this->number_of_ksimplices( k ));
            }
            gLog.message( InfoLevel::Detailed, "" );
        }

//-----------------------------------------------------------------------------

        // Cohomology variant of the CCR algorithm
        void
        SimplicialComplex::coreduce_complexCCR()
        {
            Timer tTimer;

            gLog.message( InfoLevel::Detailed, "    coreducing the complex (ccr) ... " );

            // loop over all dimensions
            for ( int p = 0; p < 3; ++p )
            {
                this->pGeneralizedCocombine(p) ;
            }

            gLog.message( InfoLevel::Detailed, "    ... coreduction time %.3f s\n", tTimer.stop() * 1e-3 );

            for ( uint k = 0; k < 4; ++k )
            {
                gLog.message( InfoLevel::Detailed, "    number of %u-cochains: %u",
                              ( unsigned int ) k,
                              ( unsigned int ) this->number_of_kcosimplices( k ));
            }
            gLog.message( InfoLevel::Detailed, "" );

        }

//-----------------------------------------------------------------------------

        // CCR algorithm as described in Computational Homology from T. Kaczynski et al. (Old implementation)
        void
        SimplicialComplex::reduce_complexCCR_old()
        {
            Timer tTimer;

            gLog.message( InfoLevel::Detailed, "    reducing the complex (ccr) ... " );

            // loop over all dimensions
            for ( int k = 3; k > 0; --k )
            {
                bool found = true;
                while (found) {
                    found = false;

                    // loop over all k-chains and (k-1)-chains of the complex
                    for (auto it = mChainsMap(k).begin(); it != mChainsMap(k).end(); ++it)
                    {
                        for (const auto &[a, tChain2]: mChainsMap(k - 1))
                        {
                            int tCoeff = it->second->getBoundary()->getCoefficient(a);
                            // Reduce if a is a boundary of b
                            if (abs(tCoeff) == 1)
                            {
                                this->reduce_pair(k, a, it, tCoeff);
                                found = true;
                                break;
                            }
                        }
                        if (found)
                        {
                            break;
                        }
                    }
                }
            }


                gLog.message( InfoLevel::Detailed, "    ... reduction time %.3f s\n", tTimer.stop() * 1e-3 );

            for ( uint k = 0; k < 4; ++k )
            {
                gLog.message( InfoLevel::Detailed, "    number of %u-chains: %u",
                              ( unsigned int ) k,
                              ( unsigned int ) this->number_of_ksimplices( k ));
            }
            gLog.message( InfoLevel::Detailed, "" );
        }

//-----------------------------------------------------------------------------

        // Cohomology variant of the CCR algorithm (old implementation)
        void
        SimplicialComplex::coreduce_complexCCR_old()
        {
            Timer tTimer;

            gLog.message( InfoLevel::Detailed, "    coreducing the complex (ccr) ... " );

            // loop over all dimensions
            //this->remove_kcochainFromMap(0, mCochainsMap(0).begin()->first);
            for ( int k = 2; k >= 0; --k )
            {
                bool found = true;
                while (found)
                {
                    found = false;

                    // loop over all k-cochains and (k+1)-cochains of the complex
                    for (const auto &[a, tCochain2]: mCochainsMap(k + 1))
                    {
                        for (auto it = mCochainsMap(k).begin(); it != mCochainsMap(k).end(); ++it)
                        {
                            int tCoeff = it->second->getCoboundary()->getCoefficient(a);
                            // Reduce if a is a coboundary of b
                            if (abs(tCoeff) == 1)
                            {
                                this->coreduce_pair(k, a, it, tCoeff);
                                found = true;
                                break;
                            }
                        }
                        if (found)
                        {
                            break;
                        }
                    }
                }
            }

            gLog.message( InfoLevel::Detailed, "    ... coreduction time %.3f s\n", tTimer.stop() * 1e-3 );

            for ( uint k = 0; k < 4; ++k )
            {
                gLog.message( InfoLevel::Detailed, "    number of %u-cochains: %u",
                              ( unsigned int ) k,
                              ( unsigned int ) this->number_of_kcosimplices( k ));
            }
            gLog.message( InfoLevel::Detailed, "" );

        }

//-----------------------------------------------------------------------------

        // Function to reduce a pair of chains (a,b), where a is a boundary of b
        void
        SimplicialComplex::reduce_pair( const uint k, const uint a, Map< index_t, Chain*>::const_iterator & bit, const int aCoeff )
        {
            Chain* bChain = bit->second;
            index_t b = bit->first;
            int val2;


            //loop over all the k-chains of the complex
            for ( auto it = bit; it!=mChainsMap( k  ).end(); ++it)
            {
                val2 = it->second->getBoundary()->getCoefficient( a );
                // Add the b to the k-chain if a is a boundary
                if ( abs( val2 ) == 1 and it->first != b)
                {
                    it->second->addChainToChain( bChain, -aCoeff*val2);
                }

            }

            // Remove a and b from the simplicial chain complex
            this->remove_kchainFromMap( k, b);
            this->remove_kchainFromMap( k - 1, a );
        }

//-----------------------------------------------------------------------------

        // Function to reduce a pair of cochains (a,b), where a is a coboundary of b
        void
        SimplicialComplex::coreduce_pair(const uint k, const uint a,
                                               Map< index_t, Cochain * >::const_iterator &bit, const int aCoeff) {
            Cochain *bCochain = bit->second;
            index_t b = bit->first;
            int val2;


            //loop over all the k-cochains of the complex
            for (auto it = bit; it != mCochainsMap(k).end(); ++it) {
                val2 = it->second->getCoboundary()->getCoefficient(a);

                // Add the b to the k-cochain if a is a coboundary
                if (abs(val2) == 1 and it->first != b) {
                    it->second->addCochainToCochain(bCochain, -aCoeff * val2);
                }

            }

            // Remove a and b from the simplicial cochain complex
            this->remove_kcochainFromMap(k + 1, a);
            this->remove_kcochainFromMap(k, b);

        }

//-----------------------------------------------------------------------------

        // pReduce from Pellikka et al.
        void
        SimplicialComplex::pReduce(const uint p)
        {
            if (p == 0)
            {
                return ;
            }
#ifdef PERFORMANCE_CHECK
            std::ofstream tFile;
            tFile.open ("Reduce.txt",std::ios_base::app);
            tFile << mNumLoopsReduce << " " << mChainsMap(0).size()+mChainsMap(1).size()+mChainsMap(2).size()+mChainsMap(3).size() << " \n";
#endif
            bool tRemoved = true;
            while (tRemoved)
            {
                tRemoved = false;
                mNumLoopsReduce+=1;

                //Loop over all the p-1 chains
                for (auto it  = mChainsMap( p-1 ).begin(), next_it = it;
                    it != mChainsMap( p-1 ).end(); it = next_it)
                {
                    ++next_it; //Work with 2 iterators because we remove entities from maps

                    //Check if the p-1 chain is a boundary of exactly one p chain
                    if (it->second->getCoboundary()->getSimplicesMap().size() == 1)
                    {
                        index_t a = it->first;
                        index_t b = it->second->getCoboundary()->getSimplicesMap().begin()->first;

                        //Update the coboundaries of the neighbors
                        Chain* tNeighbors = mChainsMap(p)(b)->getBoundary();
                        for(const auto [tID2, tCoeff2]: tNeighbors->getSimplicesMap())
                        {
                            if (tID2 != it->first )
                            {
                                mChainsMap(p-1)(tID2)->getCoboundary()->setCoefficient(b,0);
                            }
                        }
                        if (p > 1)
                        {
                            for(const auto [tID2, tCoeff2]:
                            mChainsMap(p-1)(a)->getBoundary()->getSimplicesMap())
                            {
                                mChainsMap(p-2)(tID2)->getCoboundary()->setCoefficient(a,0);
                            }
                        }

                        // Remove a and b from the complex
                        this->remove_kchainFromMap( p, b );
                        this->remove_kchainFromMap( p - 1, a );
                        tRemoved = true;
                    }
                }
#ifdef PERFORMANCE_CHECK
                tFile << mNumLoopsReduce << " " << mChainsMap(0).size()+mChainsMap(1).size()+mChainsMap(2).size()+mChainsMap(3).size() << " \n";
#endif
            }
#ifdef PERFORMANCE_CHECK
            tFile.close();
#endif
        }

//-----------------------------------------------------------------------------

        // pCombine from Pellikka et al.
        void
        SimplicialComplex::pCombine(const uint p)
        {
            index_t a;
            index_t b;
            index_t e;
            Cell< index_t > Q;

            //Loop over the p chains
            for (auto it = mChainsMap( p ).begin();  it != mChainsMap( p ).end(); ++it )
            {

                if (p == 0)
                {
                    return;
                }

                //Store the boundaries
                for (const auto & [tId, tCoeff]: it->second->getBoundary()->getSimplicesMap())
                {
                    Q.push(tId);
                }

                while (Q.size() > 0)
                {
                    a = Q.pop();
                    Chain *tBoundary = mChainsMap( p-1 )(a);
                    OrderedMap <index_t, int> & tSimplicesMap = tBoundary->getCoboundary()->getSimplicesMap();

                    // If the boundary is shared exactly by two neighbors
                    if (tSimplicesMap.size()== 2)
                    {
                        auto it2 = tSimplicesMap.begin();
                        auto it3 = ++tSimplicesMap.begin();
                        if (it2->first == it->first) //Always keep the simplex we're looking at
                        {
                            e = it2->first;
                            b = it3->first;
                        }
                        else
                        {
                            b = it2->first;
                            e = it3->first;
                        }

                        Chain* tChainAdd = mChainsMap(p)(e);
                        Chain* tChainReduce = mChainsMap(p)(b);
                        Chain* tNeighbors = tChainReduce->getBoundary();

                        int val1 = tChainAdd->getBoundary()->getCoefficient(a);
                        int val2 = tChainReduce->getBoundary()->getCoefficient(a);

                        //Add the chain
                        tChainAdd->addChainToChain(tChainReduce,-val1*val2);

                        //Update the coboundaries
                        for(const auto [tID2, tCoeff2]: tNeighbors->getSimplicesMap())
                        {
                            auto it4 = mChainsMap( p-1).find( tID2 );

                            if ( it4 != mChainsMap( p-1 ).end() )
                            {
                                it4->second->getCoboundary()->setCoefficient(b,0);
                                if (tID2 != a)
                                {
                                    it4->second->add_simplex_to_coboundary(e,-tCoeff2*val2*val1);
                                }
                            }
                        }

                        if (p > 1)
                        {
                            for(const auto [tID2, tCoeff2]: tBoundary->getBoundary()->getSimplicesMap())
                            {
                                mChainsMap(p-2)(tID2)->getCoboundary()->setCoefficient(a,0);
                            }
                        }

                        //Remove a and b from the complex
                        this->remove_kchainFromMap(p - 1, a);
                        this->remove_kchainFromMap(p, b);

                        //Add the new boundaries
                        for (const auto & [tId, tCoeff]: tChainAdd->getBoundary()->getSimplicesMap())
                        {
                            Q.push(tId);
                        }
                        unique(Q);
                    }
                }

            }
        }

//-----------------------------------------------------------------------------

        // CCR algorithm as described in Computational Homology from T. Kaczynski et al. for a given dimension p
        void
        SimplicialComplex::pGeneralizedCombine(const uint p)
        {
            // Loop over all k-chains
            for ( auto it = mChainsMap( p ).begin(), next_it = it; it != mChainsMap( p ).end(); it = next_it )
            {
                ++next_it ;
                //Get the boundary and check if it has at least one simplex
                Chain * tBoundary = it->second->getBoundary();
                OrderedMap< index_t, int > & tSimplicesMap = tBoundary->getSimplicesMap();
                index_t b = it->first;
                Chain* tChainReduce = it->second;

                if ( tSimplicesMap.size() != 0 )
                {
                    // Find the first simplex ID that correspond to a (k-1)-chain of the complex
                    auto it2 = tSimplicesMap.begin();
                    while ( !mChainsMap( p - 1 ).key_exists( it2->first ) or
                            abs( tChainReduce->getBoundary()->getCoefficient( it2->first )) != 1 )
                    {
                        ++it2;
                    }
                    if ( it2 == tSimplicesMap.end())
                    {
                        continue;
                    }
                    index_t a = it2->first;

                    // Find the neighbors of b through a
                    Chain * tNeighbors = mChainsMap( p - 1 )( a )->getCoboundary();
                    int val1 = tChainReduce->getBoundary()->getCoefficient( a );

                    //Loop over the Neighbors to update their chains
                    for ( const auto [tID, tCoeff]: tNeighbors->getSimplicesMap())
                    {
                        auto it3 = mChainsMap( p ).find( tID );

                        if ( it3 != mChainsMap( p ).end() and tID != b )
                        {
                            Chain * tChainAdd = mChainsMap( p )( tID );
                            int val2 = tChainAdd->getBoundary()->getCoefficient( a );

                            //Add the chain b to the neighbor
                            tChainAdd->addChainToChain( tChainReduce, -val1 * val2 );

                            //Update the coboundaries
                            for ( const auto [tID2, tCoeff2]: tSimplicesMap )
                            {
                                auto it4 = mChainsMap( p - 1 ).find( tID2 );

                                if ( it4 != mChainsMap( p - 1 ).end() )
                                {
                                    it4->second->getCoboundary()->setCoefficient( b,0 );
                                    if ( tID2 != a )
                                    {
                                        it4->second->add_simplex_to_coboundary( tID, -tCoeff2*val2*val1 );
                                    }
                                }
                            }
                            if ( p > 1 )
                            {
                                for ( const auto [tID2, tCoeff2]:
                                        mChainsMap( p - 1 )(a )->getBoundary()->getSimplicesMap())
                                {
                                    mChainsMap( p - 2 )( tID2 )
                                            ->getCoboundary()->setCoefficient( a,0 );
                                }
                            }
                        }
                    }
                    //Remove a and b from the complex
                    this->remove_kchainFromMap( p - 1, a );
                    this->remove_kchainFromMap( p, b );
                }

            }
        }

//-----------------------------------------------------------------------------

        // ReduceOmit from Pellikka et al.
        void
        SimplicialComplex::reduceOmit()
        {
            uint tDim = mChainsMap(3).size()==0?2:3;
            for (uint p = tDim; p >= 1; p--)
            {
                this->pReduce(p);
            }

            uint tKey;
            while (this->number_of_ksimplices(tDim) > 0)
            {
                tKey = mChainsMap(tDim).begin()->first;
                for(const auto [tID2, tCoeff2]: mChainsMap(tDim).begin()->second->getBoundary()->getSimplicesMap())
                {
                    mChainsMap(tDim-1)(tID2)->getCoboundary()->setCoefficient(tKey,0);
                }

                this->remove_kchainFromMap(tDim,tKey);
                for (uint p = tDim; p >= 1; --p)
                {
                    this->pReduce(p);
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        SimplicialComplex::reduce_complexPellikka()
        {
            Timer tTimer;

            gLog.message( InfoLevel::Detailed, "    reducing the complex (pellikka's algorithm) ... " );

            this->reduceOmit();
            for (uint p = 3; p >= 1; p--)
            {
                this->pCombine(p);
                this->pReduce(p-1);
            }

            gLog.message( InfoLevel::Detailed, "    ... reduction time %.3f s\n", tTimer.stop() * 1e-3 );

            for ( uint k = 0; k < 4; ++k )
            {
                gLog.message( InfoLevel::Detailed, "    number of %u-chains: %u",
                              ( unsigned int ) k,
                              ( unsigned int ) this->number_of_ksimplices( k ));
            }
            gLog.message( InfoLevel::Detailed, "" );
        }

//-----------------------------------------------------------------------------

        void
        SimplicialComplex::reduce_complexPellikkaGeneralized()
        {
            Timer tTimer;

            gLog.message( InfoLevel::Detailed, "    reducing the complex (generalized pellikka's algorithm) ... " );

            this->reduceOmit();
            for (uint p = 3; p >= 1; p--)
            {
                this->pGeneralizedCombine(p);
                this->pReduce(p-1);
            }

            gLog.message( InfoLevel::Detailed, "    ... reduction time %.3f s\n", tTimer.stop() * 1e-3 );

            for ( uint k = 0; k < 4; ++k )
            {
                gLog.message( InfoLevel::Detailed, "    number of %u-chains: %u",
                              ( unsigned int ) k,
                              ( unsigned int ) this->number_of_ksimplices( k ));
            }
            gLog.message( InfoLevel::Detailed, "" );
        }


//------------------------------------------------------------------------------

        // pCombine from Pellikka et al.
        void
        SimplicialComplex::pCocombine(const uint p)
        {

            index_t a;
            index_t b;
            index_t e;
            Cell< index_t > Q;

            //Loop over the p cochains
            for (auto it = mCochainsMap( p ).begin();  it != mCochainsMap( p ).end(); ++it )
            {
                if (p == 3)
                {
                    return;
                }

                //Store the boundaries
                for (const auto & [tId, tCoeff]: it->second->getCoboundary()->getSimplicesMap())
                {
                    Q.push(tId);
                }

                while (Q.size() > 0)
                {
                    a = Q.pop();
                    Cochain *tCoboundary = mCochainsMap( p+1 )(a);
                    OrderedMap <index_t, int> & tSimplicesMap = tCoboundary->getBoundary()->getSimplicesMap();

                    // If the coboundary is shared exactly by two neighbors
                    if (tSimplicesMap.size()== 2)
                    {
                        auto it2 = tSimplicesMap.begin();
                        auto it3 = ++tSimplicesMap.begin();
                        if (it2->first == it->first) //Always keep the simplex we're looking at
                        {
                            e = it2->first;
                            b = it3->first;
                        }
                        else
                        {
                            b = it2->first;
                            e = it3->first;
                        }

                        Cochain* tCochainAdd = mCochainsMap(p)(e);
                        Cochain* tCochainReduce = mCochainsMap(p)(b);
                        Cochain* tNeighbors = tCochainReduce->getCoboundary();

                        int val1 = tCochainAdd->getCoboundary()->getCoefficient(a);
                        int val2 = tCochainReduce->getCoboundary()->getCoefficient(a);

                        //Add the cochains
                        tCochainAdd->addCochainToCochain(tCochainReduce,-val1*val2);

                        //Update the boundaries
                        for(const auto [tID2, tCoeff2]: tNeighbors->getSimplicesMap())
                        {
                            auto it4 = mCochainsMap( p+1).find( tID2 );

                            if ( it4 != mCochainsMap( p+1 ).end() )
                            {
                                it4->second->getBoundary()->setCoefficient(b,0);
                                if (tID2 != a)
                                {
                                    it4->second->add_simplex_to_boundary(e,-tCoeff2*val2*val1);
                                }
                            }
                        }

                        if (p < 2)
                        {
                            for(const auto [tID2, tCoeff2]: tCoboundary->getCoboundary()->getSimplicesMap())
                            {
                                mCochainsMap(p+2)(tID2)->getBoundary()->setCoefficient(a,0);
                            }
                        }

                        //Remove a and b from the complex
                        this->remove_kcochainFromMap(p + 1, a);
                        this->remove_kcochainFromMap(p, b);

                        //Add the new coboundaries
                        for (const auto & [tId, tCoeff]: tCochainAdd->getCoboundary()->getSimplicesMap())
                        {
                            Q.push(tId);
                        }
                        unique(Q);
                    }
                }

            }
        }

//-----------------------------------------------------------------------------

        void
        SimplicialComplex::pGeneralizedCocombine(const uint p)
        {
            // reused snapshot buffers for the valid neighbors of one (a,b) elimination
            Cell< index_t > tNeighborIDs;
            Cell< Cochain * > tNeighborCochains;

            // Loop over all k-cochains
            for ( auto it = mCochainsMap( p ).begin(), next_it = it; it != mCochainsMap( p ).end(); it = next_it )
            {
                ++next_it;
                Cochain * tCoboundary = it->second->getCoboundary();
                OrderedMap< index_t, int > & tSimplicesMap = tCoboundary->getSimplicesMap();
                index_t b = it->first;
                Cochain* tCochainReduce = it->second;

                // Find the first simplex ID that correspond to a (k+1)-cochain of the complex
                if ( tSimplicesMap.size() != 0 )
                {
                    auto it2 = --tSimplicesMap.end(); //remove the lowest gradient
                    bool tOut = false;
                    while ( !mCochainsMap( p + 1 ).key_exists( it2->first ) or
                            abs( it2->second ) != 1 )
                    {
                        if (it2 == tSimplicesMap.begin())
                        {
                            tOut = true;
                            break;
                        }
                        --it2;
                    }
                    if ( tOut )
                    {
                        continue;
                    }
                    index_t a = it2->first;
                    int val1 = it2->second;

                    // Find the neighbors of b through a
                    Cochain * tCochainA = mCochainsMap( p + 1 )( a );
                    Cochain * tNeighbors = tCochainA->getBoundary();

                    //Collect the valid neighbors of b through a
                    tNeighborIDs.clear();
                    tNeighborCochains.clear();
                    for ( const auto [tID, tCoeff]: tNeighbors->getSimplicesMap())
                    {
                        auto it3 = mCochainsMap( p ).find( tID );
                        if ( it3 != mCochainsMap( p ).end() and tID != b )
                        {
                            tNeighborIDs.push( tID );
                            tNeighborCochains.push( it3->second );
                        }
                    }

                    // cleanup must stay conditional on a valid neighbor existing,
                    // as in the pre-split code
                    if ( tNeighborIDs.size() > 0 )
                    {
                        //Clear b from the boundaries of its coboundary cochains
                        for ( const auto [tID2, tCoeff2]: tSimplicesMap )
                        {
                            auto it4 = mCochainsMap( p + 1 ).find( tID2 );
                            if ( it4 != mCochainsMap( p + 1 ).end()  )
                            {
                                it4->second->getBoundary()->setCoefficient( b,0 );
                            }
                        }

                        //Clear a from the (p+2) boundaries
                        if ( p < 2 )
                        {
                            for ( const auto [tID2, tCoeff2]:
                                    tCochainA->getCoboundary()->getSimplicesMap())
                            {
                                mCochainsMap( p + 2 )( tID2 )
                                        ->getBoundary()->setCoefficient( a,0 );
                            }
                        }

                        //Loop over the Neighbors to update their cochains
                        for ( index_t k = 0; k < tNeighborIDs.size(); ++k )
                        {
                            index_t tID = tNeighborIDs( k );
                            Cochain * tCochainAdd = tNeighborCochains( k );
                            int val2 = tCochainAdd->getCoboundary()->getCoefficient( a );

                            //Add the chain b to the neighbor
                            tCochainAdd->addCochainToCochain( tCochainReduce, -val1 * val2 );

                            //Update the boundaries
                            for ( const auto [tID2, tCoeff2]: tSimplicesMap )
                            {
                                if ( tID2 != a )
                                {
                                    auto it4 = mCochainsMap( p + 1 ).find( tID2 );
                                    if ( it4 != mCochainsMap( p + 1 ).end()  )
                                    {
                                        it4->second
                                                ->add_simplex_to_boundary( tID, -tCoeff2*val2*val1 );
                                    }
                                }
                            }
                        }
                    }
                    this->remove_kcochainFromMap( p + 1, a );
                    this->remove_kcochainFromMap( p, b );
                }

            }
        }



//-----------------------------------------------------------------------------

// ReduceOmit from Pellikka et al.
        void
        SimplicialComplex::coreduceOmit()
        {
            for (uint p = 0; p <= 2; p++)
            {
                this->pCoreduce(p);
            }

            uint tKey;
            while (this->number_of_kcosimplices(0) > 0)
            {
                tKey = mCochainsMap(0).begin()->first;
                for(const auto [tID2, tCoeff2]: mCochainsMap(0).begin()->second->getCoboundary()->getSimplicesMap())
                {
                    mCochainsMap(1)(tID2)->getBoundary()->setCoefficient(tKey,0);
                }

                this->remove_kcochainFromMap(0,tKey);
                for (uint p = 0; p <= 2; p++)
                {
                    this->pCoreduce(p);
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        SimplicialComplex::coreduce_complexPellikka()
        {
            Timer tTimer;

            gLog.message( InfoLevel::Detailed, "    coreducing the complex (pellikka's algorithm) ... " );

            this->coreduceOmit();
            for (uint p = 0; p <= 2; p++)
            {
                this->pCocombine(p);
                this->pCoreduce(p+1);
            }

            gLog.message( InfoLevel::Detailed, "    ... coreduction time %.3f s\n", tTimer.stop() * 1e-3 );

            for ( uint k = 0; k < 4; ++k )
            {
                gLog.message( InfoLevel::Detailed, "    number of %u-cochains: %u",
                              ( unsigned int ) k,
                              ( unsigned int ) this->number_of_kcosimplices( k ));
            }
            gLog.message( InfoLevel::Detailed, "" );

        }

//-----------------------------------------------------------------------------

        void
        SimplicialComplex::coreduce_complexPellikkaGeneralized()
        {
            Timer tTimer;

            gLog.message( InfoLevel::Detailed, "    coreducing the complex (generalized pellikka's algorithm) ... " );

            this->coreduceOmit();
            for (uint p = 0; p <= 2; p++)
            {
                this->pGeneralizedCocombine(p);
                this->pCoreduce(p+1);
            }

            gLog.message( InfoLevel::Detailed, "    ... coreduction time %.3f s\n", tTimer.stop() * 1e-3 );

            for ( uint k = 0; k < 4; ++k )
            {
                gLog.message( InfoLevel::Detailed, "    number of %u-cochains: %u",
                              ( unsigned int ) k,
                              ( unsigned int ) this->number_of_kcosimplices( k ));
            }
            gLog.message( InfoLevel::Detailed, "" );

        }

//-----------------------------------------------------------------------------

        void
        SimplicialComplex::remove_kchainFromMap(const uint k, const index_t aID) {
            if (k <= 3) {
                delete mChainsMap(k)[aID];
                mChainsMap(k).erase_key(aID);
            }
        }

//------------------------------------------------------------------------------

        Chain *
        SimplicialComplex::get_kchain(const uint k, const index_t aID) {
            if (k <= 3) {
                return mChainsMap(k)[aID];
            } else {
                return nullptr;
            }
        }

//------------------------------------------------------------------------------

        Cochain *
        SimplicialComplex::get_kcochain(const uint k, const index_t aID) {
            if (k <= 3) {
                return mCochainsMap(k)[aID];
            } else {
                return nullptr;
            }
        }

//------------------------------------------------------------------------------

        Map<index_t, Chain *>
        SimplicialComplex::get_kchainMap(const uint k) {
            return mChainsMap(k);
        }

//------------------------------------------------------------------------------

        Map<index_t, Cochain *>
        SimplicialComplex::get_kcochainMap(const uint k) {
            return mCochainsMap(k);
        }

//------------------------------------------------------------------------------

        /*Chain *
        mesh::SimplicialComplex::boundary_of_kchain( const uint k, const index_t aID )
        {
            if ( k <= 3 )
            {
                return mChainsMap( k )[ aID ]->getBoundary();
            }
            else
            {
                return nullptr;
            }
        }*/

//------------------------------------------------------------------------------

// Function that creates the boundary matrices (with 1, -1 and 0)
// from the simplicial complex (reduced or not)
        Cell<Matrix<int> >
        SimplicialComplex::createMatrixFromBoundaryMap() {
            Cell<Matrix<int> > tBoundaryMat;
            tBoundaryMat.set_size(4, Matrix<int>());
            tBoundaryMat(0).set_size(1, this->number_of_0simplices(), 0);

    //loop over all dimensions
            for (int k = 1; k < 4; k++) {
                tBoundaryMat(k).set_size(this->number_of_ksimplices(k - 1), this->number_of_ksimplices(k), 0);
                uint tCount = 0;

        //loop over all the k-chains and (k-1)-chains to populate the matrix
                for (const auto &[tKey, tChain]: mChainsMap(k)) {
                    uint tCount2 = 0;
                    for (const auto &[tKey2, tChain2]: mChainsMap(k - 1)) {
                        tBoundaryMat(k)(tCount2, tCount) = tChain->getBoundary()->getCoefficient(tKey2);
                        tCount2++;
                    }
                    tCount++;
                }
            }
            return tBoundaryMat;
        }

//------------------------------------------------------------------------------

// Function that creates the coboundary matrices (with 1, -1 and 0)
// from the simplicial complex (coreduced or not)
        Cell<Matrix<int> >
        SimplicialComplex::createMatrixFromCoboundaryMap() {
            Cell<Matrix<int> > tCoboundaryMat;
            tCoboundaryMat.set_size(4, Matrix<int>());

            int n = this->number_of_kcosimplices(3);
            tCoboundaryMat(3).set_size(1, std::max(n, 1), 0);

    //loop over all dimensions
            for (int k = 0; k < 3; k++) {
                int m = this->number_of_kcosimplices(k + 1);
                n = this->number_of_kcosimplices(k);

                tCoboundaryMat(k).set_size(std::max(m, 1), std::max(n, 1), 0);
                uint tCount = 0;

        //loop over all the k-cochains and (k+1)-cochains to populate the matrix
                for (const auto &[tKey, tChain]: mCochainsMap(k)) {
                    uint tCount2 = 0;
                    for (const auto &[tKey2, tChain2]: mCochainsMap(k + 1)) {
                        tCoboundaryMat(k)(tCount2, tCount) = tChain->getCoboundary()->getCoefficient(tKey2);
                        tCount2++;
                    }
                    tCount++;
                }
            }
            return tCoboundaryMat;
        }

//------------------------------------------------------------------------------

        void
        SimplicialComplex::print_kchains(const uint k) {
            for (const auto &[tKey, tChain]: mChainsMap(k)) {
                std::cout << k << "-chain #" << tKey << " :" << std::endl;
                tChain->print();
            }
        }

//------------------------------------------------------------------------------

        void
        SimplicialComplex::create_kComplexField( const uint k, Mesh * aMesh, string aFieldName )
        {
            aMesh->create_field( aFieldName, EntityType::ELEMENT);
            Cell < Element* > tElements = aMesh->elements();
            for (const auto &[tKey, tChain]: this->get_kchainMap(k))
            {
                OrderedMap<index_t, int> & tSimplicesMap = tChain->getSimplicesMap();
                for( const auto& [tID, tCoeff] :  tSimplicesMap)
                {
                    aMesh->field_data(aFieldName)(tElements(tID)->index()) = tCoeff;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        SimplicialComplex::cocreate_kComplexField( const uint k, Mesh * aMesh, string aFieldName )
        {
            aMesh->create_field( aFieldName, EntityType::ELEMENT);
            Cell < Element* > tElements = aMesh->elements();
            for (const auto &[tKey, tCochain]: this->get_kcochainMap(k))
            {
                OrderedMap<index_t, int> & tSimplicesMap = tCochain->getSimplicesMap();
                for( const auto& [tID, tCoeff] :  tSimplicesMap)
                {
                    aMesh->field_data(aFieldName)(tElements(tID)->index()) = tCoeff;
                }
            }
        }

//------------------------------------------------------------------------------

    }
}
