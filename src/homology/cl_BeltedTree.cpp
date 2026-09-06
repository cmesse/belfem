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

#include "cl_BeltedTree.hpp"

#include "cl_Timer.hpp"
#include "cl_Logger.hpp"

namespace belfem
{
//------------------------------------------------------------------------------
    namespace mesh
    {
//------------------------------------------------------------------------------

        BeltedTree::BeltedTree(Mesh* aMesh, SimplicialComplex * aSimplicialComplex, Cell< Chain * > a1HomologyGenerators ) :
                m1HomologyGenerators(a1HomologyGenerators),
                mSimplicialComplex(aSimplicialComplex),
                mMesh(aMesh)
        {
            mBeltFasteners.set_size(m1HomologyGenerators.size(), nullptr);
            m1CohomologyGenerators.set_size(m1HomologyGenerators.size(), nullptr);
            for (uint i = 0; i < m1CohomologyGenerators.size(); i++)
            {
                m1CohomologyGenerators(i) = new Cochain(1,mMesh, false, false) ;
            }

            this->select_belt_fasteners();
            this->create_tree();
        }

//------------------------------------------------------------------------------

        BeltedTree::~BeltedTree()
        {
            for (uint i = 0; i < m1CohomologyGenerators.size(); i++)
            {
                delete m1CohomologyGenerators(i) ;
            }
            m1CohomologyGenerators.clear();
        }

//------------------------------------------------------------------------------

        void
        BeltedTree::select_belt_fasteners()
        {
            Cell < Edge* > & tEdges = mMesh->edges();
            for (uint i = 0 ; i < mBeltFasteners.size() ; i++)
            {
                auto it = m1HomologyGenerators(i)->getSimplicesMap().begin() ;
                for (uint j = 0; j < i; j++)
                {
                    while (tEdges(it->first)->index() == mBeltFasteners(j)->index())
                    {
                        it++ ;
                    }
                }
                mBeltFasteners(i) = tEdges(it->first);
            }
        }

//------------------------------------------------------------------------------

        void
        BeltedTree::create_tree()
        {

            Timer tTimer;

            gLog.message( InfoLevel::Detailed, "    creating the belted tree ... " );
            Cell < index_t > tTreeV ;
            Cell < Cell < index_t >> tHomologyV;
            Cell < Edge * > & tEdges = mMesh->edges();
            //Initializing the tree with the homology edges, except the belt fasteners

            for (uint i = 0 ; i < mBeltFasteners.size() ; i++)
            {
                tHomologyV.push(Cell<index_t>());
                for ( const auto& [tID, tCoeff] : m1HomologyGenerators(i)->getSimplicesMap())
                {
                    if (tID != mBeltFasteners(i)->index())
                    {
                        tHomologyV(i).push(tEdges(tID)->node(0)->index());
                        tHomologyV(i).push(tEdges(tID)->node(1)->index());

                    }
                }
                unique(tHomologyV(i));
            }

            for ( const auto& [tID, tCoeff] : m1HomologyGenerators(0)->getSimplicesMap())
            {
                if (tID != mBeltFasteners(0)->index())
                {
                    mTree.push(tID) ;
                    tTreeV.push(tEdges(tID)->node(0)->index()) ;
                    tTreeV.push(tEdges(tID)->node(1)->index()) ;
                }
            }
            unique(tTreeV);

            Cell < bool > tvisited = Cell < bool >(m1HomologyGenerators.size(),false);
            tvisited(0) = true;
            Map< index_t, Cochain * > t0CochainMap = mSimplicialComplex->get_kcochainMap(0);
            uint it = 0;
            while (it < tTreeV.size())
            {
                for (const auto& [tE,tCoeff] : t0CochainMap(tTreeV(it))->getCoboundary()->getSimplicesMap())
                {
                    if ( std::find( mTree.begin(), mTree.end(), tE) != mTree.end())
                    {
                        continue ;
                    }
                    else if ((std::find( tTreeV.begin(), tTreeV.end(), tEdges(tE)->node(0)->index()) == tTreeV.end()) !=
                             (std::find( tTreeV.begin(), tTreeV.end(), tEdges(tE)->node(1)->index()) == tTreeV.end()))
                    {
                        mTree.push(tE);
                        index_t tAddV = tEdges(tE)->node(0)->index() == tTreeV(it)?tEdges(tE)->node(1)->index():tEdges(tE)->node(0)->index();
                        tTreeV.push(tAddV);

                        for (uint i = 0 ; i < mBeltFasteners.size() ; ++i)
                        {
                            if (!tvisited(i))
                            {
                                if (std::find( tHomologyV(i).begin(), tHomologyV(i).end(), tAddV) != tHomologyV(i).end())
                                {
                                    for ( const auto& [tID, tCoeff2] : m1HomologyGenerators(i)->getSimplicesMap())
                                    {
                                        if (tID != mBeltFasteners(i)->index())
                                        {
                                            mTree.push(tID) ;
                                            tTreeV.push(tEdges(tID)->node(0)->index()) ;
                                            tTreeV.push(tEdges(tID)->node(1)->index()) ;
                                        }
                                        //unique(tTreeV);
                                    }
                                    tvisited(i) = true;
                                }
                            }
                        }
                    }
                }
                ++it;
            }

            // Add the belt fasteners
            for (uint i = 0 ; i < mBeltFasteners.size() ; i++)
            {
                mTree.push(mBeltFasteners(i)->index()) ;
            }
            //unique(mTree);

            gLog.message( InfoLevel::Detailed, "    ... belted tree created in %.3f s\n", tTimer.stop() * 1e-3 );

        }

        void
        BeltedTree::compute_cohomology()
        {

            Timer tTimer;

            gLog.message( InfoLevel::Detailed, "    compute cohomology from belted tree ... " );
            //Init
            Map<index_t, Basis * > tL ;
            int tCoeff;
            uint tNum;
            Cell < int > tOrient;
            tOrient.set_size(3,0);
            Map<index_t,bool> tSimplexExist;

            //Put all triangles in a container
            //tL.set_size(mSimplicialComplex->number_of_2simplices(), nullptr) ;

            //Loop over all belt fasteners (number of cohomology cuts)
            Cell < Element * > & tElements = mMesh->elements();
            Cell < Face * > & tFaces = mMesh->faces();
            for (uint i = 0; i < mBeltFasteners.size(); i++)
            {

                tL.clear();
                for ( const auto& [tI, tChain] : mSimplicialComplex->get_kcochainMap(2) )
                {
                    if (mMesh->number_of_dimensions() == 2)
                    {
                        tL[tI] = tElements(tI);
                    }
                    else
                    {
                        tL[tI] = tFaces(tI);
                    }
                }

                //Init the existence map to 0 for
                for ( const auto& [tI, tChain] : mSimplicialComplex->get_kcochainMap(1) )
                {
                    tSimplexExist[tI] = false;
                }

                //Impose 1 on belt fasteners and 0 on other edges of the tree
                for (index_t tE : mTree)
                {
                    if (tE != mBeltFasteners(i)->index())
                    {
                        tSimplexExist(tE) = true ;
                    }
                }

                auto it = tL.begin();

                // First, remove all triangles with 2 edges imposed
                uint tCount = 1;
                while (tCount!= 0)
                {
                    tCount = 0;
                    it = tL.begin();

                    while (it != tL.end())
                    {
                        tNum = tSimplexExist(it->second->edge(0)->index()) +
                               tSimplexExist(it->second->edge(1)->index()) +
                               tSimplexExist(it->second->edge(2)->index());


                        // Enforce the 0-circulation by computing the third coefficient on triangles with 2 edges imposed
                        if (tNum == 2 )
                        {
                            index_t tID = it->first;
                            for (uint j = 0; j < 3; ++j)
                            {
                                if (tSimplexExist(it->second->edge(j)->index()) == 0)
                                {
                                    tSimplexExist(it->second->edge(j)->index()) = true;
                                }
                            }
                            tL.erase_key(tID);
                            tCount++;
                        }
                        else
                        {
                            ++it;
                        }
                    }
                }

                // Add the belt fastener and restart with the triangles left
                m1CohomologyGenerators(i)->addSimplexToCochain(mBeltFasteners(i)->index(),m1HomologyGenerators(i)->getSimplicesMap()[mBeltFasteners(i)->index()]) ;
                tSimplexExist(mBeltFasteners(i)->index()) = true ;

                tCount = 1;
                while (tCount!= 0)
                {
                    tCount = 0;
                    it = tL.begin();
                    while (it!= tL.end())
                    {
                        tNum = tSimplexExist(it->second->edge(0)->index()) +
                               tSimplexExist(it->second->edge(1)->index()) +
                               tSimplexExist(it->second->edge(2)->index());


                        // Enforce the 0-circulation by computing the third coefficient on triangles with 2 edges imposed
                        if (tNum == 2 )
                        {
                            index_t tID = it->first;
                            for (uint j = 0; j < 3; j++)
                            {
                                tOrient(j) = mSimplicialComplex->get_kcochainMap(2)[it->second->index()]->getBoundary()->getCoefficient(it->second->edge(j)->index()) ;
                            }

                            //Compute the coefficient on the new edge
                            for (uint j = 0; j < 3; ++j)
                            {
                                if (tSimplexExist(it->second->edge(j)->index())==0)
                                {

                                    tCoeff = -1*tOrient(j)*(m1CohomologyGenerators(i)->getCoefficient(it->second->edge((j+1)%3)->index())*tOrient((j+1)%3)
                                                            +m1CohomologyGenerators(i)->getCoefficient(it->second->edge((j+2)%3)->index())*tOrient((j+2)%3)) ;
                                    m1CohomologyGenerators(i)->addSimplexToCochain(it->second->edge(j)->index(), tCoeff);

                                    tSimplexExist(it->second->edge(j)->index()) = true ;
                                    break ;
                                }
                            }
                            tL.erase_key(tID);
                            tCount++;
                        }
                        else
                        {
                            ++it;
                        }
                    }
                }
            }

            gLog.message( InfoLevel::Detailed, "    ... cohomology computed in %.3f s\n", tTimer.stop() * 1e-3 );
        }

//------------------------------------------------------------------------------

        void
        BeltedTree::create_TreeField(Mesh* aEdgeMesh , string tFieldName)
        {
            Cell < Element * > tElements = aEdgeMesh->elements();
            aEdgeMesh->create_field( tFieldName, EntityType::ELEMENT);
            for( index_t tE :  mTree)
            {
                aEdgeMesh->field_data(tFieldName)(tElements(tE)->index()) = 1;
                //tMesh->field_data(fieldName)(tMesh->edge(tID)->node(1)->index()) = 1;
            }
        }

//-----------------------------------------------------------------------------

        void
        BeltedTree::create_cohomologyField( Mesh* tEdgeMesh )
        {
            /*uint tCount = 0;
            for(const auto& tCochain : m1CohomologyGenerators)
            {
                tCount++;
                char fieldName [50];
                sprintf (fieldName, "1CohomologyGenerator(BT)%d", tCount);
                mMesh->create_field( fieldName, EntityType::NODE);
                Map< index_t, int > tSimplicesMap = tCochain->getSimplicesMap();
                for( const auto& [tID, tCoeff] :  tSimplicesMap)
                {
                    mMesh->field_data(fieldName)(mMesh->edge(tID)->node(0)->index()) = 1;
                    mMesh->field_data(fieldName)(mMesh->edge(tID)->node(1)->index()) = 1;
                }
            }*/
            uint tCount = 0;
            for(const auto& tCochain : m1CohomologyGenerators)
            {
                tCount++;
                char fieldName [50];
                /*if (mflagProp)
                {
                    sprintf (fieldName, "1CohomologyGenerator%d", tCount);
                }
                else
                {
                    sprintf (fieldName, "1CohomologyGeneratorComp%d", tCount);
                }*/
                sprintf (fieldName, "1CohomologyGenerator(BT)%d", tCount);
                tEdgeMesh->create_field( fieldName, EntityType::ELEMENT);
                OrderedMap< index_t, int > & tSimplicesMap = tCochain->getSimplicesMap();
                for( const auto& [tID, tCoeff] :  tSimplicesMap)
                {
                    tEdgeMesh->field_data(fieldName)(tEdgeMesh->element(tID)->index()) = tCoeff;
                    //tMesh->field_data(fieldName)(tMesh->edge(tID)->node(1)->index()) = 1;
                }
            }
        }

//-----------------------------------------------------------------------------

        Cell< Cochain * > &
        BeltedTree::get_cohomology()
        {
            return m1CohomologyGenerators;
        }

//-----------------------------------------------------------------------------
    }
}
