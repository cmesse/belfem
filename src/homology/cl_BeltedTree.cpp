//
// Created by grgia@ge.polymtl.ca on 31/01/24.
//

#include "cl_BeltedTree.hpp"

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
                m1CohomologyGenerators(i) = new Cochain(1,mMesh) ;
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
            for (uint i = 0 ; i < mBeltFasteners.size() ; i++)
            {
                auto it = m1HomologyGenerators(i)->getSimplicesMap().begin() ;
                for (uint j = 0; j < i; j++)
                {
                    while (mMesh->edge(it->first)->id() == mBeltFasteners(j)->id())
                    {
                        it++ ;
                    }
                }
                mBeltFasteners(i) = mMesh->edge(it->first);
            }
        }

        //------------------------------------------------------------------------------

        //todo:: Clean up and optimize this... creating the belted tree is extremely slow (probably because of the depth-first search algorithm)

        void
        BeltedTree::create_tree()
        {
            Cell < id_t > tTreeV ;
            //Initializing the tree with the homology edges, except the belt fasteners
            for (uint i = 0 ; i < mBeltFasteners.size() ; i++)
            {
                for ( const auto& [tID, tCoeff] : m1HomologyGenerators(i)->getSimplicesMap() )
                {
                    if (tID != mBeltFasteners(i)->id())
                    {
                        mTree.push(tID) ;
                        tTreeV.push(mMesh->edge(tID)->node(0)->id()) ;
                        tTreeV.push(mMesh->edge(tID)->node(1)->id()) ;
                    }
                }
            }
            this->create_TreeField("BeltedTree0") ;

            //Build the maximal possible tree
            Map< id_t, Chain * > t1ChainMap = mSimplicialComplex->get_kchainMap(1);
            Map <id_t, bool > tvisited ;
            Map < id_t, bool > tfinished ;

            //Loop over all edges
            for (const auto& [tID, tChain] : t1ChainMap)
            {
                Edge * tE = mMesh->edge(tID) ;

                //Skip if the edge is already in the tree
                if ( std::find( mTree.begin(), mTree.end(), tE->id()) != mTree.end())
                {
                    continue ;
                }

                //Add directly if only one node is shared (it will not close a loop of the tree)
                else if ( std::find( tTreeV.begin(), tTreeV.end(), tE->node(0)->id()) == tTreeV.end() || std::find( tTreeV.begin(), tTreeV.end(), tE->node(1)->id()) == tTreeV.end())
                {
                    mTree.push(tE->id());
                    tTreeV.push(tE->node(0)->id());
                    tTreeV.push(tE->node(1)->id());
                }

                // Otherwise, check if adding the edge closes a loop (using dfs algorithm)
                else
                {

                    for (id_t tNodeID : tTreeV)
                    {
                        tvisited[tNodeID] = false ;
                        tfinished[tNodeID] = false ;
                    }

                    mCycleFound = false;
                    id_t tNodeID = tE->node(0)->id();
                    mNodeToFind = tE->node(1)->id();
                    tvisited[tNodeID] = false ;
                    tvisited[mNodeToFind] = false ;
                    tfinished[tNodeID] = false ;
                    tfinished[mNodeToFind] = false ;
                    this->dfs(tNodeID, tE->id(), tfinished, tvisited) ;

                    //If no cycles are found, add the edge to the tree
                    if (!mCycleFound)
                    {
                        mTree.push(tE->id());
                        tTreeV.push(tE->node(0)->id());
                        tTreeV.push(tE->node(1)->id());
                    }

                }

            }

            // Add the belt fasteners
            for (uint i = 0 ; i < mBeltFasteners.size() ; i++)
            {
                mTree.push(mBeltFasteners(i)->id()) ;
            }
            //unique(mTree);

        }

        //------------------------------------------------------------------------------

        // dfs algorithm to check if we add a cycle to the graph
        // (as described in https://en.wikipedia.org/wiki/Cycle_(graph_theory))
        void
        BeltedTree::dfs(uint v, uint e, Map <id_t, bool > & tfinished, Map <id_t, bool > & tvisited)
        {
            //No check needed if we already found that it would create a cycle
            if (tfinished(v) || mCycleFound)
            {
                return ;
            }

            //Check if the current vertex is the node that creates a cycle
            if (v == mNodeToFind)
            {
                mCycleFound = true;
                return ;
            }

            //Update the visited edge and redo dfs on each neighbor node
            tvisited(v) = true;
            for (uint j = 0; j < mMesh->node(v)->number_of_edges(); j++)
            {
                for (uint k = 0; k < mTree.size(); k++)
                {
                    if (mTree(k) == mMesh->node(v)->edge(j)->id() && mTree(k)!= e)
                    {
                        if (mMesh->edge(mTree(k))->node(0)->id() == v)
                        {
                            this->dfs(mMesh->edge(mTree(k))->node(1)->id(), mTree(k),tfinished,tvisited);
                        }
                        else
                        {
                            this->dfs(mMesh->edge(mTree(k))->node(0)->id(),mTree(k),tfinished,tvisited);
                        }
                    }
                }
            }
            tfinished(v) = true ;
            //return false;
        }

        //------------------------------------------------------------------------------

        void
        BeltedTree::compute_cohomology()
        {

            //Init
            Cell< Element * > tL ;
            uint tInd;
            uint tCount;
            int tCoeff;
            uint tNumImp;
            uint tNum;
            Cell < int > tOrient;
            tOrient.set_size(3,0);
            Map<id_t,bool> tSimplexExist;

            //Put all triangles in a container
            tL.set_size(mSimplicialComplex->number_of_2simplices(), nullptr) ;
            tCount = 0;
            for ( const auto& [tI, tChain] : mSimplicialComplex->get_kchainMap(2) )
            {
                tL(tCount) = mMesh->element(tI);
                tCount++;
            }

            //Loop over all belt fasteners (number of cohomology cuts)
            for (uint i = 0; i < mBeltFasteners.size(); i++)
            {

                //Init the existence map to 0 for
                for ( const auto& [tI, tChain] : mSimplicialComplex->get_kchainMap(1) )
                {
                    tSimplexExist[tI] = false;
                }

                //Impose 1 on belt fasteners and 0 on other edges of the tree
                for (id_t tE : mTree)
                {
                    if (tE == mBeltFasteners(i)->id())
                    {
                        m1CohomologyGenerators(i)->addSimplexToCochain(tE,m1HomologyGenerators(i)->getSimplicesMap()[tE]) ;
                        tSimplexExist(tE) = true ;
                    }
                    else
                    {
                        tSimplexExist(tE) = true ;
                    }
                }

                tInd = 0 ;
                tNumImp = mTree.size();

                // Check all the triangles, until all the edges are imposed
                while (tNumImp < mSimplicialComplex->number_of_1simplices())
                {
                    tNum = tSimplexExist(tL(tInd)->edge(0)->id()) +
                            tSimplexExist(tL(tInd)->edge(1)->id()) +
                            tSimplexExist(tL(tInd)->edge(2)->id());


                    // Enforce the 0-circulation by computing the third coefficient on triangles with 2 edges imposed
                    if (tNum == 2 )
                    {
                        //Find the orientation of each of the triangle's edges
                        for (uint j = 0; j < 3; j++)
                        {
                            tOrient(j) = mSimplicialComplex->get_kchainMap(2)[tL(tInd)->id()]->getBoundary()->getCoefficient(tL(tInd)->edge(j)->id()) ;
                        }

                        //Compute the coefficient on the new edge
                        for (uint j = 0; j < 3; ++j)
                        {
                            if (tSimplexExist(tL(tInd)->edge(j)->id())==0)
                            {

                                tCoeff = -1*tOrient(j)*(m1CohomologyGenerators(i)->getCoefficient(tL(tInd)->edge((j+1)%3)->id())*tOrient((j+1)%3)
                                                        +m1CohomologyGenerators(i)->getCoefficient(tL(tInd)->edge((j+2)%3)->id())*tOrient((j+2)%3)) ;
                                m1CohomologyGenerators(i)->addSimplexToCochain(tL(tInd)->edge(j)->id(), tCoeff);

                                tSimplexExist(tL(tInd)->edge(j)->id()) = true ;
                                tNumImp+=1 ;
                                break ;
                            }
                        }
                    }
                    else
                    {
                        tInd += 1;
                    }

                    if (tInd >= mSimplicialComplex->number_of_2simplices())
                    {
                        //restart from zero until all the edges are imposed
                        tInd = 0;
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        BeltedTree::create_TreeField(string tFieldName)
        {
            mMesh->create_field( tFieldName, EntityType::NODE);
            for( id_t tE :  mTree)
            {
                mMesh->field_data(tFieldName)(mMesh->edge(tE)->node(0)->index()) += 1;
                mMesh->field_data(tFieldName)(mMesh->edge(tE)->node(1)->index()) += 1;
            }
        }

        //-----------------------------------------------------------------------------

        void
        BeltedTree::create_cohomologyField( )
        {
            uint tCount = 0;
            for(const auto& tCochain : m1CohomologyGenerators)
            {
                tCount++;
                char fieldName [50];
                sprintf (fieldName, "1CohomologyGenerator(BT)%d", tCount);
                mMesh->create_field( fieldName, EntityType::NODE);
                Map< id_t, int > tSimplicesMap = tCochain->getSimplicesMap();
                for( const auto& [tInd, tCoeff] :  tSimplicesMap)
                {
                    mMesh->field_data(fieldName)(mMesh->edge(tInd)->node(0)->index()) = 1;
                    mMesh->field_data(fieldName)(mMesh->edge(tInd)->node(1)->index()) = 1;
                }
            }
        }

        //-----------------------------------------------------------------------------
    }
}