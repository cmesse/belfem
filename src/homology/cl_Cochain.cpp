//
// Created by grgia@ge.polymtl.ca on 06/12/23.
//
#include <iostream>
#include "cl_Cochain.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        Cochain::Cochain( const uint aDim, Mesh * aMesh ) :
                mDim( aDim ),
                mMesh( aMesh )
        {
            if (aDim != 3)
            {
                mCoboundary = new Cochain(aDim+1, aMesh);
            }
            else
            {
                mCoboundary = nullptr;
            }
        }

        //------------------------------------------------------------------------------

        Cochain::~Cochain()
        {
            delete mCoboundary;
        }

        //------------------------------------------------------------------------------

        void
        Cochain::operator+( Cochain * aCochain )
        {
            this->addCochainToCochain( aCochain, 1 );
        }

        //-----------------------------------------------------------------------------

        void
        Cochain::operator-( Cochain * aCochain )
        {
            this->removeCochainFromCochain( aCochain );
        }

        //-----------------------------------------------------------------------------

        int
        Cochain::operator()( Chain * aChain )
        {
            int tVal = 0;
            Map< id_t, int > tChainMap = aChain->getSimplicesMap();
            for( const auto& [tID, tCoeff] : mSimplicesMap )
            {
                if (tChainMap.key_exists(tID))
                {
                    tVal += tCoeff*tChainMap(tID);
                }
            }
            return tVal ;
        }

        //-----------------------------------------------------------------------------

        void
        Cochain::addSimplexToCochain( const id_t aID, const int aCoeff )
        {

            // Add to the existing simplex, or create a new one
            if (mSimplicesMap.key_exists(aID))
            {
                mSimplicesMap[aID] += aCoeff;
            }
            else
            {
                mSimplicesMap[aID] = aCoeff;
            }

            // Update the coboundary of the cochain from what's just been added
            // 0-cochain
            if (mDim == 0)
            {
                Node * tNode = mMesh->node(aID);
                for (uint i = 0; i < tNode->number_of_edges(); i++)
                {
                    if (tNode->edge(i)->node(0)->id() == aID)
                    {
                        mCoboundary->addSimplexToCochain(tNode->edge(i)->id(),-1*aCoeff);
                    }
                    else
                    {
                        mCoboundary->addSimplexToCochain(tNode->edge(i)->id(),1*aCoeff);
                    }
                }
            }

            // 1-cochain
            else if (mDim == 1)
            {
                Edge * tEdge = mMesh->edge(aID);
                if (mMesh->number_of_dimensions() == 2)
                {
                    for (uint i = 0; i < tEdge->number_of_elements(); i++)
                    {
                        int tMul = 0 ;
                        Element* tElement = tEdge->element(i);
                        // Probably can be optimized (how to know the orientation of a given edge in a triangle)
                        for(uint j = 0; j < tElement->number_of_edges(); j++)
                        {
                            if (tElement->edge(j)->id() == tEdge->id())
                            {
                                tMul = (tElement->edge_orientation(j)?1:-1);
                                break;
                            }
                        }
                        mCoboundary->addSimplexToCochain(tElement->id(),tMul*aCoeff);
                    }
                }
                else
                {
                    for (uint i = 0; i < tEdge->number_of_faces(); i++)
                    {
                        int tMul = 0 ;
                        Face* tFace = tEdge->face(i);
                        // Probably can be optimized (how to know the orientation of a given edge in a triangle)
                        for(uint j = 0; j < tFace->number_of_edges(); j++)
                        {
                            if (tFace->edge(j)->id() == tEdge->id())
                            {
                                tMul = (tFace->edge_orientation(j)?1:-1);
                                break;
                            }
                        }
                        mCoboundary->addSimplexToCochain(tFace->id(),tMul*aCoeff);
                    }
                }
            }

            else if (mMesh->number_of_dimensions() == 3 && mDim == 2)
            {
                Face * tFace = mMesh->face(aID);
                mCoboundary->addSimplexToCochain(tFace->master()->id(),1*aCoeff);
                if (tFace->slave())
                {
                    mCoboundary->addSimplexToCochain(tFace->slave()->id(),-1*aCoeff);
                }
            }

            //remove the simplex from the map if the coefficient becomes 0
            if (mSimplicesMap[aID] == 0)
            {
                mSimplicesMap.erase_key(aID);
            }
        }

        //------------------------------------------------------------------------------

        void
        Cochain::removeSimplexToCochain(const id_t aID)
        {
            //Removing is simply adding a -1 coefficient
            this->addSimplexToCochain(aID, -1);
        }

        //-----------------------------------------------------------------------------

        void
        Cochain::addCochainToCochain(Cochain * aCochain, int aCoeff)
        {
            if (mDim == aCochain->getDim())
            {
                // Add all the simplices from the input cochain
                for( const auto& [tID, tCoeff] : aCochain->getSimplicesMap() )
                {
                    this->addSimplexToCochain(tID,tCoeff*aCoeff);
                }
            }
            else
            {
                std::cout << "Dimensions don't agree, not adding the cochain" << std::endl;
            }
        }

        //-----------------------------------------------------------------------------

        void
        Cochain::removeCochainFromCochain(Cochain * aCochain)
        {
            if (mDim == aCochain->getDim())
            {
                for( const auto& [tID, tCoeff] : aCochain->getSimplicesMap() )
                {
                    this->addSimplexToCochain(tID,-1*tCoeff);
                }
            }
            else
            {
                std::cout << "Dimensions don't agree, not adding the cochain" << std::endl;
            }
        }

        //-----------------------------------------------------------------------------

        Map< id_t, int >
        Cochain::getSimplicesMap()
        {
            return mSimplicesMap;
        }

        //------------------------------------------------------------------------------

        int
        Cochain::getCoefficient( const id_t aID )
        {
            auto it = mSimplicesMap.find(aID);
            if (it != mSimplicesMap.end()) {
                return it->second;
            }
            else {
                return 0;
            }
        }

        //------------------------------------------------------------------------------

        Cochain*
        Cochain::getCoboundary()
        {
            return mCoboundary;
        }

        //-----------------------------------------------------------------------------

        int
        Cochain::getDim()
        {
            return mDim;
        }

        //-----------------------------------------------------------------------------

        void
        Cochain::setCoefficient(const id_t aID, const int tCoeff)
        {
            if (tCoeff == 0)
            {
                mSimplicesMap.erase_key(aID);
            }
            else
            {
                mSimplicesMap[aID] = tCoeff;
            }
        }

        //-----------------------------------------------------------------------------

        void
        Cochain::print()
        {
            std::cout << mDim << "-cochain:" << std::endl;
            for(const auto& [tID, tCoeff] : mSimplicesMap)
            {
                if (tCoeff != 0)
                {
                    std::cout << "Simplex ID " << tID << " Coeff: " << tCoeff ;
                    if (mDim == 1)
                    {
                        std::cout << " Nodes: " << mMesh->edge(tID)->node(0)->id() << " to "  << mMesh->edge(tID)->node(1)->id() << std::endl;
                    }
                    else
                    {
                        std::cout << std::endl;
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------
    }
}
