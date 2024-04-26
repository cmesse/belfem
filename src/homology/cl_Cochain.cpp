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
            for( const auto& [tInd, tCoeff] : mSimplicesMap )
            {
                if (tChainMap.key_exists(tInd))
                {
                    tVal += tCoeff*tChainMap(tInd);
                }
            }
            return tVal ;
        }

        //-----------------------------------------------------------------------------

        void
        Cochain::addSimplexToCochain( const id_t aInd, const int aCoeff )
        {

            // Add to the existing simplex, or create a new one
            if (mSimplicesMap.key_exists(aInd))
            {
                mSimplicesMap[aInd] += aCoeff;
            }
            else
            {
                mSimplicesMap[aInd] = aCoeff;
            }

            // Update the coboundary of the cochain from what's just been added
            // 0-cochain
            if (mDim == 0)
            {
                Node * tNode = mMesh->node(aInd);
                for (uint i = 0; i < tNode->number_of_edges(); i++)
                {
                    if (tNode->node(0)->id() == aInd)
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
                Edge * tEdge = mMesh->edge(aInd);
                if (mMesh->number_of_dimensions() == 2)
                {
                    for (uint i = 0; i < tEdge->number_of_elements(); i++)
                    {
                        int tMul = 0 ;
                        Element* tElement = tEdge->element(i);
                        // Probably can be optimized (how to know the orientation of a given edge in a triangle)
                        if (tElement->edge(0)->id() == tEdge->id())
                        {
                            tMul = (tElement->edge_orientation(0)?1:-1);
                        }
                        else
                        {
                            tMul = (tElement->edge_orientation(1)?1:-1);
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
                        if (tFace->edge(0)->id() == tEdge->id())
                        {
                            tMul = (tFace->edge_orientation(0)?1:-1);
                        }
                        else
                        {
                            tMul = (tFace->edge_orientation(1)?1:-1);
                        }
                        mCoboundary->addSimplexToCochain(tFace->id(),tMul*aCoeff);
                    }
                }
            }

            else if (mMesh->number_of_dimensions() == 3 && mDim == 2)
            {
                Face * tFace = mMesh->face(aInd);
                mCoboundary->addSimplexToCochain(tFace->master()->id(),1*aCoeff);
                if (tFace->slave())
                {
                    mCoboundary->addSimplexToCochain(tFace->slave()->id(),-1*aCoeff);
                }
            }

            //remove the simplex from the map if the coefficient becomes 0
            if (mSimplicesMap[aInd] == 0)
            {
                mSimplicesMap.erase_key(aInd);
            }
        }

        //------------------------------------------------------------------------------

        void
        Cochain::removeSimplexToCochain(const id_t aInd)
        {
            //Removing is simply adding a -1 coefficient
            this->addSimplexToCochain(aInd, -1);
        }

        //-----------------------------------------------------------------------------

        void
        Cochain::addCochainToCochain(Cochain * aCochain, int aCoeff)
        {
            if (mDim == aCochain->getDim())
            {
                // Add all the simplices from the input cochain
                for( const auto& [tInd, tCoeff] : aCochain->getSimplicesMap() )
                {
                    this->addSimplexToCochain(tInd,tCoeff*aCoeff);
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
                for( const auto& [tInd, tCoeff] : aCochain->getSimplicesMap() )
                {
                    this->addSimplexToCochain(tInd,-1*tCoeff);
                }
            }
            else
            {
                std::cout << "Dimensions don't agree, not adding the chains" << std::endl;
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
        Cochain::getCoefficient( const id_t aInd )
        {
            if (mSimplicesMap.key_exists( aInd )) {
                return mSimplicesMap[aInd];
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
        Cochain::setCoefficient(const id_t aInd, const int tCoeff)
        {
            if (tCoeff == 0)
            {
                mSimplicesMap.erase_key(aInd);
            }
            else
            {
                mSimplicesMap[aInd] = tCoeff;
            }
        }

        //-----------------------------------------------------------------------------

        void
        Cochain::print()
        {
            std::cout << mDim << "-cochain:" << std::endl;
            for(const auto& [tInd, tCoeff] : mSimplicesMap)
            {
                if (tCoeff != 0)
                {
                    std::cout << "Simplex ID " << tInd << " Coeff: " << tCoeff ;
                    if (mDim == 1)
                    {
                        std::cout << " Nodes: " << mMesh->edge(tInd)->node(0)->id() << " to "  << mMesh->edge(tInd)->node(1)->id() << std::endl;
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
