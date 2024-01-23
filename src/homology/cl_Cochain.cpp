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
                for (uint i = 0; i < mMesh->node(aInd)->number_of_edges(); i++)
                {
                    if (mMesh->node(aInd)->edge(i)->node(0)->id() == aInd)
                    {
                        mCoboundary->addSimplexToCochain(mMesh->node(aInd)->edge(i)->id(),-1*aCoeff);
                    }
                    else
                    {
                        mCoboundary->addSimplexToCochain(mMesh->node(aInd)->edge(i)->id(),1*aCoeff);
                    }
                }
            }

                // 1-cochain
            else if (mDim == 1)
            {
                for (uint i = 0; i < mMesh->edge(aInd)->number_of_elements(); i++)
                {

                    // Probably can be optimized (how to know the orientation of a given edge in a triangle)
                    int tMul = 0;
                    for (uint j = 0; j < 3; j++)
                    {
                        if (mMesh->edge(aInd)->element(i)->node(j)->id() == mMesh->edge(aInd)->node(0)->id())
                        {
                            if (mMesh->edge(aInd)->element(i)->node((j+1)%3)->id() == mMesh->edge(aInd)->node(1)->id())
                            {
                                tMul = 1;
                                break;
                            }
                            else
                            {
                                tMul = -1;
                                break;
                            }
                        }
                    }
                    mCoboundary->addSimplexToCochain(mMesh->edge(aInd)->element(i)->id(),tMul*aCoeff);
                }
            }

            // todo: do the same with 2-cochain (elements/volumes) when we have 3-D

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
