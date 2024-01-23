//
// Created by grgia@ge.polymtl.ca on 06/12/23.
//

//todo: Maybe create a parent class for chain and cochain, as they have very similar methods

#include <iostream>
#include "cl_Chain.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        Chain::Chain( const uint aDim, Mesh * aMesh ) :
                mDim( aDim ),
                mMesh( aMesh )
        {
            if (aDim != 0)
            {
                mBoundary = new Chain(aDim-1, aMesh);
            }
            else
            {
                mBoundary = nullptr;
            }
        }

        //------------------------------------------------------------------------------

        Chain::~Chain()
        {
            delete mBoundary;
        }

        //------------------------------------------------------------------------------

        void
        Chain::operator+( Chain * aChain )
        {
            this->addChainToChain( aChain, 1 );
        }

        //-----------------------------------------------------------------------------

        void
        Chain::operator-( Chain * aChain )
        {
            this->removeChainFromChain( aChain );
        }

        //-----------------------------------------------------------------------------

        void
        Chain::addSimplexToChain( const id_t aInd, const int aCoeff )
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


            // Update the boundary of the chain from what's just been added
            // 1-chain
            if (mDim == 1)
            {
                mBoundary->addSimplexToChain(mMesh->edge(aInd)->node(0)->id(),-1*aCoeff);
                mBoundary->addSimplexToChain(mMesh->edge(aInd)->node(1)->id(),1*aCoeff);
            }

                // 2-chain
            else if (mDim == 2)
            {
                for ( int i = 0; i < 3; ++i )
                {
                    if (mMesh->element(aInd)->edge(i)->node(0)->id() == mMesh->element(aInd)->node(i)->id())
                    {
                        mBoundary->addSimplexToChain(mMesh->element(aInd)->edge(i)->id(),1*aCoeff);
                    }
                    else
                    {
                        mBoundary->addSimplexToChain(mMesh->element(aInd)->edge(i)->id(),-1*aCoeff);
                    }
                }
            }

            // todo: do the same with 3-chain (volumes/triangles) when we have 3-D

            //remove the simplex from the map if the coefficient becomes 0
            if (mSimplicesMap[aInd] == 0)
            {
                mSimplicesMap.erase_key(aInd);
            }
        }

        //------------------------------------------------------------------------------

        void
        Chain::removeSimplexToChain(const id_t aInd)
        {
            //Removing is simply adding a -1 coefficient
            this->addSimplexToChain(aInd, -1);
        }

        //-----------------------------------------------------------------------------

        void
        Chain::addChainToChain(Chain * aChain, int aCoeff)
        {
            if (mDim == aChain->getDim())
            {
                // Add all the simplices from the input chain
                for( const auto& [tInd, tCoeff] : aChain->getSimplicesMap() )
                {
                    this->addSimplexToChain(tInd,tCoeff*aCoeff);
                }
            }
            else
            {
                std::cout << "Dimensions don't agree, not adding the chains" << std::endl;
            }
        }

        //-----------------------------------------------------------------------------

        void
        Chain::removeChainFromChain(Chain * aChain)
        {
            if (mDim == aChain->getDim())
            {
                for( const auto& [tInd, tCoeff] : aChain->getSimplicesMap() )
                {
                    this->addSimplexToChain(tInd,-1*tCoeff);
                }
            }
            else
            {
                std::cout << "Dimensions don't agree, not adding the chains" << std::endl;
            }
        }

        //-----------------------------------------------------------------------------

        Map< id_t, int >
        Chain::getSimplicesMap()
        {
            return mSimplicesMap;
        }

        //------------------------------------------------------------------------------

        int
        Chain::getCoefficient( const id_t aInd )
        {
            if (mSimplicesMap.key_exists( aInd )) {
                return mSimplicesMap[aInd];
            }
            else {
                return 0;
            }
        }

        //------------------------------------------------------------------------------

        Chain*
        Chain::getBoundary()
        {
            return mBoundary;
        }

        //-----------------------------------------------------------------------------

        int
        Chain::getDim()
        {
            return mDim;
        }

        //-----------------------------------------------------------------------------

        void
        Chain::setCoefficient(const id_t aInd, const int tCoeff)
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
        Chain::print()
        {
            std::cout << mDim << "-chain:" << std::endl;
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
