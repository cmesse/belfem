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
        Chain::addSimplexToChain( const id_t aID, const int aCoeff )
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


            // Update the boundary of the chain from what's just been added
            // 1-chain
            if (mDim == 1)
            {
                mBoundary->addSimplexToChain(mMesh->edge(aID)->node(0)->id(),-1*aCoeff);
                mBoundary->addSimplexToChain(mMesh->edge(aID)->node(1)->id(),1*aCoeff);
            }

            // 2-chain
            else if (mDim == 2)
            {
                if (mMesh->number_of_dimensions() == 2)
                {
                    Element * tElement = mMesh->element(aID);
                    for ( uint i = 0; i < tElement->number_of_edges(); ++i )
                    {
                        mBoundary->addSimplexToChain(tElement->edge(i)->id(),
                                                     (tElement->edge_orientation(i)?1:-1)*aCoeff);
                    }
                }
                else
                {
                    Face * tFace = mMesh->face(aID);
                    for ( uint i = 0; i < tFace->number_of_edges(); ++i )
                    {
                        mBoundary->addSimplexToChain(tFace->edge(i)->id(),
                                                     (tFace->edge_orientation(i)?1:-1)*aCoeff);
                    }
                }
            }

            // 3-chain
            if (mMesh->number_of_dimensions() == 3 && mDim == 3)
            {
                Element * tElement = mMesh->element(aID);
                for ( uint i = 0; i < tElement->number_of_faces(); ++i )
                {
                    mBoundary->addSimplexToChain(tElement->face(i)->id(),
                                                 (tElement->face(i)->master()->id() == tElement->id()?1:-1)*aCoeff);
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
        Chain::removeSimplexToChain(const id_t aID)
        {
            //Removing is simply adding a -1 coefficient
            this->addSimplexToChain(aID, -1);
        }

        //-----------------------------------------------------------------------------

        void
        Chain::addChainToChain(Chain * aChain, int aCoeff)
        {
            if (mDim == aChain->getDim())
            {
                // Add all the simplices from the input chain
                for( const auto& [tID, tCoeff] : aChain->getSimplicesMap() )
                {
                    this->addSimplexToChain(tID,tCoeff*aCoeff);
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
                for( const auto& [tID, tCoeff] : aChain->getSimplicesMap() )
                {
                    this->addSimplexToChain(tID,-1*tCoeff);
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
        Chain::getCoefficient( const id_t aID )
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
        Chain::setCoefficient(const id_t aID, const int tCoeff)
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
        Chain::print()
        {
            std::cout << mDim << "-chain:" << std::endl;
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
