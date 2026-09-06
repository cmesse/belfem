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
#include "cl_Cochain.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        Cochain::Cochain( const uint aDim, Mesh * aMesh, const bool aIsBound, const bool aIsCobound) :
                mDim( aDim ),
                mMesh( aMesh ),
                mIsBound(aIsBound),
                mIsCobound(aIsCobound)
        {
            if (aDim != 0 and !aIsBound)
            {
                mBoundary = new Cochain(aDim-1, aMesh, true, aIsCobound);
            }
            else
            {
                mBoundary = nullptr;
            }

            if (aDim != 3 and !aIsCobound)
            {
                mCoboundary = new Cochain(aDim+1, aMesh, aIsBound, true);
            }
            else
            {
                mCoboundary = nullptr;
            }
        }

//------------------------------------------------------------------------------

        Cochain::~Cochain()
        {
            delete mBoundary;
            delete mCoboundary;
        }

//------------------------------------------------------------------------------

        void
        Cochain::operator+( Cochain * aCochain )
        {
            this->addCochainToCochain( aCochain, 1);
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
            OrderedMap< index_t, int > & tChainMap = aChain->getSimplicesMap();
            for( const auto& [tID, tCoeff] : mSimplicesMap )
            {
                auto it = tChainMap.find(tID);
                if ( it != tChainMap.end() )
                {
                    tVal += tCoeff*it->second;
                }
            }
            return tVal ;
        }

//------------------------------------------------------------------------------

        void
        Cochain::removeSimplexToCochain(const index_t aID)
        {
            //Removing is simply adding a -1 coefficient
            this->addSimplexToCochain(aID, -1);
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

        void
        Cochain::print()
        {
            Cell< Edge * > & tEdges = mMesh->edges() ;

            std::cout << mDim << "-cochain:" << std::endl;
            for(const auto& [tID, tCoeff] : mSimplicesMap)
            {
                if (tCoeff != 0)
                {
                    std::cout << "Simplex ID " << tID << " Coeff: " << tCoeff ;
                    if (mDim == 1)
                    {
                        std::cout << " Nodes: " << tEdges(tID)->node(0)->index() << " to "  << tEdges(tID)->node(1)->index() << std::endl;
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
