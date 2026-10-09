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
#include "cl_Chain.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        Chain::Chain( const uint aDim, Mesh * aMesh, const bool aIsBound, const bool aIsCobound) :
                mDim( aDim ),
                mMesh( aMesh ),
                mIsBound(aIsBound),
                mIsCobound(aIsCobound)
        {
            if (aDim != 0 and !aIsBound)
            {
                mBoundary = new Chain(aDim-1, aMesh, true, aIsCobound);
            }
            else
            {
                mBoundary = nullptr;
            }

            if (aDim != 3 and !aIsCobound)
            {
                mCoboundary = new Chain(aDim+1, aMesh, aIsBound, true);
            }
            else
            {
                mCoboundary = nullptr;
            }
        }

//------------------------------------------------------------------------------

        Chain::~Chain()
        {
            delete mBoundary;
            delete mCoboundary;
        }

//------------------------------------------------------------------------------

        void
        Chain::operator+( Chain * aChain )
        {
            this->add_chain_to_chain( aChain, 1 );
        }

//-----------------------------------------------------------------------------

        void
        Chain::operator-( Chain * aChain )
        {
            this->remove_chain_from_chain( aChain );
        }

//-----------------------------------------------------------------------------

        void
        Chain::operator*( const int aMult )
        {
            for (auto [tID, tCoeff] : mSimplicesMap)
            {
                mSimplicesMap[tID] *= aMult ;
            }
        }

//-----------------------------------------------------------------------------

        int
        Chain::operator()( Chain * aChain )
        {
            int tVal = 0;
            OrderedMap< index_t, int > & tChainMap = aChain->get_simplices_map();
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
        Chain::add_simplex_to_chain( const index_t aIndex, const int aCoeff )
        {

            // Add to the existing simplex, or create a new one
            if (mSimplicesMap.key_exists( aIndex ))
            {
                mSimplicesMap( aIndex ) += aCoeff;
            }
            else
            {
                mSimplicesMap[ aIndex ] = aCoeff;
            }


            //remove the simplex from the map if the coefficient becomes 0
            if (mSimplicesMap( aIndex ) == 0)
            {
                mSimplicesMap.erase_key( aIndex );
            }
        }

//------------------------------------------------------------------------------

        void
        Chain::remove_simplex_from_chain( const index_t aIndex )
        {
            //Removing is simply adding a -1 coefficient
            this->add_simplex_to_chain( aIndex, -1);
        }

//-----------------------------------------------------------------------------

        void
        Chain::add_chain_to_chain(Chain * aChain, int aCoeff)
        {
            if (mDim == aChain->get_dim())
            {
                // Add all the simplices from the input chain
                for( const auto& [tID, tCoeff] : aChain->get_simplices_map() )
                {
                    this->add_simplex_to_chain(tID,tCoeff*aCoeff);
                }
                if(mDim != 0 and !mIsBound)
                {
                    mBoundary->add_chain_to_chain(aChain->get_boundary(),aCoeff);
                }

                if (mDim != 3 and !mIsCobound)
                {
                    mCoboundary->add_chain_to_chain(aChain->get_coboundary(),aCoeff);
                }

            }
            else
            {
                std::cout << "Dimensions don't agree, not adding the chains" << std::endl;
            }
        }

//-----------------------------------------------------------------------------

        void
        Chain::remove_chain_from_chain( Chain * aChain )
        {
            if (mDim == aChain->get_dim())
            {
                for( const auto& [tID, tCoeff] : aChain->get_simplices_map() )
                {
                    this->add_simplex_to_chain(tID,-1*tCoeff);
                }
            }
            else
            {
                std::cout << "Dimensions don't agree, not adding the chains" << std::endl;
            }
        }

//-----------------------------------------------------------------------------


        void
        Chain::set_coefficient( const index_t aIndex, const int aCoeff )
        {
            if ( aCoeff == 0 )
            {
                mSimplicesMap.erase_key( aIndex );
            }
            else
            {
                mSimplicesMap[ aIndex ] = aCoeff;
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
                        std::cout << " Nodes: " << mMesh->edges()(tID)->node(0)->index() << " to "  << mMesh->edges()(tID)->node(1)->index() << std::endl;
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
