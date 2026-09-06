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
#include "typedefs.hpp"
#include "assert.hpp"
#include "cl_Map.hpp"
#include "cl_Cell.hpp"
#include "cl_Vertex.hpp"
#include "cl_Mesh.hpp"
#include "cl_Chain.hpp"

#ifndef BELFEM_CL_COCHAIN_HPP
#define BELFEM_CL_COCHAIN_HPP

namespace belfem
{
    namespace mesh
    {
//-----------------------------------------------------------------------------
        /**
         * @brief Cochain: a functional on chains (cohomology).
         *
         * @ingroup grp_homology
         * @see @ref homology_cohomology_theory_and_implementation
         */
        class Cochain
        {
            //Dimensions of the cochain
            int mDim ;

            //Mesh
            Mesh * mMesh ;

            //Map of the k-simplices coefficients
            OrderedMap< index_t, int > mSimplicesMap ;

            //Coboundary cochain of the given cochain
            Cochain * mCoboundary = nullptr;

            //Boundary
            Cochain * mBoundary = nullptr;

            bool mIsBound = false;
            bool mIsCobound = false;

//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            Cochain(const uint aDim, Mesh* aMesh,
                const bool aIsBound, const bool aIsCobound) ;

//-----------------------------------------------------------------------------

            ~Cochain();

//-----------------------------------------------------------------------------

            void
            operator+( Cochain * aCochain ) ;

//-----------------------------------------------------------------------------

            void
            operator-( Cochain * aCochain ) ;

//-----------------------------------------------------------------------------

            int
            operator()( Chain * aChain ) ;

//-----------------------------------------------------------------------------

            void
            addSimplexToCochain(const index_t aID, const int aCoeff ) ;

//-----------------------------------------------------------------------------

            void
            removeSimplexToCochain(const index_t aID) ;

//-----------------------------------------------------------------------------

            void
            addCochainToCochain(Cochain * aCochain, const int aCoeff) ;

//-----------------------------------------------------------------------------

            void
            removeCochainFromCochain(Cochain * aChain) ;

//-----------------------------------------------------------------------------

            OrderedMap< index_t, int > &
            getSimplicesMap() ;

//-----------------------------------------------------------------------------

            int
            getCoefficient(const index_t aID) ;

//-----------------------------------------------------------------------------

            Cochain*
            getCoboundary() ;

//-----------------------------------------------------------------------------

            Cochain*
            getBoundary() ;

//-----------------------------------------------------------------------------

            void
            add_simplex_to_coboundary(index_t aID, const int aCoeff) ;

//-----------------------------------------------------------------------------

            void
            add_simplex_to_boundary(index_t aID, const int aCoeff) ;

//-----------------------------------------------------------------------------

            int
            getDim() ;

//-----------------------------------------------------------------------------

            bool
            isBound() ;

//-----------------------------------------------------------------------------

            bool
            isCobound() ;

//-----------------------------------------------------------------------------

            void
            setCoefficient(const index_t aID, const int tCoeff) ;

//-----------------------------------------------------------------------------

            void
            print() ;

//-----------------------------------------------------------------------------

        };

//-----------------------------------------------------------------------------

        inline OrderedMap< index_t, int > &
        Cochain::getSimplicesMap()
        {
            return mSimplicesMap;
        }

//------------------------------------------------------------------------------

        inline int
        Cochain::getCoefficient( const index_t aID )
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

        inline Cochain *
        Cochain::getCoboundary()
        {
            return mCoboundary;
        }

//------------------------------------------------------------------------------

        inline Cochain *
        Cochain::getBoundary()
        {
            return mBoundary;
        }

//------------------------------------------------------------------------------

        inline void
        Cochain::add_simplex_to_coboundary(index_t aID, const int aCoeff)
        {
            mCoboundary->addSimplexToCochain(aID,aCoeff);
        }

//------------------------------------------------------------------------------

        inline void
        Cochain::add_simplex_to_boundary(index_t aID, const int aCoeff)
        {
            mBoundary->addSimplexToCochain(aID,aCoeff);
        }

//------------------------------------------------------------------------------

        inline int
        Cochain::getDim()
        {
            return mDim;
        }

//------------------------------------------------------------------------------

        inline bool
        Cochain::isBound()
        {
            return mIsBound;
        }

//------------------------------------------------------------------------------

        inline bool
        Cochain::isCobound()
        {
            return mIsCobound;
        }

//------------------------------------------------------------------------------

        inline void
        Cochain::setCoefficient(const index_t aID, const int tCoeff)
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

        inline void
        Cochain::addSimplexToCochain( const index_t aID, const int aCoeff )
        {
            if ( aCoeff == 0 ) return;

            auto & tMap = mSimplicesMap.map_data() ;
            auto tIterator = tMap.find(aID) ;

            if ( tIterator != tMap.end() )
            {
                // Key exists - update in place
                tIterator->second += aCoeff;

                // Remove if coefficient becomes 0
                if ( tIterator->second == 0 )
                {
                    tMap.erase( tIterator ) ;
                }
            }
            else // we know that aCoeff != 0
            {
                // only add if non-zero coefficient
                tMap[ aID ] = aCoeff ;
            }

        }


//-----------------------------------------------------------------------------

        inline void
        Cochain::addCochainToCochain(Cochain * aCochain, int aCoeff)
        {
            if ( aCoeff == 0 ) return;  // Early exit - nothing to add

            if (mDim == aCochain->getDim())
            {
                // Add all the simplices from the input chain
                for( const auto& [tID, tCoeff] : aCochain->getSimplicesMap() )
                {
                    this->addSimplexToCochain(tID,tCoeff*aCoeff);
                }
                if(mDim != 0 and !mIsBound)
                {
                    mBoundary->addCochainToCochain(aCochain->getBoundary(),aCoeff);
                }

                if (mDim != 3 and !mIsCobound)
                {
                    mCoboundary->addCochainToCochain(aCochain->getCoboundary(),aCoeff);
                }

            }
            else
            {
                std::cout << "Dimensions don't agree, not adding the cochain" << std::endl;
            }
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_COCHAIN_HPP
