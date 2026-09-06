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
#include "cl_OrderedMap.hpp"
#include "cl_Cell.hpp"
#include "cl_Vertex.hpp"
#include "cl_Mesh.hpp"

#ifndef BELFEM_CL_CHAIN_HPP
#define BELFEM_CL_CHAIN_HPP

namespace belfem
{
    namespace mesh
    {
//-----------------------------------------------------------------------------
        /**
         * @brief Formal sum of k-simplices (homology).
         *
         * @ingroup grp_homology
         * @see @ref homology_cohomology_theory_and_implementation
         */
        class Chain
        {
            // Dimensions of the chain
            int mDim ;

            // Mesh
            Mesh * mMesh ;

            // Map of the k-simplices coefficients
            OrderedMap< index_t, int > mSimplicesMap ;

            // Boundary chain of the given chain
            Chain * mBoundary = nullptr;

            // Coboundary chain of the given chain
            Chain * mCoboundary = nullptr;

            bool mIsBound = false;
            bool mIsCobound = false;

//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            Chain( const uint aDim, Mesh * aMesh, const bool aIsBound, const bool aIsCobound ) ;

//-----------------------------------------------------------------------------

            ~Chain();

//-----------------------------------------------------------------------------

            void
            operator+( Chain * aChain ) ;

//-----------------------------------------------------------------------------

            void
            operator-( Chain * aChain ) ;

//-----------------------------------------------------------------------------

            void
            operator*( const int aMult ) ;

//-----------------------------------------------------------------------------

            int
            operator()( Chain * aChain ) ;

//-----------------------------------------------------------------------------

            void
            addSimplexToChain(const id_t aID, const int aCoeff ) ;

//-----------------------------------------------------------------------------

            void
            removeSimplexToChain(const id_t aID) ;

//-----------------------------------------------------------------------------

            void
            addChainToChain(Chain * aChain, int aCoeff) ;

//-----------------------------------------------------------------------------

            void
            removeChainFromChain(Chain * aChain) ;

//-----------------------------------------------------------------------------

            OrderedMap< id_t, int > &
            getSimplicesMap() ;

//-----------------------------------------------------------------------------

            int
            getCoefficient(const id_t aID) ;

//-----------------------------------------------------------------------------

            Chain *
            getBoundary() ;

//-----------------------------------------------------------------------------

            Chain *
            getCoboundary() ;

//-----------------------------------------------------------------------------

            void
            add_simplex_to_boundary(id_t aID, const int aCoeff) ;

//-----------------------------------------------------------------------------

            void
            add_simplex_to_coboundary(id_t aID, const int aCoeff) ;

//-----------------------------------------------------------------------------

            int
            getDim() ;

//-----------------------------------------------------------------------------

            void
            setCoefficient(const id_t aID, const int tCoeff) ;

//-----------------------------------------------------------------------------

            void
            print() ;

//-----------------------------------------------------------------------------

        };


//-----------------------------------------------------------------------------

        inline OrderedMap< id_t, int > &
        Chain::getSimplicesMap()
        {
            return mSimplicesMap;
        }

//------------------------------------------------------------------------------

        inline int
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

        inline Chain *
        Chain::getBoundary()
        {
            return mBoundary;
        }

//-----------------------------------------------------------------------------

        inline Chain *
        Chain::getCoboundary()
        {
            return mCoboundary;
        }

//-----------------------------------------------------------------------------

        inline void
        Chain::add_simplex_to_boundary( id_t aID, const int aCoeff )
        {
            mBoundary->addSimplexToChain(aID,aCoeff);
        }

//-----------------------------------------------------------------------------

        inline void
        Chain::add_simplex_to_coboundary(id_t aID, const int aCoeff)
        {
            mCoboundary->addSimplexToChain(aID,aCoeff);
        }

//-----------------------------------------------------------------------------

        inline int
        Chain::getDim()
        {
            return mDim;
        }

//-----------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_CHAIN_HPP
