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
            add_simplex_to_chain( const index_t aIndex, const int aCoeff ) ;

//-----------------------------------------------------------------------------

            void
            remove_simplex_from_chain( const index_t aIndex ) ;

//-----------------------------------------------------------------------------

            void
            add_chain_to_chain(Chain * aChain, int aCoeff) ;

//-----------------------------------------------------------------------------

            void
            remove_chain_from_chain(Chain * aChain) ;

//-----------------------------------------------------------------------------

            OrderedMap< index_t, int > &
            get_simplices_map() ;

//-----------------------------------------------------------------------------

            int
            get_coefficient( const index_t aIndex ) ;

//-----------------------------------------------------------------------------

            Chain *
            get_boundary() ;

//-----------------------------------------------------------------------------

            Chain *
            get_coboundary() ;

//-----------------------------------------------------------------------------

            void
            add_simplex_to_boundary( index_t aIndex, const int aCoeff ) ;

//-----------------------------------------------------------------------------

            void
            add_simplex_to_coboundary( index_t aIndex, const int aCoeff ) ;

//-----------------------------------------------------------------------------

            int
            get_dim() ;

//-----------------------------------------------------------------------------

            void
            set_coefficient( const index_t aIndex, const int aCoeff ) ;

//-----------------------------------------------------------------------------

            void
            print() ;

//-----------------------------------------------------------------------------

        };


//-----------------------------------------------------------------------------

        inline OrderedMap< index_t, int > &
        Chain::get_simplices_map()
        {
            return mSimplicesMap;
        }

//------------------------------------------------------------------------------

        inline int
        Chain::get_coefficient( const index_t aIndex )
        {
            auto it = mSimplicesMap.find( aIndex );
            if (it != mSimplicesMap.end()) {
                return it->second;
            }
            else {
                return 0;
            }
        }

//------------------------------------------------------------------------------

        inline Chain *
        Chain::get_boundary()
        {
            return mBoundary;
        }

//-----------------------------------------------------------------------------

        inline Chain *
        Chain::get_coboundary()
        {
            return mCoboundary;
        }

//-----------------------------------------------------------------------------

        inline void
        Chain::add_simplex_to_boundary( index_t aIndex, const int aCoeff )
        {
            mBoundary->add_simplex_to_chain( aIndex, aCoeff );
        }

//-----------------------------------------------------------------------------

        inline void
        Chain::add_simplex_to_coboundary( index_t aIndex, const int aCoeff)
        {
            mCoboundary->add_simplex_to_chain( aIndex, aCoeff );
        }

//-----------------------------------------------------------------------------

        inline int
        Chain::get_dim()
        {
            return mDim;
        }

//-----------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_CHAIN_HPP
