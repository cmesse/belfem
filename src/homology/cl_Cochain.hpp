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
            add_simplex_to_cochain(const index_t aID, const int aCoeff ) ;

//-----------------------------------------------------------------------------

            void
            remove_simplex_from_cochain(const index_t aID) ;

//-----------------------------------------------------------------------------

            void
            add_cochain_to_cochain(Cochain * aCochain, const int aCoeff) ;

//-----------------------------------------------------------------------------

            void
            remove_cochain_from_cochain(Cochain * aChain) ;

//-----------------------------------------------------------------------------

            OrderedMap< index_t, int > &
            get_simplices_map() ;

//-----------------------------------------------------------------------------

            int
            get_coefficient(const index_t aID) ;

//-----------------------------------------------------------------------------

            Cochain*
            get_coboundary() ;

//-----------------------------------------------------------------------------

            Cochain*
            get_boundary() ;

//-----------------------------------------------------------------------------

            void
            add_simplex_to_coboundary(index_t aID, const int aCoeff) ;

//-----------------------------------------------------------------------------

            void
            add_simplex_to_boundary(index_t aID, const int aCoeff) ;

//-----------------------------------------------------------------------------

            int
            get_dim() ;

//-----------------------------------------------------------------------------

            bool
            isBound() ;

//-----------------------------------------------------------------------------

            bool
            isCobound() ;

//-----------------------------------------------------------------------------

            void
            set_coefficient(const index_t aID, const int tCoeff) ;

//-----------------------------------------------------------------------------

            void
            print() ;

//-----------------------------------------------------------------------------

        };

//-----------------------------------------------------------------------------

        inline OrderedMap< index_t, int > &
        Cochain::get_simplices_map()
        {
            return mSimplicesMap;
        }

//------------------------------------------------------------------------------

        inline int
        Cochain::get_coefficient( const index_t aID )
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
        Cochain::get_coboundary()
        {
            return mCoboundary;
        }

//------------------------------------------------------------------------------

        inline Cochain *
        Cochain::get_boundary()
        {
            return mBoundary;
        }

//------------------------------------------------------------------------------

        inline void
        Cochain::add_simplex_to_coboundary(index_t aID, const int aCoeff)
        {
            mCoboundary->add_simplex_to_cochain(aID,aCoeff);
        }

//------------------------------------------------------------------------------

        inline void
        Cochain::add_simplex_to_boundary(index_t aID, const int aCoeff)
        {
            mBoundary->add_simplex_to_cochain(aID,aCoeff);
        }

//------------------------------------------------------------------------------

        inline int
        Cochain::get_dim()
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
        Cochain::set_coefficient(const index_t aID, const int tCoeff)
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
        Cochain::add_simplex_to_cochain( const index_t aID, const int aCoeff )
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
        Cochain::add_cochain_to_cochain(Cochain * aCochain, int aCoeff)
        {
            if ( aCoeff == 0 ) return;  // Early exit - nothing to add

            if (mDim == aCochain->get_dim())
            {
                // Add all the simplices from the input chain
                for( const auto& [tID, tCoeff] : aCochain->get_simplices_map() )
                {
                    this->add_simplex_to_cochain(tID,tCoeff*aCoeff);
                }
                if(mDim != 0 and !mIsBound)
                {
                    mBoundary->add_cochain_to_cochain(aCochain->get_boundary(),aCoeff);
                }

                if (mDim != 3 and !mIsCobound)
                {
                    mCoboundary->add_cochain_to_cochain(aCochain->get_coboundary(),aCoeff);
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
