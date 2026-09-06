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

#include "cl_Mesh.hpp"
#include "cl_Chain.hpp"
#include "cl_Cochain.hpp"
#include "cl_Map.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"
#include "cl_DynamicBitset.hpp"
#include "cl_EdgeFactory.hpp"
#include "cl_Element_Factory.hpp"

#ifndef BELFEM_CL_SIMPLICIALCOMPLEX_HPP
#define BELFEM_CL_SIMPLICIALCOMPLEX_HPP

namespace belfem
{
    namespace mesh
    {
//-----------------------------------------------------------------------------
        /**
         * @brief Simplicial complex reduction engine underlying the homology and cohomology computation.
         *
         * @ingroup grp_homology
         * @see @ref homology_cohomology_algorithms
         */
        class SimplicialComplex
        {
            Mesh * mMesh ;

            Cell< Map< index_t, Chain * > > mChainsMap;

            Cell< Map< index_t, Cochain * > > mCochainsMap;

            // Edge indices of the original complex (before reduction). Usefull for the clean() function
            DynamicBitset * mOriginalEdges = nullptr ;

            uint mNumLoopsReduce = 0;
            uint mNumLoopsCoreduce = 0;



//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            SimplicialComplex( Mesh * aMesh, bool aPeriodicity = false );

//-----------------------------------------------------------------------------

            ~SimplicialComplex();

//-----------------------------------------------------------------------------

            void
            reset();

//-----------------------------------------------------------------------------

            void
            create_complex( Mesh * aMesh, bool aPeriodicity );

//-----------------------------------------------------------------------------

            void
            reduce_complexCCR();

//-----------------------------------------------------------------------------

            void
            coreduce_complexCCR();

//-----------------------------------------------------------------------------

            void
            reduce_complexCCR_old();

//-----------------------------------------------------------------------------

            void
            coreduce_complexCCR_old();

//-----------------------------------------------------------------------------

            void
            reduce_pair( const uint k, const uint a, Map<index_t, Chain*>::const_iterator & bit, const int aCoeff );

//-----------------------------------------------------------------------------

            void
            coreduce_pair( const uint k, const uint a, Map<index_t, Cochain*>::const_iterator & bit, const int aCoeff );

//-----------------------------------------------------------------------------

            void
            pReduce(const uint p);

//-----------------------------------------------------------------------------

            void
            pCombine(const uint p);

//-----------------------------------------------------------------------------

            void
            pGeneralizedCombine(const uint p);

//-----------------------------------------------------------------------------

            void
            reduceOmit();

//-----------------------------------------------------------------------------

            void
            reduce_complexPellikka();

//-----------------------------------------------------------------------------

            void
            reduce_complexPellikkaGeneralized();

//-----------------------------------------------------------------------------

            void
            pCoreduce(const uint p);

//-----------------------------------------------------------------------------

            void
            pCocombine(const uint p);

//-----------------------------------------------------------------------------

            void
            pGeneralizedCocombine(const uint p);

//-----------------------------------------------------------------------------

            void
            coreduceOmit();

//-----------------------------------------------------------------------------

            void
            coreduce_complexPellikka();

//-----------------------------------------------------------------------------

            void
            coreduce_complexPellikkaGeneralized();

//-----------------------------------------------------------------------------

            void
            remove_kchainFromMap( const uint k, const index_t aID );

//------------------------------------------------------------------------------

            void
            remove_kcochainFromMap( const uint k, const index_t aID );

//------------------------------------------------------------------------------

            uint
            number_of_0simplices() const;

//------------------------------------------------------------------------------

            uint
            number_of_1simplices() const;

//------------------------------------------------------------------------------

            uint
            number_of_2simplices() const;

//------------------------------------------------------------------------------

            uint
            number_of_ksimplices( const uint k ) const;

//------------------------------------------------------------------------------

            uint
            number_of_kcosimplices( const uint k ) const;

//------------------------------------------------------------------------------

            Chain *
            get_kchain( const uint k, const index_t aID );

//------------------------------------------------------------------------------

            Cochain *
            get_kcochain( const uint k, const index_t aID );

//------------------------------------------------------------------------------

            Map< index_t, Chain * >
            get_kchainMap( const uint k );

//------------------------------------------------------------------------------

            Map< index_t, Cochain * >
            get_kcochainMap( const uint k );

//------------------------------------------------------------------------------

            Chain *
            boundary_of_kchain( const uint k, const index_t aID );

//------------------------------------------------------------------------------

            Cell< Matrix< int > >
            createMatrixFromBoundaryMap();

//------------------------------------------------------------------------------

            Cell< Matrix< int > >
            createMatrixFromCoboundaryMap();

//------------------------------------------------------------------------------

            void
            print_kchains( const uint k );

//------------------------------------------------------------------------------

            void
            create_kComplexField(const uint k, Mesh * aMesh, string aFieldName ) ;

//------------------------------------------------------------------------------

            void
            cocreate_kComplexField( const uint k, Mesh * aMesh, string aFieldName );

//------------------------------------------------------------------------------

            const DynamicBitset &
            original_edges() const ;

//------------------------------------------------------------------------------

        };

//------------------------------------------------------------------------------

        inline const DynamicBitset &
        SimplicialComplex::original_edges() const
        {
            return * mOriginalEdges ;
        }

        inline uint
        SimplicialComplex::number_of_0simplices() const
        {
            return mChainsMap( 0 ).size();
        }

//------------------------------------------------------------------------------

        inline uint
        SimplicialComplex::number_of_1simplices() const
        {
            return mChainsMap( 1 ).size();
        }

//------------------------------------------------------------------------------

        inline uint
        SimplicialComplex::number_of_2simplices() const
        {
            return mChainsMap( 2 ).size();
        }

//------------------------------------------------------------------------------

        inline uint
        SimplicialComplex::number_of_ksimplices( const uint k ) const
        {
            if ( k < 0 )
            {
                return 0;
            }
            else
            {
                return mChainsMap( k ).size();
            }
        }

//------------------------------------------------------------------------------

        inline uint
        SimplicialComplex::number_of_kcosimplices( const uint k ) const
        {
            if ( k < 0 )
            {
                return 0;
            }
            else
            {
                return mCochainsMap( k ).size();
            }
        }


//------------------------------------------------------------------------------

        inline void
        SimplicialComplex::remove_kcochainFromMap(const uint k, const index_t aID)
        {
            if (k <= 3)
            {
                auto & tMap = mCochainsMap( k ).map_data();
                auto it = tMap.find( aID );
                if ( it != tMap.end() )
                {
                    delete it->second;
                    tMap.erase( it );
                }
            }
        }

//------------------------------------------------------------------------------

        // pReduce from Pellikka et al.
        inline void
        SimplicialComplex::pCoreduce(const uint p)
        {
            if (p == 3)
            {
                return ;
            }

#ifdef PERFORMANCE_CHECK
            std::ofstream tFile;
            tFile.open ("Coreduce.txt",std::ios_base::app);
            tFile << mNumLoopsCoreduce << " " << mCochainsMap(0).size()+mCochainsMap(1).size()+mCochainsMap(2).size()+mCochainsMap(3).size() << " \n";
#endif
            bool tRemoved = true;
            while (tRemoved)
            {
                tRemoved = false;
                mNumLoopsCoreduce+=1;

                //Loop over all the p+1 cochains
                Map< index_t, Cochain * > & tMap1 = mCochainsMap(p+1);

                for (auto it = tMap1.begin(), next_it = it; it != tMap1.end(); it = next_it)
                {
                    ++next_it;

                    //Check if the p+1 cochain is a coboundary of exactly one p cochain
                    if (it->second->getBoundary()->getSimplicesMap().size() == 1)
                    {
                        index_t a = it->first;
                        index_t b = it->second->getBoundary()->getSimplicesMap().begin()->first;

                        OrderedMap< index_t, int > & tSimpMap0 = mCochainsMap(p)(b)->getCoboundary()->getSimplicesMap() ;
                        //Update the boundaries
                        for(const auto [tID2, tCoeff2]: tSimpMap0 )
                        {
                            if (tID2 != it->first )
                            {
                                tMap1(tID2)->getBoundary()->setCoefficient(b,0);
                            }
                        }
                        if (p < 2)
                        {
                            OrderedMap< index_t, int > & tSimpMap1 = tMap1(a)->getCoboundary()->getSimplicesMap() ;
                            for(const auto [tID2, tCoeff2]: tSimpMap1 )
                            {
                                mCochainsMap(p+2)(tID2)->getBoundary()->setCoefficient(a,0);
                            }
                        }

                        //Remove a and b from the complex
                        this->remove_kcochainFromMap( p, b );
                        this->remove_kcochainFromMap( p + 1, a );
                        tRemoved = true;
                    }
                }
#ifdef PERFORMANCE_CHECK
                tFile << mNumLoopsCoreduce << " " << mCochainsMap(0).size()+mCochainsMap(1).size()+mCochainsMap(2).size()+mCochainsMap(3).size() << " \n";
#endif
            }
#ifdef PERFORMANCE_CHECK
            tFile.close();
#endif
        }
//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_SIMPLICIALCOMPLEX_HPP
