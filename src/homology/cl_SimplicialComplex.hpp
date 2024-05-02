//
// Created by grgia@ge.polymtl.ca on 07/12/23.
//

#include "cl_Mesh.hpp"
#include "cl_Chain.hpp"
#include "cl_Cochain.hpp"
#include "cl_Map.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"
#include "cl_EdgeFactory.hpp"

#ifndef BELFEM_CL_SIMPLICIALCOMPLEX_HPP
#define BELFEM_CL_SIMPLICIALCOMPLEX_HPP

namespace belfem
{
    namespace mesh
    {
//-----------------------------------------------------------------------------
        class SimplicialComplex
        {

            Cell< Map< id_t, Chain * > > mChainsMap;

            Cell< Map< id_t, Cochain * > > mCochainsMap;



            //-----------------------------------------------------------------------------
        public:
            //-----------------------------------------------------------------------------

            SimplicialComplex( Mesh * aMesh );

            //-----------------------------------------------------------------------------

            ~SimplicialComplex();

            //-----------------------------------------------------------------------------

            void
            reset();

            //-----------------------------------------------------------------------------

            void
            create_complex( Mesh * aMesh );

            //-----------------------------------------------------------------------------

            void
            reduce_complexCCR();

            //-----------------------------------------------------------------------------

            void
            coreduce_complexCCR();

            //-----------------------------------------------------------------------------

            void
            reduce_pair( const uint k, const uint a, const uint b, Chain* aChain );

            //-----------------------------------------------------------------------------

            void
            coreduce_pair( const uint k, const uint a, const uint b, Cochain* aCochain );

            //-----------------------------------------------------------------------------

            void
            pReduce(const uint p);

            //-----------------------------------------------------------------------------

            void
            pCombine(const uint p);

            //-----------------------------------------------------------------------------

            void
            reduceOmit();

            //-----------------------------------------------------------------------------

            void
            reduce_complexPellikka();

            //-----------------------------------------------------------------------------

            void
            pCoreduce(const uint p);

            //-----------------------------------------------------------------------------

            void
            pCocombine(const uint p);

            //-----------------------------------------------------------------------------

            void
            coreduceOmit();

            //-----------------------------------------------------------------------------

            void
            coreduce_complexPellikka();

            //-----------------------------------------------------------------------------

            void
            remove_kchainFromMap( const uint k, const id_t aInd );

            //------------------------------------------------------------------------------

            void
            remove_kcochainFromMap( const uint k, const id_t aInd );

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
            get_kchain( const uint k, const id_t aInd );

            //------------------------------------------------------------------------------

            Map< id_t, Chain * >
            get_kchainMap( const uint k );

            //------------------------------------------------------------------------------

            Map< id_t, Cochain * >
            get_kcochainMap( const uint k );

            //------------------------------------------------------------------------------

            Chain *
            boundary_of_kchain( const uint k, const id_t aInd );

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
            create_kComplexField( const uint k, Mesh * aMesh );

            //------------------------------------------------------------------------------

            void
            cocreate_kComplexField( const uint k, Mesh * aMesh );

            //------------------------------------------------------------------------------

        };

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
    }
}

#endif //BELFEM_CL_SIMPLICIALCOMPLEX_HPP
