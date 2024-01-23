//
// Created by grgia@ge.polymtl.ca on 06/12/23.
//

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
        class Cochain
        {
            //Dimensions of the cochain
            int mDim ;

            //Mesh
            Mesh * mMesh ;

            //Map of the k-simplices coefficients
            Map< id_t, int > mSimplicesMap ;

            //Coboundary cochain of the given cochain
            Cochain * mCoboundary = nullptr;

            //-----------------------------------------------------------------------------
        public:
            //-----------------------------------------------------------------------------

            Cochain(const uint aDim, Mesh* aMesh) ;

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
            addSimplexToCochain(const id_t aInd, const int aCoeff ) ;

            //-----------------------------------------------------------------------------

            void
            removeSimplexToCochain(const id_t aInd) ;

            //-----------------------------------------------------------------------------

            void
            addCochainToCochain(Cochain * aCochain, int aCoeff) ;

            //-----------------------------------------------------------------------------

            void
            removeCochainFromCochain(Cochain * aChain) ;

            //-----------------------------------------------------------------------------

            Map< id_t, int >
            getSimplicesMap() ;

            //-----------------------------------------------------------------------------

            int
            getCoefficient(const id_t aInd) ;

            //-----------------------------------------------------------------------------

            Cochain*
            getCoboundary() ;

            //-----------------------------------------------------------------------------

            int
            getDim() ;

            //-----------------------------------------------------------------------------

            void
            setCoefficient(const id_t aInd, const int tCoeff) ;

            //-----------------------------------------------------------------------------

            void
            print() ;

            //-----------------------------------------------------------------------------

        };
    }
}
#endif //BELFEM_CL_COCHAIN_HPP
