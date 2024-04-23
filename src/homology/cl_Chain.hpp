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

#ifndef BELFEM_CL_CHAIN_HPP
#define BELFEM_CL_CHAIN_HPP

namespace belfem
{
    namespace mesh
    {
//-----------------------------------------------------------------------------
        class Chain
        {
            //Dimensions of the chain
            int mDim ;

            //Mesh
            Mesh * mMesh ;

            //Map of the k-simplices coefficients
            Map< id_t, int > mSimplicesMap;

            //Boundary chain of the given chain
            Chain * mBoundary = nullptr;

            //-----------------------------------------------------------------------------
        public:
            //-----------------------------------------------------------------------------

            Chain(const uint aDim, Mesh* aMesh) ;

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
            addSimplexToChain(const id_t aInd, const int aCoeff ) ;

            //-----------------------------------------------------------------------------

            void
            removeSimplexToChain(const id_t aInd) ;

            //-----------------------------------------------------------------------------

            void
            addChainToChain(Chain * aChain, int aCoeff) ;

            //-----------------------------------------------------------------------------

            void
            removeChainFromChain(Chain * aChain) ;

            //-----------------------------------------------------------------------------

            Map< id_t, int >
            getSimplicesMap() ;

            //-----------------------------------------------------------------------------

            int
            getCoefficient(const id_t aInd) ;

            //-----------------------------------------------------------------------------

            Chain*
            getBoundary() ;

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
#endif //BELFEM_CL_CHAIN_HPP
