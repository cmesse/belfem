//
// Created by grgia@ge.polymtl.ca on 18/12/23.
//

#include "cl_Matrix.hpp"
#include <tuple>

#ifndef BELFEM_FN_SMITH_HPP
#define BELFEM_FN_SMITH_HPP

namespace belfem
{
//------------------------------------------------------------------------------

    void
    rowExchange(Matrix< int > &aMat, const uint i, const uint j);

//------------------------------------------------------------------------------

    void
    columnExchange(Matrix< int > &aMat, const uint i, const uint j);

//------------------------------------------------------------------------------

    void
    rowMultiply(Matrix< int > &aMat, const uint i);

//------------------------------------------------------------------------------

    void
    columnMultiply(Matrix< int > &aMat, const uint i);

//------------------------------------------------------------------------------

    void
    rowAdd(Matrix< int > &aMat, const uint i, const uint j, const int q);

//------------------------------------------------------------------------------

    void
    columnAdd(Matrix< int > &aMat, const uint i, const uint j, const int q);

//------------------------------------------------------------------------------

    void
    rowExchangeOperation(Matrix< int > &aMat, Matrix< int > &aQ, Matrix< int > &aQ_, const uint i, const uint j);

//------------------------------------------------------------------------------

    void
    rowMultiplyOperation(Matrix< int > &aMat, Matrix< int > &aQ, Matrix< int > &aQ_, const uint i);

//------------------------------------------------------------------------------

    void
    rowAddOperation(Matrix< int > &aMat, Matrix< int > &aQ, Matrix< int > &aQ_, const uint i, const uint j, const int q);

//------------------------------------------------------------------------------

    void
    columnExchangeOperation(Matrix< int > &aMat, Matrix< int > &aR, Matrix< int > &aR_, const uint i, const uint j);

//------------------------------------------------------------------------------

    void
    columnMultiplyOperation(Matrix< int > &aMat, Matrix< int > &aR, Matrix< int > &aR_, const uint i);

//------------------------------------------------------------------------------

    void
    columnAddOperation(Matrix< int > &aMat, Matrix< int > &aR, Matrix< int > &aR_, const uint i, const uint j, const int q);

//------------------------------------------------------------------------------

    void
    partRowReduce(Matrix< int > &aMat, Matrix< int > &aQ, Matrix< int > &aQ_, const uint k, const uint l);

//------------------------------------------------------------------------------

    void
    partColumnReduce(Matrix< int > &aMat, Matrix< int > &aR, Matrix< int > &aR_, const uint k, const uint l);

//------------------------------------------------------------------------------

    std::pair<uint, uint>
    smallestNonzero(Vector< int > &v, const uint k);

//------------------------------------------------------------------------------

    void
    rowPrepare(Matrix< int > &aMat, Matrix< int > &aQ, Matrix< int > &aQ_, const uint k, const uint l);

//------------------------------------------------------------------------------

    void
    rowReduce(Matrix< int > &aMat, Matrix< int > &aQ, Matrix< int > &aQ_, const uint k, const uint l);

//------------------------------------------------------------------------------

    std::tuple< Matrix< int >, Matrix< int >, uint >
    rowEchelon(Matrix< int > &aMat);

//------------------------------------------------------------------------------

    Matrix < int >
    transpose(Matrix< int > &aMat);

//------------------------------------------------------------------------------

    std::tuple< Matrix< int >, Matrix< int > >
    kernelImage(Matrix< int > &aMat);

//------------------------------------------------------------------------------

    std::pair<uint, uint>
    minNonzero(Matrix< int > &aMat, const uint k);

//------------------------------------------------------------------------------

    void
    moveMinNonzero(Matrix< int > &aMat, Matrix< int > &aQ, Matrix< int > &aQ_, Matrix< int > &aR, Matrix< int > &aR_, const uint k);

//------------------------------------------------------------------------------

    std::tuple< bool, uint, uint, int >
    checkForDivisibility(Matrix< int > &aMat, const uint k);

//------------------------------------------------------------------------------

    void
    partSmithForm(Matrix< int > &aMat, Matrix< int > &aQ, Matrix< int > &aQ_, Matrix< int > &aR, Matrix< int > &aR_, const uint k);

//------------------------------------------------------------------------------

    std::tuple< Matrix< int >, Matrix< int >,Matrix< int >,Matrix< int >, uint, uint >
    smithForm(Matrix< int > &aMat);

//------------------------------------------------------------------------------

    Matrix< int >
    SolveInt(Matrix< int > aMat, Matrix< int > &aVec);

//------------------------------------------------------------------------------
}

#endif //BELFEM_FN_SMITH_HPP
