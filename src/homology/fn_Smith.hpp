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

#include "cl_Matrix.hpp"
#include <tuple>
#include "fn_trans.hpp"

#ifndef BELFEM_FN_SMITH_HPP
#define BELFEM_FN_SMITH_HPP

namespace belfem
{
//------------------------------------------------------------------------------

    template<typename T>
    void
    rowExchange(Matrix< T > &aMat, const uint i, const uint j)
    {
        const uint n = aMat.n_cols();
        int v;
        for (uint k = 0; k < n; k++)
        {
            v = aMat(i-1,k);
            aMat(i-1,k) = aMat(j-1,k);
            aMat(j-1,k)=v;
        }
    }


//------------------------------------------------------------------------------

    template<typename T>
    void
    columnExchange(Matrix< T > &aMat, const uint i, const uint j)
    {
        const uint m = aMat.n_rows();
        int v;
        for (uint k = 0; k < m; k++)
        {
            v = aMat(k,i-1);
            aMat(k,i-1) = aMat(k,j-1);
            aMat(k,j-1)=v;
        }
    }

//------------------------------------------------------------------------------

    template<typename T>
    void
    rowMultiply(Matrix< T > &aMat, const uint i)
    {
        const uint n = aMat.n_cols();
        for (uint k = 0; k < n; k++)
        {
            aMat(i-1,k)*=-1 ;
        }
    }

//------------------------------------------------------------------------------

    template<typename T>
    void
    columnMultiply(Matrix< T > &aMat, const uint i)
    {
        const uint m = aMat.n_rows();
        for (uint k = 0; k < m; k++)
        {
            aMat(k,i-1)*=-1 ;
        }
    }

//------------------------------------------------------------------------------

    template<typename T>
    void
    rowAdd(Matrix< T > &aMat, const uint i, const uint j, const int q)
    {
        const uint n = aMat.n_cols();
        for (uint k = 0; k < n; k++)
        {
            aMat(i-1,k) += q*aMat(j-1,k) ;
        }
    }

//------------------------------------------------------------------------------

    template<typename T>
    void
    columnAdd(Matrix< T > &aMat, const uint i, const uint j, const int q)
    {
        const uint m = aMat.n_rows();
        for (uint k = 0; k < m; k++)
        {
            aMat(k,j-1) += q*aMat(k,i-1) ;
        }
    }

//------------------------------------------------------------------------------

    template<typename T>
    void
    rowExchangeOperation(Matrix< T > &aMat, Matrix< T > &aQ, Matrix< T > &aQ_, const uint i, const uint j)
    {
        rowExchange(aMat,i,j) ;
        rowExchange(aQ_,i,j) ;
        columnExchange(aQ,i,j) ;
    }

//------------------------------------------------------------------------------

    template<typename T>
    void
    rowMultiplyOperation(Matrix< T > &aMat, Matrix< T > &aQ, Matrix< T > &aQ_, const uint i)
    {
        rowMultiply(aMat,i) ;
        rowMultiply(aQ_,i) ;
        columnMultiply(aQ,i) ;
    }

//------------------------------------------------------------------------------

    template<typename T>
    void
    rowAddOperation(Matrix< T > &aMat, Matrix< T > &aQ, Matrix< T > &aQ_, const uint i, const uint j, const int q)
    {
        rowAdd(aMat,i,j,q) ;
        rowAdd(aQ_,i,j,q) ;
        columnAdd(aQ,i,j,-q) ;
    }

//------------------------------------------------------------------------------

    template<typename T>
    void
    columnExchangeOperation(Matrix< T > &aMat, Matrix< T > &aR, Matrix< T > &aR_, const uint i, const uint j)
    {
        columnExchange(aMat,i,j) ;
        rowExchange(aR_,i,j) ;
        columnExchange(aR,i,j) ;
    }

//------------------------------------------------------------------------------

    template<typename T>
    void
    columnMultiplyOperation(Matrix< T > &aMat, Matrix< T > &aR, Matrix< T > &aR_, const uint i)
    {
        columnMultiply(aMat,i) ;
        rowMultiply(aR_,i) ;
        columnMultiply(aR,i) ;
    }

//------------------------------------------------------------------------------

    template<typename T>
    void
    columnAddOperation(Matrix< T > &aMat, Matrix< T > &aR, Matrix< T > &aR_, const uint i, const uint j, const int q)
    {
        columnAdd(aMat,i,j,q) ;
        rowAdd(aR_,i,j,-q) ;
        columnAdd(aR,i,j,q) ;
    }

//------------------------------------------------------------------------------

    template<typename T>
    void
    partRowReduce(Matrix< T > &aMat, Matrix< T > &aQ, Matrix< T > &aQ_, const uint k, const uint l)
    {
        const uint m = aMat.n_rows();
        for(uint i = k+1; i < m+1; ++i)
        {
            const int q = floor(aMat(i-1,l-1)/aMat(k-1,l-1));
            rowAddOperation(aMat,aQ,aQ_,i,k,-q);
        }
    }

//------------------------------------------------------------------------------

    template<typename T>
    void
    partColumnReduce(Matrix< T > &aMat, Matrix< T > &aR, Matrix< T > &aR_, const uint k, const uint l)
    {
        const uint n = aMat.n_cols();
        for(uint i = l+1; i < n+1; ++i)
        {
            const int q = floor(aMat(k-1,i-1)/aMat(k-1,l-1));
            columnAddOperation(aMat,aR,aR_,l,i,-q);
        }
    }

//------------------------------------------------------------------------------

    template<typename T>
    std::pair<uint, uint>
    smallestNonzero(Vector< T > &v, uint k)
    {
        uint alpha = abs(v(k-1));
        uint i0 = k;
        while (alpha == 0 && k < v.length())
        {
            k+=1;
            alpha = abs(v(k-1));
            i0 = k;
        }
        for(uint i = k; i < v.length(); ++i)
        {
            if (abs(v(i)) < alpha && abs(v(i)) != 0)
            {
                alpha = abs(v(i));
                i0 = i+1;
            }
        }
        return std::pair(alpha,i0);
    }

//------------------------------------------------------------------------------

    template<typename T>
    void
    rowPrepare(Matrix< T > &aMat, Matrix< T > &aQ, Matrix< T > &aQ_, const uint k, const uint l)
    {
        Vector< T > v = aMat.col(l-1);
        std::pair<uint, uint> tPair  = smallestNonzero(v, k);
        rowExchangeOperation(aMat, aQ, aQ_, k, tPair.second);
    }

//------------------------------------------------------------------------------

    template<typename T>
    void
    rowReduce(Matrix< T > &aMat, Matrix< T > &aQ, Matrix< T > &aQ_, const uint k, const uint l)
    {
        const uint m = aMat.n_rows();

        // todo: Optimize this to put directly in the while (cont) statement...
        bool cont = false;
        for (uint i = k; i < m; ++i)
        {
            if (aMat(i,l-1) != 0)
            {
                cont = true;
                break;
            }
        }
        while (cont)
        {
            rowPrepare(aMat, aQ, aQ_,k,l);
            partRowReduce(aMat, aQ, aQ_,k,l);

            // todo: Optimize this to put directly in the while (cont) statement...
            cont = false;
            for (uint i = k; i < m; ++i)
            {
                if (aMat(i,l-1) != 0)
                {
                    cont = true;
                    break;
                }
            }
        }
    }

//------------------------------------------------------------------------------

    template<typename T>
    std::tuple< Matrix< T >, Matrix< T >, uint >
    rowEchelon(Matrix< T > &aMat)
    {
        const uint m = aMat.n_rows();
        const uint n = aMat.n_cols();
        Matrix< T > tQ = Matrix< T >(m,m,0);
        Matrix< T > tQ_ = Matrix< T >(m,m,0);
        for(uint i = 0; i < m; ++i)
        {
            tQ(i,i) = 1;
            tQ_(i,i) = 1;
        }

        uint k = 0;
        uint l = 1;
        bool cont;

        while (k < m and n > 0)
        {
            // todo: Optimize this to put directly in the while (cont) statement...
            cont = true;
            for (uint i = k; i < m; ++i)
            {
                if (aMat(i,l-1) != 0)
                {
                    cont = false;
                    break;
                }
            }
            while ( cont )
            {
                l+=1;
                if (l == n+1)
                {
                    break;
                }
                // todo: Optimize this to put directly in the while (cont) statement...
                cont = true;
                for ( uint i = k; i < m; i++ )
                {
                    if ( aMat( i, l-1 ) != 0 )
                    {
                        cont = false;
                        break;
                    }
                }
            }
            if (l == n+1)
            {
                break;
            }
            k+=1;
            rowReduce(aMat,tQ,tQ_,k,l);
        }
        return std::tuple< Matrix< T >, Matrix< T >, uint > (tQ, tQ_, k);
    }

//------------------------------------------------------------------------------

    template<typename T>
    std::tuple< Matrix< T >, Matrix< T > >
    kernelImage(Matrix< T > &aMat)
    {
        const uint n = aMat.n_cols();
        Matrix< T > tMatT = trans(aMat);
        auto [tP, tP_, k] = rowEchelon(tMatT);
        tMatT = trans(tMatT);
        Matrix< T > tPT = trans(tP_);

        Matrix< T > tKer = Matrix< T >(tPT.n_rows(),n-k,0);

        uint tCount;
        tCount = 0;
        for(uint i = k; i < n; ++i)
        {
            for(uint j = 0; j < tPT.n_rows(); ++j)
            {
                tKer(j,tCount) = tPT(j,i);
            }
            tCount++;
        }

        Matrix< T > tIm = Matrix< T >(tMatT.n_rows(),k,0);
        for(uint i = 0; i < k; ++i)
        {
            for(uint j = 0; j < tMatT.n_rows(); ++j)
            {
                tIm(j,i) = tMatT(j,i);
            }
        }

        return std::tuple< Matrix< T >, Matrix< T > > (tKer, tIm);

    }

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
