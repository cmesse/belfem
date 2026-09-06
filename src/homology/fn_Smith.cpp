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

/* ----------------------------------------------------------------------------
 * Functions on integer matrices in order to obtain the Smith normal form
 * as described in Computational Homology from T. Kaczynski et al.
------------------------------------------------------------------------------- */

#include "fn_Smith.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    // Find the entry with the smallest non-zero value in
    // submatrix B[k:end,k:end]
    std::pair<uint, uint>
    minNonzero(Matrix< int > &aMat, const uint k)
    {
        const uint m = aMat.n_rows();
        Vector< int > tv = Vector< int >(m,0);
        Vector< int > tq = Vector< int >(m,0);

        for(uint i = 0; i < m; ++i)
        {
            if(i+1 >= k)
            {
                Vector< int > v = aMat.row(i);
                std::pair<uint, uint> tPair  = smallestNonzero(v, k);
                tv(i) = tPair.first;
                tq(i) = tPair.second;
            }
        }
        std::pair<uint, uint> tPair2 = smallestNonzero(tv, k);
        return std::pair<uint, uint>(tPair2.second, tq(tPair2.second-1));
    }

//------------------------------------------------------------------------------

    // Move the minimum non zero entry to position (k,k)
    void
    moveMinNonzero(Matrix< int > &aMat, Matrix< int > &aQ, Matrix< int > &aQ_, Matrix< int > &aR, Matrix< int > &aR_, const uint k)
    {
        auto [i,j] = minNonzero(aMat,k);
        rowExchangeOperation(aMat,aQ,aQ_,k,i);
        columnExchangeOperation(aMat,aR,aR_,k,j);
    }

//------------------------------------------------------------------------------

    // Check if B[k,k] divides the entries in submatrix B[k+1:end,k+1:end]
    std::tuple< bool, uint, uint, int >
    checkForDivisibility(Matrix< int > &aMat, const uint k)
    {
        const uint m = aMat.n_rows();
        const uint n = aMat.n_cols();
        int q;
        for(uint i = k; i < m; ++i)
        {
            for(uint j = k; j < n; ++j)
            {
                q = floor(aMat(i,j)/aMat(k-1,k-1));
                if (q*aMat(k-1,k-1) != aMat(i,j))
                {
                    return std::tuple< bool, uint, uint, int >(false, i+1, j+1, q);
                }
            }
        }
        return std::tuple< bool, uint, uint, int >(true, 0, 0, 0);
    }

//------------------------------------------------------------------------------

    // Partial smith normal form up to the kth entry
    void
    partSmithForm(Matrix< int > &aMat, Matrix< int > &aQ, Matrix< int > &aQ_, Matrix< int > &aR, Matrix< int > &aR_, const uint k)
    {
        const uint m = aMat.n_rows();
        const uint n = aMat.n_cols();
        bool divisible = false;
        bool cont;

        while (!divisible)
        {
            moveMinNonzero(aMat,aQ,aQ_,aR,aR_,k);
            partRowReduce(aMat,aQ,aQ_,k,k);

            // todo: Optimize this to put directly in the if (cont) statement...
            cont = false;
            for(uint i = k; i < m; ++i)
            {
                if (aMat(i,k-1) != 0)
                {
                    cont = true;
                    break;
                }
            }
            if (cont)
            {
                continue;
            }
            cont = false;
            partColumnReduce(aMat,aR,aR_,k,k);

            // todo: Optimize this to put directly in the if (cont) statement...
            for(uint i = k; i < n; ++i)
            {
                if (aMat(k-1,i) != 0)
                {
                    cont = true;
                    break;
                }
            }
            if (cont)
            {
                continue;
            }

            auto [div, i, j, q] = checkForDivisibility(aMat,k);
            divisible = div;

            if (!divisible)
            {
                rowAddOperation(aMat,aQ,aQ_,i,k,1);
                columnAddOperation(aMat,aR,aR_,k,j,-q);
            }

        }
    }

//------------------------------------------------------------------------------

    // Smith normal form algorithm
    std::tuple< Matrix< int >, Matrix< int >,Matrix< int >,Matrix< int >, uint, uint >
    smithForm(Matrix< int > &aMat)
    {
        const uint m = aMat.n_rows();
        const uint n = aMat.n_cols();

        Matrix< int > tQ = Matrix< int >(m,m,0);
        Matrix< int > tQ_ = Matrix< int >(m,m,0);
        for(uint i = 0; i < m; ++i)
        {
            tQ(i,i) = 1;
            tQ_(i,i) = 1;
        }

        Matrix< int > tR = Matrix< int >(n,n,0);
        Matrix< int > tR_ = Matrix< int >(n,n,0);
        for(uint i = 0; i < n; ++i)
        {
            tR(i,i) = 1;
            tR_(i,i) = 1;
        }

        uint s = 0;
        uint t = 0;

        // todo: Optimize this to put directly in the while (cont) statement...
        bool cont = false;
        for (uint i = t; i < m; ++i)
        {
            for (uint j = t; j < n; ++j)
            {
                if (aMat(i,j) != 0)
                {
                    cont = true;
                    break;
                }
            }
            if(cont)
            {
                break;
            }
        }

        while (cont)
        {
            t+=1;
            partSmithForm(aMat,tQ,tQ_,tR,tR_,t);
            if (aMat(t-1,t-1) < 0 )
            {
                rowMultiplyOperation(aMat,tQ,tQ_,t);
            }
            if (aMat(t-1,t-1) == 1)
            {
                s+=1;
            }

            // todo: Optimize this to put directly in the while (cont) statement...
            cont = false;
            for (uint i = t; i < m; ++i)
            {
                for (uint j = t; j < n; ++j)
                {
                    if (aMat(i,j) != 0)
                    {
                        cont = true;
                        break;
                    }
                }
                if(cont)
                {
                    break;
                }
            }
        }

        return std::tuple< Matrix< int >, Matrix< int >,Matrix< int >,Matrix< int >, uint, uint >(tQ, tQ_, tR, tR_, s, t);
    }

//------------------------------------------------------------------------------

    // Solve integer linear system Ax = b (if possible)
    Matrix< int >
    SolveInt(Matrix< int > aMat, Matrix < int > &aVec)
    {
        const uint n = aMat.n_cols();
        auto [tQ, tQ_, tR, tR_, s, t] = smithForm(aMat);
        tQ_*=aVec;
        Matrix< int > tu = Matrix< int >(n,1,0);

        for(uint i = 0; i < t ; ++i)
        {
            if (tQ_(i,0) % aMat(i,i) == 0)
            {
                tu(i,0) = tQ_(i,0)/aMat(i,i);
            }
            else
            {
                std::cout << "Failure" << std::endl;
                return Matrix< int >();
            }
        }
        for(uint i = t; i < n; ++i)
        {
            if(tQ_(i,0) == 0)
            {
                std::cout << "Failure" << std::endl;
                return Matrix< int >();
            }
            else
            {
                tu(i,0) = 0;
            }
        }
        return tR*=tu;
    }

//------------------------------------------------------------------------------
}