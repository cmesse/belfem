//
// Created by grgia@ge.polymtl.ca on 18/12/23.
//

/* ----------------------------------------------------------------------------
 * Functions on integer matrices in order to obtain the Smith normal form
 * as described in Computational Homology from T. Kaczynski et al.
------------------------------------------------------------------------------- */

#include "fn_Smith.hpp"


namespace belfem
{
//------------------------------------------------------------------------------

    // Exchange rows i and j
    void
    rowExchange(Matrix< int > &aMat, const uint i, const uint j)
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

    // Exchange columns i and j
    void
    columnExchange(Matrix< int > &aMat, const uint i, const uint j)
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

    // Multiply row i by -1
    void
    rowMultiply(Matrix< int > &aMat, const uint i)
    {
        const uint n = aMat.n_cols();
        for (uint k = 0; k < n; k++)
        {
            aMat(i-1,k)*=-1 ;
        }
    }

//------------------------------------------------------------------------------

    // Multiply column i by -1
    void
    columnMultiply(Matrix< int > &aMat, const uint i)
    {
        const uint m = aMat.n_rows();
        for (uint k = 0; k < m; k++)
        {
            aMat(k,i-1)*=-1 ;
        }
    }

//------------------------------------------------------------------------------

    // Add q times row j to row i
    void
    rowAdd(Matrix< int > &aMat, const uint i, const uint j, const int q)
    {
        const uint n = aMat.n_cols();
        for (uint k = 0; k < n; k++)
        {
            aMat(i-1,k) += q*aMat(j-1,k) ;
        }
    }

//------------------------------------------------------------------------------

    // Add q times column i to column j
    void
    columnAdd(Matrix< int > &aMat, const uint i, const uint j, const int q)
    {
        const uint m = aMat.n_rows();
        for (uint k = 0; k < m; k++)
        {
            aMat(k,j-1) += q*aMat(k,i-1) ;
        }
    }

//------------------------------------------------------------------------------

    // Exchange rows i and j keeping track of basis change
    void
    rowExchangeOperation(Matrix< int > &aMat, Matrix< int > &aQ, Matrix< int > &aQ_, const uint i, const uint j)
    {
        rowExchange(aMat,i,j) ;
        rowExchange(aQ_,i,j) ;
        columnExchange(aQ,i,j) ;
    }

//------------------------------------------------------------------------------

    // Multiply row i by -1 keeping track of basis change
    void
    rowMultiplyOperation(Matrix< int > &aMat, Matrix< int > &aQ, Matrix< int > &aQ_, const uint i)
    {
        rowMultiply(aMat,i) ;
        rowMultiply(aQ_,i) ;
        columnMultiply(aQ,i) ;
    }

//------------------------------------------------------------------------------

    // Add q times row j to row i keeping track of basis change
    void
    rowAddOperation(Matrix< int > &aMat, Matrix< int > &aQ, Matrix< int > &aQ_, const uint i, const uint j, const int q)
    {
        rowAdd(aMat,i,j,q) ;
        rowAdd(aQ_,i,j,q) ;
        columnAdd(aQ,i,j,-q) ;
    }

//------------------------------------------------------------------------------

    // Exchange columns i and j keeping track of basis change
    void
    columnExchangeOperation(Matrix< int > &aMat, Matrix< int > &aR, Matrix< int > &aR_, const uint i, const uint j)
    {
        columnExchange(aMat,i,j) ;
        rowExchange(aR_,i,j) ;
        columnExchange(aR,i,j) ;
    }

//------------------------------------------------------------------------------

    // Multiply column i by -1 keeping track of basis change
    void
    columnMultiplyOperation(Matrix< int > &aMat, Matrix< int > &aR, Matrix< int > &aR_, const uint i)
    {
        columnMultiply(aMat,i) ;
        rowMultiply(aR_,i) ;
        columnMultiply(aR,i) ;
    }

//------------------------------------------------------------------------------

    // Add q times column i to column i keeping track of basis change
    void
    columnAddOperation(Matrix< int > &aMat, Matrix< int > &aR, Matrix< int > &aR_, const uint i, const uint j, const int q)
    {
        columnAdd(aMat,i,j,q) ;
        rowAdd(aR_,i,j,-q) ;
        columnAdd(aR,i,j,q) ;
    }

//------------------------------------------------------------------------------

    // Partial row reduction
    void
    partRowReduce(Matrix< int > &aMat, Matrix< int > &aQ, Matrix< int > &aQ_, const uint k, const uint l)
    {
        const uint m = aMat.n_rows();
        for(uint i = k+1; i < m+1; i++)
        {
            const int q = floor(aMat(i-1,l-1)/aMat(k-1,l-1));
            rowAddOperation(aMat,aQ,aQ_,i,k,-q);
        }
    }

//------------------------------------------------------------------------------

    // Partial column reduction
    void
    partColumnReduce(Matrix< int > &aMat, Matrix< int > &aR, Matrix< int > &aR_, const uint k, const uint l)
    {
        const uint n = aMat.n_cols();
        for(uint i = l+1; i < n+1; i++)
        {
            const int q = floor(aMat(k-1,i-1)/aMat(k-1,l-1));
            columnAddOperation(aMat,aR,aR_,l,i,-q);
        }
    }

//------------------------------------------------------------------------------

    //Find the smallest non zero entry of a vector
    //starting from position k
    std::pair<uint, uint>
    smallestNonzero(Vector< int > &v, uint k)
    {
        uint alpha = abs(v(k-1));
        uint i0 = k;
        while (alpha == 0 && k < v.length())
        {
            k+=1;
            alpha = abs(v(k-1));
            i0 = k;
        }
        for(uint i = k; i < v.length(); i++)
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

    // Prepare the rows by ordering from smallest non zero
    void
    rowPrepare(Matrix< int > &aMat, Matrix< int > &aQ, Matrix< int > &aQ_, const uint k, const uint l)
    {
        Vector< int > v = aMat.col(l-1);
        std::pair<uint, uint> tPair  = smallestNonzero(v, k);
        rowExchangeOperation(aMat, aQ, aQ_, k, tPair.second);
    }

//------------------------------------------------------------------------------

    // Compute the matrix that satisfies the (k,l) criterion of row echelon
    void
    rowReduce(Matrix< int > &aMat, Matrix< int > &aQ, Matrix< int > &aQ_, const uint k, const uint l)
    {
        const uint m = aMat.n_rows();

        // todo: Optimize this to put directly in the while (cont) statement...
        bool cont = false;
        for (uint i = k; i < m; i++)
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
            for (uint i = k; i < m; i++)
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

    // Compute the row echelon form of a Matrix
    std::tuple< Matrix< int >, Matrix< int >, uint >
    rowEchelon(Matrix< int > &aMat)
    {
        const uint m = aMat.n_rows();
        const uint n = aMat.n_cols();
        Matrix< int > tQ = Matrix< int >(m,m,0);
        Matrix< int > tQ_ = Matrix< int >(m,m,0);
        for(uint i = 0; i < m; i++)
        {
            tQ(i,i) = 1;
            tQ_(i,i) = 1;
        }

        uint k = 0;
        uint l = 1;
        bool cont;

        while (k < m)
        {
            // todo: Optimize this to put directly in the while (cont) statement...
            cont = true;
            for (uint i = k; i < m; i++)
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
        return std::tuple< Matrix< int >, Matrix< int >, uint > (tQ, tQ_, k);
    }

//------------------------------------------------------------------------------

    // Transpose of a matrix
    Matrix < int >
    transpose(Matrix< int > &aMat)
    {
        uint m = aMat.n_rows();
        uint n = aMat.n_cols();
        Matrix < int > tMatT = Matrix < int >(n,m,0);
        for(uint i = 0; i < n; i++)
        {
            for(uint j = 0; j < m; j++)
            {
                tMatT(i,j) = aMat(j,i);
            }
        }
        return tMatT;
    }

//------------------------------------------------------------------------------

    // Find the kernel and the image of the columns of a matrix
    std::tuple< Matrix< int >, Matrix< int > >
    kernelImage(Matrix< int > &aMat)
    {
        const uint n = aMat.n_cols();
        Matrix< int > tMatT = transpose(aMat);
        auto [tP, tP_, k] = rowEchelon(tMatT);
        tMatT = transpose(tMatT);
        Matrix< int > tPT = transpose(tP_);

        Matrix< int > tKer = Matrix< int >(tPT.n_rows(),n-k,0);

        uint tCount;
        tCount = 0;
        for(uint i = k; i < n; i++)
        {
            for(uint j = 0; j < tPT.n_rows(); j++)
            {
                tKer(j,tCount) = tPT(j,i);
            }
            tCount++;
        }

        Matrix< int > tIm = Matrix< int >(tMatT.n_rows(),k,0);
        for(uint i = 0; i < k; i++)
        {
            for(uint j = 0; j < tMatT.n_rows(); j++)
            {
                tIm(j,i) = tMatT(j,i);
            }
        }

        return std::tuple< Matrix< int >, Matrix< int > > (tKer, tIm);

    }

//------------------------------------------------------------------------------

    // Find the entry with the smallest non-zero value in
    // submatrix B[k:end,k:end]
    std::pair<uint, uint>
    minNonzero(Matrix< int > &aMat, const uint k)
    {
        const uint m = aMat.n_rows();
        Vector< int > tv = Vector< int >(m,0);
        Vector< int > tq = Vector< int >(m,0);

        for(uint i = 0; i < m; i++)
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
        for(uint i = k; i < m; i++)
        {
            for(uint j = k; j < n; j++)
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
            for(uint i = k; i < m; i++)
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
            for(uint i = k; i < n; i++)
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
        for(uint i = 0; i < m; i++)
        {
            tQ(i,i) = 1;
            tQ_(i,i) = 1;
        }

        Matrix< int > tR = Matrix< int >(n,n,0);
        Matrix< int > tR_ = Matrix< int >(n,n,0);
        for(uint i = 0; i < n; i++)
        {
            tR(i,i) = 1;
            tR_(i,i) = 1;
        }

        uint s = 0;
        uint t = 0;

        // todo: Optimize this to put directly in the while (cont) statement...
        bool cont = false;
        for (uint i = t; i < m; i++)
        {
            for (uint j = t; j < n; j++)
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
            for (uint i = t; i < m; i++)
            {
                for (uint j = t; j < n; j++)
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

        for(uint i = 0; i < t ; i++)
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
        for(uint i = t; i < n; i++)
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