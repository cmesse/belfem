//
// Created by Christian Messe on 07.09.22.
//
#ifdef BELFEM_SUITESPARSE
#include <cholmod.h>
#endif
#include "fn_rcond.hpp"
#include "assert.hpp"

namespace belfem
{
    real
    rcond( SpMatrix & aMatrix )
    {
#ifdef OMP
#ifdef BELFEM_SUITESPARSE
        // make sure that matrix uses zero-based indexing
        aMatrix.set_indexing_base( SpMatrixIndexingBase::Cpp );

        // create a cholmod data object
        cholmod_sparse tMatrix ;
        tMatrix.nrow = aMatrix.n_rows();
        tMatrix.ncol = aMatrix.n_cols();
        tMatrix.nzmax = aMatrix.number_of_nonzeros();
        tMatrix.p = aMatrix.pointers();
        tMatrix.i = aMatrix.indices();
        tMatrix.x = aMatrix.data();
        tMatrix.stype = 0 ;
        tMatrix.itype = CHOLMOD_INT ;
        tMatrix.xtype = CHOLMOD_REAL ;
        tMatrix.dtype = CHOLMOD_DOUBLE ;
        tMatrix.sorted = 1;
        tMatrix.packed = 1;

        // the parameter object
        cholmod_common tCommon ;

        // launch cholmod
        cholmod_start (&tCommon) ;

        // create the LU-factorization
        cholmod_factor * tFactor = cholmod_analyze( &tMatrix, &tCommon) ;
        cholmod_factorize( &tMatrix, tFactor, &tCommon );

        // compute the approximate conditioning
        real aCond = cholmod_rcond( tFactor, &tCommon );

        // tidy up
        cholmod_free_factor ( &tFactor, &tCommon) ;
        cholmod_finish (&tCommon);

        // return the result
        return aCond ;
#else
        BELFEM_ERROR( false, "We are not linked against CHOLMOD." );
        return BELFEM_QUIET_NAN ;
#endif
#else
        BELFEM_ERROR( false, "We are not linked against OPENMP." );
        return BELFEM_QUIET_NAN ;
#endif
    }
}
#include "fn_rcond.hpp"