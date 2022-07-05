//
// Created by Christian Messe on 2019-01-22.
//

#ifndef SSF_FSPBLAS_HPP
#define SSF_FSPBLAS_HPP

/**
 * C = alpha*A*B + beta*C
 */
#ifdef BELFEM_NETLIB
#ifdef __cplusplus
extern"C" {
#endif
    void
    dcscmm_(
            const int    * aTransposedFlag,
            const int    * aRowsA,
            const int    * aColsC,
            const int    * aColsA,
            const double * aAlpha,
            const int 	 * aParameters,
            const double * aValues,
            const int    * aIndices,
            const int    * aPointersBegin,
            const int    * aPointersEnd,
            const double * aB,
            int          * aRowsB,
            const double * aBeta,
            double       * aC,
            const int    * aRowsC,
            double       * aSwap,
            const int    * aRowsSwap
    );

    void
    dcsrmm_(
            const int    * aTransposedFlag,
            const int    * aRowsA,
            const int    * aColsC,
            const int    * aColsA,
            const double * aAlpha,
            const int 	 * aParameters,
            const double * aValues,
            const int    * aIndices,
            const int    * aPointersBegin,
            const int    * aPointersEnd,
            const double * aB,
            int          * aRowsB,
            const double * aBeta,
            double       * aC,
            const int    * aRowsC,
            double       * aSwap,
            const int    * aRowsSwap
    );

#ifdef __cplusplus
}
#endif
#elif BELFEM_MKL
#ifdef __cplusplus
extern"C" {
#endif
    void
    mkl_dcscmm_(
            const char   * aTransposedFlag,
            const int    * aRowsA,
            const int    * aColsC,
            const int    * aColsA,
            const double * aAlpha,
            const char   * aParameters,
            const double * aValues,
            const int    * aIndices,
            const int    * aPointersBegin,
            const int    * aPointersEnd,
            const double * aB,
            int          * aRowsB,
            const double * aBeta,
            double       * aC,
            const int    * aRowsC
            );

    void
    mkl_dcsrmm_(
            const char   * aTransposedFlag,
            const int    * aRowsA,
            const int    * aColsC,
            const int    * aColsA,
            const double * aAlpha,
            const char   * aParameters,
            const double * aValues,
            const int    * aIndices,
            const int    * aPointersBegin,
            const int    * aPointersEnd,
            const double * aB,
            int          * aRowsB,
            const double * aBeta,
            double       * aC,
            const int    * aRowsC
    );
#ifdef __cplusplus
}
#endif
#endif

//#ifdef BELFEM_ACCELLERATE
//#include "/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Headers/Sparse/BLAS.h"
//#endif

#endif //SSF_FSPBLAS_HPP
