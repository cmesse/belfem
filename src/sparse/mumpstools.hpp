//
// Created by Christian Messe on 10.07.20.
//

#ifndef BELFEM_MUMPSTOOLS_HPP
#define BELFEM_MUMPSTOOLS_HPP

#ifdef __cplusplus
extern"C" {
#endif

//------------------------------------------------------------------------------

    void
    mumpstools_solve(
        const int    * aIParameters,
        const double * aRParameters,
        const int    & aN,
        const int    & aNNZ,
        const int    & aNHRS,
        const int    * aRowIndices,
        const int    * aColIndices,
        const double * aValues,
        double       * aX,
        const double * aY,
        int          * aInfo ) ;

//------------------------------------------------------------------------------

    double
    mumpstools_get_determinant() ;

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif

#endif //BELFEM_MUMPSTOOLS_HPP
