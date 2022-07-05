//
// Created by Christian Messe on 11.07.20.
//

#ifndef BELFEM_PARDISOTOOLS_HPP
#define BELFEM_PARDISOTOOLS_HPP

#ifdef __cplusplus
extern"C" {
#endif

    int
    pardisotools_initialize_parameters( const int * aParameters );

    int
    pardisotools_symbolic_factorization(
            const int    & aN,
            const int    & aNNZ,
            const int    & aNRHS,
            const int    * aPointers,
            const int    * aIndices,
            const double * aValues );

    int
    pardisotools_solve(  const int    &  aN,
                         const int    & aNNZ,
                         const int    & aNRHS,
                         const int    * aPointers,
                         const int    * aIndices,
                         const double * aValues,
                         double       * aLHS,
                         const double * aRHS,
                         int          * aInfo   );

    int
    pardisotools_free() ;

    double
    pardisotools_get_determinant() ;

#ifdef __cplusplus
}
#endif

#endif //BELFEM_PARDISOTOOLS_HPP
