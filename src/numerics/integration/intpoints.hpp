//
// Created by Christian Messe on 04.10.19.
//

#ifndef BELFEM_INTPOINTS_HPP
#define BELFEM_INTPOINTS_HPP

#ifdef __cplusplus
extern"C"
{
#endif
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    intpoints_lobatto( int * n, double * w, double * xi );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    void
    intpoints_gauss( int * n, double * w, double * xi );

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef __cplusplus
}
#endif

#endif //BELFEM_INTPOINTS_HPP
