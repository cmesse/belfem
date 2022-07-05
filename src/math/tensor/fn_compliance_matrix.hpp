//
// Created by Christian Messe on 23.12.20.
//

#ifndef BELFEM_FN_COMPLIANCE_MATRIX_HPP
#define BELFEM_FN_COMPLIANCE_MATRIX_HPP

namespace belfem
{
    template< typename T >
    inline void
    compliance_matrix( const T & aE1,
                       const T & aE2,
                       const T & aE3,
                       const T & aNu23,
                       const T & aNu13,
                       const T & aNu12,
                       const T & aG23,
                       const T & aG31,
                       const T & aG12,
                       Matrix <T> & aS )
    {
        BELFEM_ASSERT( aS.n_rows() == 6
                      && aS.n_cols() == 6,
                      "target matrix must be allocated as 6x6" );

        aS.fill( 0.0 );

        aS( 0, 0 ) = 1.0    / aE1 ;
        aS( 1, 0 ) = -aNu12 / aE2 ;
        aS( 2, 0 ) = -aNu13 / aE3 ;

        aS( 0, 1 ) = aS( 1, 0 ) ;
        aS( 1, 1 ) = 1.0 / aE2 ;
        aS( 2, 1 ) = -aNu23 / aE2 ;

        aS( 0, 2 ) = aS( 2, 0 );
        aS( 1, 2 ) = aS( 2, 1 ) ;
        aS( 2, 2 ) = 1.0 / aE3 ;

        aS( 3, 3 ) = 1.0 / aG23 ;
        aS( 4, 4 ) = 1.0 / aG31 ;
        aS( 5, 5 ) = 1.0 / aG12 ;
    }
}
#endif //BELFEM_FN_COMPLIANCE_MATRIX_HPP
