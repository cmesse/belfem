//
// Created by Christian Messe on 25.11.20.
//
#include "fn_Mesh_ratio.hpp"
#include "assert.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

    void
    ratio_ar2(
            const real     & aDeltaX0,
            const real     & aLength,
            const index_t  & aNumCells,
            Vector< real > & aX )
    {
        real tR0 = 0.1 ;
        real tF0 = _check_ratio( aDeltaX0, aLength, aNumCells, tR0 );

        real tR1 = 10.0 ;
        real tF1 = _check_ratio( aDeltaX0, aLength, aNumCells, tR1 );

        real tRatio = BELFEM_QUIET_NAN ;
        real tF = BELFEM_REAL_MAX ;

        index_t tCount = 0;

        while ( abs( tF ) > 1e-12 )
        {
            tRatio = tR0 - tF0 * ( tR1 - tR0 ) / ( tF1 - tF0 );
            tF = _check_ratio( aDeltaX0, aLength, aNumCells, tRatio );

            if ( tF * tF0  > 0.0 )
            {
                tR0 = tRatio ;
                tF0 = tF ;
            }
            else
            {
                tR1 = tRatio ;
                tF1 = tF ;
            }
            BELFEM_ERROR( tCount++ < 100, "Too many iterations" );
        }

        aX.set_size( 2 * aNumCells + 1 );

        aX( 0 ) = 0.0 ;
        tCount = 0 ;

        real tDx = aDeltaX0 ;

        for ( index_t k=0; k< aNumCells; ++k )
        {
            aX( tCount + 1 ) = aX( tCount ) + 0.5 * tDx ;
            aX( tCount + 2 ) = aX( tCount ) + tDx ;
            tDx *= tRatio ;
            tCount+=2 ;
        }

        aX( 2 * aNumCells ) = aLength ;

    }

//------------------------------------------------------------------------------

    void
    ratio_adx2( const real     & aRatio,
                const real     & aLength,
                const index_t  & aNumCells,
                Vector< real > & aX )
    {
        real tDx = 1.0 ;
        real tX = 0.0 ;

        for( index_t k=0; k< aNumCells; ++k )
        {
            tX += tDx ;
            tDx *= aRatio ;
        }

        tDx = aLength / tX ;

        aX.set_size( 2 * aNumCells + 1 );

        aX( 0 ) = 0.0 ;
        index_t tCount = 0 ;

        for ( index_t k=0; k< aNumCells; ++k )
        {
            aX( tCount + 1 ) = aX( tCount ) + 0.5 * tDx ;
            aX( tCount + 2 ) = aX( tCount ) + tDx ;
            tDx *= aRatio ;
            tCount+=2 ;
        }

        aX( 2 * aNumCells ) = aLength ;

    }

//------------------------------------------------------------------------------

    // test function for ratio_ar2
    real
    _check_ratio(  const real     & aDeltaX0,
                   const real     & aLength,
                   const index_t  & aNumCells,
                   const real     & aRatio )
    {
        real tDx = aDeltaX0 ;
        real tX = 0.0 ;
        for ( index_t k=0; k<aNumCells; ++k )
        {
            tX+=tDx ;
            tDx *= aRatio ;
        }
        return tX - aLength ;
    }

//------------------------------------------------------------------------------

}