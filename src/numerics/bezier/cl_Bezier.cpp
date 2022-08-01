//
// Created by Christian Messe on 30.11.20.
//

#include "cl_Bezier.hpp"
#include "assert.h"
#include "fn_dot.hpp"
#include "intpoints.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

    Bezier::Bezier(
            const real aX0, const real aY0, const real adYdx0,
            const real aX1, const real aY1, const real adYdx1,
            const BezierType aType )
    {
        // allocate memory
        mX.set_size( 4 );
        mY.set_size( 4 );
        mWork.set_size( 4 );

        // compute the coordiantes for the basis
        switch( aType )
        {
            case( BezierType::Horizontal ) :
            {
                this->compute_basis_xwise( aX0, aY0, adYdx0,
                                           aX1, aY1, adYdx1 );
                break ;
            }
            case( BezierType::Vertical ) :
            {
                this->compute_basis_ywise( aX0, aY0, adYdx0,
                                           aX1, aY1, adYdx1 );
                break ;
            }
            default:
            {
                BELFEM_ERROR( false, "Invalid Bezier Type" );
            }
        }
    }

//------------------------------------------------------------------------------

    real
    Bezier::x_by_xi( const real aXi )
    {
        // compute the shape function
        this->compute_N( aXi );

        // return the interpolated value
        return dot( mWork, mX );
    }

//------------------------------------------------------------------------------

    real
    Bezier::y_by_xi( const real aXi )
    {
        // compute the shape function
        this->compute_N( aXi );

        // return the interpolated value
        return dot( mWork, mY );
    }

//------------------------------------------------------------------------------

    real
    Bezier::xi_by_x( const real aX )
    {
        real tXi0 = -1.0;
        real tF0 = this->x_by_xi( tXi0 ) - aX ;

        real tXi1 =  1.0;
        //real tF1 = this->x_by_xi( tXi1 ) - aX ;

        real aXi = BELFEM_QUIET_NAN ;
        real tF = BELFEM_REAL_MAX ;

        uint tCount = 0 ;

        // start regula falsi
        while( std::abs( tF ) > 1e-12 )
        {
            // aXi = tXi0 - tF0 * ( tXi1 - tXi0 ) / ( tF1 - tF0 );
            aXi = 0.5 * ( tXi0 + tXi1 );

            tF = this->x_by_xi( aXi ) - aX ;

            if ( tF * tF0 > 0.0 )
            {
                tXi0 = aXi ;
                tF = tF0 ;
            }
            else
            {
                tXi1 = aXi ;
                //tF = tF1 ;
            }

            BELFEM_ERROR( tCount++ < 100, "Too many iterations.");
        }

        return aXi ;
    }

//------------------------------------------------------------------------------

    real
    Bezier::xi_by_y( const real aY )
    {
        real tXi0 = -1.0;
        real tF0 = this->y_by_xi( tXi0 ) - aY ;

        real tXi1 =  1.0;
        //real tF1 = this->y_by_xi( tXi1 ) - aY ;

        real aXi = BELFEM_QUIET_NAN ;
        real tF = BELFEM_REAL_MAX ;

        uint tCount = 0 ;

        // start regula falsi
        while( std::abs( tF ) > 1e-12 )
        {
            //aXi = tXi0 - tF0 * ( tXi1 - tXi0 ) / ( tF1 - tF0 );
            aXi = 0.5 * ( tXi1 - tXi0 );
            tF = this->y_by_xi( aXi ) - aY ;

            if ( tF * tF0 > 0.0 )
            {
                tXi0 = aXi ;
                tF = tF0 ;
            }
            else
            {
                tXi1 = aXi ;
                // tF = tF1 ;
            }

            BELFEM_ERROR( tCount++ < 100, "Too many iterations.");
        }

        return aXi ;
    }

//------------------------------------------------------------------------------

    real
    Bezier::x( const real aY )
    {
        return this->x_by_xi( this->xi_by_y( aY ) );
    }

//------------------------------------------------------------------------------

    real
    Bezier::y( const real aX )
    {
        return this->y_by_xi( this->xi_by_x( aX ) );
    }

//------------------------------------------------------------------------------

    void
    Bezier::point( const real aXi, real & aX, real & aY )
    {
        this->compute_N( aXi );

        aX = dot( mWork, mX );
        aY = dot( mWork, mY );
    }

//------------------------------------------------------------------------------

    void
    Bezier::dpoint( const real aXi, real & adXdXi, real & adYdXi )
    {
        this->compute_dNdXi( aXi );

        adXdXi = dot( mWork, mX );
        adYdXi = dot( mWork, mY );
    }

//------------------------------------------------------------------------------

    void
    Bezier::ddpoint( const real aXi, real & ad2XdXi2, real & ad2YdXi2 )
    {
        this->compute_d2NdXi2( aXi );

        ad2XdXi2 = dot( mWork, mX );
        ad2YdXi2 = dot( mWork, mY );
    }

//------------------------------------------------------------------------------

    real
    Bezier::dydx( const real aX )
    {
        real tdXdXi ;
        real tdYdXi ;
        this->dpoint( this->xi_by_x( aX ), tdXdXi, tdYdXi ) ;

        return tdYdXi / tdYdXi ;
    }

//------------------------------------------------------------------------------

    real
    Bezier::dxdy( const real aY )
    {
        real tdXdXi ;
        real tdYdXi ;
        this->dpoint( this->xi_by_y( aY ), tdXdXi, tdYdXi ) ;

        return tdXdXi / tdXdXi ;
    }

//------------------------------------------------------------------------------

    real
    Bezier::d2ydx2( const real aX )
    {
        real tXi = this->xi_by_x( aX );
        real tdXdXi ;
        real tdYdXi ;
        real td2XdXi2 ;
        real td2YdXi2 ;

        this->dpoint( tXi, tdXdXi, tdYdXi ) ;
        this->ddpoint( tXi, td2XdXi2, td2YdXi2 );

        return td2XdXi2 * tdYdXi * tdYdXi + tdXdXi * td2YdXi2 ;
    }

//------------------------------------------------------------------------------

    real
    Bezier::d2xdy2( const real aY )
    {
        real tXi = this->xi_by_y( aY );
        real tdXdXi ;
        real tdYdXi ;
        real td2XdXi2 ;
        real td2YdXi2 ;

        this->dpoint( tXi, tdXdXi, tdYdXi ) ;
        this->ddpoint( tXi, td2XdXi2, td2YdXi2 );

        return td2YdXi2 * tdXdXi * tdXdXi + tdYdXi * td2XdXi2 ;
    }

//------------------------------------------------------------------------------

    void
    Bezier::compute_basis_xwise(
            const real aX0, const real aY0,const real adYdx0,
            const real aX1, const real aY1, const real adYdx1 )
    {
        real tDX = ( aX1 - aX0 ) / 3. ;

        mX( 0 ) = aX0 ;
        mX( 1 ) = aX0 + tDX ;
        mX( 2 ) = aX1 - tDX ;
        mX( 3 ) = aX1 ;

        mY( 0 ) = aY0 ;
        mY( 1 ) = aY0 + tDX * adYdx0 ;
        mY( 2 ) = aY1 - tDX * adYdx1 ;
        mY( 3 ) = aY1 ;
    }

//------------------------------------------------------------------------------

    void
    Bezier::compute_basis_ywise(
            const real aX0, const real aY0,const real adYdx0,
            const real aX1, const real aY1, const real adYdx1 )
    {
        real tDY = ( aY1 - aY0 ) / 3. ;

        mX( 0 ) = aX0 ;
        mX( 1 ) = aX0 - tDY * adYdx0 ;
        mX( 2 ) = aX1 + tDY * adYdx1 ;
        mX( 3 ) = aX1 ;

        mY( 0 ) = aY0 ;
        mY( 1 ) = aY0 + tDY ;
        mY( 2 ) = aY1 - tDY ;
        mY( 3 ) = aY1 ;
    }

//------------------------------------------------------------------------------

    void
    Bezier::compute_N( const real aXi )
    {
        mWork( 0 ) = std::pow(1 - aXi, 3);
        mWork( 1 ) = 3. * ((aXi * (aXi - 1.) - 1.) * aXi + 1.);
        mWork( 2 ) = 3. * ((1. - (aXi + 1.) * aXi) * aXi + 1.);
        mWork( 3 ) = ((aXi * (3. + aXi) + 3.) * aXi + 1.);

        mWork *= 0.125;
    }

//------------------------------------------------------------------------------

    void
    Bezier::compute_dNdXi( const real aXi )
    {
        mWork( 0 ) = - (1. - aXi) * (1. - aXi);
        mWork( 1 ) =  ((3. * aXi - 2.) * aXi - 1.);
        mWork( 2 ) =    1. - aXi * (2. + 3. * aXi);
        mWork( 3 ) =   (1. + aXi) * (1. + aXi);

        mWork *= 0.375;
    }

//------------------------------------------------------------------------------

    void
    Bezier::compute_d2NdXi2( const real aXi )
    {
        mWork( 0 ) =  1. - aXi;
        mWork( 1 ) =  3. * aXi - 1.;
        mWork( 2 ) = -3. * aXi - 1.;
        mWork( 3 ) =  aXi + 1.;

        mWork *= 0.75;
    }

//------------------------------------------------------------------------------

    /**
     * computes the length of the curve
     */
    real
    Bezier::compute_length( const uint aNumIntegrationPoints )
    {
        // number of points
        int tN = aNumIntegrationPoints ;

        // integration weights
        Vector< double > tW( aNumIntegrationPoints );

        // integration points
        Vector< double > tXi( aNumIntegrationPoints );

        // compute the integration points
        intpoints_gauss(
                &tN,
                tW.data(),
                tXi.data() );

        return this->compute_length( tW, tXi );
    }

//------------------------------------------------------------------------------

    /**
     * computes the length of the curve, but provide points and weights
     */
    real
    Bezier::compute_length(
            const Vector< double > & aW,
            const Vector< double > & aXi )
    {
        // check input
        BELFEM_ASSERT( aW.length() == aXi.length(),
                      "Length of weights and points does not match ( %u vs %u )",
                      ( unsigned int ) aW.length(),
                      ( unsigned int ) aXi.length() );

        uint tN = aW.length() ;

        // integrate
        real aLength = 0.0 ;
        real tdXdXi ;
        real tdYdXi ;
        for( uint k=0; k<tN; ++k )
        {
            this->dpoint( aXi( k ), tdXdXi, tdYdXi ) ;

            aLength += std::sqrt( tdXdXi * tdXdXi + tdYdXi * tdYdXi );
        }

        return aLength ;
    }
//------------------------------------------------------------------------------
}