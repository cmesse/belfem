//
// Created by Christian Messe on 13.11.19.
//

#include "cl_OrthotropicMaterial.hpp"

namespace belfem
{
//----------------------------------------------------------------------------

    OrthotropicMaterial::OrthotropicMaterial( const MaterialType aType ) :
            Material( aType )
    {

    }

//----------------------------------------------------------------------------

    void
    OrthotropicMaterial::C( Matrix< real > & aC, const real aT ) const
    {

        this->compute_C3D( this->E1( aT ),
                           this->E2( aT ),
                           this->E3( aT ),
                           this->nu23( aT ),
                           this->nu12( aT ),
                           this->nu31( aT ),
                           this->G23( aT ),
                           this->G12( aT ),
                           this->G31( aT ),
                           aC );
    }

//----------------------------------------------------------------------------


    void
    OrthotropicMaterial::C_ps( Matrix< real > & aC, const real aT ) const
    {
        this->compute_C_ps(
                this->E1( aT ),
                this->E2( aT ),
                this->nu12( aT ),
                this->G12( aT ),
                aC );
    }

//----------------------------------------------------------------------------

    void
    OrthotropicMaterial::C_rot( Matrix< real > & aC, const real aT ) const
    {
        this->compute_C_rot(
                this->E1( aT ),
                this->E2( aT ),
                this->E3( aT ),
                this->nu23( aT ),
                this->nu12( aT ),
                this->nu31( aT ),
                this->G12( aT ),
                aC );
    }

//----------------------------------------------------------------------------

    void
    OrthotropicMaterial::compute_C3D(
            const real & aE1,
            const real & aE2,
            const real & aE3,
            const real & aNu23,
            const real & aNu12,
            const real & aNu31,
            const real & aG23,
            const real & aG12,
            const real & aG31,
            Matrix< real > & aC
    ) const
    {
        BELFEM_ASSERT( aC.n_cols() == 6 && aC.n_rows() == 6,
                      "Elasticity matrix must be allocated as 6x6" );

        // help parameter
        real tX = 1.0 / ( aE2*aE2*aE3*aNu12*aNu12 + aE1*aE1*aE2*aNu31*aNu31
                          + aE1*aE3*(aE3*aNu23*aNu23 + aE2*(2.0*aNu12*aNu23*aNu31 - 1.0)) );

        aC( 0, 0 ) =  aE1*aE1*aE3*(aE3*aNu23*aNu23-aE2)*tX;
        aC( 1, 0 ) = -aE1*aE2*aE3*(aE2*aNu12+aE1*aNu23*aNu31)*tX;
        aC( 2, 0 ) = -aE1*aE2*aE3*(aE1*aNu31+aE3*aNu12*aNu23)*tX;
        aC( 3, 0 ) = 0.0;
        aC( 4, 0 ) = 0.0;
        aC( 5, 0 ) = 0.0;

        aC( 0, 1 ) = aC( 1, 0 );
        aC( 1, 1 ) = aE1*aE2*aE2*(aE1*aNu31*aNu31-aE3)*tX;
        aC( 2, 1 ) = -aE1*aE2*aE3*(aE3*aNu23+aE2*aNu12*aNu31)*tX;
        aC( 3, 1 ) = 0.0;
        aC( 4, 1 ) = 0.0;
        aC( 5, 1 ) = 0.0;

        aC( 0, 2 ) = aC( 2, 0 );
        aC( 1, 2 ) = aC( 2, 1 );
        aC( 2, 2 ) = aE2*aE3*aE3*(aE2*aNu12*aNu12-aE1)*tX;
        aC( 3, 2 ) = 0.0;
        aC( 4, 2 ) = 0.0;
        aC( 5, 2 ) = 0.0;

        aC( 0, 3 ) = 0.0;
        aC( 1, 3 ) = 0.0;
        aC( 2, 3 ) = 0.0;
        aC( 3, 3 ) = aG23;
        aC( 4, 3 ) = 0.0;
        aC( 5, 3 ) = 0.0;

        aC( 0, 4 ) = 0.0;
        aC( 1, 4 ) = 0.0;
        aC( 2, 4 ) = 0.0;
        aC( 3, 4 ) = 0.0;
        aC( 4, 4 ) = aG31;
        aC( 5, 4 ) = 0.0;

        aC( 0, 5 ) = 0.0;
        aC( 1, 5 ) = 0.0;
        aC( 2, 5 ) = 0.0;
        aC( 3, 5 ) = 0.0;
        aC( 4, 5 ) = 0.0;
        aC( 5, 5 ) = aG12;
    }

//----------------------------------------------------------------------------
    void
    OrthotropicMaterial::compute_C_ps(
            const real & aE1,
            const real & aE2,
            const real & aNu12,
            const real & aG12,
            Matrix< real > & aC
    ) const
    {
        BELFEM_ASSERT( aC.n_cols() == 3 && aC.n_rows() == 3,
                      "Elasticity matrix must be allocated as 3x3" );

        real tX = aE1/(aE1-aE2*aNu12*aNu12);
        aC( 0, 0 ) = aE1*tX;
        aC( 1, 0 ) = aNu12*aE2*tX;
        aC( 2, 0 ) = 0.0;

        aC( 0, 1 ) = aC( 1, 0 );
        aC( 1, 1 ) = aE2*tX;
        aC( 2, 1 ) = 0.0;

        aC( 0, 2 ) = 0.0;
        aC( 1, 2 ) = 0.0;
        aC( 2, 2 ) = aG12;
    }

//----------------------------------------------------------------------------

    void
    OrthotropicMaterial::compute_C_rot(
            const real & aE1,
            const real & aE2,
            const real & aE3,
            const real & aNu23,
            const real & aNu12,
            const real & aNu31,
            const real & aG12,
            Matrix< real > & aC
    ) const
    {
        BELFEM_ASSERT( aC.n_cols() == 4 && aC.n_rows() == 4,
                      "Elasticity matrix must be allocated as 4x4" );

        real tX = 1.0/(aE2*aE2*aE3*aNu12*aNu12+aE1*aE1*aE2*aNu31*aNu31
                       +aE1*aE3*(aE3*aNu23*aNu23+aE2*(2.0*aNu12*aNu23*aNu31-1.0)));

        aC( 0, 0 ) =  aE1*aE1*aE3*(aE3*aNu23*aNu23-aE2)*tX;
        aC( 1, 0 ) = -aE1*aE2*aE3*(aE2*aNu12+aE1*aNu23*aNu31)*tX;
        aC( 2, 0 ) = -aE1*aE2*aE3*(aE1*aNu31+aE3*aNu12*aNu23)*tX;
        aC( 3, 0 ) =  0.0;

        aC( 0, 1 ) =  aC( 1, 0 );
        aC( 1, 1 ) =  aE1*aE2*aE2*(aE1*aNu31*aNu31-aE3)*tX;
        aC( 2, 1 ) = -aE1*aE2*aE3*(aE3*aNu23+aE2*aNu12*aNu31)*tX;
        aC( 3, 1 ) =  0.0;

        aC( 0, 2 ) =  aC( 2, 0 );
        aC( 1, 2 ) =  aC( 2, 1 );
        aC( 2, 2 ) =  aE2*aE3*aE3*(aE2*aNu12*aNu12-aE1)*tX;
        aC( 3, 2 ) =  0.0;

        aC( 0, 3 ) = 0.0;
        aC( 1, 3 ) = 0.0;
        aC( 2, 3 ) = 0.0;
        aC( 3, 3 ) = aG12;
    }

//----------------------------------------------------------------------------

    real
    OrthotropicMaterial::E1(  const real aT ) const
    {
        BELFEM_ERROR( false, "Function E1 not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

//----------------------------------------------------------------------------

    real
    OrthotropicMaterial::E2(  const real aT ) const
    {
        BELFEM_ERROR( false, "Function E2 not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

//----------------------------------------------------------------------------

    real
    OrthotropicMaterial::E3(  const real aT ) const
    {
        BELFEM_ERROR( false, "Function E3 not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

//----------------------------------------------------------------------------

    real
    OrthotropicMaterial::nu23(  const real aT ) const
    {
        BELFEM_ERROR( false, "Function nu23 not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

//----------------------------------------------------------------------------

    real
    OrthotropicMaterial::nu12(  const real aT ) const
    {
        BELFEM_ERROR( false, "Function nu12 not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

//----------------------------------------------------------------------------

    real
    OrthotropicMaterial::nu31(  const real aT ) const
    {
        BELFEM_ERROR( false, "Function nu31 not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

//----------------------------------------------------------------------------

    real
    OrthotropicMaterial::G23(  const real aT ) const
    {
        BELFEM_ERROR( false, "Function G23 not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

//----------------------------------------------------------------------------

    real
    OrthotropicMaterial::G12(  const real aT ) const
    {
        BELFEM_ERROR( false, "Function G12 not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

//----------------------------------------------------------------------------

    real
    OrthotropicMaterial::G31(  const real aT ) const
    {
        BELFEM_ERROR( false, "Function G31 not implemented for material %s", mLabel.c_str() );
        return BELFEM_QUIET_NAN ;
    }

//----------------------------------------------------------------------------

}