//
// Created by Christian Messe on 13.11.19.
//

#include "assert.hpp"
#include "cl_IsotropicMaterial.hpp"
#include "fn_polyval.hpp"
#include "fn_dpolyval.hpp"
#include "fn_linspace.hpp"
#include "fn_polyfit.hpp"

namespace belfem
{
//----------------------------------------------------------------------------

    IsotropicMaterial::IsotropicMaterial( const MaterialType aType ) :
            Material( aType )
    {
    }

//----------------------------------------------------------------------------

    real
    IsotropicMaterial::E( const real aT ) const
    {
        BELFEM_ASSERT( mYoungPoly.length() > 0,
                      "Polynomial for Young's modulus is not set for %s",
                      mLabel.c_str());
        return polyval( mYoungPoly, aT );
    }

//----------------------------------------------------------------------------

    real
    IsotropicMaterial::nu( const real aT ) const
    {
        return 0.5 * ( this->E( aT ) / this->G( aT ) )  - 1.0;
    }

//----------------------------------------------------------------------------

    real
    IsotropicMaterial::G( const real aT ) const
    {
        BELFEM_ASSERT( mShearPoly.length() > 0,
                      "Polynomial for Shear modulus is not set for %s",
                      mLabel.c_str());
        return polyval( mShearPoly, aT );
    }

//----------------------------------------------------------------------------

    void
    IsotropicMaterial::C( Matrix< real > & aC, const real aT ) const
    {
        BELFEM_ASSERT( aC.n_cols() == 6 && aC.n_rows() == 6,
                      "Elasticity matrix must be allocated as 6x6" );

        // elasticity modulus
        real tE = this->E( aT );

        // shear modulus
        real tG = this->G( aT );

        // Poisson number
        real tNu = 0.5 * tE / tG - 1.0;

        // help values
        real tD = 1.0 - tNu * ( 1.0 + 2.0 * tNu );
        real tC0 = tE * ( tNu - 1.0 ) / tD;
        real tC1 = tE * tNu / tD;

        aC( 0, 0 ) = tC0;
        aC( 1, 0 ) = tC1;
        aC( 2, 0 ) = tC1;
        aC( 3, 0 ) = 0.0;
        aC( 4, 0 ) = 0.0;
        aC( 5, 0 ) = 0.0;

        aC( 0, 1 ) = tC1;
        aC( 1, 1 ) = tC0;
        aC( 2, 1 ) = tC1;
        aC( 3, 1 ) = 0.0;
        aC( 4, 1 ) = 0.0;
        aC( 5, 1 ) = 0.0;

        aC( 0, 2 ) = tC1;
        aC( 1, 2 ) = tC1;
        aC( 2, 2 ) = tC0;
        aC( 3, 2 ) = 0.0;
        aC( 4, 2 ) = 0.0;
        aC( 5, 2 ) = 0.0;

        aC( 0, 3 ) = 0.0;
        aC( 1, 3 ) = 0.0;
        aC( 2, 3 ) = 0.0;
        aC( 3, 3 ) = tG;
        aC( 4, 3 ) = 0.0;
        aC( 5, 3 ) = 0.0;

        aC( 0, 4 ) = 0.0;
        aC( 1, 4 ) = 0.0;
        aC( 2, 4 ) = 0.0;
        aC( 3, 4 ) = 0.0;
        aC( 4, 4 ) = tG;
        aC( 5, 4 ) = 0.0;

        aC( 0, 5 ) = 0.0;
        aC( 1, 5 ) = 0.0;
        aC( 2, 5 ) = 0.0;
        aC( 3, 5 ) = 0.0;
        aC( 4, 5 ) = 0.0;
        aC( 5, 5 ) = tG;
    }

//----------------------------------------------------------------------------

    /**
     * Elasticity matrix in plane stress
     */
    void
    IsotropicMaterial::C_ps( Matrix< real > & aC, const real aT ) const
    {
        BELFEM_ASSERT( aC.n_cols() == 3 && aC.n_rows() == 3,
                      "Elasticity matrix must be allocated as 3x3" );

        // elasticity modulus
        real tE = this->E( aT );

        // shear modulus
        real tG = this->G( aT );

        // Poisson number
        real tNu = 0.5 * tE / tG - 1.0;

        // help values
        real tD = 1.0 - tNu * tNu;
        real tC0 = tE / tD;
        real tC1 = tC0 * tNu;

        aC( 0, 0 ) = tC0;
        aC( 1, 0 ) = tC1;
        aC( 2, 0 ) = 0.0;

        aC( 0, 1 ) = tC1;
        aC( 1, 1 ) = tC0;
        aC( 2, 1 ) = 0.0;

        aC( 0, 2 ) = 0.0;
        aC( 1, 2 ) = 0.0;
        aC( 2, 2 ) = tG;
    }

//----------------------------------------------------------------------------

    void
    IsotropicMaterial::C_rot( Matrix< real > & aC, const real aT ) const
    {
        BELFEM_ASSERT( aC.n_cols() == 4 && aC.n_rows() == 4,
                      "Elasticity matrix must be allocated as 4x4" );

        // elasticity modulus
        real tE = this->E( aT );

        // shear modulus
        real tG = this->G( aT );

        // Poisson number
        real tNu = 0.5 * tE / tG - 1.0;

        // help values
        real tD = (( 2.0 * tNu + 1.0 ) * tNu - 1.0 );
        real tC0 = tE * ( 1.0 - tNu ) / tD;
        real tC1 = tE * tNu;

        aC( 0, 0 ) = tC0;
        aC( 1, 0 ) = tC1;
        aC( 2, 0 ) = tC1;
        aC( 3, 0 ) = 0.0;

        aC( 0, 0 ) = tC1;
        aC( 1, 0 ) = tC0;
        aC( 2, 0 ) = tC1;
        aC( 3, 0 ) = 0.0;

        aC( 0, 0 ) = tC1;
        aC( 1, 0 ) = tC1;
        aC( 2, 0 ) = tC0;
        aC( 3, 0 ) = 0.0;

        aC( 0, 0 ) = 0.0;
        aC( 1, 0 ) = 0.0;
        aC( 2, 0 ) = 0.0;
        aC( 3, 0 ) = tG;
    }

//----------------------------------------------------------------------------

    real
    IsotropicMaterial::rho( const real aT ) const
    {
        BELFEM_ASSERT( mDensityPoly.length() > 0,
                      "Polynomial for density is not set for %s",
                      mLabel.c_str());
        return polyval( mDensityPoly, aT );
    }

//----------------------------------------------------------------------------
//   THERMAL PROPERTIES
//----------------------------------------------------------------------------

    /**
     * thermal conductivity in W/(m*K)
     */
    real
    IsotropicMaterial::lambda( const real aT ) const
    {
        BELFEM_ASSERT( "thermal conductivity is not implemented for %s",
                      mLabel.c_str());
        return BELFEM_QUIET_NAN ;
    }

//----------------------------------------------------------------------------

    void
    IsotropicMaterial::lambda_p( Matrix< real > & aLambda, const real aT ) const
    {
        BELFEM_ASSERT( aLambda.n_cols() == 2 && aLambda.n_rows() == 2,
                      "Conductivity matrix must be allocated as 2x2" );

        real tLambda = this->lambda( aT );

        aLambda( 0, 0 ) = tLambda;
        aLambda( 1, 0 ) = 0.0;

        aLambda( 0, 1 ) = 0.0;
        aLambda( 1, 1 ) = tLambda;
    }

//----------------------------------------------------------------------------

    void
    IsotropicMaterial::lambda( Matrix< real > & aLambda, const real aT ) const
    {
        BELFEM_ASSERT( aLambda.n_cols() == 3 && aLambda.n_rows() == 3,
                      "Conductivity matrix must be allocated as 3x3" );

        real tLambda = this->lambda( aT );

        aLambda( 0, 0 ) = tLambda;
        aLambda( 1, 0 ) = 0.0;
        aLambda( 2, 0 ) = 0.0;

        aLambda( 0, 1 ) = 0.0;
        aLambda( 1, 1 ) = tLambda;
        aLambda( 2, 1 ) = 0.0;

        aLambda( 0, 2 ) = 0.0;
        aLambda( 1, 2 ) = 0.0;
        aLambda( 2, 2 ) = tLambda;
    }

//----------------------------------------------------------------------------

    real
    IsotropicMaterial::c( const real aT ) const
    {
        BELFEM_ASSERT( "specific heat capacity is not implemented for %s",
                       mLabel.c_str());
        return BELFEM_QUIET_NAN ;;
    }

//----------------------------------------------------------------------------
// Thermal Expansion
//----------------------------------------------------------------------------

    real
    IsotropicMaterial::mu( const real aT, const real aTref  ) const
    {
        BELFEM_ASSERT( mThermalExpansionPoly.length() > 0,
                      "Polynomial for thermal expansion is not set for %s",
                      mLabel.c_str() );

        return std::exp(
                this->int_alpha_dT( aT )
                - this->int_alpha_dT( aTref ) ) - 1.0;
    }

//----------------------------------------------------------------------------

    void
    IsotropicMaterial::mu_ps( Matrix< real > & aMu, const real aT, const real aTref ) const
    {
        BELFEM_ASSERT( aMu.n_cols() == 2 && aMu.n_rows() == 2,
                      "Expansion matrix must be allocated as 2x2" );

        real tMu = this->mu( aT, aTref );

        aMu( 0, 0 ) = tMu;
        aMu( 1, 0 ) = 0.0;

        aMu( 0, 1 ) = 0.0;
        aMu( 1, 1 ) = tMu;
    }

//----------------------------------------------------------------------------

    void
    IsotropicMaterial::mu( Matrix< real > & aMu, const real aT, const real aTref ) const
    {
        BELFEM_ASSERT( aMu.n_cols() == 3 && aMu.n_rows() == 3,
                      "Expansion matrix must be allocated as 3x3" );

        real tMu = this->mu( aT, aTref );

        aMu( 0, 0 ) = tMu;
        aMu( 1, 0 ) = 0.0;
        aMu( 2, 0 ) = 0.0;

        aMu( 0, 1 ) = 0.0;
        aMu( 1, 1 ) = tMu;
        aMu( 2, 1 ) = 0.0;

        aMu( 0, 2 ) = 0.0;
        aMu( 1, 2 ) = 0.0;
        aMu( 2, 2 ) = tMu;
    }

//----------------------------------------------------------------------------
// Optical Properties
//----------------------------------------------------------------------------

    real
    IsotropicMaterial::epsilon( const real aT ) const
    {
        BELFEM_ASSERT( mEpsilonPoly.length() > 0,
                      "Polynomial for thermal emmisivity is not set for %s",
                      mLabel.c_str() );
        return polyval( mEpsilonPoly, aT );
    }

//----------------------------------------------------------------------------

    real
    IsotropicMaterial::alpha( const real  aT ) const
    {
        BELFEM_ASSERT( mThermalExpansionPoly.length() > 0,
                      "Polynomial for thermal expansion is not set for %s",
                      mLabel.c_str() );

        return dpolyval( mThermalExpansionPoly, aT );
    }

//----------------------------------------------------------------------------

    real
    IsotropicMaterial::int_alpha_dT( const real aT ) const
    {
        BELFEM_ASSERT( mThermalExpansionPoly.length() > 0,
                      "Polynomial for thermal expansion is not set for %s",
                      mLabel.c_str() );
        return polyval( mThermalExpansionPoly, aT );
    }

//----------------------------------------------------------------------------

    bool
    IsotropicMaterial::has_thermal() const
    {
        return mHasThermal ;
    }

//----------------------------------------------------------------------------

    bool IsotropicMaterial::has_electric_resistivity() const
    {
        return mHasResistivity ;
    }

//----------------------------------------------------------------------------

    bool
    IsotropicMaterial::has_mechanical() const
    {
        return mHasMechanical || mYoungPoly.length() > 0 ;
    }

//----------------------------------------------------------------------------

    bool
    IsotropicMaterial::has_thermal_expansion() const
    {
        return mHasExpansion || mThermalExpansionPoly.length() > 0 ;
    }

//----------------------------------------------------------------------------

    void
    IsotropicMaterial::create_density_poly( const real aRhoRef, const real aTref )
    {
        BELFEM_ASSERT( mThermalExpansionPoly.length() > 1,
                       "Thermal expansion poly must be set and deg > 1 before computing density poly for %s",
                       this->label().c_str() );

        BELFEM_ASSERT( mTmax < BELFEM_REAL_MAX,
                       "maximum temperature must be set before computing density poly for %s",
                       this->label().c_str() );

        if( mThermalExpansionPoly.length() == 2 )
        {
            mDensityPoly.set_size( 1, aRhoRef );
            return;
        }

        // get interpolation degree
        uint tDegree = mThermalExpansionPoly.length() - 2 ;

        uint tN = 101 ;

        // temperature
        Vector< real > tT = linspace( 0.0, mTmax, tN );
        Vector< real > tRho( tN );

        // reference value
        real tRef = polyval( mThermalExpansionPoly, aTref );
        for ( uint k=0; k<tN; ++k )
        {
            // compute the length
            real tL = std::exp(  polyval( mThermalExpansionPoly, tT( k ) ) - tRef );

            // scale density
            tRho( k ) = aRhoRef / ( tL * tL * tL );
        }

        // compute density polynomial
        polyfit( tT, tRho, tDegree, mDensityPoly );

        mRhoRef = polyval( mDensityPoly, BELFEM_TREF );

        // make sure that switch is flipped
        mHasExpansion = true ;
    }

//----------------------------------------------------------------------------
}