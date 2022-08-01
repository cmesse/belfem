//
// Created by Christian Messe on 13.11.19.
//

#ifndef BELFEM_CL_ISOTROPICMATERIAL_HPP
#define BELFEM_CL_ISOTROPICMATERIAL_HPP

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Material.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    class MaterialFactory;

    class IsotropicMaterial : public Material
    {
//----------------------------------------------------------------------------
    protected:
//----------------------------------------------------------------------------

        // polynomial for Young's Modulus
        Vector< real > mYoungPoly;

        // polynomial for Shear Modulus
        Vector< real > mShearPoly;

        // polynomial for density
        Vector< real > mDensityPoly;

        // density at reference temperature
        real mRhoRef = BELFEM_QUIET_NAN ;

        // polynomial for specific heat
        // Vector< real > mSpecificHeatPoly;

        // polynomial for thermal conductivity
        // Vector< real > mThermalConductivityPoly;

        // polynomial for thermal expansion
        Vector< real > mThermalExpansionPoly;

        // polynomial for thermal radiation
        Vector< real > mEpsilonPoly;

        bool mHasThermal      = false ;
        bool mHasMechanical   = false ;
        bool mHasExpansion    = false ;
        bool mHasResistivity  = false ;
        friend MaterialFactory;

//----------------------------------------------------------------------------
    public:
//----------------------------------------------------------------------------

        IsotropicMaterial( const MaterialType aType ) ;

//----------------------------------------------------------------------------

        virtual ~IsotropicMaterial() = default;

//----------------------------------------------------------------------------

        bool
        is_isotropic() const ;

//----------------------------------------------------------------------------
//      ELASTIC PROPERTIES
//----------------------------------------------------------------------------

        /**
         * Young's Modulus in Pa
         */
        virtual real
        E( const real aT=BELFEM_TREF ) const;

//----------------------------------------------------------------------------

        /**
         * Poisson Number
         */
        virtual real
        nu( const real aT=BELFEM_TREF ) const;

//----------------------------------------------------------------------------

        /**
         * Shear Modulus in Pa
         */
        virtual real
        G( const real aT=BELFEM_TREF ) const;

//----------------------------------------------------------------------------

        /**
          * Elasticity Matrix in 3D
          */
        virtual void
        C( Matrix< real > & aC, const real aT=BELFEM_TREF ) const;

//----------------------------------------------------------------------------

        /**
         * Elasticity matrix in plane stress
         */
        virtual void
        C_ps( Matrix< real > & aC, const real aT=BELFEM_TREF  ) const;

//----------------------------------------------------------------------------

        /**
         * Elasticity matrix in rotation symmetry
         */
        virtual void
        C_rot( Matrix< real > & aC, const real aT=BELFEM_TREF  ) const;

//----------------------------------------------------------------------------

        /**
          * Density in kg/m^3
          */
        virtual real
        rho( const real aT=BELFEM_TREF ) const;

//--------------------------------------------------------------------------

        /**
         * reference density
         */
         real
         rho0() const ;

//----------------------------------------------------------------------------
//   THERMAL PROPERTIES
//----------------------------------------------------------------------------

        /**
         * specific heat capacity in J/(kg*K)
         */
        virtual real
        c( const real aT=BELFEM_TREF ) const;

//----------------------------------------------------------------------------

        /**
         * thermal conductivity in W/(m*K)
         */
        virtual real
        lambda( const real aT=BELFEM_TREF ) const;

//----------------------------------------------------------------------------

        /**
          * thermal conductivity in 2D
          */
        virtual void
        lambda_p( Matrix< real > & aLambda, const real aT=BELFEM_TREF ) const;

//----------------------------------------------------------------------------

        /**
          * thermal conductivity in 3D
          */
        virtual void
        lambda3d(  Matrix< real > & aLambda, const real aT=BELFEM_TREF ) const;

//----------------------------------------------------------------------------
// Thermal Expansion
//----------------------------------------------------------------------------

        virtual real
        alpha( const real  aT=BELFEM_TREF ) const ;

//----------------------------------------------------------------------------

        /**
         * thermal strain parameter
         *
         * int( alpha, T, aTref, aT ) )
         */
        virtual real
        mu( const real aT=BELFEM_TREF, const real aTref=BELFEM_TREF ) const;

//----------------------------------------------------------------------------

        /**
         * thermal strain matrix ( plane stress )
         */
        virtual void
        mu_ps(  Matrix< real > & aMu,
                const real aT=BELFEM_TREF,
                const real aTref=BELFEM_TREF ) const ;

//----------------------------------------------------------------------------

        /**
         * thermal strain matrix ( 3D )
         */
        virtual void
        mu(  Matrix< real > & aMu,
                const real aT=BELFEM_TREF,
                const real aTref=BELFEM_TREF ) const ;

//----------------------------------------------------------------------------
// Optical Properties
//----------------------------------------------------------------------------

        /**
         * thermal emmisivity
         */
        virtual real
        epsilon( const real aT=BELFEM_TREF ) const;

//----------------------------------------------------------------------------
// Data Information
//----------------------------------------------------------------------------

        virtual bool
        has_thermal() const ;

//----------------------------------------------------------------------------

        virtual bool
        has_mechanical() const ;

//----------------------------------------------------------------------------

        virtual bool
        has_thermal_expansion() const ;

//----------------------------------------------------------------------------

        virtual bool
        has_electric_resistivity() const ;

//----------------------------------------------------------------------------
    protected:
//----------------------------------------------------------------------------

        /**
         * integrand of alpha over T
         */
        virtual real
        int_alpha_dT( const real aT ) const;

//----------------------------------------------------------------------------

        void
        create_density_poly( const real aRhoRef, const real aTref = BELFEM_TREF );

//----------------------------------------------------------------------------
    };
//----------------------------------------------------------------------------

    inline bool
    IsotropicMaterial::is_isotropic() const
    {
        return true ;
    }

//----------------------------------------------------------------------------

    inline real
    IsotropicMaterial::rho0() const
    {
        return mRhoRef ;
    }

//----------------------------------------------------------------------------
}

#endif //BELFEM_CL_ISOTROPICMATERIAL_HPP
