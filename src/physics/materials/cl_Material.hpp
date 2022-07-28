//
// Created by Christian Messe on 13.11.19.
//

#ifndef BELFEM_CL_MATERIAL_HPP
#define BELFEM_CL_MATERIAL_HPP

#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_Matrix.hpp"
#include "cl_Material.hpp"
#include "en_Materials.hpp"

namespace belfem
{
    class Material
    {

//----------------------------------------------------------------------------
    protected:
//----------------------------------------------------------------------------
        const MaterialType mType ;

        string mLabel = "Undefined";

        // maximum temperture
        real mTmax = BELFEM_REAL_MAX ;

        // flag telling if this material is used
        bool mFlag = false ;

        PermeabilityLaw mPermeabilityLaw = PermeabilityLaw::Constant ;
        ResistivityLaw mResistivityLaw = ResistivityLaw::Constant ;

//----------------------------------------------------------------------------
    public:
//----------------------------------------------------------------------------

        Material( const MaterialType aType ) ;

//----------------------------------------------------------------------------

        virtual ~Material() = default;

//----------------------------------------------------------------------------

        const string &
        label() const;

//----------------------------------------------------------------------------

        /**
         * returns true if this is an isotropic material
         */
         virtual bool
         is_isotropic() const ;

//----------------------------------------------------------------------------

        /**
         * returns the material enum
         */
        MaterialType
        type() const ;

//----------------------------------------------------------------------------

        /**
         * returns true if thermal properties exist
         */
         virtual bool
         has_thermal() const ;

//----------------------------------------------------------------------------

        /**
         * returns true if mechanical properties exist
         */
        virtual bool
        has_mechanical() const ;

//----------------------------------------------------------------------------

        /**
         * returns true if thermal expansion data exist
         */
        virtual bool
        has_thermal_expansion() const ;

//----------------------------------------------------------------------------
        /**
         * returns if electric resistivity properties exist
         */
        virtual bool
        has_electric_resistivity() const ;

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
        C_ps( Matrix< real > & aC, const real aT=BELFEM_TREF ) const;

//----------------------------------------------------------------------------

        /**
         * Elasticity matrix in rotation symmetry
         */
        virtual void
        C_rot(  Matrix< real > & aC, const real aT=BELFEM_TREF ) const;

//----------------------------------------------------------------------------

        /**
          * Density in kg/m^3
          */
        virtual real
        rho( const real aT=BELFEM_TREF ) const;

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
         * thermal conductivity in W/(m*K)
         */
        virtual real
        lambda( const real aT, const real aB ) const;

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
        lambda( Matrix< real > & aLambda, const real aT=BELFEM_TREF ) const;

//----------------------------------------------------------------------------
// Thermal Expansion
//----------------------------------------------------------------------------

        virtual real
        alpha( const real aT=BELFEM_TREF ) const ;

//----------------------------------------------------------------------------

        /**
         * thermal strain parameter
         *
         * = exp( int( a, T, aTref, aT ) ) - 1
         */
        virtual real
        mu( const real aT=BELFEM_TREF, const real aTref=BELFEM_TREF ) const ;

//----------------------------------------------------------------------------

        /**
         * thermal strain matrix ( plane stress )
         */
        virtual void
        mu_ps( Matrix< real > & aMu, const real aT=BELFEM_TREF, const real aTref=BELFEM_TREF ) const ;

//----------------------------------------------------------------------------

        /**
         * thermal strain matrix ( 3D )
         */
        virtual void
        mu(  Matrix< real > & aMu, const real aT=BELFEM_TREF, const real aTref=BELFEM_TREF ) const ;

//----------------------------------------------------------------------------
// Optical Properties
//----------------------------------------------------------------------------

        /**
         * thermal emmisivity
         */
        virtual real
        epsilon( const real aT=BELFEM_TREF ) const ;

//----------------------------------------------------------------------------

        real
        T_max() const ;

//----------------------------------------------------------------------------
// electromagnetic properties
//----------------------------------------------------------------------------

        // permeability
        virtual real
        mu_r ( const real aH=0, const real aT=BELFEM_TREF ) const ;

        // electric resistance
        virtual
        real rho_el ( const real aJ=0, const real aT=0, const real aB=0 ) const ;


//----------------------------------------------------------------------------

        /**
         * Secant Permeability
         */
        virtual real
        nu_s( const real aH=0, const real aT=BELFEM_TREF ) const;


//----------------------------------------------------------------------------

        /**
         * return the resisivity law
         */
        ResistivityLaw
        resistivity_law() const ;

//----------------------------------------------------------------------------

        /**
         * return the permeability law
         */
        PermeabilityLaw
        permeability_law() const ;

//----------------------------------------------------------------------------

        /**
         * this flag tells if the material has been linked to a block
         */
        void
        flag();

//----------------------------------------------------------------------------

        /**
         * this flag tells if the material has been linked to a block
         */
        bool
        is_flagged() const ;

//----------------------------------------------------------------------------

        /**
         * this is only relevant for copper and silver
         */
        virtual void
        set_rrr( const real aRRR ) ;

//----------------------------------------------------------------------------

        /**
         * this is only relevant if we use a spline underneath
         */
        virtual void
        use_splines( const bool aSwitch ) ;

//----------------------------------------------------------------------------
    };

//----------------------------------------------------------------------------

    inline bool
    Material::is_isotropic() const
    {
        return false ;
    }

//----------------------------------------------------------------------------

    inline real
    Material::T_max() const
    {
        return mTmax ;
    }

//----------------------------------------------------------------------------

    inline MaterialType
    Material::type() const
    {
        return mType ;
    }

//----------------------------------------------------------------------------

    inline ResistivityLaw
    Material::resistivity_law() const
    {
        return mResistivityLaw ;
    }

//----------------------------------------------------------------------------

    inline PermeabilityLaw
    Material::permeability_law() const
    {
        return mPermeabilityLaw ;
    }

//----------------------------------------------------------------------------

    inline void
    Material::flag()
    {
        mFlag = true ;
    }

//----------------------------------------------------------------------------

    inline bool
    Material::is_flagged() const
    {
        return mFlag ;
    }

//----------------------------------------------------------------------------

    inline real
    Material::lambda( const real aT, const real aB ) const
    {
        // ignore b-field unless implemented otherwise
        return this->lambda( aT );
    }

//----------------------------------------------------------------------------
}
#endif //BELFEM_CL_MATERIAL_HPP
