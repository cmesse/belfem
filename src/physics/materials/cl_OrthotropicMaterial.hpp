//
// Created by Christian Messe on 13.11.19.
//

#ifndef BELFEM_CL_ORTHOTROPICMATERIAL_HPP
#define BELFEM_CL_ORTHOTROPICMATERIAL_HPP


#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Material.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    class OrthotropicMaterial : public Material
    {

//----------------------------------------------------------------------------
    public:
//----------------------------------------------------------------------------

        OrthotropicMaterial( const MaterialType aType ) ;

//----------------------------------------------------------------------------

        virtual ~OrthotropicMaterial() = default;

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
        C_rot( Matrix< real > & aC, const real aT=BELFEM_TREF ) const;

//----------------------------------------------------------------------------
    protected:
//----------------------------------------------------------------------------

        void
        compute_C3D(
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
                ) const;

//----------------------------------------------------------------------------

        void
        compute_C_ps(
                const real & aE1,
                const real & aE2,
                const real & aNu12,
                const real & aG12,
                Matrix< real > & aC
        ) const;

//----------------------------------------------------------------------------

        void
        compute_C_rot(
                const real & aE1,
                const real & aE2,
                const real & aE3,
                const real & aNu23,
                const real & aNu12,
                const real & aNu31,
                const real & aG12,
                Matrix< real > & aC
        ) const;

//----------------------------------------------------------------------------

        virtual real
        E1(  const real aT=BELFEM_TREF ) const;

        virtual real
        E2(  const real aT=BELFEM_TREF ) const;

        virtual real
        E3(  const real aT=BELFEM_TREF ) const;

        virtual real
        nu23(  const real aT=BELFEM_TREF ) const;

        virtual real
        nu12(  const real aT=BELFEM_TREF ) const;

        virtual real
        nu31(  const real aT=BELFEM_TREF ) const;

        virtual real
        G23(  const real aT=BELFEM_TREF ) const;

        virtual real
        G12(  const real aT=BELFEM_TREF ) const;

        virtual real
        G31(  const real aT=BELFEM_TREF ) const;

//----------------------------------------------------------------------------
    };
}

#endif //BELFEM_CL_ORTHOTROPICMATERIAL_HPP
