/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California, through
 * Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 * 
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_CL_MATERIAL_IRON_HPP
#define BELFEM_CL_MATERIAL_IRON_HPP

#include "cl_Material_Ferromagnetic.hpp"
#include "cl_Bezier.hpp"

namespace belfem
{
    namespace material
    {
        class Iron : public Ferromagnetic
        {
            Bezier * mReducedMagnetization = nullptr ;
            Bezier * mDebyeBezier = nullptr ;

            Bezier * mThermalExpansion = nullptr ;
            Vector< real > mThermalExpansionCryo ;

            Cell< Vector< real > > mDebyePolys ;
            Vector< real > mTDebyeSwitch ;

            Vector< real > mBSKohlerSwitch ;
            Cell< Vector< real > > mKohlerLongPolys ;
            Cell< Vector< real > > mKohlerTransPolys ;
            Bezier * mKohlerTransBezier = nullptr ;
        public:

            Iron( const real RRR = BELFEM_QUIET_NAN,
                  const bool aBuildTables = true );

            ~Iron() override;

        protected :

            real
            compute_mred( real T ) const override;

            real
            debye_custom( const real T ) const override ;

            real
            alpha_custom( const real T ) const override ;

            real
            kohler(const real B, const real S, const real beta) const override;

        private:

            void
            create_magnetization_curve();

            void
            create_debye_and_rho();

            void
            create_cp();

            void
            create_alpha();
            
            void
            set_constants();

            void
            create_kohler();


        };

        inline
        real Iron::compute_mred( real T ) const
        {
            real Tc = this->constant_property( MaterialProperty::Tcurie );
            if ( T  > Tc + BELFEM_EPSILON ) return 0.0 ;
            return mReducedMagnetization->y( T / Tc );
        }

    }
}

#endif //BELFEM_CL_MATERIAL_IRON_HPP
