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

#ifndef BELFEM_CL_MATERIAL_NICKEL_HPP
#define BELFEM_CL_MATERIAL_NICKEL_HPP


#include "cl_Material_Ferromagnetic.hpp"
#include "cl_Bezier.hpp"

namespace belfem
{
    namespace material
    {
        /**
         * @brief Nickel ( Ni ), face centered cubic, ferromagnetic
         *
         * ELASTIC MODULI - WHAT IS AND IS NOT REPRESENTED
         *
         * The moduli follow the same quasi-harmonic closure as the other metals
         * ( Metal::create_mech ): K and G soften with the volumetric thermal strain,
         * fitted to the single-crystal constants of Alers, Neighbours and Sato 1960
         * ( 10.1016/0022-3697(60)90125-6, Table 2 ), measured in a saturating field
         * of 10 kOe, Hill-averaged to the isotropic polycrystal; room-temperature
         * cross-check against Ledbetter and Reed 1973 ( 10.1063/1.3253127, Table 6 ).
         *
         * Nickel's Young's modulus in the DEMAGNETIZED state is lower and shows a deep
         * minimum below the Curie point ( Blanke 1989: 183.7 GPa at 293 K, 124 GPa at
         * 487 K, 195.7 GPa at 637 K; the fit carried here until 2026-09 was already
         * 18 % below its room-temperature value at 400 K ). That is the Delta-E
         * effect: domain walls move under stress and add strain; at magnetic
         * saturation the walls are immobile ( Ledbetter and Reed 1973, section 15 ).
         * The effect depends on field, stress and annealing state ( 190.5 GPa at H = 0
         * against 225.6 GPa at 6.2 kOe on one sample, Giebe and Blechschmidt 1931,
         * ibid. Table 9 ), so it is a property of the magnetic state, not of the
         * lattice. BELFEM's nickel parts sit in tesla-level fields and are saturated;
         * the saturated moduli are therefore the ones served, and the Delta-E dip is
         * deliberately not represented. A smaller intrinsic magnetic contribution
         * remains below the Curie point even at saturation ( Alers et al. 1960 ); the
         * constants are fitted over 0-300 K, so above 300 K the served curve is the
         * quasi-harmonic extrapolation and that remainder is not represented either.
         * Until 2026-09 the class carried Blanke's demagnetized
         * curve as a two-branch Wachtman fit with a Bezier bridge; it was retired
         * because it describes a state the solver does not simulate, its magnitude is
         * uncertain by tens of percent, and its sign change of dE/dT inside 465-631 K
         * would be a hazard for the Newton tangent of a thermal-stress run.
         * Consequence: an unmagnetized nickel part near 500 K is served an E about
         * 40 % too stiff, and the room-temperature E rose from 184 to about 223 GPa.
         * The thermal expansion Bezier ends at 600 K and is clamped up to T_max.
         * Decision: Christian Messe, 2026-09-15.
         */
        class Nickel : public Ferromagnetic
        {
            Bezier * mReducedMagnetization = nullptr ;
            Bezier * mThermalExpansion = nullptr ;
            Vector< real > mThermalExpansionCryo ;

            // debye temperature
            Vector< real > mDebyePoly ;

            Vector< real > mBSKohlerSwitch ;
            Bezier * mKohlerLongBezier = nullptr ;
            Cell< Vector< real > > mKohlerLongPolys ;
            Cell< Vector< real > > mKohlerTransPolys ;

        public:

            Nickel( const real RRR = BELFEM_QUIET_NAN,
                  const bool aBuildTables = true );

            ~Nickel() override;

        protected :

            real
            compute_mred( real T ) const override;

            real
            alpha_custom( const real T ) const override ;

            real
            debye_custom(const real T) const override ;

            real
            kohler( const real B, const real S, const real beta ) const override ;

        private:

            void
            create_magnetization_curve();

            void
            set_constants();

            void
            create_alpha();

            void
            create_cp();

            void
            create_debye_and_rho();

            void
            create_kohler();

        };

        inline real Nickel::compute_mred( real T ) const
        {
            real Tc = this->constant_property( MaterialProperty::Tcurie );
            if ( T  > Tc + BELFEM_EPSILON ) return 0.0 ;
            return mReducedMagnetization->y( T / Tc );
        }

        inline real Nickel::debye_custom( const real T) const
        {
            return polyval( mDebyePoly, T );
        }

    }
}

#endif //BELFEM_CL_MATERIAL_NICKEL_HPP
