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

#ifndef BELFEM_CL_MATERIAL_CHROMIUM_HPP
#define BELFEM_CL_MATERIAL_CHROMIUM_HPP

#include "cl_Material_Metal.hpp"
#include "cl_Bezier.hpp"

namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        /**
         * @brief Chromium ( Cr ), body centered cubic, antiferromagnetic
         *
         * Registered in cl_MaterialFactory as "chromium", "cr". Thermal
         * expansion, heat capacity, the Debye curve, the resistivity anchor, the
         * thermal conductivity and the Kohler magnetoresistance are fitted; the
         * elastic moduli are a Wachtman curve against Armstrong and Brown 1964
         * with the Poisson ratio derived from a constant Grueneisen parameter
         * ( Metal::create_mech ). The elastic anomaly at the Neel point is not
         * represented, see below.
         *
         * MAGNETIC TRANSITIONS - WHAT IS AND IS NOT REPRESENTED
         *
         * Chromium carries a spin density wave with a Neel transition at 311 K
         * and a spin flip transition at 123 K. Both sit inside the temperature
         * range BELFEM operates in, and both produce anomalies in alpha, cp and
         * the elastic constants. None of those anomalies is represented here:
         * neither the log-log Bezier form behind Metal::create_cp() nor a smooth
         * expansion Bezier has a channel for a lambda peak. Treat properties
         * within roughly 20 K of either transition as smoothed, and say so
         * wherever they are consumed.
         *
         * Chromium uses the Grueneisen coupling like every other metal, despite
         * the transitions. The reasoning:
         *
         *   - dln(C)/dT runs smoothly and monotonically through both 123 K and
         *     311 K with no feature at either, so the magnetic channel does not
         *     disturb the alpha/cp ratio the branch is anchored on.
         *   - Without the coupling, the expansion Bezier alone overstates alpha
         *     by 2.1x at 100 K, 4.1x at 77 K and 12.4x at 50 K measured against
         *     the fitted cp - far larger than the few percent Neel anomaly a
         *     smooth model smears out.
         *   - The split lands above the dln(C)/dT sign change at 261.8 K and
         *     below the Neel point at 311 K, so the anchor sits inside a single
         *     coherent antiferromagnetic phase.
         *
         * Two caveats. The margin above the sign change is 11.4 K, the thinnest
         * in the roster, so a refit of either cp or the expansion curve can flip
         * it - the BELFEM_ERROR in create_low_temperature_alpha() will say so.
         * And the implied Grueneisen parameter reads about 0.91 against a
         * literature 1.3 - 1.5 ; chromium genuinely has a low gamma because
         * magnetostriction partly cancels lattice expansion, but this is the one
         * cross check that does not corroborate.
         *
         * Chromium is needed in the database as an alloy constituent -
         * Alloy::create_tables() homogenizes E and nu from its components, so Cr
         * bearing alloys such as stainless steels and the nickel superalloys
         * cannot be built through that path without it. The Neel anomaly is
         * suppressed by alloying and does not propagate into them.
         */
        class Chromium : public Metal
        {
            Bezier * mThermalExpansion = nullptr ;
            Vector< real > mThermalExpansionCryo ;

            Cell< Bezier * > mDebyeBeziers ;

            Vector< real > mBSKohlerSwitch ;
            Cell< Vector< real > > mKohlerLongPolys ;
            Cell< Vector< real > > mKohlerTransPolys ;

        public:

            Chromium( const real RRR = BELFEM_QUIET_NAN,
                      const bool aBuildTables = true );

            ~Chromium() override ;

            //! Neel temperature ( spin density wave ) [K]
            static constexpr real gTNeel = 311.0 ;

            //! spin flip transition temperature [K]
            static constexpr real gTSpinFlip = 123.0 ;

        protected:

            real
            alpha_custom( const real T ) const override ;

            real debye_custom(const real T) const override ;

            real
            kohler( const real B, const real S, const real beta ) const override ;

        private:

            void
            set_constants();

            void
            create_alpha();

            void
            create_cp();

            void
            create_debye();

            void
            create_kohler();
        };

//------------------------------------------------------------------------------

    }
}
#endif //BELFEM_CL_MATERIAL_CHROMIUM_HPP
