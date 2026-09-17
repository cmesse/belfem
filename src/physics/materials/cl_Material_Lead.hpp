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

#ifndef BELFEM_CL_MATERIAL_LEAD_HPP
#define BELFEM_CL_MATERIAL_LEAD_HPP

#include "cl_Material_Metal.hpp"
#include "cl_Bezier.hpp"
#include "fn_polyval.hpp"
namespace belfem
{
    namespace material
    {
        /**
         * @brief Lead ( Pb ), face centered cubic
         *
         * ELASTIC MODULI - STATIC LEVEL, WHY
         *
         * The moduli follow the quasi-harmonic closure of Metal::create_mech. The bulk
         * modulus and its softening constant are those of the single crystal of Waldorf
         * and Alers 1962 ( 10.1063/1.1931149, Table I ); a bulk modulus is not relaxed
         * by anelasticity, so its static and dynamic values coincide. Young's modulus,
         * and with it the shear modulus and the Poisson ratio, is anchored on the
         * static-type curve of Blanke 1989 rather than on the dynamic Hill average of
         * the same crystal: lead's Zener anisotropy is 4.1 and its shear modulus relaxes
         * strongly at low frequency, so the ultrasonic Hill average ( E = 24 GPa at
         * 300 K ) is far above what a slow, elastostatic loading sees ( tensile handbooks:
         * E 16, G 5.6 GPa, nu 0.44 ). With the crystal's K and Blanke's E the served
         * values at 300 K are E 16.25, G 5.66 GPa, nu 0.435, and nu rises monotonically
         * from 0.43 at 95 K. Formula alloys that contain lead inherit this through
         * cl_Material_Alloy. For the other metals of the roster the static and dynamic
         * moduli differ by only a few percent and the dynamic ( ultrasonic ) data are
         * served.
         */
        class Lead : public Metal
        {
            // debye temperature
            Cell< Vector< real > > mDebyePolys ;
            Vector< real > mTDebyeSwitch ;

            Bezier * mThermalExpansion = nullptr ;
            Vector< real > mThermalExpansionCryo ;
            Bezier * mDebyeTemperature = nullptr ;

            Cell< Vector< real > > mKohlerLongPolys ;
            Cell< Vector< real > > mKohlerTransPolys ;
            Vector< real > mBSKohlerLongSwitch ;
            Vector< real > mBSKohlerTransSwitch ;
        public:
            Lead( const real RRR = BELFEM_QUIET_NAN,
                  const bool aBuildTables = true );

            ~Lead() override;

        protected:

            real
            alpha_custom( const real T ) const override ;

            real
            debye_custom(const real T) const override;

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

        inline real
        Lead::kohler( const real B, const real S, const real beta ) const
        {
            if ( B < BELFEM_EPSILON ) return 0.0 ;

            real BS = B * S ;

            real Along ;

            if ( BS < mBSKohlerLongSwitch( 0 ) )
            {
                Along = polyval( mKohlerLongPolys(0), BS );
            }
            else if ( BS < mBSKohlerLongSwitch( 1 ) )
            {
                Along = std::exp( polyval( mKohlerLongPolys(1), std::log(BS) ) );
            }
            else
            {
                Along = polyval( mKohlerLongPolys(2), BS );
            }

            real Atrans ;
            if ( BS < mBSKohlerTransSwitch( 0 ) )
            {
                Atrans = polyval( mKohlerTransPolys(0), BS );
            }
            else
            {
                Atrans = std::exp( polyval( mKohlerTransPolys(1), std::log(BS) ) );
            }

            real c = std::cos( beta );
            real c2 = c*c ;
            real s2 = 1.0 - c2 ;

            // Pippard's angular interpolation formula
            return Along * c2 + Atrans * s2 ;
        }

    }
}
#endif //BELFEM_CL_MATERIAL_LEAD_HPP