/*
 * BELFEM -- The Berkeley Lab Finite Element Framework
 * Copyright (c) 2026, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of any required
 * approvals from the U.S. Dept. of Energy).  All rights reserved.
 *
 * Developers: Christian Messe, Gregory Giard
 *
 * See the top-level LICENSE file for the complete license and disclaimer.
 */

#ifndef BELFEM_CL_MATERIAL_COPPER_HPP
#define BELFEM_CL_MATERIAL_COPPER_HPP

#include "cl_Material_Metal.hpp"
#include "cl_Bezier.hpp"


namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        /**
         * @brief Copper (Cu) material implementation
         *
         * Implements temperature and field-dependent properties for copper
         * using piecewise polynomial fits and empirical models.
         *
         * All properties use experimental data from standard sources
         * (Touloukian, Matula, copper.org, etc.) with smooth interpolation.
         *
         * Magnetoresistance is implemented via Kohler's rule with Bezier
         * curves for longitudinal field and polynomials for transverse field.
         */
        class Copper : public Metal
        {
            // Debye temperature piecewise polynomials
            Vector< real > mTDebyeSwitch ;
            Cell< Vector< real > > mDebyePolys ;

            // Kohler magnetoresistance curves
            Vector< real > mBSKohlerSwitch ;

            Bezier * mKohlerLongBezier = nullptr ;  // Bezier curve for longitudinal field
            Cell< Vector< real > > mKohlerLongPolys ;
            Cell< Vector< real > > mKohlerTransPolys ;

            Bezier * mThermalExpansion = nullptr ;
            Vector< real > mThermalExpansionCryo ;
        public:

            Copper( const real RRR = BELFEM_QUIET_NAN,
                    const bool aBuildTables = true );

            ~Copper() override;

        protected:

            real
            alpha_custom( const real T ) const override ;

            real
            debye_custom( const real T ) const override ;

            real
            kohler(const real B, const real S, const real beta) const override;

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

    }
}
#endif //BELFEM_CL_MATERIAL_COPPER_HPP