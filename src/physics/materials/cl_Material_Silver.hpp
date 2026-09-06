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

#ifndef BELFEM_CL_MATERIAL_SILVER_HPP
#define BELFEM_CL_MATERIAL_SILVER_HPP
#include "typedefs.hpp"
#include "cl_Material_Metal.hpp"
#include "fn_polyval.hpp"
#include "cl_Bezier.hpp"
namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        /**
         * @brief Silver (Ag) material implementation
         *
         * Implements temperature and field-dependent properties for silver
         * using polynomial fits and empirical models.
         *
         * Note: High-quality thermal conductivity data for silver with
         * well-defined RRR values is surprisingly sparse in the literature.
         *
         * Magnetoresistance is implemented via Kohler's rule with polynomial
         * fits for both longitudinal and transverse field orientations.
         */
        class Silver : public Metal
        {
            // dL/L Bezier for the thermal expansion above the split temperature
            Bezier * mThermalExpansion = nullptr ;
            Vector< real > mThermalExpansionCryo ;

            // Debye temperature piecewise polynomials
            Vector< real > mTDebyeSwitch ;
            Cell< Vector< real > > mDebyePolys ;

            // Kohler magnetoresistance curves
            Vector< real >         mBSKohlerSwitch ;
            Cell< Vector< real > > mKohlerLongPolys ;
            Cell< Vector< real > > mKohlerTransPolys ;

        public:

            Silver( const real RRR = BELFEM_QUIET_NAN,
                    const bool aBuildTables = true );
            ~Silver() override ;

        protected:

            real
            alpha_custom( const real T ) const override ;

            real
            debye_custom(const real T) const override;

            real
            kohler( const real B, const real S, const real beta ) const override;

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

#endif //BELFEM_CL_MATERIAL_SILVER_HPP