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

#ifndef BELFEM_CL_MATERIAL_ALUMINUM_HPP
#define BELFEM_CL_MATERIAL_ALUMINUM_HPP

#include "cl_Material_Metal.hpp"
#include "cl_Bezier.hpp"
#include "fn_polyval.hpp"

namespace belfem
{
    namespace material
    {
//------------------------------------------------------------------------------

        /**
         * @brief Aluminum ( Al ), face centered cubic
         *
         * Registered in cl_MaterialFactory as "aluminum", "aluminium", "al".
         * Heat capacity and thermal expansion are fitted against the Touloukian
         * dataset, the effective Debye temperature curve against Desai 1984 and
         * Cook 1975, the resistivity anchor against Hust 1984, the Kohler
         * magnetoresistance against Lüthi 1960 and Fickett 1972, and the
         * elastic moduli as a Wachtman curve against Blanke 1989 with the
         * Poisson ratio derived from a constant Grueneisen parameter
         * ( Metal::create_mech ). The implied Grueneisen parameter of the
         * cryogenic expansion branch reads 2.19 at the split temperature,
         * against a literature 2.1 - 2.2 for aluminum.
         */
        class Aluminum : public Metal
        {
            Bezier * mThermalExpansion = nullptr ;
            Vector< real > mThermalExpansionCryo ;
            Cell< Bezier * > mDebyeBeziers ;

            Vector< real > mBSKohlerSwitch ;
            Cell< Vector< real > > mKohlerLongPolys ;
            Cell< Vector< real > > mKohlerTransPolys ;


        public:

            Aluminum( const real RRR = BELFEM_QUIET_NAN,
                      const bool aBuildTables = true );

            ~Aluminum() override ;

        protected:

            real
            alpha_custom( const real T ) const override ;

            real
            debye_custom(const real T) const override ;

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

//------------------------------------------------------------------------------

    }
}
#endif //BELFEM_CL_MATERIAL_ALUMINUM_HPP
