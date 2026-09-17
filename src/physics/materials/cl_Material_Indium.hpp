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

#ifndef BELFEM_CL_MATERIAL_INDIUM_HPP
#define BELFEM_CL_MATERIAL_INDIUM_HPP


#include "cl_Material_Metal.hpp"
#include "cl_Bezier.hpp"
#include "fn_polyval.hpp"

namespace belfem
{
    namespace material
    {
        class Indium : public Metal
        {
            Bezier * mThermalExpansion = nullptr ;
            Vector< real > mThermalExpansionCryo ;

            Vector< real > mTDebyeSwitch ;
            Cell< Vector< real > > mDebyePolys ;
            Bezier * mDebyeBezier = nullptr ;

            Vector< real > mBSKohlerSwitch ;

            Bezier * mKohlerLongBezier = nullptr ;
            Cell< Vector< real > > mKohlerLongPolys ;
            Cell< Vector< real > > mKohlerTransCoeffs ;

        public:

            Indium( const real RRR = BELFEM_QUIET_NAN,
            const bool aBuildTables = true );

            ~Indium() override;

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
            create_lambda();

            void
            create_kohler();
        };

    }
}

#endif //BELFEM_CL_MATERIAL_INDIUM_HPP
