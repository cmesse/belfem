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

#ifndef BELFEM_CL_MATERIAL_TIN_HPP
#define BELFEM_CL_MATERIAL_TIN_HPP
#include "cl_Material_Metal.hpp"
#include "cl_Bezier.hpp"
#include "fn_polyval.hpp"
namespace belfem
{
    namespace material
    {
        class WhiteTin : public Metal
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
            WhiteTin( const real RRR = BELFEM_QUIET_NAN,
                      const bool aBuildTables = true );
            ~WhiteTin() override;

        protected:

            real
            alpha_custom( const real T ) const override ;

            real
            debye_custom( const real T ) const override ;

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
        WhiteTin::kohler( const real B, const real S, const real beta ) const
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
#endif //BELFEM_CL_MATERIAL_TIN_HPP