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

            Cell< Vector< real > > mYoungData ;
            Bezier * mYoungBezier = nullptr ;

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

            real
            E_custom( const real T ) const override;

            real
            dEdT_custom( const real T ) const override;
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

            void
            create_young();

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

        inline real
        Nickel::E_custom( const real T ) const
        {
            if ( T < BELFEM_EPSILON ) return mYoungData( 0 )( 0 );
            if ( T < mYoungBezier->basis_x()( 0 ) )
            {
                const Vector< real > & p = mYoungData( 0 );
                return p( 0 ) - p( 1 ) * T * std::exp( -p( 2 ) / T );
            }
            if ( T < mYoungBezier->basis_x()( 3 ) )
            {
                return mYoungBezier->y( T );
            }
            const Vector< real > & p = mYoungData( 1 );
            return p( 0 ) - p( 1 ) * T * std::exp( -p( 2 ) / T );
        }

        inline real
        Nickel::dEdT_custom( const real T ) const
        {
            if ( T < BELFEM_EPSILON ) return 0.0 ;
            if ( T < mYoungBezier->basis_x()( 0 ) )
            {
                const Vector< real > & p = mYoungData( 0 );
                return - p( 1 ) * std::exp( -p( 2 ) / T ) * ( T + p( 2 ) ) / T ;
            }
            if ( T < mYoungBezier->basis_x()( 3 ) )
            {
                return mYoungBezier->dydx( T );
            }
            const Vector< real > & p = mYoungData( 1 );
            return - p( 1 ) * std::exp( -p( 2 ) / T ) * ( T + p( 2 ) ) / T ;
        }

    }
}

#endif //BELFEM_CL_MATERIAL_NICKEL_HPP
