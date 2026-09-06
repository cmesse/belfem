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

#ifndef BELFEM_CL_MATERIAL_MAGNESIA_HPP
#define BELFEM_CL_MATERIAL_MAGNESIA_HPP

#include "cl_Material_Metal.hpp"
#include "cl_Bezier.hpp"
#include "fn_polyval.hpp"

namespace belfem
{
    namespace material
    {
        class Magnesia : public SplineLookupTable
        {

            // Young's modulus and Poisson's ratio polynomials
            Cell< Vector< real > > mYoungPolys ;
            Cell< Vector< real > > mPoissonPolys ;
            Cell< Vector< real > > mCpPolys ;
            Cell< Vector< real > > mLambdaPolys ;
            real mTYoungSwitch ;
            real mTPoissonSwitch ;
            Vector< real > mTCpSwitch ;
            Vector< real > mTLambdaSwitch ;
            Bezier * mThermalExpansion = nullptr ;
            Vector< real > mThermalExpansionCryo ;

        public:
            Magnesia();
             ~Magnesia() override;
        
        protected:

            real
            E_custom( const real T ) const override ;

            real
            nu_custom( const real T ) const override ;

            real
            alpha_custom( const real T ) const override ;

            real
            cp_custom( const real T ) const override ;

            real
            lambda_custom( const real T ) const override ;

        private:

            void
            set_constants();

            void
            create_alpha();

            void
            create_mech();

            void
            create_cp();

            void
            create_lambda();

        };

        inline real
        Magnesia::E_custom( const real T ) const
        {
            if ( T < mTYoungSwitch )
            {
                return polyval( mYoungPolys( 0 ), T*T );
            }
            else
            {
                return polyval( mYoungPolys( 1 ), T );
            }
        }

        inline real
        Magnesia::nu_custom( const real T ) const
        {
            if ( T < mTPoissonSwitch )
            {
                return polyval( mPoissonPolys( 0 ), T*T );
            }
            else
            {
                return polyval( mPoissonPolys( 1 ), T );
            }
        }

    }
}
#endif //BELFEM_CL_MATERIAL_MAGNESIA_HPP
