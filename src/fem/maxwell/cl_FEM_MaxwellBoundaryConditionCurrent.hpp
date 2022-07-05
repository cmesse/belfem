//
// Created by christian on 9/29/21.
//

#ifndef BELFEM_CL_FEM_MAXWELLBOUNDARYCONDITIONCURRENT_HPP
#define BELFEM_CL_FEM_MAXWELLBOUNDARYCONDITIONCURRENT_HPP

#include "cl_FEM_BoundaryCondition.hpp"

namespace belfem
{
    namespace fem
    {

//-----------------------------------------------------------------------------

        class MaxwellBoundaryConditionCurrent : public BoundaryCondition
        {
            // cross section for coil, must be computed later
            real mCoilCrossSection = 0.0 ;

            // current in ampere
            real mAmplitude = BELFEM_QUIET_NAN ;

            Cell< string > mCurrentFields ;

            // the function that is used to compute the current
            real
            ( MaxwellBoundaryConditionCurrent::*mFunCurrent )
            (  const real aTime ) ;


//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            MaxwellBoundaryConditionCurrent(
                    const BoundaryConditionImposing aImposing,
                    const string          & aLabel );

//-----------------------------------------------------------------------------

            ~MaxwellBoundaryConditionCurrent() = default ;

//-----------------------------------------------------------------------------

            void
            compute( const real aTime );

//-----------------------------------------------------------------------------

            void
            set(
                    const BoundaryConditionScaling aType,
                    const real        aAmplitude,
                    const real        aPeriod=BELFEM_REAL_MAX,
                    const real        aPhase=0.0 );

//-----------------------------------------------------------------------------

            void
            set_dof_manager( DofManager * aDofManager );

//-----------------------------------------------------------------------------
        private:
//-----------------------------------------------------------------------------

            void
            link_blocks( DofManager * aDofManager , const Vector< id_t > & aIDs );

//-----------------------------------------------------------------------------

            void
            compute_cross_section( DofManager * aDofManager );

//-----------------------------------------------------------------------------

            void
            compute_current_densities_2d( const real aCurrent );

//-----------------------------------------------------------------------------

            void
            compute_cut_currents( const real aCurrent );

//-----------------------------------------------------------------------------
        };

//-----------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace belfem */
#endif //BELFEM_CL_FEM_MAXWELLBOUNDARYCONDITIONCURRENT_HPP
