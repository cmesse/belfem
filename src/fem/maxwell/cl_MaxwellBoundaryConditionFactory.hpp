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

#ifndef BELFEM_CL_MAXWELLBOUNDARYCONDITIONFACTORY_HPP
#define BELFEM_CL_MAXWELLBOUNDARYCONDITIONFACTORY_HPP


#include "cl_Map.hpp"

#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_PhysicalBoundaryCondition.hpp"
#include "cl_IWG_Maxwell.hpp"
#include "cl_SourceFunction.hpp"
#include "stringtools.hpp"
#include "typedefs.hpp"

namespace belfem
{
    namespace fem
    {
        /**
         * @brief Builds the physical boundary conditions of a Maxwell problem.
         *
         * @ingroup grp_fem_maxwell
         * @see @ref fem_maxwell_maxwell_usage_guide
         */
        class MaxwellBoundaryConditionFactory
        {
            const uint mNumDimensions ;
            Cell < PhysicalBoundaryCondition * > mPhysicalBoundaryConditions ; //the boundary conditions are actually owned by the Kernel

            //! total occurrences of each global-variable base name in the
            //! deck ( pre-pass ) and how many the walk has seen so far: a
            //! base that repeats carries its section ordinal on EVERY
            //! occurrence, so section and group suffixes cannot collide
            Map< string, uint > mGlobalNameCount ;
            Map< string, uint > mGlobalNameSeen ;


        public:

            MaxwellBoundaryConditionFactory( const input::Section * aSection, const uint aNumberOfDimensions ) ;

            ~MaxwellBoundaryConditionFactory() = default ;

//-------------------------------------------------------------------------------

            Cell < PhysicalBoundaryCondition * > &
            boundary_conditions() ;

//-------------------------------------------------------------------------------

            Cell < PhysicalBoundaryCondition * >
            current_boundary_conditions() ;

//-------------------------------------------------------------------------------

            Cell < PhysicalBoundaryCondition * >
            voltage_boundary_conditions() ;

//-------------------------------------------------------------------------------

            void
            set_fields(DofManager * aField) ;

//-------------------------------------------------------------------------------

        } ;
    }
}

#endif //BELFEM_CL_MAXWELLBOUNDARYCONDITIONFACTORY_HPP
