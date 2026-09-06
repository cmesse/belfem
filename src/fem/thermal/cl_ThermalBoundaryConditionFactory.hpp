//
// Created by gregorygiard on 10/24/25.
//

#ifndef BELFEM_CL_THERMALBOUNDARYCONDITIONFACTORY_HPP
#define BELFEM_CL_THERMALBOUNDARYCONDITIONFACTORY_HPP

#include "cl_Map.hpp"
#include "cl_FEM_PhysicalBoundaryCondition.hpp"
#include "cl_SourceFunction.hpp"
#include "stringtools.hpp"
#include "typedefs.hpp"

namespace belfem
{
    namespace fem
    {
        /**
         * @brief Builds the thermal boundary conditions from the input deck.
         *
         * @ingroup grp_fem_thermal
         * @see @ref fem_thermal_index
         */
        class ThermalBoundaryConditionFactory
        {
            //! total and seen occurrences of each global-variable base
            //! name, cf. the Maxwell twin
            Map< string, uint > mGlobalNameCount ;
            Map< string, uint > mGlobalNameSeen ;


            Cell < PhysicalBoundaryCondition * > mPhysicalBoundaryConditions ; //the boundary conditions are actually owned by the Kernel

        public:

            ThermalBoundaryConditionFactory() ;

            ThermalBoundaryConditionFactory( const input::Section * aSection ) ;

            ~ThermalBoundaryConditionFactory() = default ;

//-------------------------------------------------------------------------------

            Cell < PhysicalBoundaryCondition * > &
            boundary_conditions() ;

//-------------------------------------------------------------------------------

            void
            set_fields(DofManager * aField) ;

//-------------------------------------------------------------------------------

        } ;
    }
}

#endif //BELFEM_CL_THERMALBOUNDARYCONDITIONFACTORY_HPP