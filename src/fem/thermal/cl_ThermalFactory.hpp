//
// Created by gregorygiard on 10/24/25.
//

#ifndef BELFEM_CL_THERMALFACTORY_HPP
#define BELFEM_CL_THERMALFACTORY_HPP

#include "typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Mesh.hpp"
#include "en_DomainType.hpp"
#include "cl_Vector.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_IWG_MaxwellThermal.hpp"
#include "cl_ThermalBoundaryConditionFactory.hpp"
#include "cl_InputFile.hpp"
#include "cl_FEM_Domain.hpp"

namespace belfem
{
    namespace fem
    {
        /**
         * @brief Builds the thermal kernel from an input deck, standalone or coupled to a magnetic kernel.
         *
         * @ingroup grp_fem_thermal
         * @see @ref fem_thermal_index
         */
        class ThermalFactory
        {
            const proc_t mCommRank ;

            const InputFile * mInputFile  ;

            ThermalBoundaryConditionFactory * mBoundaryConditionFactory ;

            // deleted by kernel
            bool mOwnKernelParameters = true ;
            KernelParameters * mKernelParameters = nullptr ;

            Mesh * mMesh = nullptr ;

            ModelDimensionality  mDimensionality = ModelDimensionality::UNDEFINED ;

            std::shared_ptr< Kernel > mThermalKernel = nullptr ;

            IWG_MaxwellThermal * mThermalEquation = nullptr ;
            bool mOwnThermalEquation = true ;

            DofManager * mThermalField = nullptr ;

            Cell< Domain * > mDomains ;



            bool mUseEnrichment = false ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            ThermalFactory( const string & aInputFile, Mesh * aMesh );

            ThermalFactory( const string & aInputFile, Kernel * aMagneticKernel );

            ~ThermalFactory();

//------------------------------------------------------------------------------

            std::shared_ptr< Kernel >
            create_thermal_kernel();

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            configure_solver( const input::Section * aSection, DofManager * aField );

            ThermalBoundaryConditionFactory *
            create_bc_factory();

//------------------------------------------------------------------------------

        };

    }
}

#endif //BELFEM_CL_THERMALFACTORY_HPP