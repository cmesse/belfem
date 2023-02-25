//
// Created by christian on 2/9/23.
//

#ifndef BELFEM_CL_MAXWELLMESHSYNCH_HPP
#define BELFEM_CL_MAXWELLMESHSYNCH_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "cl_FEM_KernelParameters.hpp"
#include "cl_FEM_Kernel.hpp"

namespace belfem
{
    namespace fem
    {
//----------------------------------------------------------------------------

        class MaxwellMeshSynch
        {
            const proc_t mRank;
            Kernel & mMagneticKernel ;
            Kernel & mThermalKernel ;
            Mesh   & mMagneticMesh ;
            Mesh   & mThermalMesh ;

            Vector< index_t > mTable;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            MaxwellMeshSynch( Kernel * aMagneticKernel, Kernel * aThermalKernel );

//----------------------------------------------------------------------------

            ~MaxwellMeshSynch() = default;

//----------------------------------------------------------------------------

            // pass ej and b
            void
            magnetic_to_thermal_b_and_ej();

//----------------------------------------------------------------------------

            // pass ej and b
            void
            magnetic_to_thermal_T();

//----------------------------------------------------------------------------

            // pass T
            void
            thermal_to_magnetic_T();

//--------------------------------------------------------------------------

            void
            initialize_tables();

//----------------------------------------------------------------------------
        };
    }
}
#endif //BELFEM_CL_MAXWELLMESHSYNCH_HPP
