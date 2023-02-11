//
// Created by christian on 2/9/23.
//

#ifndef BELFEM_CL_MAXWELLMESHSYNCH_HPP
#define BELFEM_CL_MAXWELLMESHSYNCH_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"

namespace belfem
{
    class MaxwellMeshSynch
    {
        const proc_t mRank ;

        Mesh & mMaxwellMesh ;
        Mesh & mThermalMesh ;
        Vector< index_t > mTable ;

//----------------------------------------------------------------------------
    public:
//----------------------------------------------------------------------------

        MaxwellMeshSynch( Mesh * aMaxwellMesh, Mesh * aThermalMesh );

//----------------------------------------------------------------------------

        ~MaxwellMeshSynch() = default ;

//----------------------------------------------------------------------------

        // pass ej and b
        void
        maxwell_to_thermal();

//----------------------------------------------------------------------------

        // pass T
        void
        thermal_to_maxwell();

//--------------------------------------------------------------------------

        void
        initialize_tables();
//----------------------------------------------------------------------------
    private:
//----------------------------------------------------------------------------

        // defines proc ownerships
        void
        set_ownerships();

//----------------------------------------------------------------------------
    };
}
#endif //BELFEM_CL_MAXWELLMESHSYNCH_HPP
