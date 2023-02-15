//
// Created by christian on 2/9/23.
//

#ifndef BELFEM_CL_IWG_MAXWELL_THERMAL2D_HPP
#define BELFEM_CL_IWG_MAXWELL_THERMAL2D_HPP

#include "cl_IWG_Timestep.hpp"
#include "cl_Material.hpp"

namespace belfem
{
    namespace fem
    {
        class IWG_Maxwell_Thermal2D : public IWG_Timestep
        {
            Cell< Material * > mMaterials ;
            Vector< real >     mTapeThickness ;
            Vector< real >     mLength ;
            Vector< uint >     mLayer ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IWG_Maxwell_Thermal2D();

//------------------------------------------------------------------------------

            ~IWG_Maxwell_Thermal2D();

//------------------------------------------------------------------------------

            void
            set_material_table( Cell< Material * > & aTapeMaterials );

//------------------------------------------------------------------------------

            void
            compute_geometry_data( Mesh * aMaxwellMesh,
                                   Mesh * aThermalMesh,
                                   const Vector< real > & aTapeThicknesses );

//------------------------------------------------------------------------------

            void
            compute_jacobian_and_rhs(
                    Element        * aElement,
                    Matrix< real > & aJacobian,
                    Vector< real > & aRHS );

//------------------------------------------------------------------------------

            void
            allocate_work_matrices( Group * aGroup ) ;

//------------------------------------------------------------------------------
        };
    }
}

#endif //BELFEM_CL_IWG_MAXWELL_THERMAL2D_HPP
