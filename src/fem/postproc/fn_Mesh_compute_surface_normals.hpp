//
// Created by Christian Messe on 14.06.20.
//
#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "Mesh_Enums.hpp"

#ifndef BELFEM_FN_MESH_COMPUTE_SURFACE_NORMALS_HPP
#define BELFEM_FN_MESH_COMPUTE_SURFACE_NORMALS_HPP

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        void
        compute_surface_normals(
                Mesh * aMesh,
                const Vector< id_t > & aGroupIDs,
                const GroupType aGroupType=GroupType::SIDESET
                );

//------------------------------------------------------------------------------

        namespace normals
        {
//------------------------------------------------------------------------------

            void
            compute_surface_normals_2d(
                    Mesh * aMesh,
                    const Vector< id_t > & aGroupIDs,
                    const GroupType aGroupType
            );

//------------------------------------------------------------------------------
            void
            compute_surface_normals_3d(
                    Mesh * aMesh,
                    const Vector< id_t > & aGroupIDs,
                    const GroupType aGroupType
            );

//------------------------------------------------------------------------------

            void
            create_normal_fields( Mesh * aMesh );

//------------------------------------------------------------------------------

            void
            collect_elements(
                    Mesh * aMesh,
                    const id_t aGroupID,
                    const GroupType aGroupType,
                    Cell< Element * > & aElements );

//------------------------------------------------------------------------------

            void
            send_surface_normals_2d( Mesh * aMesh  ) ;

//------------------------------------------------------------------------------

            void
            send_surface_normals_3d( Mesh * aMesh  ) ;

//------------------------------------------------------------------------------

            void
            receive_surface_normals_2d( Mesh * aMesh,
                                        const Vector< id_t > & aGroupIDs,
                                        const GroupType aGroupType ) ;

//------------------------------------------------------------------------------

            void
            receive_surface_normals_3d( Mesh * aMesh,
                        const Vector< id_t > & aGroupIDs,
                        const GroupType aGroupType ) ;

//------------------------------------------------------------------------------

            void
            create_commlist( Vector< proc_t > & aCommList );

//------------------------------------------------------------------------------

            void
            get_node_ids( Mesh * aMesh,
                          const Vector< id_t > & aGroupIDs,
                          const GroupType aGroupType,
                          Vector< id_t > & aNodeIDs);

//------------------------------------------------------------------------------
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_MESH_COMPUTE_SURFACE_NORMALS_HPP
