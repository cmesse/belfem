//
// Created by christian on 7/14/21.
//

#ifndef BELFEM_CL_FEM_DOFMGR_FIELDDATA_HPP
#define BELFEM_CL_FEM_DOFMGR_FIELDDATA_HPP

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"
#include "cl_Vector.hpp"
#include "cl_IWG.hpp"
#include "cl_FEM_Dof.hpp"
namespace belfem
{
    class Mesh;

    namespace fem
    {
        class Bearing;

        class Kernel;

        class DofManager;

        namespace dofmgr
        {
//------------------------------------------------------------------------------

            class FieldData
            {
                //! the parent object
                DofManager * mParent;

                //! the kernel
                Kernel * mKernel;

                //! the mesh this problem runs on
                Mesh * mMesh;

                // my rank
                const proc_t mMyRank;

                // number of procs
                const proc_t mCommSize ;

                /**
                 * this list contains the nodes as owned by each proc
                 * exists only on master
                 */
                Cell< Vector< index_t > > mNodeOwnerList ;
                Cell< Vector< index_t > > mGhostElementOwnerList ;

                index_t mMyNumberOfOwnedNodes = 0 ;
                index_t mMyNumberOfOwnedGhostElements = 0 ;

                // special purpose data for linear to higher projection
                Cell< Vector< index_t > > mAllCornerNodeIndices ;
                Vector< index_t > mMyCornerNodeIndices ;


                Cell< Vector< index_t > > mAllNonCornerNodeIndices ;
                Cell< mesh::Node * > mMyNonCornerNodes ;

                Vector< index_t > mMyNonCornerNodeIndices ;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                FieldData(
                        DofManager * aParent ) ;

//------------------------------------------------------------------------------

                ~FieldData();

//------------------------------------------------------------------------------

                void
                collect_node_owners();

//------------------------------------------------------------------------------

                void
                collect_ghost_element_owners();

//------------------------------------------------------------------------------

                /**
                 * field must be either
                 * @param aLabel
                 */
                void
                collect( const string & aLabel );

//-----------------------------------------------------------------------------

                void
                collect( const Cell< string > & aLabels ) ;

//-----------------------------------------------------------------------------

                void
                initialize_linear_projection_lists();

//------------------------------------------------------------------------------

                void
                project_linear_field_to_higher_mesh(
                        const Cell< string > & aFieldLabels ) ;

//------------------------------------------------------------------------------

                /**
                 * distributes the field data from the master to the others
                 */
                void
                distribute( const Cell< string > & aFieldLabels );

//------------------------------------------------------------------------------

                void
                reset();

//------------------------------------------------------------------------------

                void
                update_field_indices( Cell< Dof * > & aDOFs );

//------------------------------------------------------------------------------
            private:
//------------------------------------------------------------------------------

                void
                communicate_corner_node_data(
                        const Cell< string > & aFieldLabels );

//------------------------------------------------------------------------------

                void
                communicate_noncorner_node_data(
                        const Cell< string > & aFieldLabels,
                        Matrix< real > & aData ) ;


//------------------------------------------------------------------------------

                const Vector< index_t > &
                field_indices( const EntityType aType,
                               const proc_t     aTarget,
                               uint & aMultiplicity );

//------------------------------------------------------------------------------

            };
//------------------------------------------------------------------------------
        } /* end namespace dofmgr */
    } /* end namespace fem */
} /* end namespace belfem */


#endif //BELFEM_CL_FEM_DOFMGR_FIELDDATA_HPP
