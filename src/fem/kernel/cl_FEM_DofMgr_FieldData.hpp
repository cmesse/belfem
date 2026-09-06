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
                const proc_t mCommRank;

                // number of procs
                const proc_t mCommSize ;

                /**
                 * this list contains the nodes as owned by each proc
                 * exists only on master
                 */
                Cell< Vector< index_t > > mNodeOwnerList ;
                Cell< Vector< index_t > > mElementOwnerList ;

                index_t mMyNumberOfOwnedNodes = 0 ;
                index_t mMyNumberOfOwnedElements = 0 ;

                // special purpose data for linear to higher projection
                Cell< Vector< index_t > > mAllCornerNodeIndices ;
                Vector< index_t > mMyCornerNodeIndices ;


                Cell< Vector< index_t > > mAllNonCornerNodeIndices ;
                Cell< mesh::Node * > mMyNonCornerNodes ;

                Vector< index_t > mMyNonCornerNodeIndices ;

//------------------------------------------------------------------------------
            public:
//------------------------------------------------------------------------------

                FieldData( DofManager * aParent ) ;

//------------------------------------------------------------------------------

                ~FieldData();

//------------------------------------------------------------------------------

                void
                collect_node_owners();

//------------------------------------------------------------------------------

                void
                collect_element_owners();

//------------------------------------------------------------------------------

                void
                collect( const string & aLabel );

//-----------------------------------------------------------------------------

                void
                collect( const Cell< string > & aLabels ) ;

//-----------------------------------------------------------------------------

                // todo: we might not need this anymore
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

                const Cell< index_t > &
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
