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

#ifndef BELFEM_CL_FEM_KERNEL_HPP
#define BELFEM_CL_FEM_KERNEL_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Mesh.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"
#include "cl_Material.hpp"

#include "cl_FEM_Dof.hpp"
#include "cl_IF_InterpolationFunction.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_FEM_KernelParameters.hpp"
#include "cl_FEM_PhysicalBoundaryCondition.hpp"

#include "cl_FEM_DofManager.hpp"
#include "cl_CommTable.hpp"
#include "en_IWGs.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    namespace fem
    {
        class IwgFactory ;

        class Controller ;

        /**
         * @brief Top-level orchestrator; owns the mesh, materials, boundary conditions and DOF managers.
         *
         * @ingroup grp_fem_kernel
         * @see @ref fem_kernel_index
         */
        class Kernel
        {
            // rank of this proc
            const proc_t     mCommRank;

            // how many procs contribute to this kernel
            const proc_t     mCommSize;

            // pointer to parameter list
            KernelParameters * mParams;

            // pointer to mesh (points to empty mesh unless root kernel)
            Mesh             * mMesh;

            // offset for fields, in case fields already exist on mesh
            // when Kernel is generated
            const uint mFieldOffset ;

            // submesh ( points to local mesh unless this is root)
            Mesh * mSubMesh = nullptr;

            Controller       * mController = nullptr ;


            // flag telling if the kernel destroys the parameters on exit
            bool mOwnParameters = false ;

            proc_t           mMyCommIndex = 0 ;

            Vector< index_t > mNumElementsPerBlock;

            Cell< DofManager * > mDofManagers ;
            Cell< IWG * >        mIWGs ;

            Cell< Material * > mMaterials;
            Map< string, Material * > mMaterialMap;

            Cell< PhysicalBoundaryCondition * > mBoundaryConditions ;

            Cell< const mesh::CommTable * > mCommTables ;

            bool mOwnCommTables = false ;
            bool mOwnMesh = false ;
            bool mOwnSubmesh = false ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Kernel( KernelParameters * aKernelParameters );

//------------------------------------------------------------------------------

            ~Kernel();

//------------------------------------------------------------------------------

            /**
             * expose parameter object
             */
             const  KernelParameters *
             params();

//------------------------------------------------------------------------------

            /**
             * expose a field
             */
            DofManager *
            dofmgr( const uint aIndex=0 );

//------------------------------------------------------------------------------

            /**
             * returns the material object registered under aLabel
             * ( add_material ). Hard error if no such material exists.
             * Called from Group::set_material( label )
             */
             Material *
             material( const string & aLabel );

//------------------------------------------------------------------------------

            /**
             * manually add a material to the kernel. The material is now owned
             * and destroyed by the kernel
             */
            void
            add_material( const string aLabel, Material * aMaterial );

//------------------------------------------------------------------------------

            /**
             * manually add the BCs to the kernel. The BC is now owned
             * and destroyed by the kernel
             */
            void
            add_boundary_condition( PhysicalBoundaryCondition * aBC );

//------------------------------------------------------------------------------

            /**
             *expose the boundary condition container
             */
            Cell< PhysicalBoundaryCondition * > &
            boundary_conditions();

//------------------------------------------------------------------------------

            /**
             *expose a boundary condition for a specific sideset ID
             */
            PhysicalBoundaryCondition *
            boundary_condition( id_t aID );

//------------------------------------------------------------------------------

            /**
             * Compute the boundary condition and update the Dofs at a given time
             */
            void
            compute_boundary_conditions( real aTime );

//------------------------------------------------------------------------------


             const mesh::CommTable *
             comm_table( const uint aProc ) const ;

             const Cell< const mesh::CommTable * > &
             comm_tables() const ;

//------------------------------------------------------------------------------

             const proc_t &
             number_of_procs() const ;

//------------------------------------------------------------------------------

            /**
             * returns the number of fields that existed before Kernel was created
             */
            const uint &
            field_offset() const ;

//------------------------------------------------------------------------------

            /**
             * returns the mesh on the master and the submesh on other
             *
             * @return the mesh (master) or the submesh (all other ranks)
             */
            Mesh *
            mesh();

//------------------------------------------------------------------------------

            IWG *
            create_equation( const IwgType aEquationType,
                const        ModelDimensionality aModelDimensionality = ModelDimensionality::UNDEFINED,
                             const Vector< id_t > aBlocks={},
                             const Vector< id_t > aSideSets={} );

//------------------------------------------------------------------------------

            /**
             * special function to add an equation that has already
             * been created. The equation is now owned and destroyed
             * by the kernel
             */
            void
            add_equation( IWG * aEquation );

//------------------------------------------------------------------------------

            DofManager *
            create_field( IWG * aEquation );

//------------------------------------------------------------------------------

            bool
            is_master() const ;

//------------------------------------------------------------------------------

            /**
             * if this flag is set, the kernel will destroy the parameters
             */
             void
             claim_parameter_ownership( const bool aFlag = true );

//------------------------------------------------------------------------------

            void
            set_controller( Controller * aController );

//------------------------------------------------------------------------------

            Controller *
            controller() ;

//------------------------------------------------------------------------------

            /**
             * true once set_controller() has attached a Controller;
             * callers that may legitimately run without one ask this
             * before calling controller()
             */
            bool
            has_controller() const ;

//------------------------------------------------------------------------------

            // make sure that all elements have positive volume
            void
            compute_element_volumes();

//------------------------------------------------------------------------------

            uint
            number_of_dof_managers() const ;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            /**
             * count the number of elements on each block plus aura
             */
            void
            distribute_mesh();

//------------------------------------------------------------------------------

            /**
             * creates the dof connectivity graph that is used for the
             * dof numeration
             * @param aGraph
             */
            void
            create_dof_graph( Cell< Dof * > & aGraph );

//------------------------------------------------------------------------------

            /**
             * help function for parallel edge and face creation.
             * this function grabs elements based on gived blocks and
             * sidesets
             *
             * @param aBlockIDs
             * @param aSideSetIDs
             * @param aElements
             */
            void
            collect_elements(
                    const Vector< id_t >    & aBlockIDs,
                    const Vector< id_t >    & aSideSetIDs,
                    Cell< mesh::Element * > & aElements );


//------------------------------------------------------------------------------

            mesh::Basis *
            get_entity( const id_t aID, const EntityType aType ) ;

//------------------------------------------------------------------------------

            void
            partition_mesh();

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        inline const mesh::CommTable *
        Kernel::comm_table( const uint aProc ) const
        {
            BELFEM_ASSERT( mCommRank == 0, "Only master may call the comm table." );
            return mCommTables( aProc );
        }

        inline const Cell< const mesh::CommTable * > &
        Kernel::comm_tables() const
        {
            return mCommTables ;
        }

        inline const proc_t &
        Kernel::number_of_procs() const
        {
            return mCommSize ;
        }

//------------------------------------------------------------------------------

        inline const uint &
        Kernel::field_offset() const
        {
            return mFieldOffset ;
        }

//------------------------------------------------------------------------------

        inline Cell< PhysicalBoundaryCondition * > &
        Kernel::boundary_conditions()
        {
            return mBoundaryConditions ;
        }

//------------------------------------------------------------------------------

        inline PhysicalBoundaryCondition *
        Kernel::boundary_condition( id_t aID )
        {
            for ( PhysicalBoundaryCondition* tBC : mBoundaryConditions )
            {
                for (id_t tID : tBC->domains())
                {
                    if (tID == aID)
                    {
                        return tBC ;
                    }
                }
            }
            BELFEM_ERROR(false, "Sideset is not associated to a boundary contition") ;
            return nullptr ;
        }

//------------------------------------------------------------------------------

        inline Mesh *
        Kernel::mesh()
        {
            return mCommRank == 0 ? mMesh : mSubMesh ;
        }

//------------------------------------------------------------------------------

        inline uint
        Kernel::number_of_dof_managers() const
        {
            return mDofManagers.size();
        }

//------------------------------------------------------------------------------

        inline bool
        Kernel::is_master() const
        {
            return mCommRank == 0 ;
        }

//------------------------------------------------------------------------------

        inline bool
        Kernel::has_controller() const
        {
            return mController != nullptr ;
        }

//------------------------------------------------------------------------------

        inline mesh::Basis *
        Kernel::get_entity( const id_t aID, const EntityType aType )
        {
            switch( aType )
            {
                case( EntityType::NODE ) :
                {
                    return this->mesh()->node( aID );
                }
                case( EntityType::EDGE ) :
                {
                    return this->mesh()->edge( aID );
                }
                case( EntityType::FACE ) :
                {
                    return this->mesh()->face( aID );
                }
                case( EntityType::FACET ) :
                {
                    return this->mesh()->facet( aID );
                }
                case( EntityType::ELEMENT ) :
                {
                    return this->mesh()->element( aID );
                }
                default:
                {
                    BELFEM_ERROR( false, "Invalid entity type");
                    return nullptr ;
                }
            }
        }

    }
//------------------------------------------------------------------------------
}
#endif //BELFEM_CL_FEM_KERNEL_HPP
