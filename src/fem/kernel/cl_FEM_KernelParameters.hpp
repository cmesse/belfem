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

#ifndef BELFEM_CL_FEM_KERNELPARAMETERS_HPP
#define BELFEM_CL_FEM_KERNELPARAMETERS_HPP

#include "typedefs.hpp"
#include "cl_Mesh.hpp"
#include "en_IntegrationScheme.hpp"
#include "en_SolverEnums.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class Kernel ;

        class KernelParameters
        {
            // rank of this proc
            const proc_t mCommRank;
            const proc_t mCommSize;

            Kernel * mKernel ;

            // the mesh that is associated with this kernel
            Mesh * mMesh  ;

            // List of Blocks that contribute to this Kernel ( Block Indices ) ( default: all )
            Vector< index_t > mBlockIndices;

            // List of Sidesets that contribute to this Kernel ( Side Indices ) ( default: none )
            Vector< index_t > mSidesetIndices ;

            // ids of blocks that contribute to this kernel ( default : none -- set by select_blocks(); only the indices above default to all )
            Vector< id_t > mBlockIDs ;

            Vector< id_t > mSideSetIDs ;

            const uint mZero = 0 ;

            // default integration scheme for quad and hex elements
            IntegrationScheme mIntegrationScheme = IntegrationScheme::GAUSS ;

            // must be set to SCOTCH or METIS
            ReorderingMethod mReorderingMethod = ReorderingMethod::METIS ;

            // flag telling if we use metis to create the partitioning
            bool mAutoPartition = true ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            KernelParameters( Mesh & aMesh );

//------------------------------------------------------------------------------

            KernelParameters( Mesh * aMesh );

//------------------------------------------------------------------------------

            KernelParameters( Kernel * aKernel );

//------------------------------------------------------------------------------

            ~KernelParameters() = default;

//------------------------------------------------------------------------------

            /**
             * return the pointer of the associated mesh
             */
            Mesh *
            mesh();

//------------------------------------------------------------------------------

            /**
             * return the number of procs that contribute to this kernel
             */
            proc_t
            number_of_procs() const;

//------------------------------------------------------------------------------

            /**
             * return the id of the master proc
             */
            const proc_t &
            master();

//------------------------------------------------------------------------------

            /**
             * Select the IDs of the blocks that are to be used.
             * Only the Master proc that owns the mesh is responsible for this
             *
             * in this case, we use these blocks for the partitioning
             */
            void
            select_blocks( const Vector< id_t > & aBlockIDs );

            void
            select_blocks( const Cell< id_t > & aBlockIDs );

            void
            select_sidesets( const Vector< id_t > & aSidesetIDs );

            void
            select_sidesets( const Cell< id_t > & aSidesetIDs );

//------------------------------------------------------------------------------

            /**
             * return the indices of the selected blocks
             */
            const Vector< index_t > &
            block_indices() const ;

//------------------------------------------------------------------------------

            /**
             * return the ids of the selected blocks
             */
            const Vector< id_t > &
            selected_blocks() const ;

 //------------------------------------------------------------------------------

            /**
             * return the indices of the selected sidesets
             */
            const Vector< index_t > &
            sideset_indices() const ;

//------------------------------------------------------------------------------

            /**
             * return the ids of the selected sidesets
             */
            const Vector< id_t > &
            selected_sidesets() const ;

//------------------------------------------------------------------------------

            void
            set_integration_scheme( const IntegrationScheme & aIntegrationScheme ) ;

//------------------------------------------------------------------------------

            const IntegrationScheme &
            integration_scheme() const ;

//------------------------------------------------------------------------------

            bool
            auto_partition() const ;

//------------------------------------------------------------------------------

            void
            set_auto_partition( const bool aFlag );

//------------------------------------------------------------------------------

            Kernel *
            kernel() ;

//------------------------------------------------------------------------------

            ReorderingMethod
            reordering_method() const ;

            void
            set_reordering_method( const ReorderingMethod aMethod );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------
            
            /**
             * if the user does not modify anything, some default values
             * are initialized that will work for many cases
             */
            void
            init_defaults();
            
//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------


        inline proc_t
        KernelParameters::number_of_procs() const
        {
            return mCommSize;
        }

//------------------------------------------------------------------------------

        inline const proc_t &
        KernelParameters::master()
        {
            return mMesh->master();
        }

//------------------------------------------------------------------------------

        inline const Vector< index_t > &
        KernelParameters::block_indices() const
        {
            return mBlockIndices ;
        }

//------------------------------------------------------------------------------

        inline const Vector< id_t > &
        KernelParameters::selected_blocks() const
        {
            return mBlockIDs ;
        }

        inline const Vector< index_t > &
        KernelParameters::sideset_indices() const
        {
            return mSidesetIndices ;
        }

//------------------------------------------------------------------------------

        inline const Vector< id_t > &
        KernelParameters::selected_sidesets() const
        {
            return mSideSetIDs ;
        }

//------------------------------------------------------------------------------

        inline const IntegrationScheme &
        KernelParameters::integration_scheme() const
        {
            return mIntegrationScheme ;
        }

//------------------------------------------------------------------------------

        inline bool
        KernelParameters::auto_partition() const
        {
            return mAutoPartition ;
        }

//------------------------------------------------------------------------------

        inline ReorderingMethod
        KernelParameters::reordering_method() const
        {
            return mReorderingMethod ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_FEM_KERNELPARAMETERS_HPP
