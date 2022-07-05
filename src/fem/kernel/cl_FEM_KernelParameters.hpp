//
// Created by Christian Messe on 29.10.19.
//

#ifndef BELFEM_CL_FEM_KERNELPARAMETERS_HPP
#define BELFEM_CL_FEM_KERNELPARAMETERS_HPP

#include "typedefs.hpp"
#include "cl_Mesh.hpp"
#include "en_IntegrationScheme.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class KernelParameters
        {
            // rank of this proc
            const proc_t mMyRank;

            // the mesh that is associated with this kernel
            Mesh * mMesh;

            // List of Procs that contribute to this Kernel ( proc IDs, default: all )
            Vector< proc_t >  mSelectedProcs;

            proc_t mNumberOfProcs;

            // number of selected procs, legacy!
            // List of Blocks that contribute to this Kernel ( Block Indices ) ( default: all )
            Vector< index_t > mBlockIndices;

            // dimensions for each field, legacy
            Vector< uint > mNumDofsPerNode = { 1 };
            Vector< uint > mNumDofsPerEdge = { 0 };

            // integration orders for each field ( 0: auto )
            Vector< uint > mBlockIntegrationOrders = { 0 };
            Vector< uint > mSideSetIntegrationOrders = { 0 };
            const uint mZero = 0 ;

            // default integration scheme for quad and hex elements
            IntegrationScheme mIntegrationScheme = IntegrationScheme::GAUSS ;

            // flags telling if one field is enforced as linear
            Vector< uint > mEnforceLinear = { 0 };

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            KernelParameters( Mesh & aMesh );

//------------------------------------------------------------------------------

            KernelParameters( Mesh * aMesh );

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
             * select the procs that contribute to this kernel
             * EXPERIMENTAL!
             */
            void
            select_procs( const Vector < proc_t > & aProcList );

//------------------------------------------------------------------------------

            /**
             * return the number of DOFs specified per node, legacy
             */
            const Vector< uint >&
            num_dofs_per_node() const;

//------------------------------------------------------------------------------

            /**
             * return the number of DOFs specified per edge, legacy
             */
            const Vector< uint >&
            num_dofs_per_edge() const ;

//------------------------------------------------------------------------------

            /**
             * return the flags for linear enforcement if they were given
             */
            const Vector< uint >&
            linear_enforcement_flags() const;

//------------------------------------------------------------------------------

            /**
             * return the flags for linear enforcement if they were given
             */
            const uint &
            linear_enforcement_flag( const uint aIndex ) const;


//------------------------------------------------------------------------------

            /**
             * return the procs that contribute to this kernel
             */
            const Vector< proc_t > &
            selected_procs() const;

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
             */
            void
            select_blocks( const Vector< id_t > & aBlockIDs );

//------------------------------------------------------------------------------

            /**
             * return the indices of the selected blocks
             */
            const Vector< uint > &
            block_indices();

//------------------------------------------------------------------------------

            void
            set_field_dimensions( const uint aNumDofsPerNode,
                                  const uint aNumDofsPerEdge= 0 );

//------------------------------------------------------------------------------

            void
            set_field_dimensions( const Vector< uint > & aNumDofsPerNode );

//------------------------------------------------------------------------------

            void
            set_field_dimensions(
                    const Vector< uint > & aNumDofsPerNode,
                    const Vector< uint > & aNumDofsPerEdge );

//------------------------------------------------------------------------------

            void
            enforce_linear( const Vector< uint > & aEnforceLinear );

//------------------------------------------------------------------------------

            void
            set_block_integration_orders( const Vector< uint > & aBlockIntegrationOrders );

//------------------------------------------------------------------------------

            void
            set_sideset_integration_orders( const Vector< uint > & aSideSetIntegrationOrders );

//------------------------------------------------------------------------------

            void
            set_block_integration_orders( const uint & aBlockIntegrationOrder );

//------------------------------------------------------------------------------

            void
            set_sideset_integration_orders( const uint & aSideSetIntegrationOrder );

//------------------------------------------------------------------------------

            void
            set_integration_scheme( const IntegrationScheme & aIntegrationScheme ) ;

//------------------------------------------------------------------------------

            const uint &
            block_integration_order( const index_t & aIndex ) const;

//------------------------------------------------------------------------------

            const uint &
            sideset_integration_order( const index_t & aIndex ) const;

//------------------------------------------------------------------------------

            const IntegrationScheme &
            integration_scheme() const ;

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

        /**
         * return the pointer of the associated mesh
         */
        inline Mesh *
        KernelParameters::mesh()
        {
            return mMesh;
        }

//------------------------------------------------------------------------------

        inline proc_t
        KernelParameters::number_of_procs() const
        {
            return mNumberOfProcs;
        }

//------------------------------------------------------------------------------

        inline const proc_t &
        KernelParameters::master()
        {
            return mMesh->master();
        }

//------------------------------------------------------------------------------

        inline const Vector< proc_t > &
        KernelParameters::selected_procs() const
        {
            return mSelectedProcs;
        }

//------------------------------------------------------------------------------

        inline const Vector< uint > &
        KernelParameters::num_dofs_per_node() const
        {
            return mNumDofsPerNode;
        }

//------------------------------------------------------------------------------

        inline const Vector< uint > &
        KernelParameters::num_dofs_per_edge() const
        {
            return mNumDofsPerEdge;
        }

//------------------------------------------------------------------------------

        inline const Vector< uint >&
        KernelParameters::linear_enforcement_flags() const
        {
            return mEnforceLinear ;
        }

//------------------------------------------------------------------------------

        inline const uint &
        KernelParameters::linear_enforcement_flag( const uint aIndex ) const
        {
            return aIndex < mEnforceLinear.length() ? mEnforceLinear( aIndex ) : mZero ;
        }

//------------------------------------------------------------------------------

        inline const Vector< uint > &
        KernelParameters::block_indices()
        {
            return mBlockIndices;
        }

//------------------------------------------------------------------------------

        inline const uint &
        KernelParameters::block_integration_order( const index_t & aIndex ) const
        {

            return aIndex < mBlockIntegrationOrders.length() ? mBlockIntegrationOrders( aIndex ) : mZero ;
        }

//------------------------------------------------------------------------------

        inline const uint &
        KernelParameters::sideset_integration_order( const index_t & aIndex ) const
        {
            return aIndex < mSideSetIntegrationOrders.length() ? mSideSetIntegrationOrders( aIndex ) : mZero ;
        }

//------------------------------------------------------------------------------

        inline const IntegrationScheme &
        KernelParameters::integration_scheme() const
        {
            return mIntegrationScheme ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_FEM_KERNELPARAMETERS_HPP
