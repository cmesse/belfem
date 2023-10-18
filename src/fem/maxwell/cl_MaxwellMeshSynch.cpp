//
// Created by christian on 2/9/23.
//
#include "commtools.hpp"
#include "cl_IWG_TimestepOld.hpp"
#include "cl_MaxwellMeshSynch.hpp"

namespace belfem
{
    namespace fem
    {
//----------------------------------------------------------------------------

        MaxwellMeshSynch::MaxwellMeshSynch( Kernel * aMagneticKernel, Kernel * aThermalKernel ) :
                mRank( comm_rank()),
                mMagneticKernel( *aMagneticKernel ),
                mThermalKernel( *aThermalKernel ),
                mMagneticMesh( * aMagneticKernel->mesh() ),
                mThermalMesh( * aThermalKernel->mesh() )
        {

        }

//----------------------------------------------------------------------------

        void
        MaxwellMeshSynch::initialize_tables()
        {
            // initialize table
            mTable.set_size(
                    mThermalMesh.number_of_nodes(), gNoIndex );

            // loop over all nodes on thermal mesh, knowing that the IDs are
            // identical to the IDs of the ghost elements on the maxwell mess
            for ( mesh::Node * tNode: mThermalMesh.nodes())
            {
                mTable( tNode->index())
                        = mMagneticMesh.element( tNode->id())->index();
            }

            comm_barrier();
        }

//----------------------------------------------------------------------------

        void
        MaxwellMeshSynch::magnetic_to_thermal_b_and_ej()
        {
            // grab field on maxwell mesh
            Vector< real > & mElementEJ = mMagneticMesh.field_data( "elementEJ" );

            // grab field on thermal mesh
            Vector< real > & mNodeEJ = mThermalMesh.field_data( "ej" );


            index_t tNodeIndex = 0;

            for ( index_t tElementIndex: mTable )
            {
                mNodeEJ( tNodeIndex++ ) = mElementEJ( tElementIndex );
            }

            // reset index
            tNodeIndex = 0;

            // grab field on maxwell mesh
            Vector< real > & mElementB = mMagneticMesh.field_data( "elementB" );

            // grab field on thermal mesh
            Vector< real > & mNodeB = mThermalMesh.field_data( "b" );

            for ( index_t tElementIndex: mTable )
            {
                mNodeB( tNodeIndex++ ) = mElementB( tElementIndex );
            }

            // unflag all nodes on thermal mesh
            mThermalKernel.mesh()->unflag_all_nodes();
            comm_barrier();

            // synch time step
            IWG_TimestepOld * tMaxwell = reinterpret_cast< IWG_TimestepOld * >( mMagneticKernel.dofmgr( 0 )->iwg() );
            IWG_TimestepOld * tFourier = reinterpret_cast< IWG_TimestepOld * >( mThermalKernel.dofmgr( 0 )->iwg() );
            tFourier->set_timestep( tMaxwell->timestep() );

        }

//---------------------------------------------------------------------------

        void
        MaxwellMeshSynch::thermal_to_magnetic_T()
        {
            // grab field on maxwell mesh
            Vector< real > & mElementT = mMagneticMesh.field_data( "elementT" );

            // grab field on thermal mesh
            Vector< real > & mNodeT = mThermalMesh.field_data( "T" );

            index_t tNodeIndex = 0;

            for ( index_t tElementIndex: mTable )
            {
                // std::cout << "copy " << tElementIndex << " " << mNodeT( tNodeIndex ) << std::endl ;

                mElementT( tElementIndex ) = mNodeT( tNodeIndex++ );
            }

            comm_barrier();
        }

//---------------------------------------------------------------------------

        void
        MaxwellMeshSynch::magnetic_to_thermal_T()
        {
            // grab field on maxwell mesh
            Vector< real > & mElementT = mMagneticMesh.field_data( "elementT" );

            // grab field on thermal mesh
            Vector< real > & mNodeT = mThermalMesh.field_data( "T" );

            index_t tNodeIndex = 0;

            for ( index_t tElementIndex: mTable )
            {
                mNodeT( tNodeIndex++ ) = mElementT( tElementIndex ) ;
            }

            comm_barrier();
        }

//----------------------------------------------------------------------------
    }
}