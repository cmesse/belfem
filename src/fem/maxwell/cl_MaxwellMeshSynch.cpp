//
// Created by christian on 2/9/23.
//
#include "commtools.hpp"
#include "cl_MaxwellMeshSynch.hpp"

namespace belfem
{
//----------------------------------------------------------------------------

    MaxwellMeshSynch::MaxwellMeshSynch( Mesh * aMaxwellMesh, Mesh * aThermalMesh ) :
        mRank( comm_rank() ),
        mMaxwellMesh( * aMaxwellMesh ),
        mThermalMesh( * aThermalMesh )
    {
        if( mRank == 0 )
        {
            this->set_ownerships() ;
        }
        comm_barrier() ;
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
        for ( mesh::Node * tNode: mThermalMesh.nodes() )
        {
            mTable( tNode->index() )
                    = mMaxwellMesh.element( tNode->id() )->index();
        }

        comm_barrier() ;
    }

//----------------------------------------------------------------------------

    void
    MaxwellMeshSynch::set_ownerships()
    {


        Cell < mesh::Node * >    & tNodes = mThermalMesh.nodes() ;

        // get node owners from other mesh
        for( mesh::Node * tNode : tNodes )
        {
            tNode->set_owner(
                    mMaxwellMesh.element( tNode->id())->owner() );
        }

        // generate element owners
        Cell < mesh::Element * > & tElements = mThermalMesh.elements() ;
        proc_t tNumProcs = comm_size() ;
        for( mesh::Element * tElement : tElements )
        {
            proc_t tOwner = tNumProcs ;
            for( uint k=0; k<tElement->number_of_nodes(); ++k )
            {
                if( tElement->node( k )->owner() < tOwner )
                {
                    tOwner = tElement->node( k )->owner();
                }
                tElement->set_owner( tOwner );
            }
        }

    }

//----------------------------------------------------------------------------

    void
    MaxwellMeshSynch::maxwell_to_thermal()
    {
        // grab field on maxwell mesh
        Vector< real > & mElementEJ = mMaxwellMesh.field_data("elementEJ") ;

        // grab field on thermal mesh
        Vector< real > & mNodeEJ    = mThermalMesh.field_data("ej");


        index_t tNodeIndex = 0 ;

        for( index_t tElementIndex : mTable )
        {
            mNodeEJ( tNodeIndex++ ) = mElementEJ( tElementIndex );
        }

        // reset index
        tNodeIndex = 0 ;

        // grab field on maxwell mesh
        Vector< real > & mElementB = mMaxwellMesh.field_data("elementB") ;

        // grab field on thermal mesh
        Vector< real > & mNodeB    = mThermalMesh.field_data("b");

        for( index_t tElementIndex : mTable )
        {
            mNodeB( tNodeIndex++ ) = mElementB( tElementIndex );
        }

        comm_barrier() ;
    }

//---------------------------------------------------------------------------

    void
    MaxwellMeshSynch::thermal_to_maxwell()
    {
        // grab field on maxwell mesh
        Vector< real > & mElementT = mMaxwellMesh.field_data("elementT") ;

        // grab field on thermal mesh
        Vector< real > & mNodeT   = mThermalMesh.field_data("T");

        index_t tNodeIndex = 0 ;

        for( index_t tElementIndex : mTable )
        {
            mElementT( tElementIndex ) = mNodeT( tNodeIndex++) ;
        }

        comm_barrier() ;
    }


//----------------------------------------------------------------------------
}