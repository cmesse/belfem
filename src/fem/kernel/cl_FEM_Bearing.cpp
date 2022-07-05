//
// Created by Christian Messe on 15.11.19.
//

#include "cl_FEM_Bearing.hpp"
#include "cl_FEM_DofManagerBase.hpp"
#include "cl_FEM_Field.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_IWG.hpp"
#include "commtools.hpp"
namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Bearing::Bearing( DofManagerBase * aParent ) :
            mParent( aParent ),
            mID( 0 )
        {

        }

//------------------------------------------------------------------------------

        Bearing::Bearing( Field * aParent, const id_t aID, mesh::Node * aNode ) :
            mParent( aParent ),
            mID( aID ),
            mNode( aNode )
        {
            if ( aNode->owner() == comm_rank() )
            {
                // get number of dofs per node
                uint tNumDOFs = mParent->number_of_dofs_per_node();

                // allocate container
                mDOFs.set_size( tNumDOFs, nullptr );

                // populate container
                for ( uint k = 0; k < tNumDOFs; ++k )
                {
                    mDOFs( k ) = mParent->dof( mParent->calculate_dof_id( aNode, k ) );
                }
             }
        }

//------------------------------------------------------------------------------

        Bearing::Bearing( DofManager * aParent,
                          const id_t aID, mesh::Node * aNode ) :
                mParent( aParent ),
                mID( aID ),
                mNode( aNode )
        {
            if ( aNode->owner() == comm_rank() )
            {
                // get the entity types
                const Vector< index_t > & tTypes = aParent->iwg()->dof_entity_types();

                // get the total number of types
                uint tNumTypes = tTypes.length() ;

                // count node dof types
                index_t tCount = 0 ;
                for( uint k=0; k<tNumTypes; ++k )
                {
                    if( static_cast< EntityType >( tTypes( k ) ) == EntityType::NODE )
                    {
                        ++tCount ;
                    }
                }

                // allocate container
                mDOFs.set_size( tCount, nullptr );

                // reset counter
                tCount = 0 ;

                // loop over all types
                for( uint k=0; k<tNumTypes; ++k )
                {
                    // check if this is a node dof
                    if( static_cast< EntityType >( tTypes( k ) ) == EntityType::NODE )
                    {
                        // compute the DOF id
                        id_t tDofID = aParent->calculate_dof_id( aNode, tTypes( k ) );

                        // get the dof. Note: in very rare cases, this could cause a problem
                        // if the dof does not exist on the current proc.
                        // This must be fixed with a piece of code that changes the owner of

                        // test if dof exists, otherwise this contains a nullptr which might throw an error
                        if( aParent->dof_exists( tDofID ) )
                        {
                            mDOFs( tCount++ ) = aParent->dof( tDofID );
                        }
                    }
                }
            }

        }

//------------------------------------------------------------------------------

        void
        Bearing::impose_dirichlet( const real & aValue, const uint aDofType )
        {
            // need this if statement, since empty bearing contains no dofs
            if ( mDOFs.size() > aDofType )
            {
                // catch possible error
                BELFEM_ERROR( mDOFs( aDofType ) != nullptr,
                         "Tried to impose a boundary condition at dof %u on bearing at node %lu, but it does not exist",
                         ( unsigned int ) aDofType,
                         ( long unsigned int ) mNode->id() );

                mDOFs( aDofType )->fix( aValue );
            }
        }

//------------------------------------------------------------------------------

        void
        Bearing::free()
        {
            for ( Dof * tDOF : mDOFs )
            {
                tDOF->free();
            }
        }

//------------------------------------------------------------------------------
    }
}