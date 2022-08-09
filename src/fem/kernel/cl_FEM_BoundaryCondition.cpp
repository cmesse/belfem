//
// Created by christian on 9/29/21.
//

#include "cl_FEM_BoundaryCondition.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_Dof.hpp"
#include "cl_FEM_Kernel.hpp"

namespace belfem
{
    namespace fem
    {
//-----------------------------------------------------------------------------

        BoundaryConditionScaling
        bc_scaling_type( const string & aString )
        {
            string tString = string_to_lower( aString );

            if( tString == "constant" || tString == "const" || tString == "const." )
            {
                return BoundaryConditionScaling::Constant ;
            }
            else if ( tString == "ramp" )
            {
                return BoundaryConditionScaling::Ramp ;
            }
            else if ( tString == "sine" )
            {
                return BoundaryConditionScaling::Sine ;
            }
            else if ( tString == "sawtooth" )
            {
                return BoundaryConditionScaling::Sawtooth;
            }
            else if ( tString == "square" )
            {
                return BoundaryConditionScaling::Square ;
            }
            else if ( tString == "triangle" )
            {
                return BoundaryConditionScaling::Triangle ;
            }
            else if ( tString == "sigmoid" )
            {
                return BoundaryConditionScaling::Sigmoid ;
            }
            else
            {
                BELFEM_ERROR( false, "unknown bc scaling type : %s", aString.c_str() );
                return BoundaryConditionScaling::UNDEFINED ;
            }
        }
//-----------------------------------------------------------------------------

        BoundaryCondition::BoundaryCondition(
                const BoundaryConditionPhysics   aPhysics,
                const BoundaryConditionImposing  aImposing,
                const string                  aLabel ) :
                mPhysics( aPhysics ),
                mImposing( aImposing ),
                mLabel( aLabel )
        {

        }

//-----------------------------------------------------------------------------

        void
        BoundaryCondition::set_dof_manager( DofManager * aDofManager )
        {
            mParent = aDofManager ;

            if( mSideSetIDs.length() > 0 )
            {
                this->link_sidesets( aDofManager, mSideSetIDs );
            }
            if( mBlockIDs.length() > 0 )
            {
                this->link_blocks( aDofManager, mBlockIDs );
            }
            if( mBearingID != gNoID )
            {
                comm_barrier() ;
                if( comm_rank() == aDofManager->parent()->master() )
                {
                    mBearingID = aDofManager->mesh()->vertex( mBearingID )->node( 0 )->id() ;
                    send( aDofManager->parent()->comm_table(), mBearingID );
                }
                else
                {
                    receive( aDofManager->parent()->master(), mBearingID );
                }

                this->link_bearing( aDofManager, mBearingID );
            }
            this->collect_nodes();
        }


//-----------------------------------------------------------------------------

        void
        BoundaryCondition::set_sidesets( const Vector< id_t > & aSideSetIDs )
        {
            mSideSetIDs = aSideSetIDs ;
        }

//-----------------------------------------------------------------------------

        void
        BoundaryCondition::set_penalty( const real aValue )
        {
            mPenalty = aValue ;
        }
//-----------------------------------------------------------------------------

        void
        BoundaryCondition::set_bearing( const id_t aVertexID )
        {
            mBearingID = aVertexID ;
        }

//-----------------------------------------------------------------------------

        void
        BoundaryCondition::set_blocks( const Vector< id_t > & aBlockIDs )
        {
            mBlockIDs = aBlockIDs ;
        }

//-----------------------------------------------------------------------------

        void
        BoundaryCondition::link_bearing(  DofManager * aDofManager, const id_t aBearingID )
        {
            if( aDofManager->mesh()->node_exists( aBearingID ) )
            {
                mNodes.set_size( 1, nullptr );
                mNodes( 0 ) = aDofManager->mesh()->node( aBearingID );
            }
        }

//-----------------------------------------------------------------------------

        void
        BoundaryCondition::link_sidesets(  DofManager * aDofManager, const Vector< id_t > & aIDs )
        {
            // reset containers
            mSideSets.clear();

            // initialize counter
            index_t tCount = 0 ;

            // count existing sidesets
            for( id_t tID : aIDs )
            {
                if( aDofManager->sideset_exists( tID ) )
                {
                    // increment counter
                    ++tCount ;
                }
            }

            if( tCount > 0 )
            {
                // pick sideset from parent
                mSideSets.set_size( tCount, nullptr );

                // reset counter
                tCount = 0;

                for ( id_t tID: aIDs )
                {
                    if ( aDofManager->sideset_exists( tID ) )
                    {
                        aDofManager->sideset( tID )->set_boundary_condition( this );
                        mSideSets( tCount++ ) = aDofManager->sideset( tID );
                    }
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        BoundaryCondition::link_blocks(  DofManager * aDofManager, const Vector< id_t > & aIDs )
        {
            // sanity check
            //this->check_if_linked_to_dof_manager();

            // reset containers
            mBlocks.clear();

            // initialize counter
            index_t tCount = 0 ;

            // count existing blocks
            for( id_t tID : aIDs )
            {
                if( aDofManager->block_exists( tID ) )
                {
                    // increment counter
                    ++tCount ;
                }
            }

            if( tCount > 0 )
            {
                // allocate container
                mBlocks.set_size( tCount, nullptr );

                // reset counter
                tCount = 0 ;

                for( id_t tID : aIDs )
                {
                    if( mParent->block_exists( tID ) )
                    {
                        aDofManager->block( tID )->set_boundary_condition( this );
                        mBlocks( tCount++ ) = aDofManager->block( tID );
                    }
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        BoundaryCondition::check_if_linked_to_dof_manager()
        {
            BELFEM_ERROR( mParent != nullptr, "Boundary Condition %s is not linked to a dof manager",
                         mLabel.c_str() );
        };

//-----------------------------------------------------------------------------

        void
        BoundaryCondition::collect_nodes()
        {
            // reset container
            mNodes.clear() ;

            // get mesh
            Mesh * tMesh = mParent->mesh();

            // unflag all nodes
            tMesh->unflag_all_nodes();

            // flag all nodes on the selected groups
            if( mBlockIDs.length() > 0 )
            {
                for ( id_t tID : mBlockIDs )
                {
                    if( tMesh->block_exists( tID ) )
                    {
                        tMesh->block( tID )->flag_nodes() ;
                    }
                }
            }
            else
            {
                for ( id_t tID : mSideSetIDs )
                {
                    if( tMesh->sideset_exists( tID ) )
                    {
                        tMesh->sideset( tID )->flag_all_nodes();
                    }
                }
            }




            // count nodes
            index_t tCount = 0 ;
            Cell< mesh::Node * > & tNodes = tMesh->nodes() ;

            for( mesh::Node * tNode : tNodes )
            {
                // check if node is flagged
                if( tNode->is_flagged() )
                {
                    // increment counter
                    ++tCount ;
                }
            }

            // allocate memory
            mNodes.set_size( tCount, nullptr );

            // reset counter
            tCount = 0 ;
            for( mesh::Node * tNode : tNodes )
            {
                // check if node is flagged
                if( tNode->is_flagged() )
                {
                   mNodes( tCount++ ) = tNode ;
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        BoundaryCondition::compute( const real aTimeStamp  )
        {
            BELFEM_ERROR( false, "invalid call to base class");
        }

//-----------------------------------------------------------------------------

        void
        BoundaryCondition::update_dof_types()
        {
            // sanity check
            this->check_if_linked_to_dof_manager();

            switch( mImposing )
            {
                case ( BoundaryConditionImposing::Dirichlet ) :
                case( BoundaryConditionImposing::Neumann ) :
                {
                    const Cell< string > & tFields =
                            mImposing == BoundaryConditionImposing::Neumann ?
                            mParent->iwg()->flux_fields() :
                            mParent->iwg()->dof_fields() ;

                    // map with dof types
                    Map< string, uint > tMap ;

                    // populate map
                    for( uint f=0; f<tFields.size(); ++f )
                    {
                        tMap[ tFields( f ) ] = f ;
                    }

                    // populate dof types
                    mDofTypes.set_size( mFields.size() );

                    uint tCount = 0 ;

                    for( const string & tField : mFields )
                    {
                        mDofTypes( tCount++ ) = tMap( tField );
                    }
                    break ;
                }
                default :
                {
                    // do nothing
                }
            }

        }

//-----------------------------------------------------------------------------

        void
        BoundaryCondition::set_scaling(
                const BoundaryConditionScaling aType,
                const real        aPeriod,
                const real        aPhase,
                const real        aFuzzyness )
        {
            mScaling = aType ;

            mPeriod = aPeriod ;
            mPhase = aPhase ;
            mFuziness = aFuzzyness ;

            if( aPeriod < BELFEM_REAL_MAX )
            {
                mFrequency = 1.0 / aPeriod;
                mOmega = 2.0 * constant::pi / aPeriod;
                mTimeOffset = aPhase / mOmega;
            }
            else
            {
                mFrequency = 0.0 ;
                mOmega = 0.0 ;
                mTimeOffset = 0.0 ;
            }

            // link scaling function
            switch( aType )
            {
                case( BoundaryConditionScaling::Constant ) :
                {
                    mFunScale = & BoundaryCondition::scale_constant ;
                    break ;
                }
                case( BoundaryConditionScaling::Ramp ) :
                {
                    mFunScale = & BoundaryCondition::scale_ramp ;
                    break ;
                }
                case( BoundaryConditionScaling::Sigmoid ) :
                {
                    mFunScale = & BoundaryCondition::scale_sigmoid ;
                    break ;
                }
                case( BoundaryConditionScaling::Sine ) :
                {
                    mFunScale = & BoundaryCondition::scale_sine ;
                    break ;
                }
                case( BoundaryConditionScaling::Sawtooth ) :
                {
                    mFunScale = & BoundaryCondition::scale_sawtooth ;
                    break ;
                }
                case( BoundaryConditionScaling::Square ) :
                {
                    mFunScale = & BoundaryCondition::scale_square ;
                    break ;
                }
                case( BoundaryConditionScaling::Triangle ) :
                {
                    mFunScale = & BoundaryCondition::scale_triangle ;
                    break ;
                }
                case( BoundaryConditionScaling::UNDEFINED ) :
                {
                    // if this is chosen, the BC must have another way
                    // of computing the function
                    mFunScale = nullptr ;
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Invalid BC scaling");
                }
            }
        }

//-----------------------------------------------------------------------------
    }
}