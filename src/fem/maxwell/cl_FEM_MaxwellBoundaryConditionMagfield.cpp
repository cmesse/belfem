//
// Created by christian on 10/8/21.
//

#include "cl_FEM_MaxwellBoundaryConditionMagfield.hpp"
#include "cl_FEM_DofManager.hpp"
#include "fn_norm.hpp"

namespace belfem
{
    namespace fem
    {

//-----------------------------------------------------------------------------

        MaxwellBoundaryConditionMagfield::MaxwellBoundaryConditionMagfield(
                const BoundaryConditionImposing aImposing,
                const string             & aLabel ) :
                BoundaryCondition( BoundaryConditionPhysics::Magfield,
                                   aImposing,
                                   aLabel )
        {
        }

//-----------------------------------------------------------------------------

        void
        MaxwellBoundaryConditionMagfield::set(
                const BoundaryConditionScaling aType,
                const Vector< real > & aAmplitude,
                const real aPeriod,
                const real aPhase )
        {
            mSubType = MagfieldBcType::Wave ;
            mAmplitude = aAmplitude ;
            this->set_scaling( aType, aPeriod, aPhase );
        }

//-----------------------------------------------------------------------------

        void
        MaxwellBoundaryConditionMagfield::set(
                const MagfieldBcType aType )
        {
            BELFEM_ERROR(  aType != MagfieldBcType::Wave, "need more parameters" );
            mSubType = aType ;
        }

//-----------------------------------------------------------------------------

        void
        MaxwellBoundaryConditionMagfield::compute( const real aTime )
        {
            if( mSubType == MagfieldBcType::Wave )
            {
                // time dependent scaling
                real tScale = this->scale( aTime );

                // get the number of dimensions
                uint tNumberOfDimensions = mParent->mesh()->number_of_dimensions();

                // loop over all fields
                for ( uint i = 0; i < mBoundaryFields.size(); ++i )
                {
                    // get data from mesh
                    Vector< real > & tData = mParent->mesh()->field_data( mBoundaryFields( i ) );

                    real tAmplitude = mAmplitude( i ) * tScale;
                    mData( i ) = tAmplitude ;

                    // write values into nodes
                    for ( index_t tIndex: mNodeIndices )
                    {
                        tData( tIndex ) = tAmplitude;
                    }
                }
                comm_barrier();
                mParent->collect_fields( mBoundaryFields );
                comm_barrier();
                mParent->distribute_fields( mBoundaryFields );

                // add data to mesh
                uint tCount = 0;
                for ( string & tLabel: mBoundaryFields )
                {
                    mParent->mesh()->global_variable_data( mLabel + "_" + tLabel ) = mAmplitude( tCount++ ) * tScale;
                }


                // check if this is a phi field
                if ( mIsHPhi )
                {
                    // BELFEM_ASSERT( mDofFields.size() == 1, "only one field is allowed in mDofFields" );
                    mData *= - constant::nu0 ;

                    if( mImposing == BoundaryConditionImposing::Dirichlet )
                    {
                        Vector< real > & tPhi = mParent->mesh()->field_data( "phi" );

                        // must convert B-Field to H-Field, also flip sign since b = - nu0 * grad phi
                        tScale *= -constant::nu0;

                        if ( tNumberOfDimensions == 2 )
                        {
                            real tHx = mAmplitude( 0 ) * tScale;
                            real tHy = mAmplitude( 1 ) * tScale;

                            for ( Dof * tDof: mFixedDofs )
                            {
                                // get the node
                                mesh::Node * tNode = tDof->node();

                                // compute the value of phi
                                tDof->fix(( tHx * tNode->x() + tHy * tNode->y()));

                                // write value into field list
                                tPhi( tNode->index()) = tDof->value();
                            }
                        }
                        else if ( tNumberOfDimensions == 3 )
                        {
                            real tHx = mAmplitude( 0 ) * tScale;
                            real tHy = mAmplitude( 1 ) * tScale;
                            real tHz = mAmplitude( 2 ) * tScale;

                            for ( Dof * tDof: mFixedDofs )
                            {
                                // get the node
                                mesh::Node * tNode = tDof->node();

                                // compute the value of phi
                                tDof->fix( tHx * tNode->x()
                                           + tHy * tNode->y()
                                           + tHz * tNode->z());

                                // write valye into field list
                                tPhi( tNode->index()) = tDof->value();
                            }
                        }
                    } // end if dirichlet
                }
                else // this is a HA-Model
                {
                    if( mImposing == BoundaryConditionImposing::Dirichlet )
                    {
                        if ( tNumberOfDimensions == 2 )
                        {
                            real tBx = mAmplitude( 0 ) * tScale;
                            real tBy = mAmplitude( 1 ) * tScale;

                            Vector< real > & tAz = mParent->mesh()->field_data( "az" );

                            for ( Dof * tDof: mFixedDofs )
                            {
                                // get the node
                                mesh::Node * tNode = tDof->node();

                                // compute the value of Az
                                tDof->fix( tBx * tNode->y()
                                           - tBy * tNode->x());

                                // write valye into field list
                                tAz( tNode->index()) = tDof->value();

                            }
                        }
                        else if ( tNumberOfDimensions == 3 )
                        {
                            real tBx = mAmplitude( 0 ) * tScale;
                            real tBy = mAmplitude( 1 ) * tScale;
                            real tBz = mAmplitude( 2 ) * tScale;

                            Vector< real > & tAx = mParent->mesh()->field_data( mDofFields( 0 ));
                            Vector< real > & tAy = mParent->mesh()->field_data( mDofFields( 1 ));
                            Vector< real > & tAz = mParent->mesh()->field_data( mDofFields( 2 ));


                            // reset counter
                            tCount = 0;

                            for ( index_t k: mNodeIndices )
                            {
                                // get dofs
                                Dof * tDofX = mFixedDofs( tCount++ );
                                Dof * tDofY = mFixedDofs( tCount++ );
                                Dof * tDofZ = mFixedDofs( tCount++ );

                                // get node coords
                                real tX = tDofX->node()->x();
                                real tY = tDofY->node()->y();
                                real tZ = tDofZ->node()->z();

                                // compute dofs
                                tDofX->fix( 0.5 * ( tBy * tZ - tBz * tY ));
                                tDofY->fix( 0.5 * ( tBz * tX - tBx * tZ ));
                                tDofZ->fix( 0.5 * ( tBx * tY - tBy * tX ));

                                // write dofs into fields
                                tAx( k ) = tDofX->value();
                                tAy( k ) = tDofY->value();
                                tAz( k ) = tDofZ->value();
                            }
                        }
                    }
                }

                comm_barrier();
                mParent->collect_fields( mDofFields );
                comm_barrier();
                mParent->distribute_fields( mDofFields );

                // wait for other procs
                comm_barrier();
            }
        }

//-----------------------------------------------------------------------------

        void
        MaxwellBoundaryConditionMagfield::link_bearing(  DofManager * aDofManager, const id_t aBearingID )
        {
            BoundaryCondition::link_bearing( aDofManager, aBearingID );

            // check if this is a h-phi formulation
            mIsHPhi = is_maxwell_hphi( aDofManager->iwg()->type() ) ;

            if( mNodes.size() == 0 )
            {
                return;
            }

            if( mIsHPhi )
            {
                mDofFields = { "phi" };

                // get the type id of phi
                uint tPhi = aDofManager->iwg()->doftype( "phi" );

                if( mImposing == BoundaryConditionImposing::Dirichlet )
                {
                    // compute dof ID
                    id_t tID = aDofManager->calculate_dof_id( mNodes( 0 ), tPhi );

                    // only do something if dof exists on this proc
                    if( aDofManager->dof_exists( tID ) )
                    {
                        mFixedDofs.set_size( 1, nullptr );

                        Dof * tDof = aDofManager->dof(  tID );
                        tDof->fix( 0.0 );

                        // save dof in container
                        mFixedDofs( 0 ) = tDof ;
                    }
                }

            }
            else  // this is an a-formulation
            {
                if( mImposing == BoundaryConditionImposing::Dirichlet )
                {
                    if ( aDofManager->mesh()->number_of_dimensions() == 2 )
                    {
                        mDofFields = { "az" };

                        uint tAz = aDofManager->iwg()->doftype( "az" );


                        // compute dof ID
                        id_t tID =  aDofManager->calculate_dof_id( mNodes( 0 ), tAz ) ;

                        // only do something if dof exists on this proc
                        if( aDofManager->dof_exists( tID ) )
                        {
                            mFixedDofs.set_size( 1, nullptr );

                            Dof * tDof = aDofManager->dof(  tID );
                            tDof->fix( 0.0 );

                            // save dof in container
                            mFixedDofs( 0 ) = tDof ;
                        }
                    }
                    else if ( aDofManager->mesh()->number_of_dimensions() == 3 )
                    {
                        mDofFields = { "ax", "ay", "az" };

                        mFixedDofs.set_size( 3 , nullptr );

                        uint tAx = aDofManager->iwg()->doftype( "ax" );
                        uint tAy = aDofManager->iwg()->doftype( "ay" );
                        uint tAz = aDofManager->iwg()->doftype( "az" );

                        // compute DoF IDs
                        id_t tIDx = aDofManager->calculate_dof_id( mNodes( 0 ), tAx ) ;
                        id_t tIDy = aDofManager->calculate_dof_id( mNodes( 0 ), tAy ) ;
                        id_t tIDz = aDofManager->calculate_dof_id( mNodes( 0 ), tAz ) ;

                        uint tCount = 0 ;

                        if( aDofManager->dof_exists( tIDx ) )
                        {
                            ++tCount ;
                        }
                        if( aDofManager->dof_exists( tIDy ) )
                        {
                            ++tCount ;
                        }
                        if( aDofManager->dof_exists( tIDx ) )
                        {
                            ++tCount ;
                        }

                        if( tCount > 0 )
                        {
                            mFixedDofs.set_size( tCount, nullptr );
                            tCount = 0 ;
                        }

                        if( aDofManager->dof_exists( tIDx ) )
                        {
                            Dof * tDof = aDofManager->dof( tIDx );
                            tDof->fix( 0.0 );
                            mFixedDofs( tCount++ ) = tDof;
                        }

                        if( aDofManager->dof_exists( tIDy ) )
                        {
                            Dof * tDof = aDofManager->dof( tIDy );
                            tDof->fix( 0.0 );
                            mFixedDofs( tCount++ ) = tDof;
                        }

                        if( aDofManager->dof_exists( tIDz ) )
                        {
                            Dof * tDof = aDofManager->dof( tIDz );
                            tDof->fix( 0.0 );
                            mFixedDofs( tCount++ ) = tDof;
                        }
                    }
                }
            }
        }


//-----------------------------------------------------------------------------

        void
        MaxwellBoundaryConditionMagfield::link_sidesets(
                DofManager * aDofManager, const Vector< id_t > & aIDs )
        {
            BoundaryCondition::link_sidesets( aDofManager, aIDs );

            // get pointer to mesh
            Mesh * tMesh = aDofManager->mesh();

            // check if this is a h-phi formulation
            mIsHPhi = is_maxwell_hphi( aDofManager->iwg()->type() ) ;


            mData.set_size( tMesh->number_of_dimensions(), 0.0 );


            // unflag all nodes on the mesh
            tMesh->unflag_all_nodes();

            if( aDofManager->enforce_linear_interpolation() )
            {
                for( id_t tID : aIDs )
                {
                    if( tMesh->sideset_exists( tID ) )
                    {
                        tMesh->sideset( tID )->flag_corner_nodes();
                    }
                }
            }
            else
            {
                for( id_t tID : aIDs )
                {
                    if( tMesh->sideset_exists( tID ) )
                    {
                        tMesh->sideset( tID )->flag_all_nodes();
                    }
                }
            }

            // count nodes that belong to this boundary condition
            index_t tCount = 0 ;
            Cell< mesh::Node * > & tNodes = tMesh->nodes() ;
            for( mesh::Node * tNode : tNodes )
            {
                if( tNode->is_flagged() )
                {
                    ++tCount ;
                }
            }

            // allocate memory with local node indices
            mNodeIndices.set_size( tCount );

            // populate node indices
            tCount = 0 ;
            for( mesh::Node * tNode : tNodes )
            {
                if( tNode->is_flagged() )
                {
                    mNodeIndices( tCount++ ) = tNode->index() ;
                }
            }

            if( mIsHPhi )
            {
                mDofFields = { "phi" };

                // get the type id of phi
                uint tPhi = aDofManager->iwg()->doftype( "phi" );

                if( mImposing == BoundaryConditionImposing::Dirichlet )
                {
                    mFixedDofs.set_size( tCount, nullptr );
                    tCount = 0;
                    for ( mesh::Node * tNode: tNodes )
                    {
                        if ( tNode->is_flagged())
                        {
                            // get dof
                            Dof * tDof = aDofManager->dof( aDofManager->calculate_dof_id( tNode, tPhi ));
                            tDof->fix( 0.0 );

                            mFixedDofs( tCount++ ) = tDof;
                        }
                    }
                }

            }
            else  // this is a n a-formulation
            {
                if( mImposing == BoundaryConditionImposing::Dirichlet )
                {
                    if ( tMesh->number_of_dimensions() == 2 )
                    {
                        mDofFields = { "az" };

                        uint tAz = aDofManager->iwg()->doftype( "az" );

                        mFixedDofs.set_size( tCount, nullptr );
                        tCount = 0;

                        for ( mesh::Node * tNode: tNodes )
                        {
                            if ( tNode->is_flagged())
                            {
                                // get dof
                                Dof * tDof = aDofManager->dof( aDofManager->calculate_dof_id( tNode, tAz ));
                                tDof->fix( 0.0 );

                                mFixedDofs( tCount++ ) = tDof;
                            }
                        }

                    }
                    else if ( tMesh->number_of_dimensions() == 3 )
                    {
                        mDofFields = { "ax", "ay", "az" };

                        mFixedDofs.set_size( 3 * tCount, nullptr );
                        tCount = 0;

                        uint tAx = aDofManager->iwg()->doftype( "ax" );
                        uint tAy = aDofManager->iwg()->doftype( "ay" );
                        uint tAz = aDofManager->iwg()->doftype( "az" );

                        for ( mesh::Node * tNode: tNodes )
                        {
                            if ( tNode->is_flagged())
                            {
                                // fix ax
                                Dof * tDofX = aDofManager->dof( aDofManager->calculate_dof_id( tNode, tAx ));
                                tDofX->fix( 0.0 );
                                mFixedDofs( tCount++ ) = tDofX;

                                // fix Ay
                                Dof * tDofY = aDofManager->dof( aDofManager->calculate_dof_id( tNode, tAy ));
                                tDofY->fix( 0.0 );
                                mFixedDofs( tCount++ ) = tDofY;


                                Dof * tDofZ = aDofManager->dof( aDofManager->calculate_dof_id( tNode, tAz ));
                                tDofZ->fix( 0.0 );
                                mFixedDofs( tCount++ ) = tDofZ;
                            }
                        }
                    }
                }
            }

            if( tMesh->number_of_dimensions() == 2 )
            {
                mBoundaryFields = { "bx", "by" };
            }
            else if( tMesh->number_of_dimensions() == 3 )
            {
                mBoundaryFields = { "bx", "by", "bz" } ;
            }

            for( string tLabel : mBoundaryFields )
            {
                string tName = mLabel + "_" + tLabel ;

                // create global variables
                if( ! aDofManager->mesh()->global_variable_exists( tName ) )
                {
                    aDofManager->mesh()->create_global_variable( tName );
                }
            }

        }

//-----------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace belfem */