//
// Created by christian on 9/30/21.
//


#include "constants.hpp"
#include "assert.hpp"
#include "cl_FEM_MaxwellBoundaryConditionCurrent.hpp"
#include "stringtools.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_Kernel.hpp"
#include "fn_det.hpp"
#include "commtools.hpp"
#include "fn_sum.hpp"
#include "fn_Mesh_compute_volume.hpp"
namespace belfem
{
    namespace fem
    {


//-----------------------------------------------------------------------------

        MaxwellBoundaryConditionCurrent::MaxwellBoundaryConditionCurrent(
                const BoundaryConditionImposing aImposing,
                const string          & aLabel ) :
                    BoundaryCondition( BoundaryConditionPhysics::Current,
                                       aImposing,
                                       aLabel )
        {

        }

//-----------------------------------------------------------------------------

        void
        MaxwellBoundaryConditionCurrent::set(
                const BoundaryConditionScaling aType,
                const real        aAmplitude,
                const real        aPeriod,
                const real        aPhase )
        {
            // copy input
            mAmplitude = aAmplitude ;

            this->set_scaling( aType, aPeriod, aPhase );
        }

//-----------------------------------------------------------------------------

        void
        MaxwellBoundaryConditionCurrent::set_dof_manager( DofManager * aDofManager )
        {
            // call function from parent
            BoundaryCondition::set_dof_manager( aDofManager );

            // set bc labels
            switch( aDofManager->iwg()->type() )
            {
                case( IwgType::MAXWELL_HA_TRI3 ) :
                case( IwgType::MAXWELL_HA_TRI6 ) :
                case( IwgType::MAXWELL_HPHI_TRI3 ) :
                case( IwgType::MAXWELL_HPHI_TRI6 ) :
                {
                    mCurrentFields = { "jz" } ;
                    break ;
                }
                case( IwgType::MAXWELL_HA_TET4 ) :
                case( IwgType::MAXWELL_HA_TET10 ) :
                case( IwgType::MAXWELL_HPHI_TET4 ) :
                case( IwgType::MAXWELL_HPHI_TET10 ) :
                {
                    mCurrentFields = { "jx", "jy", "jz" } ;
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "invalid equation type");
                }
            }

            // create global variable
            if( ! aDofManager->mesh()->global_variable_exists( mLabel ) )
            {
                aDofManager->mesh()->create_global_variable( mLabel );
            }
        }

//-----------------------------------------------------------------------------

        void
        MaxwellBoundaryConditionCurrent::link_blocks(
                DofManager * aDofManager,
                const Vector< id_t > & aIDs )
        {
            BoundaryCondition::link_blocks( aDofManager, aIDs );

            if( mCoilCrossSection == 0.0 )
            {
                // for 2D, we can also compute the cross section
                this->compute_cross_section( aDofManager );
            }
        }

//-----------------------------------------------------------------------------

        void
        MaxwellBoundaryConditionCurrent::compute_cross_section( DofManager * aDofManager )
        {
            if( aDofManager->mesh()->number_of_dimensions() == 2 )
            {
                // reset cross section (only for h-a)
                if( mBlockIDs.length() > 0 )
                {
                    mCoilCrossSection = mesh::compute_volume( aDofManager->mesh(), mBlockIDs );
                }

                comm_barrier();
            }
            else
            {
                BELFEM_ERROR( false, "can't compute cross section for 3d coil" );
            }
        }

//-----------------------------------------------------------------------------

        void
        MaxwellBoundaryConditionCurrent::compute( const real aTime )
        {
            // compute the current
            real tI = mAmplitude * this->scale( aTime );

            mParent->mesh()->global_variable_data( mLabel ) = tI ;

            switch( mParent->iwg()->type() )
            {
                case( IwgType::MAXWELL_HA_TRI3 ) :
                case( IwgType::MAXWELL_HA_TRI6 ) :
                {
                    this->compute_current_densities_2d( tI );

                    mParent->synchronize_fields( mCurrentFields );

                    break ;
                }
                case( IwgType::MAXWELL_HPHI_TRI3 ) :
                case( IwgType::MAXWELL_HPHI_TRI6 ) :
                {

                    // write current info to fields on mesh
                    this->compute_cut_currents( tI );

                    // values for visualization on coils
                    // values for superconductors will be overwritten by L2
                    this->compute_current_densities_2d( tI );


                    mParent->synchronize_fields( mCurrentFields );

                    break;
                }
                case( IwgType::MAXWELL_HA_TET4 ) :
                case( IwgType::MAXWELL_HA_TET10 ) :
                {
                    BELFEM_ERROR( false, "compute() not implemented for 3d mesh");
                    break ;
                }
                case( IwgType::MAXWELL_HPHI_TET4 ) :
                case( IwgType::MAXWELL_HPHI_TET10 ) :
                {
                    BELFEM_ERROR( false, "compute() not implemented for 3d mesh");
                    break ;
                }
                default :
                {
                    // pass
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        MaxwellBoundaryConditionCurrent::compute_current_densities_2d( const real aCurrent )
        {
            real tJ = BELFEM_QUIET_NAN ;
            if( mCoilCrossSection > 0.0 )
            {
                // compute the current density
                real tJ = aCurrent / mCoilCrossSection;

                // failsafe
                if ( std::abs( tJ ) < BELFEM_EPSILON )
                {
                    tJ = BELFEM_EPSILON;
                }
            }

            // get field from mesh
            Vector< real > & tValues =
                    mParent->parent()->mesh()->field_data( "jz" );

            // write node values into field
            for ( mesh::Node * tNode: mNodes )
            {
                tValues( tNode->index() ) = tJ;
            }
        }

//-----------------------------------------------------------------------------

        void
        MaxwellBoundaryConditionCurrent::compute_cut_currents( const real aCurrent )
        {
            Vector< real > & tValues =
                    mParent->parent()->mesh()->field_data("lambda_I");

            // get mesh
            Mesh * tMesh = mParent->parent()->mesh() ;

            for( id_t tID : mSideSetIDs )
            {
                if ( tMesh->cut_exists( tID ) )
                {
                    // grab all facets on the cut
                    Cell< mesh::Facet * > & tFacets = tMesh->cut( tID )->facets() ;

                    // write into to mesh
                    for( mesh::Facet * tFacet : tFacets )
                    {
                        tValues( tFacet->index() ) = aCurrent ;
                    }
                }
            }
        }

//-----------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace belfem */