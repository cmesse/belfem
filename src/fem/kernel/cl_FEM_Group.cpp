//
// Created by Christian Messe on 08.11.19.
//

#include "cl_FEM_Group.hpp"
#include "assert.hpp"
#include "commtools.hpp"
#include "meshtools.hpp"
#include "en_Materials.hpp"
#include "cl_Material.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_DofManagerBase.hpp"
#include "cl_FEM_Block.hpp"
#include "cl_FEM_BoundaryCondition.hpp"
#include "fn_IF_initialize_integration_points_on_facet.hpp"
#include "fn_IF_initialize_shape_function.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Group::Group(
                DofManagerBase * aParent,
                const GroupType aGroupType,
                const ElementType aElementType,
                const id_t aID,
                const index_t aNumberOfElements,
                const bool aOwnElements ) :
                mParent( aParent ),
                mType( aGroupType ),
                mElementType( aElementType ),
                mID( aID ),
                mNumberOfElements( aNumberOfElements ),
                mOwnElements( aOwnElements ),
                mMyRank( comm_rank()),
                mMeshID( aID )
        {
            this->link_material_functions();
            this->create_calculator();
        }

//------------------------------------------------------------------------------

        void
        Group::delete_pointers()
        {
            if( mCalc != nullptr )
            {
                delete mCalc ;
            }

            if ( mOwnElements )
            {
                for ( auto tElement: mElements )
                {
                    delete tElement;
                }
            }
        }

//------------------------------------------------------------------------------

        Vector< real > &
        Group::field_data( const string & aLabel )
        {
            return mParent->field_data( aLabel );
        }

//------------------------------------------------------------------------------

        void
        Group::set_material( const MaterialType aMaterial )
        {
            // get material object from kernel. If it didn't exist, it is created
            mMaterial = mParent->parent()->get_material( aMaterial );

            BELFEM_ASSERT( mMaterial != nullptr, "No material created" );
        }

//------------------------------------------------------------------------------

        void
        Group::set_material( Material * aMaterial )
        {
            mMaterial = aMaterial;
        }

//------------------------------------------------------------------------------

        void
        Group::create_element_map()
        {
            // reset map
            mElementMap.clear();

            // loop over all elements
            for ( Element * tElement: mElements )
            {
                // add element to map
                mElementMap[ tElement->element()->id() ] = tElement;
            }
        }

//------------------------------------------------------------------------------

        void
        Group::link_material_functions()
        {
            if( mParent == nullptr )
            {
                return;
            }

            uint tNumDim = mParent->mesh()->number_of_dimensions();

            if ( tNumDim == 2 )
            {
                mThermalConductivity
                        = &Group::thermal_conductivity_2d;

                mLinearElasticity
                        = &Group::linear_elasticity_PlaneStress;
            }
            else if ( tNumDim == 3 )
            {
                mThermalConductivity
                        = &Group::thermal_conductivity_3d;

                mLinearElasticity
                        = &Group::linear_elasticity_3d;
            }
            else
            {
                BELFEM_ERROR( false, "Invalid number of dimensions: %u",
                             ( unsigned int ) tNumDim );

            }
        }

//------------------------------------------------------------------------------

        void
        Group::create_calculator()
        {

            if ( mParent->iwg() != nullptr && mCalc == nullptr )
            {
                // create the calculator object
                mCalc = new Calculator( this, mParent->iwg()->model_dimensionality() );

            }
        }


//------------------------------------------------------------------------------

        IntegrationData *
        Group::master_integration( const uint aSideSetIndex )
        {
            BELFEM_ERROR( false,
                         "Group::volume_integration() not implemented for this class." );
            return nullptr;
        }

//------------------------------------------------------------------------------

        IntegrationData *
        Group::slave_integration( const uint aSideSetIndex )
        {
            BELFEM_ERROR( false,
                         "Group::volume_integration() not implemented for this class." );
            return nullptr;
        }

//------------------------------------------------------------------------------

        ElementType
        Group::master_type() const
        {
            BELFEM_ERROR( false,
                         "Group::master_type() not implemented for this class." );
            return ElementType::UNDEFINED;
        }

//------------------------------------------------------------------------------

        ElementType
        Group::slave_type() const
        {
            BELFEM_ERROR( false,
                         "Group::slave_type() not implemented for this class." );
            return ElementType::UNDEFINED;
        }

//------------------------------------------------------------------------------

        void
        Group::set_boundary_condition( const BoundaryCondition * aBoundaryCondition )
        {
            mBoundaryCondition = aBoundaryCondition;
        }

//------------------------------------------------------------------------------

        const BoundaryCondition *
        Group::boundary_condition() const
        {
            BELFEM_ASSERT( mBoundaryCondition != nullptr,
                         "Boundary condition requested but not set." );
            return mBoundaryCondition ;
        }

//------------------------------------------------------------------------------

        uint
        Group::number_of_thin_shell_layers() const
        {
            BELFEM_ERROR( false, "number_of_thin_shell_layers() can only be called for sidesets or shells");
            return 0 ;
        }

//------------------------------------------------------------------------------

        uint
        Group::number_of_ghost_sidesets() const
        {
            BELFEM_ERROR( false, "number_of_ghost_sidesets() can only be called for sidesets or shells");
            return 0 ;
        }

//------------------------------------------------------------------------------

        const Material *
        Group::thin_shell_material( const uint aLayerIndex ) const
        {
            BELFEM_ERROR( false, "thin_shell_material() can only be called for sidesets or shells");
            return nullptr ;
        }

//------------------------------------------------------------------------------

        /**
         * dummy function, throws error unless sideset or shell
         * @return
         */
        real
        Group::thin_shell_thickness( const uint aLayerIndex ) const
        {
            BELFEM_ERROR( false, "thin_shell_thickness() can only be called for sidesets or shells");
            return BELFEM_QUIET_NAN ;
        }

//------------------------------------------------------------------------------

        /**
         * dummy function, throws error unless sideset or shell
         * @return
         */
        real
        Group::thin_shell_thickness() const
        {
            BELFEM_ERROR( false, "thin_shell_thickness() can only be called for sidesets or shells");
            return BELFEM_QUIET_NAN ;
        }

//------------------------------------------------------------------------------

        const IntegrationData *
        Group::thinshell_integration() const
        {
            BELFEM_ERROR( false,
                           "Function not implemented for this group type" );

            // this will never happen
            return nullptr ;

        }

//------------------------------------------------------------------------------
    }
}