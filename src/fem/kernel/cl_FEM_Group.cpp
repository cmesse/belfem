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

#include "cl_FEM_Group.hpp"
#include "assert.hpp"
#include "commtools.hpp"
#include "meshtools.hpp"
#include "cl_Material.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_FEM_DofManagerBase.hpp"
#include "cl_FEM_Block.hpp"
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
                mCommRank( comm_rank()),
                mMeshID( aID )
        {
            if ( aParent != nullptr )
            {
                this->create_calculator();
            }
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
                for ( auto tElement : mAuraElements )
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
        Group::set_material( Material * aMaterial )
        {
            mMaterial = aMaterial;
        }

 //------------------------------------------------------------------------------

        void
        Group::set_material( const string & aLabel )
        {
            this->set_material( mParent->parent()->material( aLabel ) );
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
            for ( Element * tElement: mAuraElements )
            {
                // add element to map
                mElementMap[ tElement->element()->id() ] = tElement;
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
                         "Group::master_integration() not implemented for this class." );
            return nullptr;
        }

//------------------------------------------------------------------------------

        IntegrationData *
        Group::slave_integration( const uint aSideSetIndex )
        {
            BELFEM_ERROR( false,
                         "Group::slave_integration() not implemented for this class." );
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
                           "thinshell_integration() not implemented for this group type" );

            // this will never happen
            return nullptr ;

        }

//------------------------------------------------------------------------------

        void
        Group::initialize_lookup_tables( const uint aOrder )
        {

            BELFEM_ERROR( false,
                           "initialize_lookup_tables() not implemented for this group type" );
        }

//------------------------------------------------------------------------------
    }
}
