//
// Created by christian on 7/13/21.
//

#include "cl_FEM_DofManagerBase.hpp"
#include "commtools.hpp"
#include "cl_FEM_Kernel.hpp"
#include "assert.hpp"


namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        DofManagerBase::DofManagerBase(  const DofManagerType  aType,
                                         Kernel        * aParent  ) :
            mType( aType ),
            mParent( aParent ),
            mMesh( aParent->mesh() ),
            mMyRank( comm_rank() )
        {
        }

//------------------------------------------------------------------------------

        Dof *
        DofManagerBase::dof( const id_t aID )
        {
            BELFEM_ERROR( false, "invalid call to abstract base class function." );
            return nullptr ;
        }

//------------------------------------------------------------------------------

        bool
        DofManagerBase::enforce_linear_interpolation() const
        {
            BELFEM_ERROR( false, "invalid call to abstract base class function." );
            return false ;
        }

//------------------------------------------------------------------------------

        id_t
        DofManagerBase::calculate_dof_id( const mesh::Node * aNode , const uint aDofType ) const
        {
            BELFEM_ERROR( false, "invalid call to abstract base class function." );
            return gNoID ;
        }

//------------------------------------------------------------------------------

        id_t
        DofManagerBase::calculate_dof_id( const mesh::Edge * aEdge , const uint aDofType ) const
        {
            BELFEM_ERROR( false, "invalid call to abstract base class function." );
            return gNoID ;
        }

//------------------------------------------------------------------------------

        id_t
        DofManagerBase::calculate_dof_id( const mesh::Facet * aFacet , const uint aDofType ) const
        {
            BELFEM_ERROR( false, "invalid call to abstract base class function." );
            return gNoID ;
        }

//------------------------------------------------------------------------------

        Block *
        DofManagerBase::block( const id_t aID )
        {
            BELFEM_ERROR( false, "invalid call to abstract base class function." );
            return nullptr ;
        }

//------------------------------------------------------------------------------

        SideSet *
        DofManagerBase::sideset( const id_t aID )
        {
            BELFEM_ERROR( false, "invalid call to abstract base class function." );
            return nullptr ;
        }

//------------------------------------------------------------------------------

        uint
        DofManagerBase::number_of_dofs_per_node() const
        {
            BELFEM_ERROR( false, "invalid call to abstract base class function." );
            return BELFEM_UINT_MAX ;
        }

//------------------------------------------------------------------------------

        uint
        DofManagerBase::number_of_dofs_per_edge() const
        {
            BELFEM_ERROR( false, "invalid call to abstract base class function." );
            return BELFEM_UINT_MAX ;
        }

//------------------------------------------------------------------------------


        uint
        DofManagerBase::sideset_integration_order() const
        {
            BELFEM_ERROR( false, "invalid call to abstract base class function." );
            return BELFEM_UINT_MAX ;
        }

//------------------------------------------------------------------------------

        uint
        DofManagerBase::block_integration_order() const{
            BELFEM_ERROR( false, "invalid call to abstract base class function." );
            return BELFEM_UINT_MAX ;
        }

//------------------------------------------------------------------------------

        IntegrationScheme
        DofManagerBase::integration_scheme() const
        {
            BELFEM_ERROR( false, "invalid call to abstract base class function." );
            return IntegrationScheme::UNDEFINED ;
        }

//------------------------------------------------------------------------------

        bool
        DofManagerBase::block_exists( const id_t & aID ) const
        {
            BELFEM_ERROR( false, "invalid call to abstract base class function." );
            return false ;
        }

//------------------------------------------------------------------------------

        bool
        DofManagerBase::sideset_exists( const id_t & aID ) const
        {
            BELFEM_ERROR( false, "invalid call to abstract base class function." );
            return false ;
        }

//------------------------------------------------------------------------------

        void
        DofManagerBase::collect_fields( const Cell< string > & aFieldLabels )
        {
            BELFEM_ERROR( false, "invalid call to abstract base class function." );
        }

//------------------------------------------------------------------------------

        void
        DofManagerBase::distribute_fields( const Cell< string > & aFieldLabels )
        {
            BELFEM_ERROR( false, "invalid call to abstract base class function." );
        }

//------------------------------------------------------------------------------

        void
        DofManagerBase::synchronize_fields( const Cell< string > & aFieldLabels )
        {
            BELFEM_ERROR( false, "invalid call to abstract base class function." );
        }

//------------------------------------------------------------------------------

        void
        DofManagerBase::print_worst_dof()
        {
            BELFEM_ERROR( false, "invalid call to abstract base class function." );
        }

//------------------------------------------------------------------------------

        void
        DofManagerBase::write_residuals_to_mesh()
        {
            BELFEM_ERROR( false, "invalid call to abstract base class function." );
        }

//------------------------------------------------------------------------------

        bool
        DofManagerBase::dof_exists( const id_t aID ) const
        {
            BELFEM_ERROR( false, "invalid call to abstract base class function." );
            return false ;
        }

//------------------------------------------------------------------------------

        bool
        DofManagerBase::is_master() const
        {
            return mParent->is_master() ;
        }

//------------------------------------------------------------------------------
    }
}