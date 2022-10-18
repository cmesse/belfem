//
// Created by christian on 7/13/21.
//

#ifndef BELFEM_CL_FEM_DOFMANAGERBASE_HPP
#define BELFEM_CL_FEM_DOFMANAGERBASE_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Mesh.hpp"
#include "en_IntegrationScheme.hpp"

namespace belfem
{
    class Mesh ;

    enum class DofManagerType
    {
        OLD,
        NEW,
        UNDEFINED
    };

    namespace fem
    {
        class Kernel;
        class IWG ;
        class Block ;
        class SideSet ;
        class Dof ;

        /**
         * this class is a preliminary baseclass so that both the old field and the new
         * dof manager can be used with the group class
         */
        class DofManagerBase
        {
            const DofManagerType mType ;

//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            Kernel * mParent ;

            Mesh   * mMesh ;

            IWG  * mIWG = nullptr ;

            const proc_t mMyRank ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            DofManagerBase( const DofManagerType  aType,
                                  Kernel        * aParent ) ;

//------------------------------------------------------------------------------

            virtual ~DofManagerBase() = default ;

//--------------------------------------------------------------------------

            DofManagerType
            type() const ;

//------------------------------------------------------------------------------

            /**
             * return the rank of this proc ( faster than comm_rank() )
             */
            proc_t
            rank() const ;

//------------------------------------------------------------------------------

            Mesh *
            mesh();

//------------------------------------------------------------------------------

            Kernel *
            parent();

//------------------------------------------------------------------------------

            IWG *
            iwg();

//------------------------------------------------------------------------------

            Vector< real > &
            field_data( const string & aLabel ) ;

//------------------------------------------------------------------------------

            /**
             * return a specific dof
             */
            virtual Dof *
            dof( const id_t aID );

//------------------------------------------------------------------------------

            virtual bool
            enforce_linear_interpolation() const ;

//------------------------------------------------------------------------------

            virtual id_t
            calculate_dof_id( const mesh::Node * aNode , const uint aDofType )  const;

//------------------------------------------------------------------------------

            virtual id_t
            calculate_dof_id( const mesh::Edge * aEdge , const uint aDofType )  const;

//------------------------------------------------------------------------------

            virtual id_t
            calculate_dof_id( const mesh::Facet * aFacet , const uint aDofType )  const;

//------------------------------------------------------------------------------

            virtual Block *
            block( const id_t aID );

//------------------------------------------------------------------------------

            virtual  SideSet *
            sideset( const id_t aID ) ;

//------------------------------------------------------------------------------

            virtual uint
            number_of_dofs_per_node() const ;

//------------------------------------------------------------------------------

            virtual uint
            number_of_dofs_per_edge() const ;

//------------------------------------------------------------------------------

            virtual uint
            sideset_integration_order() const ;

//------------------------------------------------------------------------------

            virtual uint
            block_integration_order() const ;

//------------------------------------------------------------------------------

            virtual IntegrationScheme
            integration_scheme() const ;

//------------------------------------------------------------------------------

            virtual bool
            block_exists( const id_t & aID ) const ;

//------------------------------------------------------------------------------

            virtual bool
            sideset_exists( const id_t & aID ) const ;

//------------------------------------------------------------------------------

            virtual void
            collect_fields( const Cell< string > & aFieldLabels ) ;

//------------------------------------------------------------------------------

            virtual void
            distribute_fields( const Cell< string > & aFieldLabels ) ;

//------------------------------------------------------------------------------

            /**
             * collect and distribute afterwards
             * @param aFieldLabels
             */
            virtual void
            synchronize_fields( const Cell< string > & aFieldLabels ) ;

//------------------------------------------------------------------------------

            bool
            is_master() const ;

//------------------------------------------------------------------------------

            virtual bool
            dof_exists( const id_t aID ) const ;

//----------------------------------------------------------------------------

            virtual void
            print_worst_dof() ;

//----------------------------------------------------------------------------

            virtual void
            write_residuals_to_mesh();

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline DofManagerType
        DofManagerBase::type() const
        {
            return mType ;
        }


//------------------------------------------------------------------------------

        inline proc_t
        DofManagerBase::rank() const
        {
            return mMyRank;
        }

//------------------------------------------------------------------------------
        inline Mesh *
        DofManagerBase::mesh()
        {
            return mMesh ;
        }

//---------------------------------------------------------------------------

        inline Kernel *
        DofManagerBase::parent()
        {
            return mParent ;
        }

//------------------------------------------------------------------------------

        inline IWG *
        DofManagerBase::iwg()
        {
            return mIWG ;
        }

//---------------------------------------------------------------------------

        inline Vector< real > &
        DofManagerBase::field_data( const string & aLabel )
        {
            return mMesh->field_data( aLabel );
        }

//---------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_FEM_DOFMANAGERBASE_HPP
