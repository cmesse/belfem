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
         * abstract interface of the dof manager as seen by the group classes
         * and the IWG; DofManager is its only implementation
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

            const proc_t mCommRank ;

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

           virtual void
           initialize();

//------------------------------------------------------------------------------

            /**
             * true once initialize() has frozen the free/fixed dof split and
             * the matrix graphs; a first fix() of a free dof after that point
             * changes a classification the solver containers are sized for
             * ( see SideSet::impose_dirichlet )
             */
            virtual bool
            is_initialized() const
            {
                return false ;
            }

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
            return mCommRank;
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
