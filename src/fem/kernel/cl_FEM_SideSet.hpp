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

#ifndef BELFEM_CL_FEM_SIDESET_HPP
#define BELFEM_CL_FEM_SIDESET_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_FEM_Group.hpp"
#include "en_FEM_BoundaryConditionImposing.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class Field ;
        class DofManager ;
        class DofManagerBase ;
        class Dof ;
        class Block ;

//------------------------------------------------------------------------------

        /**
         * shared Dirichlet pin for a single node dof, used by
         * SideSet::impose_dirichlet and the Maxwell factory pre-fix loop so
         * both sites apply ONE shape test. If the dof is condensed
         * onto a single source with unit weight, the SOURCE is fixed to
         * aValue / weight and the function returns true ( incrementing
         * aNumFirstFlips when the source was still free ); every other shape
         * returns false and the CALLER keeps its historical behaviour.
         * BELFEM_ERRORs on a non-unit single-source weight, a broken source,
         * or a first free-to-fixed flip after aParent->is_initialized().
         * Call only for non-duplicate nodes: duplicate-pair condensations
         * ( cuts, thin shells ) can carry an inhomogeneous relation that a
         * plain source pin would violate.
         */
        bool
        pin_dirichlet_dof(
                const DofManagerBase * aParent,
                      Dof            * aDof,
                const real             aValue,
                const id_t             aGroupID,
                const id_t             aNodeID,
                      index_t        & aNumFirstFlips );

//------------------------------------------------------------------------------

        class SideSet  : public Group
        {
            // cell with sideset integration information
            Cell< IntegrationData * > mSideSetIntegrationData ;


//----------------------------------------------------------------------------
        protected:
//----------------------------------------------------------------------------

            const ElementType mMasterType = ElementType::EMPTY ;
            const ElementType mSlaveType  = ElementType::EMPTY ;

            // side nodes
            Cell< mesh::Node * > mNodes;

            // container for BC values per dof type
            Vector< real > mBcValues ;

            Cell< BoundaryConditionImposing > mBcTypes ;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------
            // container for temperature for alpha BC
            real mTinf = BELFEM_QUIET_NAN ;

            // conainer with reference blocks for sideset integrations
            Cell< IntegrationData * > mMasterIntegration ;
            Cell< IntegrationData * > mSlaveIntegration ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            SideSet( DofManager * aParent,
                     const id_t aID,
                     Cell< mesh::Facet * > & aFacets,
                     const GroupType aGroupType = GroupType::SIDESET );

//------------------------------------------------------------------------------

            /*
             *  a constructor with an empty parent, just used for tests
             */
            SideSet( const ElementType aElementType,
                     const ElementType aMasterType,
                     const ElementType aSlaveType,
                     const GroupType aGroupType = GroupType::SIDESET  );

//------------------------------------------------------------------------------

            ~SideSet() override;

//------------------------------------------------------------------------------

            void
            impose_dirichlet( const real aValue, const uint aDofType=0 );

//------------------------------------------------------------------------------

            void
            impose_neumann( const real aValue, const uint aDofType=0 );

//------------------------------------------------------------------------------

            /**
             * imposing an alpha value requires setting the value in the field
             * by default, a value can be set, but the alpha and T0 field can be
             * overwritten
             */
            void
            impose_alpha( const real aAlpha = BELFEM_QUIET_NAN, const real aTinf=BELFEM_QUIET_NAN );

//------------------------------------------------------------------------------

            void
            free();

//------------------------------------------------------------------------------

            void
            set_boundary_conditions();

//------------------------------------------------------------------------------

            // get bc type
            BoundaryConditionImposing
            bc_type( const index_t & aDimension ) const ;

//------------------------------------------------------------------------------

            // get bc type
            const real &
            bc_value( const index_t & aDimension ) const ;

//------------------------------------------------------------------------------

            /**
             * expose the node container
             */
             Cell< mesh::Node * > &
             nodes() ;

//------------------------------------------------------------------------------

            uint
            number_of_boundary_conditions() const ;

//------------------------------------------------------------------------------

            /**
             * returns type of sideset master
             */
             ElementType
             master_type() const override ;

//------------------------------------------------------------------------------

             /**
              * returns type of sideset slave
              */
             ElementType
             slave_type() const override ;

//------------------------------------------------------------------------------

            /**
             * integration data on master element
             * @param aSideSetIndex
             * @return
             */
            IntegrationData *
            master_integration( const uint aSideSetIndex ) override ;

//------------------------------------------------------------------------------

            /**
             * integration data on slave element
             * @param aSideSetIndex
             * @return
             */
            IntegrationData *
            slave_integration( const uint aSideSetIndex ) override ;

//------------------------------------------------------------------------------

            void
            initialize_lookup_tables( const uint aIntegrationOrder ) override;

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            // surface elements
            virtual void
            initialize_elements( Cell< mesh::Facet * > & aFacets  );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            collect_nodes( Cell< mesh::Facet * > & aFacets );


//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------

        // get bc type
        inline BoundaryConditionImposing
        SideSet::bc_type( const index_t & aDimension ) const
        {
            return mBcTypes( aDimension );
        }

//------------------------------------------------------------------------------

        // get bc type
        inline const real &
        SideSet::bc_value( const index_t & aDimension ) const
        {
            return mBcValues( aDimension );
        }

//------------------------------------------------------------------------------

        inline
        Cell< mesh::Node * > &
        SideSet::nodes()
        {
            return mNodes ;
        }

//------------------------------------------------------------------------------

        inline uint
        SideSet::number_of_boundary_conditions() const
        {
            return mBcTypes.size() ;
        }

//------------------------------------------------------------------------------

        inline ElementType
        SideSet::master_type() const
        {
            return mMasterType ;
        }

//------------------------------------------------------------------------------

        inline ElementType
        SideSet::slave_type() const
        {
            return mSlaveType ;
        }

//------------------------------------------------------------------------------

        inline IntegrationData *
        SideSet::master_integration( const uint aSideSetIndex )
        {
            return mMasterIntegration( aSideSetIndex );
        }

//------------------------------------------------------------------------------

        inline IntegrationData *
        SideSet::slave_integration( const uint aSideSetIndex )
        {
            return mSlaveIntegration( aSideSetIndex );
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_FEM_SIDESET_HPP
