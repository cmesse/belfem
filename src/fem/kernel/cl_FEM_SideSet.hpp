//
// Created by Christian Messe on 03.11.19.
//

#ifndef BELFEM_CL_FEM_SIDESET_HPP
#define BELFEM_CL_FEM_SIDESET_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_SideSet.hpp"
#include "cl_FEM_Element.hpp"
#include "cl_IF_InterpolationFunction.hpp"
#include "cl_FEM_Group.hpp"
#include "en_FEM_BoundaryConditionImposing.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class Field ;
        class DofManager ;
        class Block ;

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

            // thin shell materials, owned by kernel
            Cell< Material * > mThinShellMaterials ;

            // thin shell thicknesses
            Vector< real > mThinShellThicknesses ;
            real mThinShellThickness = BELFEM_QUIET_NAN ;

            uint mNumberOfThinShellLayers = 0 ;
            uint mNumberOfGhostSideSets = 0 ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            // deprecated
            SideSet(                Field * aParent,
                    const id_t              aID,
                    Cell< mesh::Facet * > & aFacets );

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

            virtual ~SideSet();

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

            ElementType
            element_type() const;

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
             master_type() const ;
//------------------------------------------------------------------------------

             /**
              * returns type of sideset slave
              */
             ElementType
             slave_type() const ;

//------------------------------------------------------------------------------

            /**
             * integration data on master element
             * @param aSideSetIndex
             * @return
             */
            const IntegrationData *
            master_integration( const uint aSideSetIndex ) ;

//------------------------------------------------------------------------------

            /**
             * integration data on slave element
             * @param aSideSetIndex
             * @return
             */
            const IntegrationData *
            slave_integration( const uint aSideSetIndex ) ;

//------------------------------------------------------------------------------

            void
            set_thin_shells( Cell< Material* > & aMaterials,
                             const Vector< real > & aLayerThicknesses );

//------------------------------------------------------------------------------

            inline uint
            number_of_thin_shell_layers() const ;

//------------------------------------------------------------------------------

            inline uint
            number_of_ghost_sidesets() const ;

//------------------------------------------------------------------------------

            inline const Material *
            thin_shell_material( const uint aLayerIndex ) const ;

//------------------------------------------------------------------------------

            inline real
            thin_shell_thickness( const uint aLayerIndex ) const ;

//------------------------------------------------------------------------------

            inline real
            thin_shell_thickness() const ;

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

            void
            initialize_lookup_tables( const uint aIntegrationOrder );

//------------------------------------------------------------------------------

            void
            initialize_work_matrices();

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------

        inline ElementType
        SideSet::element_type() const
        {
            return mElementType;
        }
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

        inline const IntegrationData *
        SideSet::master_integration( const uint aSideSetIndex )
        {
            return mMasterIntegration( aSideSetIndex );
        }

//------------------------------------------------------------------------------

        inline const IntegrationData *
        SideSet::slave_integration( const uint aSideSetIndex )
        {
            return mSlaveIntegration( aSideSetIndex );
        }

//------------------------------------------------------------------------------

        inline uint
        SideSet::number_of_thin_shell_layers() const
        {
            return mNumberOfThinShellLayers ;
        }

//------------------------------------------------------------------------------

        inline uint
        SideSet::number_of_ghost_sidesets() const
        {
            return mNumberOfGhostSideSets ;
        }

//------------------------------------------------------------------------------

        inline const Material *
        SideSet::thin_shell_material( const uint aLayerIndex ) const
        {
            return mThinShellMaterials( aLayerIndex );
        }

//------------------------------------------------------------------------------

        inline real
        SideSet::thin_shell_thickness( const uint aLayerIndex ) const
        {
            return mThinShellThicknesses( aLayerIndex );
        }

//------------------------------------------------------------------------------

        inline real
        SideSet::thin_shell_thickness() const
        {
            return mThinShellThickness ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_FEM_SIDESET_HPP
