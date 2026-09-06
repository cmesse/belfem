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

#ifndef BELFEM_CL_FEM_GROUP_HPP
#define BELFEM_CL_FEM_GROUP_HPP

#include "typedefs.hpp"
#include "constants.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "Mesh_Enums.hpp"
#include "cl_Mesh.hpp"
#include "meshtools.hpp"
#include "cl_IF_InterpolationFunction.hpp"
#include "cl_Material.hpp"
#include "en_DomainType.hpp"
#include "en_FEM_GroupActivationMode.hpp"
#include "cl_IF_IntegrationData.hpp"
#include "cl_FEM_Calculator.hpp"

namespace belfem
{
    class Material;
    enum class MaterialType;

    namespace fem
    {
        class DofManagerBase;
        class Element;
        class Block ;
        class BoundaryCondition ;

//------------------------------------------------------------------------------

        class Group
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            // pointer to parent
            DofManagerBase * mParent ;

            // pointer to swap object
            Calculator * mCalc = nullptr ;

            // either block or sideset
            const GroupType mType ;

            const ElementType mElementType;

            const id_t mID;

            const index_t mNumberOfElements;

            // flag that tells if elements are destoyed by destructor
            const bool mOwnElements ;

            const proc_t  mCommRank;

            // detailed domain type, mainly used for Maxwell
            DomainType mDomainType = DomainType::Default ;

            // container for node coordinates
            Matrix< real > mNodeCoords ;

            // must be set by child
            uint mNumberOfNodesPerElement = BELFEM_UINT_MAX ;

            uint mIntegrationOrder = 0 ;

            //bool mIsIsogeometric = true ;

            //IntegrationData * mIntegrationData = nullptr ;

            // this is meant if the geometry data are different
            //IntegrationData * mGeometryIntegrationData = nullptr ;

            // shape function for element boundaries
            Cell< Cell< Matrix < real > > > mBoundaryN ;

            // Element container
            Cell< Element * > mElements;

            Cell< Element * > mAuraElements;

            // pointer to material ( owned by kernel )
            Material * mMaterial = nullptr;

            // pointer to boundary condition, if set
            const BoundaryCondition * mBoundaryCondition = nullptr ;

            Map< id_t, Element * > mElementMap ; // map to access element by ID

            // empty sidesets and blocks have the id zero.
            // the fake ID helps to access the underlying objects on the mesh
            id_t mMeshID ;

            //! flag telling if block has a right hand side
            bool mHasRHS = true ;

            //! activation mode determining how this group is used for computation
            GroupActivationMode mActivationMode = GroupActivationMode::GeometryAndDofs ;

            //! in case for example if bubble functions are used
            Cell< IntegrationData * > mEnrichmentData ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Group(  DofManagerBase      * aParent,
                    const GroupType     aGroupType,
                    const ElementType   aElementType,
                    const        id_t   aID,
                    const     index_t   aNumberOfElements,
                    const        bool   aOwnElements=true );

//------------------------------------------------------------------------------

            virtual ~Group() = default ;

//------------------------------------------------------------------------------

            /**
             * return the type of this group
             */
             GroupType
             type() const ;

//------------------------------------------------------------------------------

            /**
             * return the id of the group
             */
             id_t
             id() const;

//------------------------------------------------------------------------------

            /**
             * return the element type of this block
             */
            virtual ElementType
            element_type() const;

//------------------------------------------------------------------------------

            /**
             * return the field of this block
             */
            DofManagerBase *
            parent();

//------------------------------------------------------------------------------

            /*
             * return the number of elements of this block
             */
            index_t
            number_of_elements() const;

//------------------------------------------------------------------------------

            /*
             * return the number of nodes per element
             */
            uint
            number_of_nodes_per_element() const;

//------------------------------------------------------------------------------

            /*
             * return the number of edges per element
             */
            uint
            number_of_edges_per_element() const;

//------------------------------------------------------------------------------

            /*
             * return the number of faces per element
             */
            uint
            number_of_faces_per_element() const;

//------------------------------------------------------------------------------

            /**
             * expose the integration data
             */
             const IntegrationData *
             integration() const ;

//------------------------------------------------------------------------------

            const IntegrationData *
            enrichment_data( const uint aFacet ) const ;

//------------------------------------------------------------------------------

            virtual const IntegrationData *
            thinshell_integration() const ;

//------------------------------------------------------------------------------

            /*
             * expose element container
             */
            Cell< Element * > &
            elements();

            Cell< Element * > &
            aura_elements();

//------------------------------------------------------------------------------

            bool
            element_exists( const id_t aID );

//------------------------------------------------------------------------------

            /**
             * expose calculator object
             */
             Calculator *
             calculator() ;

//------------------------------------------------------------------------------

            /**
             * expose container for node coords
             */
             Matrix< real > &
             node_coords() ;

//------------------------------------------------------------------------------

            /**
             * access scal
             * ars from the parent
             */
             Vector< real > &
             field_data( const string & aLabel );

//------------------------------------------------------------------------------

            void
            set_material( Material * aMaterial );

//------------------------------------------------------------------------------

            void
            set_material( const string & aLabel );

//------------------------------------------------------------------------------

            Material *
            material();

//------------------------------------------------------------------------------

            const Material *
            material() const;

//------------------------------------------------------------------------------

            // needed if you want to access element by id
            void
            create_element_map() ;

//------------------------------------------------------------------------------

            // get an element using its ID instead of index
            Element *
            element( const id_t aID );

//------------------------------------------------------------------------------
// RHS flags
//------------------------------------------------------------------------------

            // flag telling field if an RHS side exists
            bool
            has_rhs() const ;

            // set or unset the rhs flag
            void
            set_rhs_flag( const bool aFlag );

//------------------------------------------------------------------------------
// Helpers
//------------------------------------------------------------------------------

            /**
             * if the sideset or block is empty, ID is zero.
             * use mesh_id if you really need the ID of the corresponding
             * mesh object
             *
             * @return the ID of the corresponding mesh object
             */
            const id_t &
            mesh_id() const ;

            void
            set_mesh_id( const id_t & aID );

//------------------------------------------------------------------------------

            /**
             * for special purpose integration, sideset only
             */
            virtual IntegrationData *
            master_integration( const uint aSideSetIndex );

            virtual IntegrationData *
            slave_integration( const uint aSideSetIndex );

            void
            delete_pointers();

//------------------------------------------------------------------------------

            const InterpolationFunction *
            interpolation_function() const ;

//------------------------------------------------------------------------------

            void
            set_domain_type( const DomainType aType );

//------------------------------------------------------------------------------

            DomainType
            domain_type() const ;

//------------------------------------------------------------------------------

            /**
             * returns type of sideset master
             */
            virtual ElementType
            master_type() const ;
//------------------------------------------------------------------------------

            /**
             * returns type of sideset slave
             */
            virtual ElementType
            slave_type() const ;

//------------------------------------------------------------------------------

            /**
             * dummy function, throws error unless sideset or shell
             * @return
             */
            virtual uint
            number_of_thin_shell_layers() const ;

//------------------------------------------------------------------------------

            /**
             * dummy function, throws error unless sideset or shell
             * @return
             */
            virtual uint
            number_of_ghost_sidesets() const ;

//------------------------------------------------------------------------------

            /**
             * dummy function, throws error unless sideset or shell
             * @return
             */
            virtual real
            thin_shell_thickness( const uint aLayerIndex ) const ;

//------------------------------------------------------------------------------

            /**
             * dummy function, throws error unless sideset or shell
             * @return
             */
            virtual real
            thin_shell_thickness() const ;

//------------------------------------------------------------------------------

            /**
             * function telling if values from this group are computed
             * ( default : true )
             * Returns true for GeometryAndDofs, false otherwise
             */
             bool
             is_active() const ;

//------------------------------------------------------------------------------

            /**
             * returns the activation mode of this group
             */
             GroupActivationMode
             activation_mode() const ;

//------------------------------------------------------------------------------

            /**
             * set the activation mode of the group
             */
            void
            set_activation_mode( const GroupActivationMode aMode );

//------------------------------------------------------------------------------

            /**
             * set the active flag of the group (backward compatibility)
             * Maps bool to GeometryAndDofs (true) or Inactive (false)
             */
            void
            activate( bool aFlag );

//------------------------------------------------------------------------------

            virtual void
            initialize_lookup_tables( const uint aOrder );

            bool
            has_enrichment();

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------
            void
            create_calculator();

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline DofManagerBase *
        Group::parent()
        {
            return mParent;
        }

//------------------------------------------------------------------------------

        inline GroupType
        Group::type() const
        {
            return mType ;
        }

//------------------------------------------------------------------------------

        inline id_t
        Group::id() const
        {
            return mID;
        }

//------------------------------------------------------------------------------

        inline index_t
        Group::number_of_elements() const
        {
            return mNumberOfElements;
        }

//------------------------------------------------------------------------------

        inline uint
        Group::number_of_nodes_per_element() const
        {
            return mNumberOfNodesPerElement ;
        }

//------------------------------------------------------------------------------

        inline uint
        Group::number_of_edges_per_element() const
        {
            return mesh::number_of_edges( this->element_type() );
        }

//------------------------------------------------------------------------------

        inline uint
        Group::number_of_faces_per_element() const
        {
            return mesh::number_of_faces( this->element_type() );
        }

//------------------------------------------------------------------------------

        inline Cell< Element * > &
        Group::elements()
        {
            return mElements;
        }

//------------------------------------------------------------------------------

        inline Cell< Element * > &
        Group::aura_elements()
        {
            return mAuraElements;
        }

//------------------------------------------------------------------------------

        inline Matrix< real > &
        Group::node_coords()
        {
            return mNodeCoords ;
        }

//------------------------------------------------------------------------------

        // get an element using its ID instead of index
        inline Element *
        Group::element( const id_t aID )
        {
            return mElementMap( aID );
        }

//------------------------------------------------------------------------------

        inline  Material *
        Group::material()
        {
            return mMaterial;
        }

//------------------------------------------------------------------------------

        inline  const Material *
        Group::material() const
        {
            return mMaterial;
        }

//------------------------------------------------------------------------------

        // flag telling field if an RHS side exists
        inline bool
        Group::has_rhs() const
        {
            return mHasRHS ;
        }

//------------------------------------------------------------------------------

        // set or unset the rhs flag
        inline void
        Group::set_rhs_flag( const bool aFlag )
        {
            mHasRHS = aFlag ;
        }

//------------------------------------------------------------------------------

        inline const id_t &
        Group::mesh_id() const
        {
            return mMeshID ;
        }

//------------------------------------------------------------------------------

        inline void
        Group::set_mesh_id( const id_t & aID )
        {
            mMeshID = aID ;
        }

//------------------------------------------------------------------------------


        inline void
        Group::set_domain_type( const DomainType aType )
        {
            mDomainType = aType ;
        }

//------------------------------------------------------------------------------

        inline DomainType
        Group::domain_type() const
        {
            return mDomainType ;
        }

//------------------------------------------------------------------------------

        inline bool
        Group::is_active() const
        {
            return mActivationMode == GroupActivationMode::GeometryAndDofs ;
        }

//------------------------------------------------------------------------------

        inline GroupActivationMode
        Group::activation_mode() const
        {
            return mActivationMode ;
        }

//------------------------------------------------------------------------------

        inline void
        Group::set_activation_mode( const GroupActivationMode aMode )
        {
            mActivationMode = aMode ;
        }

//------------------------------------------------------------------------------

        inline void
        Group::activate( bool aFlag )
        {
            mActivationMode = aFlag ? GroupActivationMode::GeometryAndDofs : GroupActivationMode::Inactive ;
        }

//------------------------------------------------------------------------------

        inline Calculator *
        Group::calculator()
        {
            return mCalc ;
        }

//------------------------------------------------------------------------------

        inline ElementType
        Group::element_type() const
        {
            return mElementType;
        }

//------------------------------------------------------------------------------

        inline const IntegrationData *
        Group::integration() const
        {
            return mCalc->integration() ;
        }

//------------------------------------------------------------------------------

        inline const IntegrationData *
        Group::enrichment_data( const uint aFacet ) const
        {
            return mEnrichmentData( aFacet );
        }

//------------------------------------------------------------------------------

        inline bool Group::has_enrichment()
        {
            return mEnrichmentData.size() > 0 ;
        }

//------------------------------------------------------------------------------

        inline bool
        Group::element_exists( const id_t aID )
        {
            return mElementMap.key_exists( aID );
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_FEM_GROUP_HPP
