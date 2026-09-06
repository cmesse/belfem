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
#ifndef BELFEM_CL_FEM_DOF_HPP
#define BELFEM_CL_FEM_DOF_HPP

#include "typedefs.hpp"
#include "cl_Graph_Vertex.hpp"
#include "cl_Node.hpp"
#include "cl_Edge.hpp"
#include "cl_Cell.hpp"
#include "cl_Facet.hpp"

namespace belfem
{
//------------------------------------------------------------------------------

    namespace fem
    {
        class Dof : public graph::Vertex
        {
            // type of the current dof, needed for vectors
            // temperature: always zero
            // displacement: 1: ux, 2: uy, 3:uz
            const uint mTypeID;

            // pointer to node or edge of this flag
            mesh::Basis * mMeshBasis ;

            // value for multiplicities, needed for parallel
            const uint mIndexOnEntity ;

            // index of dof on field
            const index_t mIndexOnField  ;

            // index of field on mesh, set by constructor
            index_t mFieldIndex = gNoIndex ;

            real mDirichletValue = 0.0 ; // BELFEM_QUIET_NAN;

            // tells if a Dirichlet boundary condition is imposed
            bool mFixedFlag = false;

            // index for DOF on current proc
            index_t mMyIndex = gNoIndex;

            uint    mNumberOfSources = 0 ;
            Dof  ** mSources = nullptr ;
            real *  mCoefficients = nullptr ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Dof( const id_t aID, const uint aType, mesh::Node * aNode );

//------------------------------------------------------------------------------

            Dof(
                const id_t aID,
                const uint aType,
                mesh::Edge * aEdge,
                const uint aIndexOnEdge,
                const index_t aDofIndexOnField );

//------------------------------------------------------------------------------

            Dof(const id_t aID,
                const uint aType,
                mesh::Face * aFace,
                uint aIndexOnFace,
                const index_t aDofIndexOnField );

//------------------------------------------------------------------------------

            Dof( const id_t aID,
                const uint aType,
                mesh::Element * aElement,
                const uint aIndexOnElement,
                const index_t aDofIndexOnField );

//------------------------------------------------------------------------------

            Dof( const id_t aID,
                const uint aType,
                mesh::Facet * aFacet ,
                const uint aIndexOnFacet,
                const index_t aDofIndexOnField );

//------------------------------------------------------------------------------

            ~Dof() override;

//------------------------------------------------------------------------------

            bool
            mesh_basis_is_flagged();

//------------------------------------------------------------------------------

            mesh::Basis *
            mesh_basis();

//------------------------------------------------------------------------------

            /**
             * expose the node that is linked to this dof
             */
            mesh::Node *
            node();

//------------------------------------------------------------------------------

            /**
             * expose the edge that is linked to this dof
             */
            mesh::Edge *
            edge();

//------------------------------------------------------------------------------

            /**
             * expose the face that is linked to this dof
             */
            mesh::Face *
            face();

//------------------------------------------------------------------------------

            /**
             * expose the element that is linked to this dof
             */
            mesh::Element *
            element();

//------------------------------------------------------------------------------

            /**
             * impose a Diriclet boundary condition on this node
             */
            void
            fix( const real aDirichletValue );

//------------------------------------------------------------------------------

            /**
             * free this dof ( remove Diriclet condition )
             */
            void
            free();

//------------------------------------------------------------------------------

            /**
             * tells if this dof is fixed
             */
            bool
            is_fixed() const;

//------------------------------------------------------------------------------

            /**
             * tells if this dof is hanging
             */
            bool
            is_hanging() const;

//------------------------------------------------------------------------------

            /**
             * tells if basis of dof is hanging, needed during initialization
             */
            bool
            basis_is_hanging() const;

//------------------------------------------------------------------------------

            /**
             * returns true if this dof is linked to a node
             */
            bool
            is_node() const;

//------------------------------------------------------------------------------

            /**
             * returns true if this dof is linked to an edge
             */
            bool
            is_edge() const;

//------------------------------------------------------------------------------

            /**
             * returns true if this dof is linked to a face
             */
            bool
            is_face() const;

//------------------------------------------------------------------------------

            /**
             * returns true if this dof is linked to a cell ( element-interior dof )
             */
            bool
            is_cell() const;

//------------------------------------------------------------------------------

            /**
             * returns true if this dof is linked to a facet
             */
            bool
            is_lambda() const;

//------------------------------------------------------------------------------

            EntityType
            entity_type() const ;

//------------------------------------------------------------------------------

            /**
             * return the Dirichlet value
             */
            real &
            value();

//------------------------------------------------------------------------------

            /**
             * access a dof
             */
            Dof *
            dof( const uint aIndex );

//------------------------------------------------------------------------------

            /**
             * return the type of this dof
             */
             uint
             type_id() const ;

//------------------------------------------------------------------------------

            /**
             * how many dofs are connected to this dof
             */
            uint
            number_of_dofs() const ;

//------------------------------------------------------------------------------

            /**
             *local index on mesh entity
             */
            uint
            index_on_entity() const ;

//------------------------------------------------------------------------------

            /**
             * set proc local index
             */
             void
             set_my_index( const index_t aIndex );
//------------------------------------------------------------------------------

            void
            set_field_index( const index_t aIndex );

//------------------------------------------------------------------------------

            /**
             * returns the corresponding index of the dof on the field
             * @return
             */
            index_t
            dof_index_on_field() const ;

//------------------------------------------------------------------------------

            /**
             * returns the index of the field on the mesh
             */
            index_t
            field_index() const ;

//------------------------------------------------------------------------------

            /**
             * get proc local index
             */
            index_t
            my_index() const;

//------------------------------------------------------------------------------

            void
            set_sources( Cell< Dof * > & aSources, Vector< real > & aWeights );

            void
            set_sources( Cell< Dof * > & aSources, Cell< real > & aWeights );

//------------------------------------------------------------------------------

            void
            set_source( Dof * aSource, const real aWeight);

//------------------------------------------------------------------------------

            void
            reset_sources();

//------------------------------------------------------------------------------

            uint
            number_of_sources() const;

//------------------------------------------------------------------------------

            Dof *
            source( const uint aIndex );

//------------------------------------------------------------------------------

            real
            weight( const uint aIndex ) const;

//------------------------------------------------------------------------------

            size_t
            memory() const;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline bool
        Dof::mesh_basis_is_flagged()
        {
            return mMeshBasis->is_flagged() ;
        }

//------------------------------------------------------------------------------

        inline mesh::Basis *
        Dof::mesh_basis()
        {
            return mMeshBasis ;
        }

//------------------------------------------------------------------------------
        inline mesh::Node *
        Dof::node()
        {
            BELFEM_ASSERT( this->is_node(),
                          "Tried to access DOF %lu as node, but it is not a node ",
                          ( long unsigned int ) mMyIndex );

            return reinterpret_cast< mesh::Node * >( mMeshBasis ) ;
        }

//------------------------------------------------------------------------------

        inline mesh::Edge *
        Dof::edge()
        {
            BELFEM_ASSERT( this->is_edge(),
                          "Tried to access DOF %lu as edge, but it is not an edge",
                          ( long unsigned int ) mMyIndex );

            return reinterpret_cast< mesh::Edge * >( mMeshBasis ) ;
        }

//------------------------------------------------------------------------------

        inline mesh::Face *
        Dof::face()
        {
            BELFEM_ASSERT( this->is_face(),
                          "Tried to access DOF %lu as face, but it is not a Face",
                          ( long unsigned int ) mMyIndex );

            return reinterpret_cast< mesh::Face * >( mMeshBasis ) ;
        }

//------------------------------------------------------------------------------

        inline mesh::Element *
        Dof::element()
        {
            BELFEM_ASSERT( this->is_cell(),
                           "Tried to access DOF %lu as element, but it is not an element ",
                           ( long unsigned int ) mMyIndex );

            return reinterpret_cast< mesh::Element * >( mMeshBasis );
        }

//------------------------------------------------------------------------------

        inline void
        Dof::fix( const real aDirichletValue )
        {
            mFixedFlag = true;
            mDirichletValue = aDirichletValue;
        }

//------------------------------------------------------------------------------

        inline void
        Dof::free()
        {
            mFixedFlag = false;
            mDirichletValue = BELFEM_QUIET_NAN;
        }

//------------------------------------------------------------------------------

        inline bool
        Dof::is_fixed() const
        {
            return mFixedFlag;
        }

//------------------------------------------------------------------------------

        inline bool
        Dof::is_hanging() const
        {
            return mNumberOfSources > 0 ;
        }

//------------------------------------------------------------------------------

        inline bool
        Dof::basis_is_hanging() const
        {
            return mMeshBasis->is_hanging() ;
        }

//------------------------------------------------------------------------------

        inline bool
        Dof::is_node() const
        {
            return mMeshBasis->entity_type() == EntityType::NODE ;
        }

//------------------------------------------------------------------------------

        inline EntityType
        Dof::entity_type() const
        {
            return mMeshBasis->entity_type() ;
        }


//------------------------------------------------------------------------------

        inline bool
        Dof::is_edge() const
        {
            return mMeshBasis->entity_type() == EntityType::EDGE ;
        }

//------------------------------------------------------------------------------

        inline bool
        Dof::is_face() const
        {
            return mMeshBasis->entity_type() == EntityType::FACE ;
        }

//------------------------------------------------------------------------------

        inline bool
        Dof::is_cell() const
        {
            return mMeshBasis->entity_type() == EntityType::CELL ;
        }

//------------------------------------------------------------------------------

        inline bool
        Dof::is_lambda() const
        {
            return mMeshBasis->entity_type() == EntityType::FACET ;
        }

//------------------------------------------------------------------------------

        inline real &
        Dof::value()
        {
            return mDirichletValue;
        }

//------------------------------------------------------------------------------

        inline Dof *
        Dof::dof( const uint aIndex )
        {
            BELFEM_ASSERT( aIndex < mVertexCounter,
                          "Index %u for DOF %lu out of bounds. ( must be less than %u )",
                          ( unsigned int ) aIndex,
                          ( long unsigned int ) this->id(),
                          ( unsigned int ) mVertexCounter );

            return reinterpret_cast< Dof * >( mVertices[ aIndex ] );
        }

//------------------------------------------------------------------------------

        inline uint
        Dof::type_id() const
        {
            return mTypeID;
        }

//------------------------------------------------------------------------------

        inline void
        Dof::set_my_index( const index_t aIndex )
        {
            mMyIndex = aIndex ;
        }

//------------------------------------------------------------------------------

        inline uint
        Dof::number_of_dofs() const
        {
            return this->number_of_vertices();
        }

//------------------------------------------------------------------------------

        inline uint
        Dof::index_on_entity() const
        {
            return mIndexOnEntity ;
        }

//------------------------------------------------------------------------------

        inline void
        Dof::set_field_index( const index_t aIndex )
        {
            mFieldIndex = aIndex ;
        }

//------------------------------------------------------------------------------

        inline index_t
        Dof::dof_index_on_field() const
        {
            return mIndexOnField ;
        }

//------------------------------------------------------------------------------

        inline index_t
        Dof::field_index() const
        {
            return mFieldIndex ;
        }

//------------------------------------------------------------------------------

        inline index_t
        Dof::my_index() const
        {
            return mMyIndex;
        }

//------------------------------------------------------------------------------

        inline uint
        Dof::number_of_sources() const
        {
            return mNumberOfSources ;
        }

//------------------------------------------------------------------------------

        inline Dof *
        Dof::source( const uint aIndex )
        {
            BELFEM_ASSERT( aIndex < mNumberOfSources,
                           "Source index %u out of range for dof %lu (must be < %u).",
                           ( unsigned int ) aIndex,
                           ( long unsigned int ) this->id(),
                           ( unsigned int ) mNumberOfSources );

            return mSources[ aIndex ];
        }

//------------------------------------------------------------------------------

        inline real
        Dof::weight( const uint aIndex ) const
        {
            BELFEM_ASSERT( aIndex < mNumberOfSources,
                           "Source index %u out of range for dof %lu (must be < %u).",
                           ( unsigned int ) aIndex,
                           ( long unsigned int ) this->id(),
                           ( unsigned int ) mNumberOfSources );
            return mCoefficients[ aIndex ];
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_FEM_DOF_HPP
