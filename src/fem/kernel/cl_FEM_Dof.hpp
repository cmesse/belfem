//
// Created by Christian Messe on 24.10.19.
//

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
            // pointer to node or edge of this flag
            mesh::Vertex * mMeshVertex;

            // type of the current dof, needed for vectors
            // temperature: always zero
            // displacement: 1: ux, 2: uy, 3:uz
            const uint mTypeID;

            // index of field on mesh, set by constructor
            index_t mFieldIndex = gNoIndex ;

            // index of dof on field
            const index_t mIndexOnField  ;

            real mDirichletValue = 0.0 ; // BELFEM_QUIET_NAN;

            // tells if a Dirichlet boundary condition is imposed
            bool mFixedFlag = false;

            // index for DOF on current proc
            index_t mMyIndex = gNoIndex;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Dof( const id_t aID, const uint & aType, mesh::Node * aNode );

//------------------------------------------------------------------------------

            Dof( const id_t aID, const uint & aType, mesh::Edge * aEdge, const index_t aDofIndexOnField );

//------------------------------------------------------------------------------

            Dof( const id_t aID, const uint & aType, mesh::Face * aFace, const index_t aDofIndexOnField );

//------------------------------------------------------------------------------

            Dof( const id_t aID, const uint & aType, mesh::Facet * aFacet );

//------------------------------------------------------------------------------

            ~Dof() = default;

//------------------------------------------------------------------------------

            bool
            mesh_vertex_is_flagged();

//------------------------------------------------------------------------------

            mesh::Vertex *
            mesh_vertex();

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
             * impose a Diriclet boundary condition on this node
             */
            void
            fix( const real & aDirichletValue );

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
            dof( const uint & aIndex );

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
             * set proc local index
             */
             void
             set_my_index( const index_t & aIndex );
//------------------------------------------------------------------------------

            void
            set_field_index( const uint & aIndex );

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
        };

//------------------------------------------------------------------------------

        inline bool
        Dof::mesh_vertex_is_flagged()
        {
            return mMeshVertex->is_flagged() ;
        }

//------------------------------------------------------------------------------

        inline mesh::Vertex *
        Dof::mesh_vertex()
        {
            return mMeshVertex ;
        }

//------------------------------------------------------------------------------
        inline mesh::Node *
        Dof::node()
        {
            BELFEM_ASSERT( this->is_node(),
                          "Tried to access DOF %lu as node, but it is linked to edge %lu ",
                          ( long unsigned int ) mMyIndex,
                          ( long unsigned int ) mMeshVertex->id() );

            return reinterpret_cast< mesh::Node * >( mMeshVertex ) ;
        }

//------------------------------------------------------------------------------

        inline mesh::Edge *
        Dof::edge()
        {
            BELFEM_ASSERT( this->is_edge(),
                          "Tried to access DOF %lu as edge, bot it is linked to node %lu ",
                          ( long unsigned int ) mMyIndex,
                          ( long unsigned int ) mMeshVertex->id() );

            return reinterpret_cast< mesh::Edge * >( mMeshVertex ) ;
        }

//------------------------------------------------------------------------------

        inline void
        Dof::fix( const real & aDirichletValue )
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
        Dof::is_node() const
        {
            return mMeshVertex->entity_type() == EntityType::NODE ;
        }

//------------------------------------------------------------------------------

        inline EntityType
        Dof::entity_type() const
        {
            return mMeshVertex->entity_type() ;
        }


//------------------------------------------------------------------------------

        inline bool
        Dof::is_edge() const
        {
            return mMeshVertex->entity_type() == EntityType::EDGE ;
        }

//------------------------------------------------------------------------------

        inline bool
        Dof::is_face() const
        {
            return mMeshVertex->entity_type() == EntityType::FACE ;
        }

//------------------------------------------------------------------------------

        inline bool
        Dof::is_lambda() const
        {
            return mMeshVertex->entity_type() == EntityType::FACET ;
        }


//------------------------------------------------------------------------------

        inline real &
        Dof::value()
        {
            return mDirichletValue;
        }

//------------------------------------------------------------------------------

        inline Dof *
        Dof::dof( const uint & aIndex )
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
        Dof::set_my_index( const index_t & aIndex )
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

        inline void
        Dof::set_field_index( const index_t & aIndex )
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
    }
}

#endif //BELFEM_CL_FEM_DOF_HPP
