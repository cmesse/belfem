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

#ifndef BELFEM_CL_MESH_BASIS_HPP
#define BELFEM_CL_MESH_BASIS_HPP

#include <limits>

#include "typedefs.hpp"
#include "cl_Graph_Vertex.hpp"
#include "Mesh_Enums.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace mesh
    {
        class Node ;
        class Edge ;
        class Face ;
        class Element ;

        class Basis : public graph::Vertex
        {
            // for hanging basis
            //! cohomology cut trunks stack one source per crossing cut plus
            //! the original, so this counter must hold hundreds; the former
            //! uint8_t wrapped silently on a 464-cut deck
            uint16_t   mNumberOfSources   = 0 ;
            Basis ** mSources           = nullptr ;
            real  *  mWeights      = nullptr ;

            uint8_t    mNumberOfDofs = 0 ;
            graph::Vertex ** mDofs ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Basis();

            ~Basis() override ;

//------------------------------------------------------------------------------

            /**
             * returns the type of this vertex
             */
            virtual EntityType
            entity_type() const ;

//------------------------------------------------------------------------------

            virtual uint
            number_of_nodes() const ;

//------------------------------------------------------------------------------

            virtual uint
            number_of_edges() const ;

//------------------------------------------------------------------------------

            virtual uint
            number_of_faces() const ;

//------------------------------------------------------------------------------

            virtual Node *
            node( const uint aIndex );

//------------------------------------------------------------------------------

            virtual const Node *
            node( const uint aIndex ) const ;

//------------------------------------------------------------------------------

            virtual Edge *
            edge( const uint aIndex );

//------------------------------------------------------------------------------

            virtual const Edge *
            edge( const uint aIndex ) const ;

//------------------------------------------------------------------------------

            virtual Face *
            face( const uint aIndex=0 );

//------------------------------------------------------------------------------

            virtual const Face *
            face( const uint aIndex=0 ) const ;

//------------------------------------------------------------------------------

            virtual uint
            number_of_elements() const ;

//------------------------------------------------------------------------------

            virtual Element *
            element( const uint aIndex );

//------------------------------------------------------------------------------

            bool
            is_hanging() const ;

//------------------------------------------------------------------------------

            uint
            number_of_sources() const ;

//------------------------------------------------------------------------------

            Basis *
            source( const uint aIndex );

//------------------------------------------------------------------------------

            const Basis *
            source( const uint aIndex ) const ;

//------------------------------------------------------------------------------

            Node *
            source_node( const uint aIndex );

//------------------------------------------------------------------------------

            const Node *
            source_node( const uint aIndex ) const ;

//------------------------------------------------------------------------------

            real
            weight( const uint aIndex ) const ;

//------------------------------------------------------------------------------

            void
            set_weight( const uint aIndex, const real aValue );

//------------------------------------------------------------------------------

            void
            set_sources( Cell< Basis * > & aSources, const Vector< real > & aCoefficients );

//------------------------------------------------------------------------------

            void
            set_sources( Cell< Basis * > & aSources, const Cell< real > & aCoefficients );

//------------------------------------------------------------------------------

            void
            set_sources( Cell< Node * > & aSources, const Vector< real > & aCoefficients );

//------------------------------------------------------------------------------

            void
            reset_source_container();

//------------------------------------------------------------------------------

            //! Allocates exactly the number of source slots that subsequent
            //! add_source() calls must append. Capacity is deliberately not
            //! stored on every Basis; callers own this invariant.
            void
            allocate_source_container( uint aNumSources );

//------------------------------------------------------------------------------

            //! Appends one source to the allocation sized by
            //! allocate_source_container().
            void
            add_source( Basis * aSource, const real aWeight=1.0 );

//------------------------------------------------------------------------------

            void
            flag_sources();

//------------------------------------------------------------------------------

            void
            increment_dof_counter();

//------------------------------------------------------------------------------

            void
            allocate_dof_container();

//------------------------------------------------------------------------------

            void
            reset_dof_container();

//------------------------------------------------------------------------------

            void
            insert_dof( graph::Vertex * aDof );

//------------------------------------------------------------------------------

            uint
            number_of_dofs() const ;

//------------------------------------------------------------------------------

            graph::Vertex *
            dof( const uint aIndex );

//------------------------------------------------------------------------------
        };

//----------------------------------------------------------------------------

        inline bool
        Basis::is_hanging() const
        {
            return mNumberOfSources > 0 ;
        }

//----------------------------------------------------------------------------

        inline uint
        Basis::number_of_sources() const
        {
            return mNumberOfSources ;
        }


//------------------------------------------------------------------------------

        inline void
        Basis::increment_dof_counter()
        {
            // a wrapped counter would make allocate_dof_container skip its
            // malloc and insert_dof write through an uninitialized pointer
            BELFEM_ERROR( mNumberOfDofs < std::numeric_limits< decltype( mNumberOfDofs ) >::max(),
                          "Dof counter of basis %lu is full ( max %lu )",
                          ( long unsigned int ) this->id(),
                          ( long unsigned int ) std::numeric_limits< decltype( mNumberOfDofs ) >::max() );

            ++mNumberOfDofs ;
        }

//------------------------------------------------------------------------------

        inline void
        Basis::allocate_dof_container()
        {
            if( mNumberOfDofs > 0 )
            {
                mDofs = ( graph::Vertex ** ) malloc( mNumberOfDofs * sizeof ( graph::Vertex * ) );
                std::fill( mDofs, mDofs + mNumberOfDofs, nullptr );
                mNumberOfDofs = 0 ;
            }
        }

//------------------------------------------------------------------------------

        inline void
        Basis::reset_dof_container()
        {
            if( mNumberOfDofs > 0 )
            {
                free( mDofs );
                mDofs = nullptr ;
                mNumberOfDofs = 0 ;
            }
        }

//------------------------------------------------------------------------------

        inline void
        Basis::insert_dof( graph::Vertex * aDof )
        {
            BELFEM_ERROR( mNumberOfDofs < std::numeric_limits< decltype( mNumberOfDofs ) >::max(),
                          "Dof counter of basis %lu is full ( max %lu )",
                          ( long unsigned int ) this->id(),
                          ( long unsigned int ) std::numeric_limits< decltype( mNumberOfDofs ) >::max() );

            mDofs[ mNumberOfDofs++ ] = aDof ;
        }

//----------------------------------------------------------------------------

        inline uint
        Basis::number_of_dofs() const
        {
            return mNumberOfDofs ;
        }

//----------------------------------------------------------------------------

        inline graph::Vertex *
        Basis::dof( const uint aIndex )
        {
            BELFEM_ASSERT( aIndex < ( uint ) mNumberOfDofs, "Index %u out of bounds for basis %lu (expect < %u).",
                           ( unsigned int ) aIndex,
                           ( long unsigned int ) this->id(),
                           ( unsigned int ) mNumberOfDofs );

            return mDofs[ aIndex ];
        }

//----------------------------------------------------------------------------

        inline Basis *
        Basis::source( const uint aIndex )
        {
            BELFEM_ASSERT( aIndex < ( uint ) mNumberOfSources, "Index %u out of bounds for basis %lu (must be < %u).",
                    ( unsigned int ) aIndex,
                    ( long unsigned int ) this->id(),
                    ( unsigned int ) mNumberOfSources );

            return mSources[ aIndex ];
        }

//----------------------------------------------------------------------------

        inline const Basis *
        Basis::source( const uint aIndex ) const
        {
            BELFEM_ASSERT( aIndex < ( uint ) mNumberOfSources, "Index %u out of bounds for basis %lu (must be < %u).",
                           ( unsigned int ) aIndex,
                           ( long unsigned int ) this->id(),
                           ( unsigned int ) mNumberOfSources );

            return mSources[ aIndex ];
        }

//----------------------------------------------------------------------------

        inline real
        Basis::weight( const uint aIndex ) const
        {
            BELFEM_ASSERT( aIndex < mNumberOfSources, "Index %u out of bounds for basis %lu (must be < %u).",
                    ( unsigned int ) aIndex,
                           ( long unsigned int ) this->id(),
                           ( unsigned int ) mNumberOfSources );

            return mWeights[ aIndex ];
        }

//----------------------------------------------------------------------------

        inline void
        Basis::set_weight( const uint aIndex, const real aValue )
        {
            BELFEM_ASSERT( aIndex < ( uint ) mNumberOfSources, "Index %u out of bounds for basis %lu (must be < %u).",
                    ( unsigned int ) aIndex,
                           ( long unsigned int ) this->id(),
                           ( unsigned int ) mNumberOfSources );

            mWeights[ aIndex ] = aValue ;
        }
//----------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_MESH_BASIS_HPP
