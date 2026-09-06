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

#include "cl_Mesh_Basis.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace mesh
    {
//----------------------------------------------------------------------------

        Basis::Basis() :
                graph::Vertex()
        {
            // mesh data always assume that owner is zero at first
            this->set_owner( 0 );
        }

//----------------------------------------------------------------------------

        Basis::~Basis()
        {
            // delete T-Matrix
            if( mNumberOfSources != 0 )
            {
                free( mSources );
                free( mWeights );
            }

            // delete dofs
            if( mNumberOfDofs != 0 )
            {
                free( mDofs );
            }
        }

//------------------------------------------------------------------------------

        uint
        Basis::number_of_elements() const
        {
            BELFEM_ERROR( false, "Invalid call to abstract basis class : number_of_elements()");
            return 0 ;
        }

//------------------------------------------------------------------------------

        Element *
        Basis::element( const uint aIndex )
        {
            BELFEM_ERROR( false, "Invalid call to abstract basis class : element()");
            return nullptr ;
        }

//------------------------------------------------------------------------------

        EntityType
        Basis::entity_type() const
        {
            BELFEM_ERROR( false, "Invalid call to abstract basis class : entity_type()");
            return EntityType::UNDEFINED ;
        }

//------------------------------------------------------------------------------

        uint
        Basis::number_of_nodes() const
        {
            BELFEM_ERROR( false, "Invalid call to abstract basis class : number_of_nodes()");
            return 0 ;
        }

//------------------------------------------------------------------------------

        uint
        Basis::number_of_edges() const
        {
            BELFEM_ERROR( false, "Invalid call to abstract basis class : number_of_edges()");
            return 0 ;
        }

//------------------------------------------------------------------------------

        uint
        Basis::number_of_faces() const
        {
            BELFEM_ERROR( false, "Invalid call to abstract basis class : number_of_faces()" );
            return 0 ;
        }

//------------------------------------------------------------------------------

        Node *
        Basis::node( const uint aIndex )
        {
            BELFEM_ERROR( false,
                          "Invalid call to abstract basis class : node()" );
            return nullptr ;
        }

//------------------------------------------------------------------------------

        const Node *
        Basis::node( const uint aIndex ) const
        {
            BELFEM_ERROR( false,
                          "Invalid call to abstract basis class : node() const" );
            return nullptr ;
        }


//------------------------------------------------------------------------------

        Edge *
        Basis::edge( const uint aIndex )
        {
            BELFEM_ERROR( false,
                          "Invalid call to abstract basis class : edge()" );
            return nullptr ;
        }


//------------------------------------------------------------------------------

        const Edge *
        Basis::edge( const uint aIndex ) const
        {
            BELFEM_ERROR( false,
                          "Invalid call to abstract basis class : edge() const" );
            return nullptr ;
        }


//------------------------------------------------------------------------------

        Face *
        Basis::face( const uint aIndex )
        {
            BELFEM_ERROR( false,
                          "Invalid call to abstract basis class : face()" );
            return nullptr ;
        }

//------------------------------------------------------------------------------

        const Face *
        Basis::face( const uint aIndex ) const
        {
            BELFEM_ERROR( false,
                          "Invalid call to abstract basis class : face() const" );
            return nullptr ;
        }

//------------------------------------------------------------------------------

        Node *
        Basis::source_node( const uint aIndex )
        {
            BELFEM_ASSERT( mSources[ aIndex ]->entity_type() == EntityType::NODE,
                           "Source basis %u of basis %lu is not a node",
                           ( unsigned int ) aIndex,
                           ( long unsigned int ) this->id() );

            return reinterpret_cast< Node * >( mSources[ aIndex ] );
        }

//------------------------------------------------------------------------------

        const Node *
        Basis::source_node( const uint aIndex ) const
        {
            BELFEM_ASSERT( mSources[ aIndex ]->entity_type() == EntityType::NODE,
                           "Source basis %u of basis %lu is not a node",
                           ( unsigned int ) aIndex,
                           ( long unsigned int ) this->id() );

            return reinterpret_cast< Node * >( mSources[ aIndex ] );
        }

//------------------------------------------------------------------------------

        void
        Basis::set_sources( Cell< Basis * >      & aSources,
                              const Vector< real > & aCoefficients )
        {
            // sanity check
            BELFEM_ASSERT( mNumberOfSources == 0,
                           "Coefficients for basis %lu already assigned",
                           ( long unsigned int ) this->id() );

            // sanity check
            BELFEM_ASSERT( aSources.size() == aCoefficients.length(),
                           "Number of source basis and coefficients for basis %lu does not match (%u vs %u)",
                           ( long unsigned int ) this->id(),
                           ( unsigned int ) aSources.size(),
                           ( unsigned int ) aCoefficients.length() );

            // the counter is deliberately narrow; a count beyond its width must
            // abort loudly instead of wrapping ( a silent wrap dropped the
            // original node from 465-entry cut-trunk source lists )
            BELFEM_ERROR( aSources.size() <= std::numeric_limits< decltype( mNumberOfSources ) >::max(),
                          "Basis %lu gets %lu sources, but the source counter holds at most %lu",
                          ( long unsigned int ) this->id(),
                          ( long unsigned int ) aSources.size(),
                          ( long unsigned int ) std::numeric_limits< decltype( mNumberOfSources ) >::max() );

            // set the number of sources
            mNumberOfSources = aSources.size() ;

            // allocate memory
            mSources = ( Basis ** ) malloc( mNumberOfSources * sizeof( Basis * ) );

            // copy sources
            for( uint k=0; k<mNumberOfSources; ++k )
            {
                mSources[ k ] = aSources( k );
            }

            // allocate memory
            mWeights = ( real * ) malloc( mNumberOfSources * sizeof( real ) );

            // copy coefficients
            for( uint k=0; k<mNumberOfSources; ++k )
            {
                mWeights[ k ] = aCoefficients( k );
            }
        }

//------------------------------------------------------------------------------

        void
        Basis::set_sources( Cell< Basis * >      & aSources,
                              const Cell< real > & aCoefficients )
        {
            // sanity check
            BELFEM_ASSERT( mNumberOfSources == 0,
                           "Coefficients for basis %lu already assigned",
                           ( long unsigned int ) this->id() );

            // sanity check
            BELFEM_ASSERT( aSources.size() == aCoefficients.size(),
                           "Number of source basis and coefficients for basis %lu does not match (%u vs %u)",
                           ( long unsigned int ) this->id(),
                           ( unsigned int ) aSources.size(),
                           ( unsigned int ) aCoefficients.size() );

            // the counter is deliberately narrow; a count beyond its width must
            // abort loudly instead of wrapping ( a silent wrap dropped the
            // original node from 465-entry cut-trunk source lists )
            BELFEM_ERROR( aSources.size() <= std::numeric_limits< decltype( mNumberOfSources ) >::max(),
                          "Basis %lu gets %lu sources, but the source counter holds at most %lu",
                          ( long unsigned int ) this->id(),
                          ( long unsigned int ) aSources.size(),
                          ( long unsigned int ) std::numeric_limits< decltype( mNumberOfSources ) >::max() );

            // set the number of sources
            mNumberOfSources = aSources.size() ;

            // allocate memory
            mSources = ( Basis ** ) malloc( mNumberOfSources * sizeof( Basis * ) );

            // copy sources
            for( uint k=0; k<mNumberOfSources; ++k )
            {
                mSources[ k ] = aSources( k );
            }

            // allocate memory
            mWeights = ( real * ) malloc( mNumberOfSources * sizeof( real ) );

            // copy coefficients
            for( uint k=0; k<mNumberOfSources; ++k )
            {
                mWeights[ k ] = aCoefficients( k );
            }
        }

 //------------------------------------------------------------------------------

        void
        Basis::reset_source_container()
        {
            if( mSources != nullptr )
            {
                free( mSources );
                mSources = nullptr ;
            }
            if( mWeights != nullptr )
            {
                free ( mWeights );
                mWeights = nullptr ;
            }
            mNumberOfSources = 0 ;
        }

//------------------------------------------------------------------------------

        void
        Basis::allocate_source_container( uint aNumSources )
        {
            // sanity check
            BELFEM_ASSERT( mNumberOfSources == 0,
                           "Coefficients for basis %lu already assigned",
                           ( long unsigned int ) this->id() );

            // Capacity is deliberately not stored on every Basis: callers
            // must append exactly aNumSources entries through add_source().
            // The only local bound is therefore the source-counter width.
            BELFEM_ERROR( aNumSources <= std::numeric_limits< decltype( mNumberOfSources ) >::max(),
                          "Basis %lu gets %lu sources, but the source counter holds at most %lu",
                          ( long unsigned int ) this->id(),
                          ( long unsigned int ) aNumSources,
                          ( long unsigned int ) std::numeric_limits< decltype( mNumberOfSources ) >::max() );

            mSources = ( Basis ** ) malloc( aNumSources * sizeof( Basis * ) );
            mWeights = ( real * ) malloc( aNumSources * sizeof( real ) );
        }

//------------------------------------------------------------------------------

        void
        Basis::add_source( Basis * aSource, const real aWeight )
        {
            BELFEM_ERROR( mNumberOfSources < std::numeric_limits< decltype( mNumberOfSources ) >::max(),
                          "Source counter of basis %lu is full ( max %lu )",
                          ( long unsigned int ) this->id(),
                          ( long unsigned int ) std::numeric_limits< decltype( mNumberOfSources ) >::max() );

            mSources[ mNumberOfSources ] = aSource ;
            mWeights[ mNumberOfSources ] = aWeight ;
            ++mNumberOfSources ;
        }

//------------------------------------------------------------------------------

        void
        Basis::set_sources(   Cell< Node * > & aSources,
                              const Vector< real > & aCoefficients )
        {
            // sanity check
            BELFEM_ASSERT( mNumberOfSources == 0,
                           "Coefficients for basis %lu already assigned",
                           ( long unsigned int ) this->id() );


            // sanity check
            BELFEM_ASSERT( aSources.size() == aCoefficients.length(),
                           "Number of source nodes and coefficients for basis %lu does not match (%u vs %u)",
                           ( long unsigned int ) this->id(),
                           ( unsigned int ) aSources.size(),
                           ( unsigned int ) aCoefficients.length() );


            // the counter is deliberately narrow; a count beyond its width must
            // abort loudly instead of wrapping ( a silent wrap dropped the
            // original node from 465-entry cut-trunk source lists )
            BELFEM_ERROR( aSources.size() <= std::numeric_limits< decltype( mNumberOfSources ) >::max(),
                          "Basis %lu gets %lu sources, but the source counter holds at most %lu",
                          ( long unsigned int ) this->id(),
                          ( long unsigned int ) aSources.size(),
                          ( long unsigned int ) std::numeric_limits< decltype( mNumberOfSources ) >::max() );

            // set the number of sources
            mNumberOfSources = aSources.size() ;

            // allocate memory
            mSources = ( Basis ** ) malloc( mNumberOfSources * sizeof( Basis * ) );

            // copy sources
            for( uint k=0; k<mNumberOfSources; ++k )
            {
                mSources[ k ] = reinterpret_cast< Basis * > ( aSources( k ) );
            }

            // allocate memory
            mWeights = ( real * ) malloc( mNumberOfSources * sizeof( real ) );

            // copy coefficients
            for( uint k=0; k<mNumberOfSources; ++k )
            {
                mWeights[ k ] = aCoefficients( k );
            }

        }

//------------------------------------------------------------------------------

        void
        Basis::flag_sources()
        {
            for( uint k=0; k<mNumberOfSources; ++k )
            {
                mSources[ k ]->flag() ;
            }
        }

//------------------------------------------------------------------------------

    }
}
