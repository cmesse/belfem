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

#include <cmath>

#include "cl_FEM_Bearing.hpp"
#include "cl_FEM_DofManagerBase.hpp"
#include "cl_FEM_DofManager.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_IWG.hpp"
#include "commtools.hpp"
#include "cl_Logger.hpp"

namespace belfem
{
    namespace fem
    {
        //! a periodic pair hangs with exactly weight 1 ( a stored literal,
        //! no arithmetic on it ); anything measurably off unit is a
        //! different condensation shape and must not be pinned
        constexpr real gBearingUnitWeightTol = BELFEM_EPSILON ;

//------------------------------------------------------------------------------

        Bearing::Bearing( DofManagerBase * aParent ) :
            mParent( aParent ),
            mID( 0 )
        {

        }

//------------------------------------------------------------------------------

        Bearing::Bearing( DofManager * aParent,
                          const id_t aID, mesh::Node * aNode ) :
                mParent( aParent ),
                mID( aID ),
                mNode( aNode )
        {

        }

//------------------------------------------------------------------------------

        Bearing::~Bearing()
        {
            if( mNumDofs != 0 )
            {
                ::free( mDOFs );
            }
        }

//------------------------------------------------------------------------------

        void
        Bearing::allocate_dof_container( const index_t aNumDofs )
        {
            BELFEM_ASSERT(
                mNumDofs == 0,
                "Dof container of bearing %lu has already been allocated",
                ( long unsigned int ) mID );

            mNumDofs = aNumDofs;
            mDOFs = ( Dof ** ) malloc( aNumDofs * sizeof( Dof * ) );
        }

//------------------------------------------------------------------------------

        void
        Bearing::insert_dof( Dof * aDof, const index_t aIndex )
        {
            BELFEM_ASSERT(
                aIndex < mNumDofs,
                "Invalid dof index for bearing %lu ( Node %lu ; is: %lu, expect < %lu )",
                ( long unsigned int ) mNode->id(),
                ( long unsigned int ) aIndex,
                ( long unsigned int ) mNumDofs,
                ( long unsigned int ) mID );

            mDOFs[ aIndex ] = aDof ;
        }

//------------------------------------------------------------------------------

        void
        Bearing::impose_dirichlet( const real aValue, const uint aDofType )
        {
            // Empty bearings are legitimate on WORKER ranks: linking is
            // master-only, and synchronize_dirichlet_bcs() redistributes the
            // authoritative fixed set before the graphs are built. On the
            // MASTER — which holds the full mesh and every dof — a deck-named
            // bearing that resolved to nothing is a usage error, and the old
            // silent return here is what let a dead bearing run for months.
            // mID == 0 is the empty-lookup sentinel: gmsh point tags
            // start at 1.
            if ( mNumDofs == 0 )
            {
                // mID == 0 is the empty-lookup sentinel: the requested id was
                // dropped by the map lookup, so name the deck list instead
                if ( mID == 0 )
                {
                    BELFEM_ERROR( comm_rank() != 0,
                        "bearing/gauge target invalid: an id in the deck's 'nodes' list does not exist as a point on the mesh" );
                }
                else
                {
                    BELFEM_ERROR( comm_rank() != 0,
                        "bearing/gauge target invalid: point %lu ( node %lu ) carries no dof of the requested field",
                        ( long unsigned int ) mID,
                        ( long unsigned int ) mNode->id() );
                }
                return ;
            }

            BELFEM_ASSERT(
                aDofType < mNumDofs,
                "Invalid dof index for bearing %lu ( Node %lu ; is: %lu, expect < %lu )",
                ( long unsigned int ) mID,
                ( long unsigned int ) mNode->id(),
                ( long unsigned int ) aDofType,
                ( long unsigned int ) mNumDofs );

            Dof * tDof = mDOFs[ aDofType ];

            // A condensed ( hanging ) dof never enters the free/fixed split,
            // so pinning it directly is a silent no-op — the original defect: a
            // bearing on the periodic TARGET face. Alias and source share ONE
            // independent dof; pin the SOURCE, and compute_hanging_dofs()
            // writes the value back onto the alias unconditionally. Only the
            // single-source unit-weight shape ( the periodic pair ) is safely
            // pinnable this way — any other condensation shape would be
            // over-constrained, and a hanging source violates the flattening
            // invariant of the T-matrix construction. BEARING LINKING is
            // master-only ( worker dofs may carry sources, but worker
            // bearings never hold dofs ), so this branch runs on rank 0 and
            // synchronize_dirichlet_bcs() distributes the flag. LIMITATION:
            // a per-timestep gauge whose FIRST imposition lands on a
            // condensed node arrives after the free/fixed graphs are frozen —
            // the reroute then flips a graph-frozen free dof to fixed, which
            // the fixed-value gather is not sized for. A gauge on a condensed
            // node therefore still requires a factory-time bearing on the
            // same point to pre-fix the source.
            if ( tDof->is_hanging() )
            {
                BELFEM_ERROR( tDof->number_of_sources() == 1,
                    "Bearing %lu ( node %lu ): the dof is condensed onto %u sources and cannot be pinned - pick an unconstrained node",
                    ( long unsigned int ) mID,
                    ( long unsigned int ) mNode->id(),
                    ( unsigned int ) tDof->number_of_sources() );

                const real tWeight = tDof->weight( 0 );

                BELFEM_ERROR( std::isfinite( tWeight )
                              && std::abs( tWeight - 1.0 ) <= gBearingUnitWeightTol,
                    "Bearing %lu ( node %lu ): the dof hangs on its source with weight %g != 1 and cannot be pinned - pick an unconstrained node",
                    ( long unsigned int ) mID,
                    ( long unsigned int ) mNode->id(),
                    ( double ) tWeight );

                Dof * tSource = tDof->source( 0 );

                BELFEM_ERROR( tSource != nullptr && ! tSource->is_hanging(),
                    "Bearing %lu ( node %lu ): the source dof is null or itself hanging - condensation chains are flattened by construction, this is an internal error",
                    ( long unsigned int ) mID,
                    ( long unsigned int ) mNode->id() );

                // announce the rerouting once — the gauge path re-imposes
                // every timestep with a fresh value, and only the first call
                // finds the source still free
                const bool tFirstImposition = ! tSource->is_fixed() ;

                // the ONE independent dof of the pair. The alias reconstructs
                // as weight * source, so divide the requested value by the
                // ( unit-tolerance ) weight to land exactly on aValue
                tSource->fix( aValue / tWeight );

                if ( tFirstImposition && comm_rank() == 0 )
                {
                    message( InfoLevel::Default,
                        "    Bearing %lu: node %lu is condensed onto node %lu; pinning that source dof instead",
                        ( long unsigned int ) mID,
                        ( long unsigned int ) mNode->id(),
                        ( long unsigned int ) tSource->mesh_basis()->id() );
                }
            }
            else
            {
                tDof->fix( aValue );
            }
        }

//------------------------------------------------------------------------------

        void
        Bearing::free()
        {
            for ( index_t k=0; k<mNumDofs; ++k )
            {
                mDOFs[ k ]->free() ;
            }
        }

//------------------------------------------------------------------------------
    }
}