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

#include "cl_ThinShell.hpp"
namespace belfem
{
    namespace mesh
    {
        ThinShell::ThinShell( SideSet * aSideSet, SideSet * aGhostSideSet ) :
            mSideSet( aSideSet ),
            mGhostSideSet( aGhostSideSet ),
            mFacets( aSideSet->facets() )
        {
            // don't write this sideset to the mesh
            mSideSet->hide( true );

            // also hide ghost sideset
            if ( aGhostSideSet != nullptr )
            {
                aGhostSideSet->hide( true );
            }
        }

        void
        ThinShell::set_label( const string & aLabel )
        {
            mSideSet->label() = aLabel;
        }

        void
        ThinShell::set_thicknesses( const Vector< real > & aThicknesses )
        {
            BELFEM_ASSERT( aThicknesses.length() == mBlocks.size(),
                "need as many thicknesses as blocks ( is %u, expect %u )" ,
                ( unsigned int ) aThicknesses.length(), ( unsigned int ) mBlocks.size()
                );

            mThicknesses = aThicknesses;

            index_t tCount = 0 ;
            for ( Block * tBlock : mBlocks )
            {
                tBlock->set_thickness( mThicknesses( tCount++ ) );
            }
        }

        void
        ThinShell::move_node_indices( Cell< index_t > & aNodeIndices )
        {
            mNodeIndices = std::move( aNodeIndices );
            aNodeIndices.clear();
        }

    }
}