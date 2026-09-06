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

#ifndef CL_FEM_THINSHELL_HPP
#define CL_FEM_THINSHELL_HPP

#include "cl_Block.hpp"
#include "cl_Facet.hpp"
#include "cl_SideSet.hpp"


namespace belfem
{
    namespace mesh
    {
        class ThinShellFactory ;

        /**
         * the thin shell class is a data container needed to link
         * element masters and slaves to thin shell blocks
         */
        class ThinShell
        {
            // facets that are connected to the elements
            // in the same order as the elements in the blocks
            // this is needed later to link the air elements
            // with the shell elements in the FEM kernel
            // sideset is destroyed by mesh
            SideSet *         mSideSet  ;

            SideSet *         mGhostSideSet  ;

            Cell< Facet * > & mFacets ;

            // blocks that contain the thin shell elements
            // each block represents one layer
            // blocks are destroyed by mesh
            Cell< Block * > mBlocks ;

            Vector< real >  mThicknesses ;

            Cell< string >  mMaterials;

            // node indices, needed for hanging dofs
            Cell< index_t > mNodeIndices ;

            // empty cell, needed as dummy output
            Cell< Facet * > mNull ;

            // Blocks for side connectors
            Cell< Block * > mSideConnectorBlocks ;

            // additional nodes
            Cell< Node * > mSideConnectorNodes ;

            // additional edges
            Cell< Edge * > mSideConnectorEdges ;

            Cell< SideSet * > mSideConnectorSideSets ;


            // map: side-connector element id -> adjacent shell facet id.
            // Currently unpopulated: it has no producer and no consumer in
            // the tree; kept for the post-processor's perpendicular-H
            // recovery path.
            Map< id_t, std::pair< id_t, id_t > > mConnectorFacetMap ;

            friend class ThinShellFactory ;

//------------------------------------------------------------------------------
        public :
//------------------------------------------------------------------------------

            ThinShell( SideSet * aSideSet, SideSet * aGhostSideSet ) ;

            ~ThinShell() = default ;

            void
            set_thicknesses( const Vector < real > & aThicknesses );

            void
            move_node_indices( Cell< index_t > & aNodeIndices ) ;

            void
            set_label( const string & aLabel ) ;

            Cell< Facet * > &
            facets() ;

            Cell< Block * > &
            blocks() ;

            void
            set_materials( const Cell< string > & aMaterials ) ;

            const Cell< string > &
            materials() const ;

            ElementType
            element_type() const ;

            const Vector< real > &
            thicknesses() const ;

            const string &
            label() const ;

            id_t
            id() const ;

            id_t
            ghost_id() const ;

            size_t
            memory() const ;

            Cell< index_t > &
            node_indices() ;

            Cell< Facet * > &
            ghost_facets();

            Cell< Block * > &
            side_connector_blocks();

            Cell< SideSet * > &
            side_connector_sidesets();

            Map< id_t, std::pair< id_t, id_t > > &
            connector_facet_map();

            const Map< id_t, std::pair< id_t, id_t > > &
            connector_facet_map() const;

            Map< id_t, real > &
            connector_thicknesses();

            const Map< id_t, real > &
            connector_thicknesses() const;

        private:


            Cell< Node * > &
            side_connector_nodes();

            Cell< Edge * > &
            side_connector_edges();

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        inline id_t
        ThinShell::ghost_id() const
        {
            if ( mGhostSideSet != nullptr )
            {
                return mGhostSideSet->id();
            }
            else
            {
                return gNoID ;
            }
        }

        inline Cell< Facet * > &
        ThinShell::facets()
        {
            return mFacets ;
        }

//------------------------------------------------------------------------------

        inline Cell< Facet * > &
        ThinShell::ghost_facets()
        {
            if ( mGhostSideSet != nullptr )
            {
                return mGhostSideSet->facets();
            }
            else
            {
                return mNull ;
            }
        }

//------------------------------------------------------------------------------

        inline Cell< Block * > &
        ThinShell::blocks()
        {
            return mBlocks ;
        }

//------------------------------------------------------------------------------

        inline void
        ThinShell::set_materials( const Cell< string > & aMaterials )
        {
            BELFEM_ERROR( aMaterials.size() == mBlocks.size(),
                "size of materials must match number of blocks (%u vs %u)",
                ( unsigned int ) aMaterials.size(), ( unsigned int ) mBlocks.size()
                ) ;

            mMaterials = aMaterials ;
        }

//------------------------------------------------------------------------------

        inline const Cell< string > &
        ThinShell::materials() const
        {
            return mMaterials ;
        }

        inline ElementType
        ThinShell::element_type() const
        {
            return mSideSet->element_type() ;
        }

        inline const Vector< real > &
        ThinShell::thicknesses() const
        {
            return mThicknesses ;
        }

        inline const string &
        ThinShell::label() const
        {
            return mSideSet->label() ;
        }

        inline id_t
        ThinShell::id() const
        {
            return mSideSet->id() ;
        }

        inline size_t
        ThinShell::memory() const
        {
            size_t aMem = sizeof( ThinShell )
                + mBlocks.size() * sizeof( Block * )
                + mThicknesses.length() * sizeof( real ) ;

            // account for material strings
            for( const string & tMaterial : mMaterials )
            {
                aMem += sizeof( string ) + tMaterial.capacity() * sizeof( char );
            }

            return aMem ;
        }

        inline Map< id_t, std::pair< id_t, id_t > > &
        ThinShell::connector_facet_map()
        {
            return mConnectorFacetMap ;
        }

        inline const Map< id_t, std::pair< id_t, id_t > > &
        ThinShell::connector_facet_map() const
        {
            return mConnectorFacetMap ;
        }

        inline Cell< index_t > &
        ThinShell::node_indices()
        {
            return mNodeIndices ;
        }

        inline Cell< Block * > &
        ThinShell::side_connector_blocks()
        {
            return mSideConnectorBlocks ;
        }

        inline Cell< SideSet * > &
        ThinShell::side_connector_sidesets()
        {
            return mSideConnectorSideSets ;
        }

        inline Cell< Node * > &
        ThinShell::side_connector_nodes()
        {
            return mSideConnectorNodes ;
        }

        inline Cell< Edge * > &
        ThinShell::side_connector_edges()
        {
            return mSideConnectorEdges ;
        }

//------------------------------------------------------------------------------

    }
}
#endif //CL_FEM_THINSHELL_HPP
