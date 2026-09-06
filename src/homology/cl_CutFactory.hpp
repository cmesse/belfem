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

#ifndef BELFEM_CL_CUTFACTORY_HPP
#define BELFEM_CL_CUTFACTORY_HPP

#include "commtools.hpp"
#include "typedefs.hpp"
#include "cl_Mesh.hpp"
#include "cl_SimplicialComplex.hpp"
#include "cl_Cohomology.hpp"
#include "cl_Homology.hpp"

#include "en_CutAlgorithm.hpp"
#include "../mesh/cl_Protoshell.hpp"
#include "cl_Topology.hpp"
#include "cl_Mesh_PeriodicityFactory.hpp"
namespace belfem
{
    namespace mesh
    {
//-----------------------------------------------------------------------------

        /**
         * @brief Main orchestrator for automatic cut generation on multiply-connected domains.
         *
         * @ingroup grp_homology
         * @see @ref homology_homology_usage_guide
         */
        class CutFactory
        {
            const proc_t mCommRank ;
            Mesh * mMesh ;
            id_t mMaxID = 0 ;

            Topology * mTopology ;
            PeriodicityFactory * mPeriodicFactory ;
            const Vector< id_t > & mThinShellSidesets ;

            Cell< Node * > & mAbstractNodes ;
            Cell< Node * > & mOrphanedNodes ;
            Cell< Protoshell * > & mProtoshells ;

            const CutAlgorithm mAlgorithm ;
            const bool mUseEnrichment ;

            Vector< id_t > mThinShellBoundaryIDs ;

            Cell< SideSet * > mThinShellBoundaries ;
            Cell< SideSet * > mCuts ;
            Cell< SideSet * > mTemporaryThinShellSidesets ;

            //Cell<Cell< id_t >> mCurrentTerminals ;
            Cell<Cell< id_t >> mTerminals ;
            Cell< id_t > mThinShellIndices; //Indices of thin shell terminals in mTerminals ;

            Vector< id_t>  mOuterBoundaries ;

            SimplicialComplex * mSimplicialComplex = nullptr ;
            Cohomology        * mCohomology = nullptr ;
            Homology          * mRelativeHomology = nullptr ;


            Cell< Node * >    mThinShellDuplicates ;
            Cell< Node * >    mThinShellMasterNodes ;
            Cell< Node * >    mThinShellSlaveNodes ;

            Map< id_t, Vector< index_t > * > mIndicesOfOriginalTerminalNodes ;

            //Map< id_t, DynamicBitset * > mNodesOnOriginalTerminalSidesets ;

            // this matrix lists the cuts based on the cohomologies
            //Vector< id_t > mTriBlocks ;
            //Vector< id_t > mTetBlocks ;
            //Vector< id_t > mPyraBlocks ;

            //Map< id_t, Element * > mElementMap ;
            //Map< id_t, Node * >    mNodeMap ;



            // future user settings
            bool mSuggestHomologies = true ; // this needs to be turned on for the incidence matrix to work!

            // poisson gives prettier cuts, but is slower
            bool mUsePoissonInsteadOfRcm = false ;


//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            CutFactory(
                Mesh * aMesh,
                Topology * aTopology,
                Cell< Protoshell * >  & aProtoshells,
                const CutAlgorithm aAlgorithm,
                const bool aUseEnrichment );

//-----------------------------------------------------------------------------

            ~CutFactory();

//-----------------------------------------------------------------------------

            void
            set_terminals( const Cell<Cell< id_t >> & aTerminals, const Cell<id_t> & aThinShellIndices );

//-----------------------------------------------------------------------------

            void
            set_periodicity( PeriodicityFactory * aPeriodicFactory );

//-----------------------------------------------------------------------------

            void
            run();

//-----------------------------------------------------------------------------

            // returns the abstract nodes that represent the dof for the cuts
            Cell< Node * > &
            abstract_nodes();

//-----------------------------------------------------------------------------

            // returns the orphaned nodes if they exist
            Cell< Node * > &
            orphaned_nodes();

//-----------------------------------------------------------------------------

            Cell< SideSet * > &
            cuts();

//-----------------------------------------------------------------------------

            bool
            create_thin_shell_cuts();

//-----------------------------------------------------------------------------

            void
            duplicate_nodes_on_face_sidesets();

            void
            relink_slave_elements_with_duplicate_nodes();

            void
            duplicate_and_relink_facets();

            void
            relink_non_thinshell_facets();

            void
            close_terminal_loops();

//-----------------------------------------------------------------------------

            void
            orient_terminal_curves();

            void
            orient_terminal_curves_sub(
                SideSet * aSideSet,
                Vector< real > & aWorkA,
                Vector< real > & aWorkB,
                Vector< real > & aN );

            void
            orient_terminal_curves_2D();

//-----------------------------------------------------------------------------

            void
            save_curve_debug_meshes();

//-----------------------------------------------------------------------------

            inline Vector< id_t > &
            thin_shell_boundaries()
            {
                return mThinShellBoundaryIDs ;
            }

//-----------------------------------------------------------------------------

            Cell< Facet * > &
            thin_shell_facets( const id_t aID )
            {
                return mTemporaryThinShellSidesets( aID )->facets();
            }

//-----------------------------------------------------------------------------

            Cell< Node * > &
            thin_shell_master_nodes()
            {
                return mThinShellMasterNodes ;
            }

//-----------------------------------------------------------------------------

            Cell< Node * > &
            thin_shell_slave_nodes()
            {
                return mThinShellSlaveNodes ;
            }

//-----------------------------------------------------------------------------
        private:
//-----------------------------------------------------------------------------


            void
            link_node_duplicates_and_originals();

//-----------------------------------------------------------------------------

            void
            collect_orphan_nodes();

//-----------------------------------------------------------------------------
            void
            collect_nodes_and_elements_on_blocks(
                const Vector< id_t > & aBlockIDs,
                Cell< Element * > & aElements,
                Cell< Node * > & aNodes );

//-----------------------------------------------------------------------------

            void
            compute_rcm_problem();

//-----------------------------------------------------------------------------

            void
            compute_cohomologies();

//-----------------------------------------------------------------------------

            void
            compute_thin_cuts_and_duplicate_interface_nodes();

//-----------------------------------------------------------------------------

            void
            restore_thin_shell_sidesets();

//-----------------------------------------------------------------------------

            void
            write_debug_cohomology( Homology * aSuggestedHomology ) ;

//-----------------------------------------------------------------------------

            void
            create_curves_for_thinshells();

            void
            create_side_curves_for_thinshells_3d();

            void
            collect_boundary_sidesets( Vector< id_t > & aIDs );

//-----------------------------------------------------------------------------

            ElementType
            check_element_types() ;

//-----------------------------------------------------------------------------

            // for debugging
            void
            save_edges( const uint aIndex );


//-----------------------------------------------------------------------------

            SideSet *
            create_cut_sideset_2d(
                    const uint aIndex,
                    const ElementType aType,
                    Cell< Element * > & aElements,
                    const Vector< int > & aCases );

//-----------------------------------------------------------------------------

            SideSet *
            create_cut_sideset_3d(
                    const uint aIndex,
                    const ElementType aType,
                    Cell< Element * > & aElements,
                    const Vector< int > & aCases );

//-----------------------------------------------------------------------------


            void
            flag_nodes_and_facets_of_tape_sidesets( Cell< Node * > & aNodes ) ;

//-----------------------------------------------------------------------------

            Facet *
            create_facet( Element * aElement, Cell< Element * > & aCandidates );

//-----------------------------------------------------------------------------

            bool
            connect_facet_to_slave( Facet * aFacet );

//-----------------------------------------------------------------------------

            void
            unflag_symmetry_sidesets();

//-----------------------------------------------------------------------------

            void
            compute_element_adjacencies();

//-----------------------------------------------------------------------------

            void
            compute_poisson_problem();

        };
    }
}
#endif //BELFEM_CL_CUTFACTORY_HPP
