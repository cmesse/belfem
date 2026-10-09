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

#ifndef CL_MAXWELLFACTORY_HPP
#define CL_MAXWELLFACTORY_HPP
#include <en_CutAlgorithm.hpp>

#include "typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Mesh.hpp"
#include "en_DomainType.hpp"
#include "cl_Vector.hpp"
#include "cl_FEM_Kernel.hpp"
#include "cl_IWG_Maxwell.hpp"
#include "cl_MaxwellBoundaryConditionFactory.hpp"
#include "cl_InputFile.hpp"
#include "cl_FEM_Domain.hpp"
#include "cl_Topology.hpp"
#include "cl_Protoshell.hpp"

namespace belfem
{
    namespace fem
    {
        /**
         * @brief High-level orchestrator that builds a Maxwell problem from an input deck.
         *
         * @ingroup grp_fem_maxwell
         * @see @ref fem_maxwell_maxwell_usage_guide
         */
        class MaxwellFactory
        {

            const proc_t mCommRank ;

            const InputFile * mInputFile  ;

            // -- user settings

            mesh::CutAlgorithm mCutAlgorithm = mesh::CutAlgorithm::PellikkaGeneralized ;

            bool mSaveElementConnectivitiesToBfm = false ;
            // --- end user settings
            // deleted by kernel
            bool mOwnKernelParameters = true ;
            KernelParameters * mKernelParameters = nullptr ;

            //! the mesh this factory creates ( the file mesh on the master,
            //! an empty one on the workers ); owned until create_magnetic_kernel()
            //! hands it to the kernel, borrowed for reading afterwards
            Mesh * mMesh = nullptr ;
            bool mOwnMesh = true ;
            string mMeshPath ;

            bool mComputeCohomologies = true ;

            maxwell::Formulation mFormulation = maxwell::Formulation::HPhi ;

            ModelDimensionality  mDimensionality = ModelDimensionality::UNDEFINED ;
            bool                 mUseEnrichment = false ;

            Vector< proc_t > mCommTable ;

            uint mMeshDimension = 0 ;
            uint mElementOrder = 0 ;

            std::shared_ptr< Kernel > mMagneticKernel ;

            IWG_Maxwell * mMagneticEquation = nullptr ;
            bool mOwnMagneticEquation = true ;

            DofManager * mMagneticField = nullptr ;
            mesh::Topology * mTopology = nullptr ;

            Cell< Domain * > mDomains ;

            Cell< mesh::Curve * > mCurves ;

            MaxwellBoundaryConditionFactory * mBoundaryConditionFactory = nullptr ;
            mesh::PeriodicityFactory * mPeriodicFactory = nullptr ;

            //Terminal information to send to the cut factory
            Cell< Cell<id_t> > mTerminals;
            Cell<id_t> mThinShellTerminalIndices ;

            Map< id_t, string > mMaterialBlockAssignment ;

            Vector< index_t > mHangingFacetIndices ;
            Matrix< real > mHangingFacetTMatrices ;

            Cell<  Protoshell * > mProtoshells ;

            Map< string, id_t > mMaterialIDs ;
            Map<string, Material *> mMaterialMap ;

            // these sidesets are the templates from which we create
            // the thin shells
            Cell< mesh::SideSet * > mTapes ;
            Map< id_t, mesh::ThinShell * > mThinShellMap ;

            Cell< Cell< index_t > > mThinShellNodeIndices ;
            Cell< mesh::Node * >    mThinShellMasterNodes ;
            Cell< mesh::Node * >    mThinShellSlaveNodes ;

            struct EdgeWorkData
            {
                Cell< mesh::Node * > NodesOnThinShell ;
                Cell< mesh::Node * > NodesOnVolume ;

                Cell< mesh::Edge * > EdgesOnThinShell ;
                Cell< mesh::Edge * > EdgesOnVolume ;
            };

            string mLabel ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            MaxwellFactory( const string & aInputFile );

            ~MaxwellFactory();

            // owns raw pointers ( the mesh until handover, input file,
            // factories, topology ); a copy would double-delete them
            MaxwellFactory( const MaxwellFactory & ) = delete ;
            MaxwellFactory & operator=( const MaxwellFactory & ) = delete ;
            MaxwellFactory( MaxwellFactory && ) = delete ;
            MaxwellFactory & operator=( MaxwellFactory && ) = delete ;

            const string &
            mesh_path() const ;

//------------------------------------------------------------------------------

            std::shared_ptr< Kernel >
            create_magnetic_kernel();

//------------------------------------------------------------------------------

            std::shared_ptr< Controller >
            create_controller() ;

//------------------------------------------------------------------------------

            Cell <PhysicalBoundaryCondition *> &
            boundary_conditions() ;

//------------------------------------------------------------------------------

            Cell <PhysicalBoundaryCondition *>
            current_BCs() ;

//------------------------------------------------------------------------------

            Mesh *
            mesh();

//------------------------------------------------------------------------------

            string
            label() const ;

//------------------------------------------------------------------------------

            /**
             * whether a constant permeability counts as mu0 for the thin-shell
             * normal-field recovery beside an h-conductor: within a relative
             * rounding allowance of 1e-6, so mu0 written to seven significant
             * digits passes and a real conductor's susceptibility does not.
             * False for NaN
             */
            static bool
            mu_is_mu0( const real aMu );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            Mesh *
            read_mesh();

//------------------------------------------------------------------------------

            void
            read_domain_types();

//------------------------------------------------------------------------------

            /**
             * print the lines in which the stored mesh configuration differs
             * from the current one, so a rebuild says WHICH setting changed
             */
            void
            report_config_difference(
                    const string & aStored,
                    const string & aCurrent );

//------------------------------------------------------------------------------

            void
            create_curves( const input::Section * aSection );

//------------------------------------------------------------------------------

            void
            create_periodic( const input::Section * aSection );

//------------------------------------------------------------------------------

            IWG_Maxwell *
            create_equation( const maxwell::Formulation aFormulation );

//------------------------------------------------------------------------------

            void
            set_block_types_in_magnetic_equation();

//------------------------------------------------------------------------------
            //! orient the facets for the cut, rank 0 only. Between two
            //! domain types, the master follows the DomainType enum order:
            //! conductor > coil > ferro > buffer > air ( see
            //! en_DomainType.hpp ). In 3D, facets with one domain type on
            //! both sides inherit the orientation of their neighbors, and
            //! then every thin-shell sideset is made consistent: each
            //! connected patch takes one winding by a vote, which may give
            //! some cross-type facets the lower domain type as master;
            //! patches with the same two blocks on their sides get the same
            //! master block; and a sideset with one block on both of its
            //! sides is refused
            void
            fix_facet_masters();

//------------------------------------------------------------------------------

            //! give every connected patch of each thin-shell sideset one
            //! consistent winding, and the patches of one sideset that have
            //! the same two blocks on their sides one master block; 3D only,
            //! needs facet-to-facet
            //! connectivity, runs at the end of fix_facet_masters and
            //! before flip_thin_shell_sidesets
            void
            align_thin_shell_facets();

//------------------------------------------------------------------------------

            //! refuse a thin-shell sideset with a facet that lacks an element
            //! on one side, or with a block that is master on some facets and
            //! slave on others: the cut opens a shell block by block. 3D,
            //! after align_thin_shell_facets
            void
            check_thin_shell_block_sides();

//------------------------------------------------------------------------------

            //! parse a thin-shell sideset list with gmsh-style signs
            //! ( sidesets : -5, -6, 7:20 ; ): absolute values go to
            //! sidesets(), signed ids additionally to flipped_sidesets()
            void
            read_signed_sidesets(
                    const input::Section * aSection,
                    const string         & aKey,
                    Protoshell           * aShell );

//------------------------------------------------------------------------------

            //! flip the orientation of every facet of the deck-signed
            //! sidesets; runs after fix_facet_masters ( which would undo
            //! it ) and before the cut / thin-shell pipeline reads the
            //! facet windings for the layer normals
            void
            flip_thin_shell_sidesets();

//------------------------------------------------------------------------------

            void
            create_cuts();

//------------------------------------------------------------------------------

            void
            create_cuts_sub_master();

//------------------------------------------------------------------------------

            void
            set_block_and_sideset_names();

//------------------------------------------------------------------------------

            void
            create_cuts_sub_slave();

//------------------------------------------------------------------------------

            void
            create_thinshells();

//------------------------------------------------------------------------------

            void
            create_terminal_list();

//------------------------------------------------------------------------------

            void
            create_edges_and_faces_on_mesh();

//------------------------------------------------------------------------------

            void
            create_hanging_edges_and_facets();

//------------------------------------------------------------------------------

            void
            create_postprocessors();

//------------------------------------------------------------------------------

            void
            init_fields();

//------------------------------------------------------------------------------

            void
            configure_solver( const input::Section * aSection, DofManager * aField );

//------------------------------------------------------------------------------

            void
            create_block_to_material_map();

//------------------------------------------------------------------------------

            void
            collect_material_labels_from_domains( Vector< id_t > & aBlockIDs, Cell< string > & aMaterialLabels );

//------------------------------------------------------------------------------

            void
            set_physical_tags_for_elements();

//------------------------------------------------------------------------------

            void
            synch_material_map();

//------------------------------------------------------------------------------

            void
            create_materials();

            void
            assign_materials();

            void
            delete_unused_materials();

//------------------------------------------------------------------------------

            void
            read_thin_shell_data();

//------------------------------------------------------------------------------

            void
            find_autopins( Cell< mesh::Node * > & aPins );

//------------------------------------------------------------------------------

            /**
             * fixes the automatic gauge pins on the master mesh, one per
             * connected phi component, except in a component whose phi is
             * already fixed by a symmetry plane or a dirichlet sideset.
             * Called when the deck has no bearing
             */
            void
            impose_autopins();

//------------------------------------------------------------------------------

            /**
             * walks the phi component of aPin and returns the first node
             * flagged on plane 1 ( a fixed phi ), or nullptr. Element flag 0
             * must mark the phi blocks. aVisited is the queue and, on return,
             * the list of nodes flagged on plane 0 by this search
             */
            mesh::Node *
            find_fixed_phi_contact(
                    mesh::Node           * aPin,
                    Cell< mesh::Node * > & aVisited );

//------------------------------------------------------------------------------

            /**
             * pushes the unvisited originals of the flagged elements of aNode
             */
            void
            push_phi_neighbors(
                    mesh::Node           * aNode,
                    Cell< mesh::Node * > & aVisited );

//------------------------------------------------------------------------------

            void
            hang_thinshell_edges_on_nodes_bottom(
                EdgeWorkData  & aWork,
                mesh::Facet   * aFacets,
                mesh::Element * aElement );

            void
            hang_thinshell_edges_on_nodes_top(
                EdgeWorkData  & aWork,
                mesh::Facet   * aFacets,
                mesh::Element * aElement );

            void
            hang_thinshell_edges_on_edges_bottom(
                EdgeWorkData  & aWork,
                mesh::Facet   * aFacets,
                mesh::Element * aElement );

            void
            hang_thinshell_edges_on_edges_top(
                EdgeWorkData  & aWork,
                mesh::Facet   * aFacets,
                mesh::Element * aElement );

//------------------------------------------------------------------------------
        };

    }
}

#endif //CL_MAXWELLFACTORY_HPP
