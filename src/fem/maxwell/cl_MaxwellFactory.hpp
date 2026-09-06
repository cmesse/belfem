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

            Mesh * mMesh = nullptr ;
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
            //! make sure that masters follow the DomainType enum order:
            //! conductor > coil > ferro > buffer > air ( see en_DomainType.hpp )
            void
            fix_facet_masters();

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
