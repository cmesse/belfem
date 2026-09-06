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

#include "cl_Mesh.hpp"
#include "cl_Cochain.hpp"
#include "cl_BeltedTree.hpp"
#include "cl_SimplicialComplex.hpp"
#include "cl_Chain.hpp"
#include "fn_Smith.hpp"
#include "cl_Map.hpp"
#include "cl_Cell.hpp"
#include "cl_Progressbar.hpp"

#ifndef BELFEM_CL_COHOMOLOGY_HPP
#define BELFEM_CL_COHOMOLOGY_HPP

namespace belfem
{
    namespace mesh
    {
//-----------------------------------------------------------------------------
        /**
         * @brief Computes the cohomology groups H^k of a mesh complex.
         *
         * @ingroup grp_homology
         * @see @ref homology_cohomology_theory_and_implementation
         */
        class Cohomology
        {

            Cell <Cell< Cochain * >> mGenerators;

            Cell <Cell< int >> mOrders;

            Map<int, Matrix< int >> mV;

            Map<int, Matrix< int >> mW;

            Map<int, Matrix< int >> mU;

            Map<int, Matrix< int >> mB;

            Map<int, uint > ms;

            Map<int, uint > mt;

            Cell< Matrix< int > > mD;

            SimplicialComplex * mSimplicialComplex = nullptr;

            bool mflagProp = false;

            Mesh * mMesh ;

            Progressbar * mProgress = nullptr ;


//-----------------------------------------------------------------------------
        public:
//-----------------------------------------------------------------------------

            Cohomology( SimplicialComplex * aSimplicialComplex, Mesh * aMesh );

//-----------------------------------------------------------------------------

            Cohomology( SimplicialComplex * aSimplicialComplex, Mesh * aMesh, BeltedTree * aBTree );

//-----------------------------------------------------------------------------

            ~Cohomology();

//-----------------------------------------------------------------------------

            void
            cohomologyGroupOfChainComplex();

//-----------------------------------------------------------------------------

            void
            quotientGroup();

//-----------------------------------------------------------------------------

            void
            generatorsOfCohomology();

//-----------------------------------------------------------------------------

            void
            clean();

//-----------------------------------------------------------------------------

            void
            check();

//-----------------------------------------------------------------------------

            void
            create_kGeneratorsField(const uint k, Mesh * aMesh, string aFieldName );

//-----------------------------------------------------------------------------

            Cell< Matrix< int > > &
            get_CoboundaryMatrix();

//-----------------------------------------------------------------------------

            Cell <Cell< Cochain * >> &
            get_Generators();

//-----------------------------------------------------------------------------

            void
            updatekGeneratorsFromHomology(Cell< Chain * > & tkGenerators, const uint k);

//-----------------------------------------------------------------------------

            Matrix< real >
            coefficient_TMatrix(Cell< Chain * > & mkGenerators, const uint k);

//-----------------------------------------------------------------------------
        private:
 //-----------------------------------------------------------------------------

            void
            clean_spfa();

            /**
             * pocket census + removal on the rectified generators
             * ( todo/cut_pocket_removal_rules.md ). Logs every zero-graph
             * component per generator; fires Tier-A pure pockets if
             * aFireTierA is set. Precondition: called from clean_spfa()
             * while the periodic slave flags are still set.
             */
            void
            remove_cut_pockets( const bool aFireTierA );

        };
    }
}

#endif //BELFEM_CL_COHOMOLOGY_HPP
