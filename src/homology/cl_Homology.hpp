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
#include "cl_Chain.hpp"
#include "cl_SimplicialComplex.hpp"
#include "cl_Element_Factory.hpp"
#include "fn_Smith.hpp"
#include "cl_Map.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"

#ifndef BELFEM_CL_HOMOLOGY_HPP
#define BELFEM_CL_HOMOLOGY_HPP

namespace belfem
{
    namespace mesh
    {
//-----------------------------------------------------------------------------
        /**
         * @brief Computes the homology groups H_k of a mesh complex.
         *
         * @ingroup grp_homology
         * @see @ref homology_homology_usage_guide
         */
        class Homology
        {

            Cell <Cell< Chain * >> mGenerators;

            Cell <Cell< int >> mOrders;

            Map<int, Matrix< int >> mV;

            Map<int, Matrix< int >> mW;

            Map<int, Matrix< int >> mU;

            Map<int, Matrix< int >> mB;

            Map<int, uint > ms;

            Map<int, uint > mt;

            Cell< Matrix< int > > mD;

            SimplicialComplex * mSimplicialComplex;

            //Mesh
            Mesh * mMesh ;


//-----------------------------------------------------------------------------
        public:

//-----------------------------------------------------------------------------

            Homology( Mesh * aMesh, Cell<Cell< id_t >> aTerminals,
                      Cell< id_t > aThinShellIndices ) ;

//-----------------------------------------------------------------------------

            Homology(SimplicialComplex * aSimplicialComplex, Mesh * aMesh);

//-----------------------------------------------------------------------------

            ~Homology();

//-----------------------------------------------------------------------------

            void
            reset();

//-----------------------------------------------------------------------------

            void
            homologyGroupOfChainComplex();

//-----------------------------------------------------------------------------

            void
            quotientGroup();

//-----------------------------------------------------------------------------

            void
            generatorsOfHomology();

//-----------------------------------------------------------------------------

            void
            create_kGeneratorsField(const uint k, Mesh * aMesh, string tFieldName );

//-----------------------------------------------------------------------------

            void
            suggest_Homology(Cell<Cell< id_t >> aTerminals,
                             Cell< id_t > aThinShellIndices );

//-----------------------------------------------------------------------------

            Cell< Matrix< int > > &
            get_BoundaryMatrix();

//-----------------------------------------------------------------------------

            Cell <Cell< Chain * >> &
            get_Generators();

//-----------------------------------------------------------------------------

            void
            create_self_intersecting(const uint tGeneratorIndex);

//-----------------------------------------------------------------------------

            Cell < int >
            generators_orientation(Cell< Vector <real> > & aDirections);

//-----------------------------------------------------------------------------

            void
            reorient_generators();

//-----------------------------------------------------------------------------

        };
    }
}

#endif //BELFEM_CL_HOMOLOGY_HPP
