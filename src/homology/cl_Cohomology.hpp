//
// Created by grgia@ge.polymtl.ca on 18/12/23.
//

#include "cl_Mesh.hpp"
#include "cl_Cochain.hpp"
#include "cl_SimplicialComplex.hpp"
#include "cl_Chain.hpp"
#include "fn_Smith.hpp"
#include "cl_Map.hpp"
#include "cl_Cell.hpp"

#ifndef BELFEM_CL_COHOMOLOGY_HPP
#define BELFEM_CL_COHOMOLOGY_HPP

namespace belfem
{
    namespace mesh
    {
//-----------------------------------------------------------------------------
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

            SimplicialComplex* mSimplicialComplex;

            bool mflagProp = false;

            //Mesh
            Mesh * mMesh ;


            //-----------------------------------------------------------------------------
        public:

            //-----------------------------------------------------------------------------

            Cohomology( Mesh * aMesh, Mesh * aEdgeMesh );

            //-----------------------------------------------------------------------------

            Cohomology(SimplicialComplex * aSimplicialComplex, Mesh * aMesh);

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
            create_kGeneratorsField(const uint k, Mesh * tMesh );

            //-----------------------------------------------------------------------------

            void
            generatorsFromField(Mesh * tEdgeMesh);

            //-----------------------------------------------------------------------------

            Cell< Matrix< int > >
            get_CoboundaryMatrix();

            //-----------------------------------------------------------------------------

            Cell <Cell< Cochain * >>
            get_Generators();

            //-----------------------------------------------------------------------------

            void
            updatekGeneratorsFromHomology(Cell< Chain * > mkGenerators, const uint k);

            //-----------------------------------------------------------------------------

        };
    }
}

#endif //BELFEM_CL_COHOMOLOGY_HPP
