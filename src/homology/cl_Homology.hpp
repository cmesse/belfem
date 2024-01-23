//
// Created by grgia@ge.polymtl.ca on 18/12/23.
//

#include "cl_Mesh.hpp"
#include "cl_Chain.hpp"
#include "cl_SimplicialComplex.hpp"
#include "fn_Smith.hpp"
#include "cl_Map.hpp"
#include "cl_Cell.hpp"

#ifndef BELFEM_CL_HOMOLOGY_HPP
#define BELFEM_CL_HOMOLOGY_HPP

namespace belfem
{
    namespace mesh
    {
//-----------------------------------------------------------------------------
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

            SimplicialComplex* mSimplicialComplex;

            bool mflagProp = false;

            //Mesh
            Mesh * mMesh ;


            //-----------------------------------------------------------------------------
        public:

            //-----------------------------------------------------------------------------

            Homology( Mesh * aMesh);

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

            /*void
            createMatrixFromBoundaryMap();

            //-----------------------------------------------------------------------------

            Vector< int >
            canonicalCoordinates(Chain* aC);*/

            //-----------------------------------------------------------------------------

            void
            create_kGeneratorsField(const uint k );

            //-----------------------------------------------------------------------------

            void
            suggest_Homology();

            //-----------------------------------------------------------------------------

            Cell< Matrix< int > >
            get_BoundaryMatrix();

            //-----------------------------------------------------------------------------

            Cell <Cell< Chain * >>
            get_Generators();

            //-----------------------------------------------------------------------------

        };
    }
}

#endif //BELFEM_CL_HOMOLOGY_HPP
