//
// Created by grgia@ge.polymtl.ca on 18/12/23.
//

#include "cl_Cohomology.hpp"
#include "fn_inv.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        Cohomology::Cohomology( Mesh * aMesh, Mesh * aEdgeMesh ) :
                mMesh(aMesh)
        {
            mGenerators.set_size(4,Cell< Cochain * >());
            mOrders.set_size(4,Cell< int >());
            this->generatorsFromField(aEdgeMesh);
        }

        //-----------------------------------------------------------------------------

        Cohomology::Cohomology( SimplicialComplex * aSimplicialComplex, Mesh * aMesh ) :
                mSimplicialComplex( aSimplicialComplex ),
                mMesh( aMesh )
        {
            mGenerators.set_size(4,Cell< Cochain * >());
            mOrders.set_size(4,Cell< int >());

            mD = mSimplicialComplex->createMatrixFromCoboundaryMap();
            this->cohomologyGroupOfChainComplex();
            this->generatorsOfCohomology();
        }

        //------------------------------------------------------------------------------

        Cohomology::~Cohomology()
        {
            for(uint i = 0; i < mGenerators.size(); i++)
            {
                for (uint j = 0; j < mGenerators(i).size(); j++)
                {
                    delete mGenerators(i)(j);
                }
            }
            mGenerators.clear();
            mOrders.clear();
            mV.clear();
            mW.clear();
            mU.clear();
            mB.clear();
            mD.clear();
        }

        //-----------------------------------------------------------------------------

        void
        Cohomology::cohomologyGroupOfChainComplex()
        {
            mV.clear();
            mW.clear();

            // loop over all dimensions
            for(uint k = 0; k < 4; k++)
            {
                // populate the matrices W and V (kernel and image of the coboundary)
                auto [w, v] = kernelImage(mD(k));
                mW[k] = w;
                mV[k+1] = v;
            }
            mV[0] = Matrix< int >(mW[0].n_rows(),1,0);

            //Compute the quotient groups
            this->quotientGroup();
        }

        //-----------------------------------------------------------------------------

        void
        Cohomology::quotientGroup()
        {
            uint n;

            //loop over all dimensions
            for(uint k = 0; k < 4; k++)
            {
                n = mV[k].n_cols();
                mB[k].set_size(mW[k].n_cols(),n,0);

                //Solve W\V to get the A matrix
                for (uint i = 0; i < n; i++)
                {
                    Matrix < int > tM = Matrix< int >(mV[k].col(i));
                    Matrix < int > tv = SolveInt(mW[k],tM);

                    for (uint j = 0; j < mW[k].n_cols(); j++)
                    {
                        mB[k](j,i) = tv(j,0);
                    }
                }

                // Compute the Smith normal form of W\V
                auto [tQ, tQ_, tR, tR_, s, t] = smithForm(mB[k]);
                ms[k] = s;
                mt[k] = t;

                //Get the U matrix, with the cohomology generators as the last columns.
                mU[k] = mW[k];
                mU[k]*=tQ;
            }
        }

        //-----------------------------------------------------------------------------

        void
        Cohomology::generatorsOfCohomology()
        {
            uint tCount;
            uint tCount2;

            //loop over all dimensions
            for (uint k = 0; k < 4; k++)
            {
                tCount = 0;

                // loop over the last columns of U
                for (uint j = ms[k]+1; j < mU[k].n_cols()+1; j++)
                {

                    //Get the order (0 means infinity... Maybe should change that)
                    if(j > mt[k])
                    {
                        mOrders(k).push(0);
                    }
                    else
                    {
                        mOrders(k).push(mB[k](j-1,j-1));
                    }

                    // Create the generator from column data
                    Cochain * tCochain = new Cochain(k,mMesh);
                    mGenerators(k).push(tCochain);
                    tCount2 = 0;
                    for(const auto& [tKey, tCochain2] : mSimplicialComplex->get_kcochainMap(k))
                    {
                        mGenerators(k)(tCount)->addCochainToCochain(tCochain2,mU[k](tCount2,j-1));
                        tCount2++;
                    }
                    tCount++;
                }
            }
        }

        //-----------------------------------------------------------------------------

        void
        Cohomology::create_kGeneratorsField(const uint k, Mesh * tMesh )
        {
            uint tCount = 0;
            for(const auto& tCochain : mGenerators(k))
            {
                tCount++;
                char fieldName [50];
                /*if (mflagProp)
                {
                    sprintf (fieldName, "1CohomologyGenerator%d", tCount);
                }
                else
                {
                    sprintf (fieldName, "1CohomologyGeneratorComp%d", tCount);
                }*/
                sprintf (fieldName, "1CohomologyGenerator%d", tCount);
                tMesh->create_field( fieldName, EntityType::ELEMENT);
                Map< id_t, int > tSimplicesMap = tCochain->getSimplicesMap();
                for( const auto& [tInd, tCoeff] :  tSimplicesMap)
                {
                    tMesh->field_data(fieldName)(tMesh->element(tInd)->index()) = tCoeff;
                    //tMesh->field_data(fieldName)(tMesh->edge(tInd)->node(1)->index()) = 1;
                }
            }
        }

        //-----------------------------------------------------------------------------

        void
        Cohomology::generatorsFromField(Mesh * tEdgeMesh)
        {
            char fieldName [50];
            uint tCount = 1 ;
            sprintf (fieldName, "1CohomologyGenerator%d", tCount);
            while ( tEdgeMesh->field_exists(fieldName) )
            {
                mGenerators(1).push(new Cochain(1, mMesh));
                for(uint i = 0; i < tEdgeMesh->number_of_elements(); i++)
                {
                    mGenerators(1)(tCount-1)->addSimplexToCochain(tEdgeMesh->elements()(i)->id(),tEdgeMesh->field_data(fieldName)(i));
                }
                tCount++;
                sprintf (fieldName, "1CohomologyGenerator%d", tCount);
            }
        }

        //-----------------------------------------------------------------------------

        Cell< Matrix< int > >
        Cohomology::get_CoboundaryMatrix()
        {
            return this->mD;
        }

        //-----------------------------------------------------------------------------

        Cell <Cell< Cochain * >>
        Cohomology::get_Generators()
        {
            return this->mGenerators;
        }

        //-----------------------------------------------------------------------------

        void
        Cohomology::updatekGeneratorsFromHomology(Cell< Chain * > tkGenerators, const uint k)
        {
            const uint n = mGenerators(k).size();
            BELFEM_ASSERT(tkGenerators.size() == n, "Error, not the same number of k-generators" );
            Matrix< real > tTransMat = Matrix< real >(mGenerators(k).size(),mGenerators(k).size(),0);

            //Populate the transformation matrix
            for(uint i = 0; i < n; i++)
            {
                Cochain * tCoGen = mGenerators(k)(i);
                for(uint j = 0; j < n; j++)
                {
                    tTransMat(i,j) = tCoGen->operator()(tkGenerators(j));
                }
            }

            Matrix< real > tTransInv = inv(tTransMat);

            Cell< Cochain * > tNewGenerators = Cell< Cochain * >();
            for(uint i = 0; i < n; i++)
            {
                Cochain * tGen = new Cochain(k, mMesh);
                for(uint j = 0; j < n; j++)
                {
                    tGen->addCochainToCochain(mGenerators(k)(j),tTransInv(i,j));
                }
                tNewGenerators.push(tGen);
            }

            for(uint i = 0; i < n; i++)
            {
                delete mGenerators(k)(i);
            }
            mGenerators(k) = tNewGenerators;

            mflagProp = true;
        }

        //-----------------------------------------------------------------------------
    }
}
//------------------------------------------------------------------------------