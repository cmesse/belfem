//
// Created by grgia@ge.polymtl.ca on 18/12/23.
//

//todo: Maybe create a parent class for homology and cohomology, as they have very similar methods

#include "cl_Homology.hpp"
#include "en_FEM_DomainType.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        Homology::Homology( Mesh * aMesh, Mesh * aEdgeMesh ):
                mMesh( aMesh )
        {
            mGenerators.set_size(4,Cell< Chain * >());
            mOrders.set_size(4,Cell< int >());
            this->generatorsFromField(aEdgeMesh);
        }

        Homology::Homology( Mesh * aMesh ) :
                mMesh( aMesh )
        {
            this->suggest_Homology();
        }

        //------------------------------------------------------------------------------

        Homology::Homology( SimplicialComplex * aSimplicialComplex, Mesh * aMesh ) :
                mSimplicialComplex( aSimplicialComplex ),
                mMesh( aMesh )
        {
            mGenerators.set_size(4,Cell< Chain * >());
            mOrders.set_size(4,Cell< int >());

            mD = mSimplicialComplex->createMatrixFromBoundaryMap();
            this->homologyGroupOfChainComplex();
            this->generatorsOfHomology();
        }

        //------------------------------------------------------------------------------

        Homology::~Homology()
        {
            this->reset();
            if (mflagProp)
            {
                delete mSimplicialComplex;
            }
        }

        //-----------------------------------------------------------------------------

        void
        Homology::reset()
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
        Homology::homologyGroupOfChainComplex()
        {
            mV.clear();
            mW.clear();

            // loop over all dimensions
            for(uint k = 0; k < 4; k++)
            {
                // populate the matrices W and V (kernel and image of the boundary)
                auto [w, v] = kernelImage(mD(k));
                mW[k] = w;
                mV[k-1] = v;
            }
            mV[3] = Matrix< int >(mW[3].n_rows(),1,0);

            //Compute the quotient groups
            this->quotientGroup();
        }

        //-----------------------------------------------------------------------------

        void
        Homology::quotientGroup()
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

                //Get the U matrix, with the homology generators as the last columns.
                mU[k] = mW[k];
                mU[k]*=tQ;
            }
        }

        //-----------------------------------------------------------------------------

        void
        Homology::generatorsOfHomology()
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
                    Chain * tChain = new Chain(k,mMesh);
                    mGenerators(k).push(tChain);
                    tCount2 = 0;
                    for(const auto& [tKey, tChain2] : mSimplicialComplex->get_kchainMap(k))
                    {
                        mGenerators(k)(tCount)->addChainToChain(tChain2,mU[k](tCount2,j-1));
                        tCount2++;
                    }
                    tCount++;
                }
            }
        }

        //-----------------------------------------------------------------------------


        // Create field in the mesh to visualize the generators
        void
        Homology::create_kGeneratorsField(const uint k, Mesh * tMesh )
        {
            uint tCount = 0;
            for(const auto& tChain : mGenerators(k))
            {
                tCount++;
                char fieldName [50];
                /*if (mflagProp)
                {
                    sprintf (fieldName, "1HomologyGenerator%d", tCount);
                }
                else
                {
                    sprintf (fieldName, "1HomologyGeneratorComp%d", tCount);
                }*/

                sprintf (fieldName, "1HomologyGenerator%d", tCount);
                tMesh->create_field( fieldName, EntityType::ELEMENT);
                Map< id_t, int > tSimplicesMap = tChain->getSimplicesMap();
                for( const auto& [tInd, tCoeff] :  tSimplicesMap)
                {
                    tMesh->field_data(fieldName)(tMesh->element(tInd)->index()) = tCoeff;
                    //mMesh->field_data(fieldName)(mMesh->edge(tInd)->node(1)->index()) = 1;
                }
            }
        }

        //-----------------------------------------------------------------------------

        void
        Homology::generatorsFromField(Mesh * tEdgeMesh)
        {
            char fieldName [50];
            uint tCount = 1 ;
            sprintf (fieldName, "1HomologyGenerator%d", tCount);
            while ( tEdgeMesh->field_exists(fieldName) )
            {
                mGenerators(1).push(new Chain(1, mMesh));
                for(uint i = 0; i < tEdgeMesh->number_of_elements(); i++)
                {
                    mGenerators(1)(tCount-1)->addSimplexToChain(tEdgeMesh->elements()(i)->id(),tEdgeMesh->field_data(fieldName)(i));
                }
                tCount++;
                sprintf (fieldName, "1HomologyGenerator%d", tCount);
            }
        }

        //-----------------------------------------------------------------------------

        void
        Homology::suggest_Homology()
        {
            //Reset attributes
            this->reset();
            mGenerators.set_size(4,Cell< Chain * >());
            mOrders.set_size(4,Cell< int >());

            //Flag the desired simplices (around the conductors)
            for( mesh::SideSet * tSideSet : mMesh->sidesets() )
            {
                // make sure that sideset is not empty
                if ( tSideSet->number_of_facets() > 0 )
                {
                    // get first facet
                    mesh::Facet * tFacet = tSideSet->facet_by_index( 0 );

                    // make sure that it has a slave
                    if ( tFacet->has_slave())
                    {
                        fem::DomainType tMasterType
                                = static_cast< fem::DomainType >( tFacet->master()->physical_tag());

                        fem::DomainType tSlaveType
                                = static_cast< fem::DomainType >( tFacet->slave()->physical_tag());

                        // check if this sideset is an interface
                        if ( (tMasterType == fem::DomainType::Coil || tMasterType == fem::DomainType::Conductor) && tSlaveType == fem::DomainType::Air )
                        {
                            // Flag edges and nodes of conductor boundaries
                            for ( uint i = 0; i < tSideSet->number_of_facets(); i++ )
                            {
                                mMesh->edge(tSideSet->facet_by_index( i )->id())->flag();
                                for ( uint j = 0; j < mMesh->edge(tSideSet->facet_by_index( i )->id())->number_of_nodes(); j++ )
                                {
                                    mMesh->edge(tSideSet->facet_by_index( i )->id())->node(j)->flag();
                                }
                            }
                        }
                    }
                }
            }

            //Create a new simplicial complex for the suggested homology
            mSimplicialComplex = new SimplicialComplex(mMesh);
            mSimplicialComplex->reduce_complexCCR();
            mMesh->unflag_everything();
            mflagProp = true;

            //Compute the suggested homology
            mD = mSimplicialComplex->createMatrixFromBoundaryMap();
            this->homologyGroupOfChainComplex();
            this->generatorsOfHomology();
        }

        Cell< Matrix< int > >
        Homology::get_BoundaryMatrix()
        {
            return this->mD;
        }

        //-----------------------------------------------------------------------------

        Cell <Cell< Chain * >>
        Homology::get_Generators()
        {
            return this->mGenerators;
        }

        //-----------------------------------------------------------------------------
    }
}
//------------------------------------------------------------------------------