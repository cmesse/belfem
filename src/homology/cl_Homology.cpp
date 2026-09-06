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

#include "cl_Homology.hpp"
#include "en_DomainType.hpp"
#include "cl_Logger.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        Homology::Homology( Mesh * aMesh, Cell<Cell< id_t >> aTerminals,
                            Cell< id_t > aThinShellIndices) :
                mMesh( aMesh )
        {
            //start by suggesting homology for the conductors that have input terminals
            this->suggest_Homology(aTerminals, aThinShellIndices );

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

//--------------------delete tHomology ;----------------------------------------------------------

        Homology::~Homology()
        {
            this->reset();
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

//------------------------------------------------------------------------------

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

                    Matrix < int > tM ( mV(k).n_rows(), 1 );
                    tM.set_col( 0, mV(k).col(i) );
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

//------------------------------------------------------------------------------

        void
        Homology::generatorsOfHomology()
        {
            uint tCount;
            uint tCount2;

            //loop over all dimensions
            for (uint k = 0; k < 4; k++)
            {
                tCount = mGenerators(k).size();

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
                    Chain * tChain = new Chain(k,mMesh, false, false);
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

//------------------------------------------------------------------------------


        // Create field in the mesh to visualize the generators
        void
        Homology::create_kGeneratorsField(const uint k, Mesh * aMesh, string tFieldName )
        {
            Cell < Element * > & tElements = aMesh->elements();
            uint tCount = 0;
            for(const auto& tChain : mGenerators(k))
            {
                tCount++;
                aMesh->create_field( tFieldName+std::to_string(tCount), EntityType::ELEMENT);
                OrderedMap< index_t, int > & tSimplicesMap = tChain->getSimplicesMap();
                for( const auto& [tID, tCoeff] :  tSimplicesMap)
                {
                    aMesh->field_data(tFieldName+std::to_string(tCount))(tElements(tID)->index()) = tCoeff;
                    /*if ( mMesh->edges()(tID)->is_periodic() )
                    {
                        std::cout << "HERE" << std::endl;
                        aMesh->field_data(tFieldName+std::to_string(tCount))(tElements(mMesh->edges()(tID)->periodic()->index())->index()) = tCoeff;
                    }*/
                }
            }
        }

//-----------------------------------------------------------------------------

        void
        Homology::suggest_Homology(Cell<Cell< id_t >> aTerminals,
                                   Cell< id_t > aThinShellIndices )
        {

            //Reset attributes
            this->reset();
            mGenerators.set_size(4,Cell< Chain * >());
            mOrders.set_size(4,Cell< int >());
            mMesh->unflag_everything() ;

            uint tCountThinshells = 0 ;

            // Loop over the list of terminals. The size of aTerminals is the total number of terminal pairs
            // (bulk and thin shells), but the thin shells are represented by empty lists in aTerminals.
            // Therefore we get the terminal information from thin shells from aThinShellTerminals
            for (uint i = 0 ; i < aTerminals.size() ; ++i)
            {
                bool aIsThinShell = false ;
                if (tCountThinshells+1 <= aThinShellIndices.size())
                {
                    if (i == aThinShellIndices(tCountThinshells)) aIsThinShell = true ;
                }
                //For bulk conductors (voltage or current terminals, doesn't matter)
                if (!aIsThinShell)
                {
                    //Only input terminals in 2-D
                    for (uint j = 0; j < aTerminals(i).size()/2; ++j)
                    {
                        if (mMesh->number_of_dimensions() == 2)
                        {
                            mMesh->block( aTerminals( i )( j ))->flag_edges();
                            mMesh->block( aTerminals( i )( j ))->flag_corner_nodes();
                            mMesh->block( aTerminals( i )( j ))->flag_elements();
                        }
                        else
                        {
                            mMesh->sideset( aTerminals( i )( j ))->flag_edges();
                            mMesh->sideset( aTerminals( i )( j ))->flag_corner_nodes();
                            for (Face * tFace: mMesh->faces())
                            {
                                uint tCount = 0;
                                for ( uint k = 0; k < tFace->number_of_edges(); ++k )
                                {
                                    if (tFace->edge( k )->is_flagged())
                                    {
                                        tCount++;
                                    }
                                }
                                if (tCount == tFace->number_of_edges())
                                {
                                    tFace->flag();
                                }
                            }
                        }
                    }

                    //Initialize the generator
                    mGenerators(1).push( new Chain(1,mMesh,true,false) ) ;

                    //Initialize the terminal 2-chain
                    Chain * tTotalChain = new Chain(2,mMesh,false,false);

                    //First, for this terminal, check if it touches a curve. If so, we add it to the suggested homology
                    for (Curve * tCurve : mMesh->curves())
                    {
                        bool tAddCurve = false ;
                        for ( Segment * tSegment : tCurve->segments() )
                        {
                            if (tSegment->edge()->is_flagged())
                            {
                                tAddCurve = true;
                                break ;
                            }
                        }
                        if (tAddCurve)
                        {
                            for ( Segment * tSegment : tCurve->segments() )
                            {
                                Edge * tEdge = tSegment->edge() ;
                                mGenerators(1)(mGenerators(1).size()-1)->addSimplexToChain(tEdge->index(),-1*tSegment->edge_direction()) ;
                            }
                        }
                    }

                    mSimplicialComplex = new SimplicialComplex(mMesh);
                    for (auto [tID, tChain] : mSimplicialComplex->get_kchainMap(2))
                    {
                        tTotalChain->addChainToChain(tChain,1) ;
                    }

                    //The suggested homology is simply the boundary of the input terminal
                    mGenerators(1)(mGenerators(1).size()-1)->addChainToChain(tTotalChain->getBoundary(),-1) ;

                    mMesh->unflag_everything() ;
                    delete mSimplicialComplex ;

                    //Also output terminal for 3D
                    if (mMesh->number_of_dimensions() == 3)
                    {
                        for (uint j = aTerminals(i).size()/2; j < aTerminals(i).size(); ++j)
                        {
                            mMesh->sideset( aTerminals( i )( j ))->flag_edges();
                            mMesh->sideset( aTerminals( i )( j ))->flag_corner_nodes();
                            for (Face * tFace: mMesh->faces())
                            {
                                uint tCount = 0;
                                for (uint k = 0; k < tFace->number_of_edges(); ++k)
                                {
                                    if (tFace->edge( k )->is_flagged())
                                    {
                                        tCount++;
                                    }
                                }
                                if (tCount == tFace->number_of_edges())
                                {
                                    tFace->flag();
                                }
                            }
                        }

                        //Initialize the generator
                        mGenerators(1).push( new Chain(1,mMesh,true,false) ) ;

                        //Initialize the terminal 2-chain
                        Chain * tTotalChain2 = new Chain(2,mMesh,false,false);

                        //First, for this terminal, check if it touches a curve. If so, we add it to the suggested homology
                        for (Curve * tCurve : mMesh->curves())
                        {
                            bool tAddCurve = false ;
                            for ( Segment * tSegment : tCurve->segments() )
                            {
                                if (tSegment->edge()->is_flagged())
                                {
                                    tAddCurve = true;
                                    break ;
                                }
                            }
                            if (tAddCurve)
                            {
                                for ( Segment * tSegment : tCurve->segments() )
                                {
                                    Edge * tEdge = tSegment->edge() ;
                                    mGenerators(1)(mGenerators(1).size()-1)->addSimplexToChain(tEdge->index(),-1*tSegment->edge_direction()) ;
                                }
                            }
                        }
                        //The suggested homology is simply the boundary of the input terminal
                        mSimplicialComplex = new SimplicialComplex(mMesh);
                        for (auto [tID, tChain] : mSimplicialComplex->get_kchainMap(2))
                        {
                            tTotalChain2->addChainToChain(tChain,1) ;
                        }

                        //If the output terminal is the same as the input terminal, the geometry loops on itself and the inbound/outbound logic is reverted for the output terminal
                        tTotalChain->removeChainFromChain(tTotalChain2);
                        real tCoeff = -1.0 ;
                        if (tTotalChain->getSimplicesMap().size()==0) tCoeff = 1.0 ;

                        //The suggested homology is simply the boundary of the input terminal
                        mGenerators(1)(mGenerators(1).size()-1)->addChainToChain(tTotalChain2->getBoundary(),tCoeff) ;

                        mMesh->unflag_everything() ;
                        delete tTotalChain2 ;
                        delete mSimplicialComplex ;
                    }
                    delete tTotalChain ;
                }

                //For thin shells
                else
                {
                    //Only input terminals in 2-D
                    mGenerators(1).push( new Chain(1,mMesh,true,false) ) ;
                    Chain * tChainInput = mGenerators(1)(mGenerators(1).size()-1) ;
                    for ( uint j = 0 ; j < aTerminals(aThinShellIndices(tCountThinshells)).size()/2 ; ++j )
                    {
                        for ( Segment * tSegment : mMesh->curve(aTerminals(aThinShellIndices(tCountThinshells))(j))->segments() )
                        {
                            Edge * tEdge = tSegment->edge() ;
                            tChainInput->addSimplexToChain(tEdge->index(),-1*tSegment->edge_direction()) ;
                        }
                    }

                    //Also output terminal in 3D
                    if (mMesh->number_of_dimensions() == 3)
                    {
                        //Output terminal
                        mGenerators(1).push( new Chain(1,mMesh,true,false) ) ;
                        Chain * tChainOutput = mGenerators(1)(mGenerators(1).size()-1) ;
                        for ( uint j = aTerminals(aThinShellIndices(tCountThinshells)).size()/2 ; j < aTerminals(aThinShellIndices(tCountThinshells)).size() ; ++j )
                        {
                            for ( Segment * tSegment : mMesh->curve(aTerminals(aThinShellIndices(tCountThinshells))(j))->segments() )
                            {
                                Edge * tEdge = tSegment->edge() ;
                                tChainOutput->addSimplexToChain(tEdge->index(),-1*tSegment->edge_direction()) ;
                            }
                        }
                    }
                    tCountThinshells++;
                }
            }


        }

//-----------------------------------------------------------------------------

        Cell< Matrix< int > > &
        Homology::get_BoundaryMatrix()
        {
            return mD;
        }

//-----------------------------------------------------------------------------

        Cell <Cell< Chain * >> &
        Homology::get_Generators()
        {
            return mGenerators;
        }

//-----------------------------------------------------------------------------

        //Function that create self-intersection for a given generator (only for test purposes)
        void
        Homology::create_self_intersecting(const uint tGeneratorIndex)
        {
            Cell < Edge* > & tEdges = mMesh->edges();
            for (auto & [ tInd, tCoeff2 ]: mGenerators(1)(tGeneratorIndex)->getSimplicesMap())
            {
                tEdges( tInd )->node( 0 )->flag();
                tEdges( tInd )->node( 1 )->flag();
                tEdges(tInd)->flag();
            }

            if( mMesh->number_of_dimensions() == 2 )
            {
                Cell < Element *> & tElements = mMesh->elements();
                for(Element * tElement : tElements)
                {
                    uint tCountEdge = 0;
                    uint tCountNode = 0;
                    for (uint j = 0; j < 3; ++j)
                    {
                        if(tElement->node(j)->is_flagged())
                        {
                            tCountNode++;
                        }
                        if(tElement->edge(j)->is_flagged())
                        {
                            tCountEdge++;
                        }
                    }

                    if (tCountNode == 1 && tCountEdge == 0)
                    {
                        for(uint j = 0 ; j < 3 ; ++j)
                        {
                            mGenerators(1)(tGeneratorIndex)->addSimplexToChain(tElement->edge(j)->index(),tElement->edge_direction(j)?-1:1);
                        }
                        break;
                    }
                }
            }
            else
            {
                Cell < Face *> & tFaces = mMesh->faces();
                for(Face * tFace : tFaces)
                {
                    uint tCountEdge = 0;
                    uint tCountNode = 0;
                    for ( uint j = 0; j < 3; ++j )
                    {
                        if(tFace->node(j)->is_flagged())
                        {
                            tCountNode++;
                        }
                        if(tFace->edge(j)->is_flagged())
                        {
                            tCountEdge++;
                        }
                    }

                    if (tCountNode == 1 && tCountEdge == 0)
                    {
                        for(uint j = 0 ; j < 3 ; ++j)
                        {
                            mGenerators(1)(tGeneratorIndex)->addSimplexToChain(tFace->edge(j)->index(),tFace->edge_direction(j)?-1:1);
                        }
                        break;
                    }
                }
            }

        }

//-----------------------------------------------------------------------------

        //Function that returns the orientation of every homology generator along an input direction for each
        Cell < int >
        Homology::generators_orientation(Cell< Vector <real> > & aDirections)
        {
            Cell < int > tOrientation = Cell < int >(aDirections.size(),0) ;
            //Loop over the given directions for each conductors;
            for(uint i = 0 ; i < aDirections.size(); i++)
            {
                Vector <real> & tDirection = aDirections(i);

                Cell< Edge * > & tMeshEdges = mMesh->edges();

                //Ordering the edges from the loop
                Cell <Edge*> tEdges( mGenerators(1)(i)->getSimplicesMap().size(), nullptr);
                Cell <int> tCoeffs( mGenerators(1)(i)->getSimplicesMap().size(), 0);
                tEdges(0) = tMeshEdges(mGenerators(1)(i)->getSimplicesMap().begin()->first);
                tCoeffs(0) = mGenerators(1)(i)->getSimplicesMap().begin()->second ;
                uint tCount = 1;


                while (tCount < tEdges.size())
                {
                    for(const auto & [tID, tCoeff]: mGenerators(1)(i)->getSimplicesMap())
                    {
                        if (tCoeff > 0 and tCoeffs(tCount-1) > 0)
                        {
                            if (tEdges(tCount-1)->node(1)->index() == tMeshEdges( tID )->node(0)->index())
                            {
                                tEdges(tCount) = tMeshEdges( tID );
                                tCoeffs(tCount++) = tCoeff;
                                break;
                            }
                        }
                        else if (tCoeff < 0 and tCoeffs(tCount-1) > 0)
                        {
                            if (tEdges(tCount-1)->node(1)->index() == tMeshEdges( tID )->node(1)->index())
                            {
                                tEdges(tCount) = tMeshEdges( tID );
                                tCoeffs(tCount++) = tCoeff;
                                break;
                            }
                        }
                        else if (tCoeff > 0 and tCoeffs(tCount-1) < 0)
                        {
                            if (tEdges(tCount-1)->node(0)->index() == tMeshEdges( tID )->node(0)->index())
                            {
                                tEdges(tCount) = tMeshEdges( tID );
                                tCoeffs(tCount++) = tCoeff;
                                break;
                            }
                        }
                        else if (tCoeff < 0 and tCoeffs(tCount-1) < 0)
                        {
                            if (tEdges(tCount-1)->node(0)->index() == tMeshEdges( tID )->node(1)->index())
                            {
                                tEdges(tCount) = tMeshEdges( tID );
                                tCoeffs(tCount++) = tCoeff;
                                break;
                            }
                        }
                    }
                }

                //Loop over the edges to compute each contribution
                real tSum = 0;
                for(uint j = 1; j < tEdges.size(); ++j)
                {
                    Edge * tEdge1 = tEdges(j-1);
                    int tCoeff1 = tCoeffs(j-1);
                    Vector<real> tv1 = Vector<real>(3,0) ;
                    if (tCoeff1 > 0)
                    {
                        tv1(0) = tEdge1->node(1)->x()-tEdge1->node(0)->x();
                        tv1(1) = tEdge1->node(1)->y()-tEdge1->node(0)->y();
                        tv1(2) = tEdge1->node(1)->z()-tEdge1->node(0)->z();
                    }
                    else
                    {
                        tv1(0) = tEdge1->node(0)->x()-tEdge1->node(1)->x();
                        tv1(1) = tEdge1->node(0)->y()-tEdge1->node(1)->y();
                        tv1(2) = tEdge1->node(0)->z()-tEdge1->node(1)->z();
                    }

                    Edge * tEdge2 = tEdges(j);
                    int tCoeff2 = tCoeffs(j);
                    Vector<real> tv2 = Vector<real>(3,0) ;
                    if (tCoeff2 > 0)
                    {
                        tv2(0) = tEdge2->node(1)->x()-tEdge2->node(0)->x();
                        tv2(1) = tEdge2->node(1)->y()-tEdge2->node(0)->y();
                        tv2(2) = tEdge2->node(1)->z()-tEdge2->node(0)->z();
                    }
                    else
                    {
                        tv2(0) = tEdge2->node(0)->x()-tEdge2->node(1)->x();
                        tv2(1) = tEdge2->node(0)->y()-tEdge2->node(1)->y();
                        tv2(2) = tEdge2->node(0)->z()-tEdge2->node(1)->z();
                    }

                    //Triple product
                    tSum += tDirection(0)*(tv1(1)*tv2(2)-tv1(2)*tv2(1))-
                            tDirection(1)*(tv1(0)*tv2(2)-tv1(2)*tv2(0))+
                            tDirection(2)*(tv1(0)*tv2(1)-tv1(1)*tv2(0)) ;


                }

                //Compute the last contribution (last edge and first edge)
                Edge * tEdge1 = tEdges(tEdges.size()-1);
                int tCoeff1 = tCoeffs(tEdges.size()-1);
                Vector<real> tv1 = Vector<real>(3,0) ;
                if (tCoeff1 > 0)
                {
                    tv1(0) = tEdge1->node(1)->x()-tEdge1->node(0)->x();
                    tv1(1) = tEdge1->node(1)->y()-tEdge1->node(0)->y();
                    tv1(2) = tEdge1->node(1)->z()-tEdge1->node(0)->z();
                }
                else
                {
                    tv1(0) = tEdge1->node(0)->x()-tEdge1->node(1)->x();
                    tv1(1) = tEdge1->node(0)->y()-tEdge1->node(1)->y();
                    tv1(2) = tEdge1->node(0)->z()-tEdge1->node(1)->z();
                }

                Edge * tEdge2 = tEdges(0);
                int tCoeff2 = tCoeffs(0);
                Vector<real> tv2 = Vector<real>(3,0) ;
                if (tCoeff2 > 0)
                {
                    tv2(0) = tEdge2->node(1)->x()-tEdge2->node(0)->x();
                    tv2(1) = tEdge2->node(1)->y()-tEdge2->node(0)->y();
                    tv2(2) = tEdge2->node(1)->z()-tEdge2->node(0)->z();
                }
                else
                {
                    tv2(0) = tEdge2->node(0)->x()-tEdge2->node(1)->x();
                    tv2(1) = tEdge2->node(0)->y()-tEdge2->node(1)->y();
                    tv2(2) = tEdge2->node(0)->z()-tEdge2->node(1)->z();
                }

                //Triple product
                tSum += tDirection(0)*(tv1(1)*tv2(2)-tv1(2)*tv2(1))-
                        tDirection(1)*(tv1(0)*tv2(2)-tv1(2)*tv2(0))+
                        tDirection(2)*(tv1(0)*tv2(1)-tv1(1)*tv2(0)) ;

                tOrientation(i) = (tSum>0)?1:-1;
            }
            return tOrientation;
        }

//-----------------------------------------------------------------------------

        //Function that reorient every homology generator along an input direction for each
        void
        Homology::reorient_generators()
        {
            // Fixes the sense of the Amperian loop each suggested generator
            // represents, and with it the sign of the imposed current.
            //
            // Downstream, updatekGeneratorsFromHomology normalises the cut
            // cochain to < c_i, gamma_i > = +1, a duplicate node ties as
            // phi_dup = phi_orig + I ( CutSet::create_duplicates ), and the air
            // uses h = -grad( phi ). A loop crossing the cut from the original
            // to the duplicate side therefore integrates to
            // oint h . dl = -I, so gamma must run CLOCKWISE for a positive
            // declared current to come out along +z by the right hand rule.
            //
            // This used to read +1 in 2-D and -1 in 3-D, and the 2-D half was
            // wrong: measured 2026-08-28 on both a bulk and a thin shell 2-D
            // model, every conductor carried the declared current MAGNITUDE
            // with the sign inverted, uniformly. 3-D measured correct and keeps
            // the -1 it already had, so this edit is a no-op there. The
            // measurements, not the Stokes argument above, are what justify
            // this sign; the run detail is in the 2026-08-28 devlog.
            //
            // The split was never a property of the machinery below; what
            // differs is how suggest_Homology BUILDS the generator, since 3-D
            // gets the difference of the input and output terminal boundaries
            // while 2-D gets the input loop alone, and the two come out with
            // opposite sense. CutProcessor::check_edge carries a further
            // dimension-dependent flip on mIs2D; it is not needed to justify
            // this repair and is deliberately left alone.
            for (uint i = 0 ; i < mGenerators(1).size() ; ++i)
            {
                mGenerators(1)(i)->operator*(-1) ;
            }
        }
    }
}
//------------------------------------------------------------------------------
