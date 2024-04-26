//
// Created by grgia@ge.polymtl.ca on 07/12/23.
//

#include "cl_SimplicialComplex.hpp"

namespace belfem
{
//------------------------------------------------------------------------------
    namespace mesh
    {
        SimplicialComplex::SimplicialComplex( Mesh * aMesh )
        {
            mChains.set_size( 4, Cell< Chain * >());
            mCochains.set_size( 4, Cell< Cochain * >());
            mChainsMap.set_size( 4, Map< id_t, Chain * >());
            mCochainsMap.set_size( 4, Map< id_t, Cochain * >());
            mBoundaryMap.set_size( 4, Map< id_t, Chain * >());
            mCoboundaryMap.set_size( 4, Map< id_t, Cochain * >());
            this->create_complex( aMesh );
        }

        //------------------------------------------------------------------------------

        SimplicialComplex::~SimplicialComplex()
        {
            this->reset();
        }

        //------------------------------------------------------------------------------

        void
        SimplicialComplex::reset()
        {
            // delete map
            mChainsMap( 0 ).clear();
            mBoundaryMap( 0 ).clear();

            mCochainsMap( 0 ).clear();
            mCoboundaryMap( 0 ).clear();

            // delete pointers
            for ( Chain * tChain: mChains( 0 ))
            {
                delete tChain;
            }

            for ( Cochain * tCochain: mCochains( 0 ))
            {
                delete tCochain;
            }

            mChains( 0 ).clear();
            mCochains( 0 ).clear();

            // delete map
            mChainsMap( 1 ).clear();
            mBoundaryMap( 1 ).clear();

            mCochainsMap( 1 ).clear();
            mCoboundaryMap( 1 ).clear();

            // delete pointers
            for ( Chain * tChain: mChains( 1 ))
            {
                delete tChain;
            }

            for ( Cochain * tCochain: mCochains( 1 ))
            {
                delete tCochain;
            }

            mChains( 1 ).clear();
            mCochains( 1 ).clear();

            // delete map
            mChainsMap( 2 ).clear();
            mBoundaryMap( 2 ).clear();

            mCochainsMap( 2 ).clear();
            mCoboundaryMap( 2 ).clear();

            // delete pointers
            for ( Chain * tChain: mChains( 2 ))
            {
                delete tChain;
            }

            for ( Cochain * tCochain: mCochains( 2 ))
            {
                delete tCochain;
            }

            mChains( 2 ).clear();
            mCochains( 2 ).clear();

            // delete map
            mChainsMap( 3 ).clear();
            mBoundaryMap( 3 ).clear();

            mCochainsMap( 3 ).clear();
            mCoboundaryMap( 3 ).clear();

            // delete pointers
            for ( Chain * tChain: mChains( 3 ))
            {
                delete tChain;
            }

            for ( Cochain * tCochain: mCochains( 3 ))
            {
                delete tCochain;
            }

            mChains( 3 ).clear();
            mCochains( 3 ).clear();
        }

        //------------------------------------------------------------------------------

        void
        SimplicialComplex::create_complex( Mesh * aMesh )
        {
            // restore factory settings
            this->reset();

            // initialize counter
            index_t tCount = 0;

            uint tDim = aMesh->number_of_dimensions();
            // allocate memory for 0simplices
            mChains( 0 ).set_size( aMesh->number_of_nodes(), nullptr );
            mCochains( 0 ).set_size( aMesh->number_of_nodes(), nullptr );

            // loop over all nodes on mesh
            for ( Node * tNode: aMesh->nodes())
            {
                if ( tNode->is_flagged())
                {
                    mChains( 0 )( tCount ) = new Chain( 0, aMesh );
                    mChains( 0 )( tCount )->addSimplexToChain( tNode->id(), 1 );

                    mCochains( 0 )( tCount ) = new Cochain( 0, aMesh );
                    mCochains( 0 )( tCount )->addSimplexToCochain( tNode->id(), 1 );

                    // add entry to map
                    mChainsMap( 0 )[ tNode->id() ] = mChains( 0 )( tCount );
                    mCochainsMap( 0 )[ tNode->id() ] = mCochains( 0 )( tCount++ );

                    //populate the 1-boundary map
                    mBoundaryMap( 0 )[ tNode->id() ] = mChainsMap( 0 )[ tNode->id() ]->getBoundary();
                    mCoboundaryMap( 0 )[ tNode->id() ] = mCochainsMap( 0 )[ tNode->id() ]->getCoboundary();
                    tNode->unflag();
                }
            }

            // initialize counter
            tCount = 0;

            // allocate memory for 1simplices
            mChains( 1 ).set_size( aMesh->number_of_edges(), nullptr );
            mCochains( 1 ).set_size( aMesh->number_of_edges(), nullptr );

            // loop over all edges on mesh
            for ( Edge * tEdge: aMesh->edges())
            {
                if ( tEdge->is_flagged())
                {
                    mChains( 1 )( tCount ) = new Chain( 1, aMesh );
                    mChains( 1 )( tCount )->addSimplexToChain( tEdge->id(), 1 );

                    mCochains( 1 )( tCount ) = new Cochain( 1, aMesh );
                    mCochains( 1 )( tCount )->addSimplexToCochain( tEdge->id(), 1 );

                    // add entry to map
                    mChainsMap( 1 )[ tEdge->id() ] = mChains( 1 )( tCount );
                    mCochainsMap( 1 )[ tEdge->id() ] = mCochains( 1 )( tCount++ );

                    //populate the 1-boundary map
                    mBoundaryMap( 1 )[ tEdge->id() ] = mChainsMap( 1 )[ tEdge->id() ]->getBoundary();
                    mCoboundaryMap( 1 )[ tEdge->id() ] = mCochainsMap( 1 )[ tEdge->id() ]->getCoboundary();
                    tEdge->unflag();
                }
            }

            // initialize counter
            tCount = 0;

            if (tDim == 3)
            {
                mChains( 2 ).set_size( aMesh->number_of_faces(), nullptr );
                mCochains( 2 ).set_size( aMesh->number_of_faces(), nullptr );

                mChains( 3 ).set_size( aMesh->number_of_elements(), nullptr );
                mCochains( 3 ).set_size( aMesh->number_of_elements(), nullptr );
                for ( Face * tFace: aMesh->faces())
                {
                    if ( tFace->is_flagged())
                    {
                        mChains( 2 )( tCount ) = new Chain( 2, aMesh );
                        mChains( 2 )( tCount )->addSimplexToChain( tFace->id(), 1 );

                        mCochains( 2 )( tCount ) = new Cochain( 2, aMesh );
                        mCochains( 2 )( tCount )->addSimplexToCochain( tFace->id(), 1 );

                        // add entry to map
                        mChainsMap( 2 )[ tFace->id() ] = mChains( 2 )( tCount );
                        mCochainsMap( 2 )[ tFace->id() ] = mCochains( 2 )( tCount++ );

                        //populate the 2-boundary map
                        mBoundaryMap( 2 )[ tFace->id() ] = mChainsMap( 2)[ tFace->id() ]->getBoundary();
                        mCoboundaryMap( 2 )[ tFace->id() ] = mCochainsMap( 2 )[ tFace->id() ]->getCoboundary();
                        tFace->unflag();
                    }
                }
            }
            else
            {
                mChains( 2 ).set_size( aMesh->number_of_elements(), nullptr );
                mCochains( 2 ).set_size( aMesh->number_of_elements(), nullptr );
            }

            // initialize counter
            tCount = 0;
            // loop over all elements on mesh
            for ( Element * tElement: aMesh->elements())
            {
                if ( tElement->is_flagged())
                {
                    mChains( tDim )( tCount ) = new Chain( tDim, aMesh );
                    mChains( tDim )( tCount )->addSimplexToChain( tElement->id(), 1 );

                    mCochains( tDim )( tCount ) = new Cochain( tDim, aMesh );
                    mCochains( tDim )( tCount )->addSimplexToCochain( tElement->id(), 1 );

                    // add entry to map
                    mChainsMap( tDim )[ tElement->id() ] = mChains( tDim )( tCount );
                    mCochainsMap( tDim )[ tElement->id() ] = mCochains( tDim )( tCount++ );

                    //populate the 2-boundary map
                    mBoundaryMap( tDim )[ tElement->id() ] = mChainsMap( tDim )[ tElement->id() ]->getBoundary();
                    mCoboundaryMap( tDim )[ tElement->id() ] = mCochainsMap( tDim )[ tElement->id() ]->getCoboundary();
                    tElement->unflag();
                }
            }

        }


        //------------------------------------------------------------------------------

        // CCR algorithm as described in Computational Homology from T. Kaczynski et al.
        void
        SimplicialComplex::reduce_complexCCR()
        {
            // loop over all dimensions
            for ( int k = 3; k > 0; --k )
            {
                bool found = true;
                while ( found )
                {
                    found = false;

                    // loop over all k-chains and (k-1)-chains of the complex
                    for ( const auto & [b, tChain1]: mChainsMap( k ))
                    {
                        for ( const auto & [a, tChain2]: mChainsMap( k - 1 ))
                        {
                            // Reduce if a is a boundary of b
                            if ( abs( mBoundaryMap( k )[ b ]->getCoefficient( a )) == 1 )
                            {
                                this->reduce_pair( k, a, b );
                                found = true;
                                break;
                            }
                        }
                        if ( found )
                        {
                            break;
                        }
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------

        // Cohomology variant of the CCR algorithm
        void
        SimplicialComplex::coreduce_complexCCR()
        {
            // loop over all dimensions
            for ( int k = 2; k >= 0; --k )
            {
                bool found = true;
                while ( found )
                {
                    found = false;

                    // loop over all k-cochains and (k+1)-cochains of the complex
                    for ( const auto & [a, tCochain2]: mCochainsMap( k + 1 ))
                    {
                        for ( const auto & [b, tCochain1]: mCochainsMap( k ))
                        {
                            // Reduce if a is a coboundary of b
                            if ( abs( mCoboundaryMap( k )[ b ]->getCoefficient( a )) == 1 )
                            {
                                this->coreduce_pair( k, a, b );
                                found = true;
                                break;
                            }
                        }
                        if ( found )
                        {
                            break;
                        }
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------

        // Function to reduce a pair of chains (a,b), where a is a boundary of b
        void
        SimplicialComplex::reduce_pair( int k, uint a, uint b )
        {
            int val1;
            int val2;

            val1 = mBoundaryMap( k )[ b ]->getCoefficient( a );

            //loop over all the k-chains of the complex
            for ( const auto & [tKey, tChain]: mChainsMap( k ))
            {
                val2 = mBoundaryMap( k )[ tKey ]->getCoefficient( a );
                // Add the b to the k-chain if a is a boundary
                if ( abs( val2 ) == 1 and tKey != b )
                {
                    tChain->addChainToChain( mChainsMap( k )[ b ], -val1 * val2 );
                }
            }

            // Remove a and b from the simplicial chain complex
            this->remove_kchainFromMap( k, b );
            this->remove_kchainFromMap( k - 1, a );
            for ( const auto & [tKey, tChain]: mBoundaryMap( k ))
            {
                mBoundaryMap( k )[ tKey ]->setCoefficient( a, 0 );
            }
        }

        //-----------------------------------------------------------------------------

        // Function to reduce a pair of cochains (a,b), where a is a coboundary of b
        void
        SimplicialComplex::coreduce_pair( int k, uint a, uint b )
        {
            int val1;
            int val2;

            val1 = mCoboundaryMap( k )[ b ]->getCoefficient( a );

            //loop over all the k-cochains of the complex
            for ( const auto & [tKey, tCochain]: mCochainsMap( k ))
            {
                val2 = mCoboundaryMap( k )[ tKey ]->getCoefficient( a );

                // Add the b to the k-cochain if a is a coboundary
                if ( abs( val2 ) == 1 and tKey != b )
                {
                    tCochain->addCochainToCochain( mCochainsMap( k )[ b ], -val1 * val2 );
                }
            }

            // Remove a and b from the simplicial cochain complex
            this->remove_kcochainFromMap( k, b );
            this->remove_kcochainFromMap( k + 1, a );
            for ( const auto & [tKey, tCochain]: mCoboundaryMap( k ))
            {
                mCoboundaryMap( k )[ tKey ]->setCoefficient( a, 0 );
            }
        }

        //-----------------------------------------------------------------------------

        // pReduce from Pellikka et al.
        bool
        SimplicialComplex::pReduce(const uint p)
        {
            bool tRemoved = false;
            if (p == 0)
            {
                return tRemoved;
            }
            uint tCount = 0;
            for ( auto it = mChainsMap( p-1 ).begin(); it != mChainsMap( p-1 ).end() ;)
            {
                tCount = 0;
                uint tpChainReduced;
                uint tp1ChainReduced;
                for (const auto & [tKey, tpChain]: mChainsMap( p ))
                {
                    if(tpChain->getBoundary()->getCoefficient(it->first) != 0)
                    {
                        tCount += 1;
                        tpChainReduced = tKey;
                        tp1ChainReduced = it->first;
                    }
                    if (tCount > 1)
                    {
                        break;
                    }
                }

                if (tCount == 1)
                {
                    --it;
                    this->remove_kchainFromMap( p, tpChainReduced );
                    this->remove_kchainFromMap( p - 1, tp1ChainReduced );
                    tRemoved = true;
                }
                else
                {
                    it++;
                }
            }
            return tRemoved;
        }

        //-----------------------------------------------------------------------------

        // pCombine from Pellikka et al.
        void
        SimplicialComplex::pCombine(const uint p)
        {
            int val1;
            int val2;
            uint tCount;
            id_t tpChainAdd;
            id_t tpChainReduced;
            id_t tCp_1;
            Cell< id_t > Q;

            for (const auto & [tKey, tCp]: mChainsMap( p ))
            {

                if (p == 0)
                {
                    return;
                }
                for (const auto & [tId, tCoeff]: tCp->getBoundary()->getSimplicesMap())
                {
                    Q.push(tId);
                }

                while (Q.size() > 0)
                {
                    tCp_1 = Q.pop();
                    tCount = 0;
                    for (const auto & [t1Key, tpChain]: mChainsMap( p ))
                    {
                        if ( tpChain->getBoundary()->getCoefficient( tCp_1 ) != 0 )
                        {
                            tCount += 1;
                            if ( tCount == 1 )
                            {
                                tpChainAdd = t1Key;
                            }
                            else if ( tCount == 2 )
                            {
                                tpChainReduced = t1Key;
                            }

                        }
                        if ( tCount > 2 )
                        {
                            break;
                        }
                    }
                    if ( tCount == 2)
                    {
                        val1 = mBoundaryMap( p )[ tpChainReduced ]->getCoefficient( tCp_1 );
                        val2 = mBoundaryMap( p )[ tpChainAdd ]->getCoefficient( tCp_1 );
                        mChainsMap( p )[tpChainAdd]->addChainToChain( mChainsMap( p )[tpChainReduced], -val1*val2);
                        this->remove_kchainFromMap( p, tpChainReduced );
                        this->remove_kchainFromMap( p - 1, tCp_1);
                        mBoundaryMap( p )[ tpChainAdd ]->setCoefficient( tCp_1, 0 );
                        for (const auto & [tId, tCoeff]: mBoundaryMap( p )[ tpChainAdd ]->getSimplicesMap())
                        {
                            Q.push(tId);
                        }
                        unique(Q);

                    }
                }
            }
        }

        //-----------------------------------------------------------------------------

        // ReduceOmit from Pellikka et al.
        void
        SimplicialComplex::reduceOmit()
        {
            for (uint p = 2; p > 1; p--)
            {
                this->pReduce(p);
            }

            /*uint tKey;
            while (this->number_of_ksimplices(2) > 0)
            {
                tKey = mChainsMap(2).begin()->first;
                this->remove_kchainFromMap(2,tKey);
                for (uint p = 2; p > 1; p--)
                {
                    this->pReduce(p);
                }
            }*/
        }

        //-----------------------------------------------------------------------------

        void
        SimplicialComplex::reduce_complexPellikka()
        {
            this->reduceOmit();
            for (uint p = 2; p >= 1; p--)
            {
                this->pCombine(p);
                this->pReduce(p-1);
            }
        }

        //-----------------------------------------------------------------------------


        // pReduce from Pellikka et al.
        bool
        SimplicialComplex::pCoreduce(const uint p)
        {
            bool tRemoved = false;
            if (p == 3)
            {
                return tRemoved;
            }
            uint tCount = 0;
            for ( auto it = mCochainsMap( p+1 ).begin(); it != mCochainsMap( p+1 ).end() ;)
            {
                tCount = 0;
                uint tpChainReduced;
                uint tp1ChainReduced;
                for (const auto & [tKey, tpChain]: mCochainsMap( p ))
                {
                    if(tpChain->getCoboundary()->getCoefficient(it->first) != 0)
                    {
                        tCount += 1;
                        tpChainReduced = tKey;
                        tp1ChainReduced = it->first;
                    }
                    if (tCount > 1)
                    {
                        break;
                    }
                }

                if (tCount == 1)
                {
                    --it;
                    this->remove_kcochainFromMap( p + 1, tp1ChainReduced );
                    this->remove_kcochainFromMap( p, tpChainReduced );
                    tRemoved = true;
                }
                else
                {
                    it++;
                }
            }
            return tRemoved;
        }

        //-----------------------------------------------------------------------------

        // pCombine from Pellikka et al.
        void
        SimplicialComplex::pCocombine(const uint p)
        {
            int val1;
            int val2;
            uint tCount;
            id_t tpChainAdd;
            id_t tpChainReduced;
            id_t tCp_1;
            Cell< id_t > Q;

            for (const auto & [tKey, tCp]: mCochainsMap( p ))
            {

                if (p == 3)
                {
                    return;
                }
                for (const auto & [tId, tCoeff]: tCp->getCoboundary()->getSimplicesMap())
                {
                    Q.push(tId);
                }

                while (Q.size() > 0)
                {
                    tCp_1 = Q.pop();
                    tCount = 0;
                    for (const auto & [t1Key, tpChain]: mCochainsMap( p ))
                    {
                        if ( tpChain->getCoboundary()->getCoefficient( tCp_1 ) != 0 )
                        {
                            tCount += 1;
                            if ( tCount == 1 )
                            {
                                tpChainAdd = t1Key;
                            }
                            else if ( tCount == 2 )
                            {
                                tpChainReduced = t1Key;
                            }

                        }
                        if ( tCount > 2 )
                        {
                            break;
                        }
                    }
                    if ( tCount == 2)
                    {
                        val1 = mCoboundaryMap( p )[ tpChainReduced ]->getCoefficient( tCp_1 );
                        val2 = mCoboundaryMap( p )[ tpChainAdd ]->getCoefficient( tCp_1 );
                        mCochainsMap( p )[tpChainAdd]->addCochainToCochain( mCochainsMap( p )[tpChainReduced], -val1*val2);
                        this->remove_kcochainFromMap( p, tpChainReduced );
                        this->remove_kcochainFromMap( p + 1, tCp_1);
                        mCoboundaryMap( p )[ tpChainAdd ]->setCoefficient( tCp_1, 0 );
                        for (const auto & [tId, tCoeff]: mCoboundaryMap( p )[ tpChainAdd ]->getSimplicesMap())
                        {
                            Q.push(tId);
                        }
                        unique(Q);
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------

        // ReduceOmit from Pellikka et al.
        void
        SimplicialComplex::coreduceOmit()
        {
            id_t tKey = mCochainsMap(0).begin()->first;
            this->remove_kcochainFromMap(0,tKey);
            for (uint p = 0; p < 2; p++)
            {
                this->pCoreduce(p);
            }

            /*uint tKey;
            while (this->number_of_ksimplices(2) > 0)
            {
                tKey = mChainsMap(2).begin()->first;
                this->remove_kchainFromMap(2,tKey);
                for (uint p = 2; p > 1; p--)
                {
                    this->pReduce(p);
                }
            }*/
        }

        //-----------------------------------------------------------------------------

        void
        SimplicialComplex::coreduce_complexPellikka()
        {
            this->coreduceOmit();
            for (uint p = 0; p <= 1; p++)
            {
                this->pCocombine(p);
                this->pCoreduce(p+1);
            }
        }

        //-----------------------------------------------------------------------------

        void
        SimplicialComplex::remove_kchainFromMap( const uint k, const id_t aInd )
        {
            if ( k <= 3 )
            {
                mChainsMap( k ).erase_key( aInd );
                mBoundaryMap( k ).erase_key( aInd );
            }
        }

        //------------------------------------------------------------------------------

        void
        SimplicialComplex::remove_kcochainFromMap( const uint k, const id_t aInd )
        {
            if ( k <= 3 )
            {
                mCochainsMap( k ).erase_key( aInd );
                mCoboundaryMap( k ).erase_key( aInd );
            }
        }

        //------------------------------------------------------------------------------

        Cell< Map< id_t, Chain * > >
        SimplicialComplex::get_boundaryMap()
        {
            return mBoundaryMap;
        }

        //------------------------------------------------------------------------------

        Cell< Map< id_t, Cochain * > >
        SimplicialComplex::get_coboundaryMap()
        {
            return mCoboundaryMap;
        }

        //------------------------------------------------------------------------------

        Chain *
        SimplicialComplex::get_kchain( const uint k, const id_t aInd )
        {
            if ( k <= 3 )
            {
                return mChainsMap( k )[ aInd ];
            }
            else
            {
                return nullptr;
            }
        }

        //------------------------------------------------------------------------------

        Map< id_t, Chain * >
        SimplicialComplex::get_kchainMap( const uint k )
        {
            return mChainsMap( k );
        }

        //------------------------------------------------------------------------------

        Map< id_t, Cochain * >
        SimplicialComplex::get_kcochainMap( const uint k )
        {
            return mCochainsMap( k );
        }

        //------------------------------------------------------------------------------

        Chain *
        SimplicialComplex::boundary_of_kchain( const uint k, const id_t aInd )
        {
            if ( k <= 3 )
            {
                return mBoundaryMap( k )[ aInd ];
            }
            else
            {
                return nullptr;
            }
        }

        //------------------------------------------------------------------------------

        // Function that creates the boundary matrices (with 1, -1 and 0)
        // from the simplicial complex (reduced or not)
        Cell< Matrix< int > >
        SimplicialComplex::createMatrixFromBoundaryMap()
        {
            Cell< Matrix< int > > tBoundaryMat;
            tBoundaryMat.set_size( 4, Matrix< int >());
            tBoundaryMat( 0 ).set_size( 1, this->number_of_0simplices(), 0 );

            //loop over all dimensions
            for ( int k = 1; k < 4; k++ )
            {
                tBoundaryMat( k ).set_size( this->number_of_ksimplices( k - 1 ), this->number_of_ksimplices( k ), 0 );
                uint tCount = 0;

                //loop over all the k-chains and (k-1)-chains to populate the matrix
                for ( const auto & [tKey, tChain]: mChainsMap( k ))
                {
                    uint tCount2 = 0;
                    for ( const auto & [tKey2, tChain2]: mChainsMap( k - 1 ))
                    {
                        tBoundaryMat( k )( tCount2, tCount ) = mBoundaryMap( k )[ tKey ]->getCoefficient( tKey2 );
                        tCount2++;
                    }
                    tCount++;
                }
            }
            return tBoundaryMat;
        }

        //------------------------------------------------------------------------------

        // Function that creates the coboundary matrices (with 1, -1 and 0)
        // from the simplicial complex (coreduced or not)
        Cell< Matrix< int > >
        SimplicialComplex::createMatrixFromCoboundaryMap()
        {
            Cell< Matrix< int > > tCoboundaryMat;
            tCoboundaryMat.set_size( 4, Matrix< int >());

            int n = this->number_of_kcosimplices( 3 );
            tCoboundaryMat( 3 ).set_size( 1, std::max( n, 1 ), 0 );

            //loop over all dimensions
            for ( int k = 0; k < 3; k++ )
            {
                int m = this->number_of_kcosimplices( k + 1 );
                n = this->number_of_kcosimplices( k );

                tCoboundaryMat( k ).set_size( std::max( m, 1 ), std::max( n, 1 ), 0 );
                uint tCount = 0;

                //loop over all the k-cochains and (k+1)-cochains to populate the matrix
                for ( const auto & [tKey, tChain]: mCochainsMap( k ))
                {
                    uint tCount2 = 0;
                    for ( const auto & [tKey2, tChain2]: mCochainsMap( k + 1 ))
                    {
                        tCoboundaryMat( k )( tCount2, tCount ) = mCoboundaryMap( k )[ tKey ]->getCoefficient( tKey2 );
                        tCount2++;
                    }
                    tCount++;
                }
            }
            return tCoboundaryMat;
        }

        //------------------------------------------------------------------------------

        void
        SimplicialComplex::print_kchains( const uint k )
        {
            for ( const auto & [tKey, tChain]: mChainsMap( k ))
            {
                std::cout << k << "-chain #" << tKey << " :" << std::endl;
                tChain->print();
            }
        }

        //------------------------------------------------------------------------------

        void
        SimplicialComplex::create_kComplexField( const uint k, Mesh * aMesh )
        {
            char fieldName[50];
            sprintf( fieldName, "1ChainComplex" );
            aMesh->create_field( fieldName, EntityType::NODE );
            for ( const auto & [tKey, tChain]: this->get_kchainMap( k ))
            {
                Map< id_t, int > tSimplicesMap = tChain->getSimplicesMap();
                for ( const auto & [tInd, tCoeff]: tSimplicesMap )
                {
                    aMesh->field_data( fieldName )( aMesh->edge( tInd )->node( 0 )->index()) = 1;
                    aMesh->field_data( fieldName )( aMesh->edge( tInd )->node( 1 )->index()) = 1;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        SimplicialComplex::cocreate_kComplexField( const uint k, Mesh * aMesh )
        {
            for ( const auto & [tKey, tChain]: this->get_kcochainMap( k ))
            {
                char fieldName[50];
                sprintf( fieldName, "1CochainComplex%d", tKey );
                aMesh->create_field( fieldName, EntityType::NODE );
                Map< id_t, int > tSimplicesMap = tChain->getSimplicesMap();
                for ( const auto & [tInd, tCoeff]: tSimplicesMap )
                {
                    aMesh->field_data( fieldName )( aMesh->edge( tInd )->node( 0 )->index()) = 1;
                    aMesh->field_data( fieldName )( aMesh->edge( tInd )->node( 1 )->index()) = 1;
                }
            }
        }

        //------------------------------------------------------------------------------

    }
}