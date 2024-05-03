//
// Created by grgia@ge.polymtl.ca on 07/12/23.
//

#include "cl_SimplicialComplex.hpp"

#include "cl_Timer.hpp"

namespace belfem
{
//------------------------------------------------------------------------------
    namespace mesh
    {
        SimplicialComplex::SimplicialComplex( Mesh * aMesh )
        {
            mChainsMap.set_size( 4, Map< id_t, Chain * >());
            mCochainsMap.set_size( 4, Map< id_t, Cochain * >());
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
            // delete pointers
            for ( const auto & [tKey, tChain]: mChainsMap( 0 ))
            {
                delete tChain;
            }

            // delete map
            mChainsMap( 0 ).clear();
            mCochainsMap( 0 ).clear();

            // delete pointers
            for ( const auto & [tKey, tChain]: mChainsMap( 1 ))
            {
                delete tChain;
            }

            // delete map
            mChainsMap( 1 ).clear();
            mCochainsMap( 1 ).clear();

            // delete pointers
            for ( const auto & [tKey, tChain]: mChainsMap( 2 ))
            {
                delete tChain;
            }

            // delete map
            mChainsMap( 2 ).clear();
            mCochainsMap( 2 ).clear();

            // delete pointers
            for ( const auto & [tKey, tChain]: mChainsMap( 3 ))
            {
                delete tChain;
            }

            // delete map
            mChainsMap( 3 ).clear();
            mCochainsMap( 3 ).clear();

        }

        //------------------------------------------------------------------------------

        void
        SimplicialComplex::create_complex( Mesh * aMesh )
        {
            // restore factory settings
            this->reset();

            uint tDim = aMesh->number_of_dimensions();

            // loop over all nodes on mesh
            for ( Node * tNode: aMesh->nodes())
            {
                if ( tNode->is_flagged())
                {
                    Chain* tChain = new Chain( 0, aMesh );
                    tChain->addSimplexToChain( tNode->id(), 1 );

                    Cochain* tCochain = new Cochain( 0, aMesh );
                    tCochain->addSimplexToCochain( tNode->id(), 1 );

                    // add entry to map
                    mChainsMap( 0 )[ tNode->id() ] = tChain;
                    mCochainsMap( 0 )[ tNode->id() ] = tCochain;

                    tNode->unflag();
                }
            }


            // loop over all edges on mesh
            for ( Edge * tEdge: aMesh->edges())
            {
                if ( tEdge->is_flagged())
                {
                    Chain* tChain = new Chain( 1, aMesh );
                    tChain->addSimplexToChain( tEdge->id(), 1 );

                    Cochain* tCochain = new Cochain( 1, aMesh );
                    tCochain->addSimplexToCochain( tEdge->id(), 1 );

                    // add entry to map
                    mChainsMap( 1 )[ tEdge->id() ] = tChain;
                    mCochainsMap( 1 )[ tEdge->id() ] = tCochain;

                    tEdge->unflag();
                }
            }

            if (tDim == 3)
            {
                for ( Face * tFace: aMesh->faces())
                {
                    if ( tFace->is_flagged())
                    {
                        Chain* tChain = new Chain( 2, aMesh );
                        tChain->addSimplexToChain( tFace->id(), 1 );

                        Cochain* tCochain = new Cochain( 2, aMesh );
                        tCochain->addSimplexToCochain( tFace->id(), 1 );

                        // add entry to map
                        mChainsMap( 2 )[ tFace->id() ] = tChain;
                        mCochainsMap( 2 )[ tFace->id() ] = tCochain;

                        tFace->unflag();
                    }
                }
            }

            // loop over all elements on mesh
            for ( Element * tElement: aMesh->elements())
            {
                if ( tElement->is_flagged())
                {
                    Chain* tChain = new Chain( tDim, aMesh );
                    tChain->addSimplexToChain( tElement->id(), 1 );

                    Cochain* tCochain = new Cochain( tDim, aMesh );
                    tCochain->addSimplexToCochain( tElement->id(), 1 );

                    // add entry to map
                    mChainsMap( tDim )[ tElement->id() ] = tChain;
                    mCochainsMap( tDim )[ tElement->id() ] = tCochain;

                    tElement->unflag();
                }
            }

        }


        //------------------------------------------------------------------------------

        // CCR algorithm as described in Computational Homology from T. Kaczynski et al.
        void
        SimplicialComplex::reduce_complexCCR()
        {
            Timer tTimer;
            std::cout << "Reducing the complex (CCR) ..." << std::endl;

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
                            if ( abs( tChain1->getBoundary()->getCoefficient( a )) == 1 )
                            {
                                this->reduce_pair( k, a, b, tChain1 );
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

            std::cout << "Complex reduced in: "<< tTimer.stop()*1e-3 << " s" << std::endl;

            std::cout << "Number of 0-chains: " << this->number_of_ksimplices(0) <<std::endl;
            std::cout << "Number of 1-chains: " << this->number_of_ksimplices(1) <<std::endl;
            std::cout << "Number of 2-chains: " << this->number_of_ksimplices(2) <<std::endl;
            std::cout << "Number of 3-chains: " << this->number_of_ksimplices(3) <<std::endl;
            std::cout << std::endl;
        }

        //-----------------------------------------------------------------------------

        // Cohomology variant of the CCR algorithm
        void
        SimplicialComplex::coreduce_complexCCR()
        {
            Timer tTimer;
            std::cout << "Coreducing the complex (CCR) ..." << std::endl;

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
                        for ( const auto & [b, tCochain1]: mCochainsMap( k  ))
                        {
                            // Reduce if a is a coboundary of b
                            if ( abs( tCochain1->getCoboundary()->getCoefficient( a )) == 1 )
                            {
                                this->coreduce_pair( k, a, b, tCochain1 );
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

            std::cout << "Complex coreduced in: "<< tTimer.stop()*1e-3 << " s" << std::endl;

            std::cout << "Number of 0-cochains: " << this->number_of_kcosimplices(0) <<std::endl;
            std::cout << "Number of 1-cochains: " << this->number_of_kcosimplices(1) <<std::endl;
            std::cout << "Number of 2-cochains: " << this->number_of_kcosimplices(2) <<std::endl;
            std::cout << "Number of 3-cochains: " << this->number_of_kcosimplices(3) <<std::endl;
            std::cout << std::endl;
        }

        //-----------------------------------------------------------------------------

        // Function to reduce a pair of chains (a,b), where a is a boundary of b
        void
        SimplicialComplex::reduce_pair( const uint k, const uint a, const uint b, Chain* aChain )
        {
            int val1;
            int val2;

            val1 = aChain->getBoundary()->getCoefficient( a );

            //loop over all the k-chains of the complex
            for ( const auto & [tKey, tChain]: mChainsMap( k ))
            {
                val2 = tChain->getBoundary()->getCoefficient( a );
                // Add the b to the k-chain if a is a boundary
                if ( abs( val2 ) == 1 and tKey != b )
                {
                    tChain->addChainToChain( aChain, -val1 * val2 );
                }
            }

            // Remove a and b from the simplicial chain complex
            this->remove_kchainFromMap( k, b );
            this->remove_kchainFromMap( k - 1, a );
            /*for ( const auto & [tKey, tChain]: mChainsMap( k ))
            {
                tChain->getBoundary()->setCoefficient( a, 0 );
            }*/
        }

        //-----------------------------------------------------------------------------

        // Function to reduce a pair of cochains (a,b), where a is a coboundary of b
        void
        SimplicialComplex::coreduce_pair( const uint k, const uint a, const uint b, Cochain* aCochain )
        {
            int val1;
            int val2;

            val1 = aCochain->getCoboundary()->getCoefficient( a );

            //loop over all the k-cochains of the complex
            for ( const auto & [tKey, tCochain]: mCochainsMap( k ))
            {
                val2 = tCochain->getCoboundary()->getCoefficient( a );

                // Add the b to the k-cochain if a is a coboundary
                if ( abs( val2 ) == 1 and tKey != b )
                {
                    tCochain->addCochainToCochain( aCochain, -val1 * val2 );
                }
            }

            // Remove a and b from the simplicial cochain complex
            this->remove_kcochainFromMap( k, b );
            this->remove_kcochainFromMap( k + 1, a );
            /*for ( const auto & [tKey, tCochain]: mCochainsMap( k ))
            {
                tCochain->getCoboundary()->setCoefficient( a, 0 );
            }*/
        }

        //-----------------------------------------------------------------------------

        // pReduce from Pellikka et al.
        void
        SimplicialComplex::pReduce(const uint p)
        {
            if (p == 0)
            {
                return ;
            }
            uint tCount = 0;
            bool tRemoved = true;
            while (tRemoved)
            {
                tRemoved = false;
                for (auto it = mChainsMap( p-1 ).begin(), next_it = it;
                    it != mChainsMap( p-1 ).end(); it = next_it)
                {
                    ++next_it;
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
                        this->remove_kchainFromMap( p - 1, tp1ChainReduced );
                        this->remove_kchainFromMap( p, tpChainReduced );
                        tRemoved = true;
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------

        // pCombine from Pellikka et al.
        void
        SimplicialComplex::pCombine(const uint p)
        {
            int val1;
            int val2;
            uint tCount;
            Chain* tpChainAdd;
            Chain* tpChainReduced;
            id_t tKeyReduced;
            id_t tCp_1;
            Cell< id_t > Q;

            for (auto it = mChainsMap( p ).begin();  it != mChainsMap( p ).end(); ++it )
            {

                if (p == 0)
                {
                    return;
                }
                for (const auto & [tId, tCoeff]: it->second->getBoundary()->getSimplicesMap())
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
                                tpChainAdd = tpChain;
                                val1 = tpChain->getBoundary()->getCoefficient( tCp_1 );
                            }
                            else if ( tCount == 2 )
                            {
                                tpChainReduced = tpChain;
                                val2 = tpChain->getBoundary()->getCoefficient( tCp_1 );
                                tKeyReduced = t1Key;
                            }

                        }
                        if ( tCount > 2 )
                        {
                            break;
                        }
                    }
                    if ( tCount == 2)
                    {
                        tpChainAdd->addChainToChain( tpChainReduced, -val1*val2);
                        this->remove_kchainFromMap( p, tKeyReduced );
                        this->remove_kchainFromMap( p - 1, tCp_1);
                        //tpChainAdd->getBoundary()->setCoefficient( tCp_1, 0 );
                        for (const auto & [tId, tCoeff]: tpChainAdd->getBoundary()->getSimplicesMap())
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
            for (uint p = 3; p > 1; p--)
            {
                this->pReduce(p);
            }

            uint tKey;
            while (this->number_of_ksimplices(2) > 0)
            {
                tKey = mChainsMap(2).begin()->first;
                this->remove_kchainFromMap(2,tKey);
                for (uint p = 3; p >= 1; p--)
                {
                    this->pReduce(p);
                }
            }
        }

        //-----------------------------------------------------------------------------

        void
        SimplicialComplex::reduce_complexPellikka()
        {
            Timer tTimer;
            std::cout << "Reducing the complex (Pellikka's algorithm) ..." << std::endl;

            this->reduceOmit();
            for (uint p = 3; p >= 1; p--)
            {
                this->pCombine(p);
                this->pReduce(p-1);
            }

            std::cout << "Complex reduced in: "<< tTimer.stop()*1e-3 << " s" << std::endl;

            std::cout << "Number of 0-chains: " << this->number_of_ksimplices(0) <<std::endl;
            std::cout << "Number of 1-chains: " << this->number_of_ksimplices(1) <<std::endl;
            std::cout << "Number of 2-chains: " << this->number_of_ksimplices(2) <<std::endl;
            std::cout << "Number of 3-chains: " << this->number_of_ksimplices(3) <<std::endl;
            std::cout << std::endl;
        }

//-----------------------------------------------------------------------------

        // pReduce from Pellikka et al.
        void
        SimplicialComplex::pCoreduce(const uint p)
        {
            if (p == 3)
            {
                return ;
            }
            uint tCount = 0;
            bool tRemoved = true;
            while (tRemoved)
            {
                tRemoved = false;
                for (auto it = mCochainsMap( p+1 ).begin(), next_it = it;
                    it != mCochainsMap( p+1 ).end(); it = next_it)
                {
                    ++next_it;
                    tCount = 0;
                    uint tpChainReduced;
                    uint tp1ChainReduced;
                    for (const auto & [tKey, tpChain]: mCochainsMap( p ))
                    {
                        if(tpChain->getCoboundary()->getCoefficient(it->first) != 0 )
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
                        this->remove_kcochainFromMap( p + 1, tp1ChainReduced );
                        this->remove_kcochainFromMap( p, tpChainReduced );
                        tRemoved = true;
                    }
                }
            }
        }

        //-----------------------------------------------------------------------------

        // pCombine from Pellikka et al.
        void
        SimplicialComplex::pCocombine(const uint p)
        {
            int val1;
            int val2;
            uint tCount;
            Cochain* tpChainAdd;
            Cochain* tpChainReduced;
            id_t tKeyReduced;
            id_t tCp_1;
            Cell< id_t > Q;

            for (auto it = mCochainsMap( p ).begin();  it != mCochainsMap( p ).end(); ++it )
            {

                if (p == 3)
                {
                    return;
                }
                for (const auto & [tId, tCoeff]: it->second->getCoboundary()->getSimplicesMap())
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
                                tpChainAdd = tpChain;
                                val1 = tpChain->getCoboundary()->getCoefficient( tCp_1 );
                            }
                            else if ( tCount == 2 )
                            {
                                tpChainReduced = tpChain;
                                val2 = tpChain->getCoboundary()->getCoefficient( tCp_1 );
                                tKeyReduced = t1Key;
                            }

                        }
                        if ( tCount > 2 )
                        {
                            break;
                        }
                    }
                    if ( tCount == 2)
                    {
                        tpChainAdd->addCochainToCochain( tpChainReduced, -val1*val2);
                        this->remove_kcochainFromMap( p, tKeyReduced );
                        this->remove_kcochainFromMap( p + 1, tCp_1);
                        //tpChainAdd->getCoboundary()->setCoefficient( tCp_1, 0 );
                        for (const auto & [tId, tCoeff]: tpChainAdd->getCoboundary()->getSimplicesMap())
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
            /*auto it = mChainsMap(0).begin();
            this->remove_kcochainFromMap(0,it->first);
            for (uint p = 0; p <= 2; p++)
            {
                this->pCoreduce(p);
            }*/

            uint tKey;
            while (this->number_of_kcosimplices(0) > 0)
            {
                tKey = mCochainsMap(0).begin()->first;
                this->remove_kcochainFromMap(0,tKey);
                for (uint p = 0; p <= 2; p++)
                {
                    this->pCoreduce(p);
                }
            }
        }

        //-----------------------------------------------------------------------------

        void
        SimplicialComplex::coreduce_complexPellikka()
        {
            Timer tTimer;
            std::cout << "Coreducing the complex (Pellikka's algorithm) ..." << std::endl;

            this->coreduceOmit();
            for (uint p = 0; p <= 2; p++)
            {
                this->pCocombine(p);
                this->pCoreduce(p+1);
            }

            std::cout << "Complex coreduced in: "<< tTimer.stop()*1e-3 << " s" << std::endl;

            std::cout << "Number of 0-cochains: " << this->number_of_kcosimplices(0) <<std::endl;
            std::cout << "Number of 1-ccohains: " << this->number_of_kcosimplices(1) <<std::endl;
            std::cout << "Number of 2-cochains: " << this->number_of_kcosimplices(2) <<std::endl;
            std::cout << "Number of 3-cochains: " << this->number_of_kcosimplices(3) <<std::endl;
            std::cout << std::endl;
        }

        //-----------------------------------------------------------------------------

        void
        SimplicialComplex::remove_kchainFromMap( const uint k, const id_t aID )
        {
            if ( k <= 3 )
            {
                delete mChainsMap( k )[aID];
                mChainsMap( k ).erase_key( aID );
            }
        }

        //------------------------------------------------------------------------------

        void
        SimplicialComplex::remove_kcochainFromMap( const uint k, const id_t aID )
        {
            if ( k <= 3 )
            {
                delete mCochainsMap( k )[aID];
                mCochainsMap( k ).erase_key( aID );
            }
        }

        //------------------------------------------------------------------------------

        Chain *
        SimplicialComplex::get_kchain( const uint k, const id_t aID )
        {
            if ( k <= 3 )
            {
                return mChainsMap( k )[ aID ];
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
        SimplicialComplex::boundary_of_kchain( const uint k, const id_t aID )
        {
            if ( k <= 3 )
            {
                return mChainsMap( k )[ aID ]->getBoundary();
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
                        tBoundaryMat( k )( tCount2, tCount ) = tChain->getBoundary()->getCoefficient( tKey2 );
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
                        tCoboundaryMat( k )( tCount2, tCount ) = tChain->getCoboundary()->getCoefficient( tKey2 );
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
                for ( const auto & [tID, tCoeff]: tSimplicesMap )
                {
                    aMesh->field_data( fieldName )( aMesh->edge( tID )->node( 0 )->index()) = 1;
                    aMesh->field_data( fieldName )( aMesh->edge( tID )->node( 1 )->index()) = 1;
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
                for ( const auto & [tID, tCoeff]: tSimplicesMap )
                {
                    aMesh->field_data( fieldName )( aMesh->edge( tID )->node( 0 )->index()) = 1;
                    aMesh->field_data( fieldName )( aMesh->edge( tID )->node( 1 )->index()) = 1;
                }
            }
        }

        //------------------------------------------------------------------------------

    }
}