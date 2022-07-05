//
// Created by Christian Messe on 28.04.20.
//
#include "assert.hpp"
#include "cl_BS_TMatrix.hpp"

namespace belfem
{
    namespace bspline
    {
//------------------------------------------------------------------------------

        TMatrix::TMatrix( const uint aNumberOfDimensions,
                const uint aBsplineOrder,
                const uint aLangrangeOrder ) :
            mNumberOfDimensions( aNumberOfDimensions ),
            mBsplineOrder( aBsplineOrder ),
            mLagrangeOrder( aLangrangeOrder == 0 ? aBsplineOrder : aLangrangeOrder ),
            mNumberOfBasis( std::pow( aBsplineOrder + 1, aNumberOfDimensions ) ),
            mNumberOfNodes( aLangrangeOrder == 0 ?
                        std::pow( aBsplineOrder + 1, aNumberOfDimensions ) :
                        std::pow( aLangrangeOrder + 1, aNumberOfDimensions ) ),
            mBsplineType( this->select_type( aNumberOfDimensions, aBsplineOrder ) ),
            mLagrangeType(  aLangrangeOrder == 0 ?
                            this->select_type( aNumberOfDimensions, aBsplineOrder ) :
                            this->select_type( aNumberOfDimensions, aLangrangeOrder ) )
        {

            // write indices for lagrange and B-Splines
            mBasisIndex.set_size( mNumberOfDimensions, mNumberOfBasis );
            mNodeIndex.set_size( mNumberOfDimensions, mNumberOfNodes );

            // populate indices
            this->set_index( mBsplineType, mBasisIndex );
            this->set_index( mLagrangeType, mNodeIndex );

            // parameter coordinates for nodes
            this->populate_node_param_coords();

            this->compute_lagrange_matrix();

        }

//------------------------------------------------------------------------------

        ElementType
        TMatrix::select_type( const uint aNumDimensions, const uint aOrder )
        {
            switch( aNumDimensions )
            {
                case( 1 ):
                {
                    switch( aOrder )
                    {
                        case( 1 ):
                        {
                            return ElementType::LINE2;
                        }
                        case( 2 ):
                        {
                            return ElementType::LINE3;
                        }
                        case( 3 ):
                        {
                            return ElementType::LINE4;
                        }
                        default:
                        {
                            BELFEM_ERROR( false,
                                    "aOrder = %u is invalid",
                                     ( unsigned int ) aOrder );
                            return ElementType::UNDEFINED ;
                        }
                    }
                }
                    case( 2 ):
                    {
                        switch( aOrder )
                        {
                            case( 1 ):
                            {
                                return ElementType::QUAD4;
                            }
                            case( 2 ):
                            {
                                return ElementType::QUAD9;
                            }
                            case( 3 ):
                            {
                                return ElementType::QUAD16;
                            }
                            default:
                            {
                                BELFEM_ERROR( false,
                                     "aOrder = %u is invalid",
                                     ( unsigned int ) aOrder );
                                return ElementType::UNDEFINED ;
                            }
                        }
                    }
                case( 3 ):
                {
                    switch( aOrder )
                    {
                        case( 1 ):
                        {
                            return ElementType::HEX8;
                        }
                        case( 2 ):
                        {
                            return ElementType::HEX27;
                        }
                        case( 3 ):
                        {
                            return ElementType::HEX64;
                        }
                        default:
                        {
                            BELFEM_ERROR( false,
                                         "aOrder = %u is invalid",
                                         ( unsigned int ) aOrder );
                            return ElementType::UNDEFINED ;
                        }
                    }
                }
                default:
                {
                    BELFEM_ERROR( false, "aNumDimensions = %u is invalid",
                                 ( unsigned int ) aNumDimensions );
                    return ElementType::UNDEFINED ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        TMatrix::set_index( const ElementType aType, Matrix< index_t > & aIndex )
        {

            switch( aType )
            {
                case( ElementType::LINE2 ) :
                {
                    aIndex( 0, 0 ) = 0;
                    aIndex( 0, 1 ) = 1;
                    break;
                }
                case( ElementType::LINE3 ) :
                {
                    aIndex( 0, 0 ) = 0;
                    aIndex( 0, 1 ) = 2;
                    aIndex( 0, 2 ) = 1;
                    break;
                }
                case( ElementType::LINE4 ) :
                {
                    aIndex( 0, 0 ) = 0;
                    aIndex( 0, 1 ) = 3;
                    aIndex( 0, 2 ) = 1;
                    aIndex( 0, 3 ) = 2;
                    break;
                }
                case( ElementType::QUAD4 ) :
                {
                    aIndex( 0, 0 ) = 0 ;
                    aIndex( 1, 0 ) = 0 ;
                    aIndex( 0, 1 ) = 1 ;
                    aIndex( 1, 1 ) = 0 ;
                    aIndex( 0, 2 ) = 1 ;
                    aIndex( 1, 2 ) = 1 ;
                    aIndex( 0, 3 ) = 0 ;
                    aIndex( 1, 3 ) = 1 ;
                    break;
                }
                case( ElementType::QUAD9 ):
                {
                    aIndex( 0, 0 ) = 0 ;
                    aIndex( 1, 0 ) = 0 ;
                    aIndex( 0, 1 ) = 2 ;
                    aIndex( 1, 1 ) = 0 ;
                    aIndex( 0, 2 ) = 2 ;
                    aIndex( 1, 2 ) = 2 ;
                    aIndex( 0, 3 ) = 0 ;
                    aIndex( 1, 3 ) = 2 ;
                    aIndex( 0, 4 ) = 1 ;
                    aIndex( 1, 4 ) = 0 ;
                    aIndex( 0, 5 ) = 2 ;
                    aIndex( 1, 5 ) = 1 ;
                    aIndex( 0, 6 ) = 1 ;
                    aIndex( 1, 6 ) = 2 ;
                    aIndex( 0, 7 ) = 0 ;
                    aIndex( 1, 7 ) = 1 ;
                    aIndex( 0, 8 ) = 1 ;
                    aIndex( 1, 8 ) = 1 ;
                    break;
                }
                case( ElementType::QUAD16 ):
                {
                    aIndex( 0,  0 ) = 0 ;
                    aIndex( 1,  0 ) = 0 ;
                    aIndex( 0,  1 ) = 3 ;
                    aIndex( 1,  1 ) = 0 ;
                    aIndex( 0,  2 ) = 3 ;
                    aIndex( 1,  2 ) = 3 ;
                    aIndex( 0,  3 ) = 0 ;
                    aIndex( 1,  3 ) = 3 ;
                    aIndex( 0,  4 ) = 1 ;
                    aIndex( 1,  4 ) = 0 ;
                    aIndex( 0,  5 ) = 2 ;
                    aIndex( 1,  5 ) = 0 ;
                    aIndex( 0,  6 ) = 3 ;
                    aIndex( 1,  6 ) = 1 ;
                    aIndex( 0,  7 ) = 3 ;
                    aIndex( 1,  7 ) = 2 ;
                    aIndex( 0,  8 ) = 2 ;
                    aIndex( 1,  8 ) = 3 ;
                    aIndex( 0,  9 ) = 1 ;
                    aIndex( 1,  9 ) = 3 ;
                    aIndex( 0, 10 ) = 0 ;
                    aIndex( 1, 10 ) = 2 ;
                    aIndex( 0, 11 ) = 0 ;
                    aIndex( 1, 11 ) = 1 ;
                    aIndex( 0, 12 ) = 1 ;
                    aIndex( 1, 12 ) = 1 ;
                    aIndex( 0, 13 ) = 2 ;
                    aIndex( 1, 13 ) = 1 ;
                    aIndex( 0, 14 ) = 2 ;
                    aIndex( 1, 14 ) = 2 ;
                    aIndex( 0, 15 ) = 1 ;
                    aIndex( 1, 15 ) = 2 ;
                    break;
                }
                case( ElementType::HEX8 ) :
                {
                    aIndex( 0, 0 ) = 0 ;
                    aIndex( 1, 0 ) = 0 ;
                    aIndex( 2, 0 ) = 0 ;
                    aIndex( 0, 1 ) = 1 ;
                    aIndex( 1, 1 ) = 0 ;
                    aIndex( 2, 1 ) = 0 ;
                    aIndex( 0, 2 ) = 1 ;
                    aIndex( 1, 2 ) = 1 ;
                    aIndex( 2, 2 ) = 0 ;
                    aIndex( 0, 3 ) = 0 ;
                    aIndex( 1, 3 ) = 1 ;
                    aIndex( 2, 3 ) = 0 ;
                    aIndex( 0, 4 ) = 0 ;
                    aIndex( 1, 4 ) = 0 ;
                    aIndex( 2, 4 ) = 1 ;
                    aIndex( 0, 5 ) = 1 ;
                    aIndex( 1, 5 ) = 0 ;
                    aIndex( 2, 5 ) = 1 ;
                    aIndex( 0, 6 ) = 1 ;
                    aIndex( 1, 6 ) = 1 ;
                    aIndex( 2, 6 ) = 1 ;
                    aIndex( 0, 7 ) = 0 ;
                    aIndex( 1, 7 ) = 1 ;
                    aIndex( 2, 7 ) = 1 ;
                    break;
                }
                case( ElementType::HEX27 ):
                {
                    aIndex( 0 ,  0 ) = 0 ;
                    aIndex( 1 ,  0 ) = 0 ;
                    aIndex( 2 ,  0 ) = 0 ;
                    aIndex( 0 ,  1 ) = 2 ;
                    aIndex( 1 ,  1 ) = 0 ;
                    aIndex( 2 ,  1 ) = 0 ;
                    aIndex( 0 ,  2 ) = 2 ;
                    aIndex( 1 ,  2 ) = 2 ;
                    aIndex( 2 ,  2 ) = 0 ;
                    aIndex( 0 ,  3 ) = 0 ;
                    aIndex( 1 ,  3 ) = 2 ;
                    aIndex( 2 ,  3 ) = 0 ;
                    aIndex( 0 ,  4 ) = 0 ;
                    aIndex( 1 ,  4 ) = 0 ;
                    aIndex( 2 ,  4 ) = 2 ;
                    aIndex( 0 ,  5 ) = 2 ;
                    aIndex( 1 ,  5 ) = 0 ;
                    aIndex( 2 ,  5 ) = 2 ;
                    aIndex( 0 ,  6 ) = 2 ;
                    aIndex( 1 ,  6 ) = 2 ;
                    aIndex( 2 ,  6 ) = 2 ;
                    aIndex( 0 ,  7 ) = 0 ;
                    aIndex( 1 ,  7 ) = 2 ;
                    aIndex( 2 ,  7 ) = 2 ;
                    aIndex( 0 ,  8 ) = 1 ;
                    aIndex( 1 ,  8 ) = 0 ;
                    aIndex( 2 ,  8 ) = 0 ;
                    aIndex( 0 ,  9 ) = 2 ;
                    aIndex( 1 ,  9 ) = 1 ;
                    aIndex( 2 ,  9 ) = 0 ;
                    aIndex( 0 , 10 ) = 1 ;
                    aIndex( 1 , 10 ) = 2 ;
                    aIndex( 2 , 10 ) = 0 ;
                    aIndex( 0 , 11 ) = 0 ;
                    aIndex( 1 , 11 ) = 1 ;
                    aIndex( 2 , 11 ) = 0 ;
                    aIndex( 0 , 12 ) = 0 ;
                    aIndex( 1 , 12 ) = 0 ;
                    aIndex( 2 , 12 ) = 1 ;
                    aIndex( 0 , 13 ) = 2 ;
                    aIndex( 1 , 13 ) = 0 ;
                    aIndex( 2 , 13 ) = 1 ;
                    aIndex( 0 , 14 ) = 2 ;
                    aIndex( 1 , 14 ) = 2 ;
                    aIndex( 2 , 14 ) = 1 ;
                    aIndex( 0 , 15 ) = 0 ;
                    aIndex( 1 , 15 ) = 2 ;
                    aIndex( 2 , 15 ) = 1 ;
                    aIndex( 0 , 16 ) = 1 ;
                    aIndex( 1 , 16 ) = 0 ;
                    aIndex( 2 , 16 ) = 2 ;
                    aIndex( 0 , 17 ) = 2 ;
                    aIndex( 1 , 17 ) = 1 ;
                    aIndex( 2 , 17 ) = 2 ;
                    aIndex( 0 , 18 ) = 1 ;
                    aIndex( 1 , 18 ) = 2 ;
                    aIndex( 2 , 18 ) = 2 ;
                    aIndex( 0 , 19 ) = 0 ;
                    aIndex( 1 , 19 ) = 1 ;
                    aIndex( 2 , 19 ) = 2 ;
                    aIndex( 0 , 20 ) = 1 ;
                    aIndex( 1 , 20 ) = 1 ;
                    aIndex( 2 , 20 ) = 1 ;
                    aIndex( 0 , 21 ) = 1 ;
                    aIndex( 1 , 21 ) = 1 ;
                    aIndex( 2 , 21 ) = 0 ;
                    aIndex( 0 , 22 ) = 1 ;
                    aIndex( 1 , 22 ) = 1 ;
                    aIndex( 2 , 22 ) = 2 ;
                    aIndex( 0 , 23 ) = 0 ;
                    aIndex( 1 , 23 ) = 1 ;
                    aIndex( 2 , 23 ) = 1 ;
                    aIndex( 0 , 24 ) = 2 ;
                    aIndex( 1 , 24 ) = 1 ;
                    aIndex( 2 , 24 ) = 1 ;
                    aIndex( 0 , 25 ) = 1 ;
                    aIndex( 1 , 25 ) = 0 ;
                    aIndex( 2 , 25 ) = 1 ;
                    aIndex( 0 , 26 ) = 1 ;
                    aIndex( 1 , 26 ) = 2 ;
                    aIndex( 2 , 26 ) = 1 ;
                    break;
                }
                case( ElementType::HEX64 ):
                {
                    aIndex( 0 ,  0 ) = 0 ;
                    aIndex( 1 ,  0 ) = 0 ;
                    aIndex( 2 ,  0 ) = 0 ;
                    aIndex( 0 ,  1 ) = 3 ;
                    aIndex( 1 ,  1 ) = 0 ;
                    aIndex( 2 ,  1 ) = 0 ;
                    aIndex( 0 ,  2 ) = 3 ;
                    aIndex( 1 ,  2 ) = 3 ;
                    aIndex( 2 ,  2 ) = 0 ;
                    aIndex( 0 ,  3 ) = 0 ;
                    aIndex( 1 ,  3 ) = 3 ;
                    aIndex( 2 ,  3 ) = 0 ;
                    aIndex( 0 ,  4 ) = 0 ;
                    aIndex( 1 ,  4 ) = 0 ;
                    aIndex( 2 ,  4 ) = 3 ;
                    aIndex( 0 ,  5 ) = 3 ;
                    aIndex( 1 ,  5 ) = 0 ;
                    aIndex( 2 ,  5 ) = 3 ;
                    aIndex( 0 ,  6 ) = 3 ;
                    aIndex( 1 ,  6 ) = 3 ;
                    aIndex( 2 ,  6 ) = 3 ;
                    aIndex( 0 ,  7 ) = 0 ;
                    aIndex( 1 ,  7 ) = 3 ;
                    aIndex( 2 ,  7 ) = 3 ;
                    aIndex( 0 ,  8 ) = 1 ;
                    aIndex( 1 ,  8 ) = 0 ;
                    aIndex( 2 ,  8 ) = 0 ;
                    aIndex( 0 ,  9 ) = 2 ;
                    aIndex( 1 ,  9 ) = 0 ;
                    aIndex( 2 ,  9 ) = 0 ;
                    aIndex( 0 , 10 ) = 0 ;
                    aIndex( 1 , 10 ) = 1 ;
                    aIndex( 2 , 10 ) = 0 ;
                    aIndex( 0 , 11 ) = 0 ;
                    aIndex( 1 , 11 ) = 2 ;
                    aIndex( 2 , 11 ) = 0 ;
                    aIndex( 0 , 12 ) = 0 ;
                    aIndex( 1 , 12 ) = 0 ;
                    aIndex( 2 , 12 ) = 1 ;
                    aIndex( 0 , 13 ) = 0 ;
                    aIndex( 1 , 13 ) = 0 ;
                    aIndex( 2 , 13 ) = 2 ;
                    aIndex( 0 , 14 ) = 3 ;
                    aIndex( 1 , 14 ) = 1 ;
                    aIndex( 2 , 14 ) = 0 ;
                    aIndex( 0 , 15 ) = 3 ;
                    aIndex( 1 , 15 ) = 2 ;
                    aIndex( 2 , 15 ) = 0 ;
                    aIndex( 0 , 16 ) = 3 ;
                    aIndex( 1 , 16 ) = 0 ;
                    aIndex( 2 , 16 ) = 1 ;
                    aIndex( 0 , 17 ) = 3 ;
                    aIndex( 1 , 17 ) = 0 ;
                    aIndex( 2 , 17 ) = 2 ;
                    aIndex( 0 , 18 ) = 2 ;
                    aIndex( 1 , 18 ) = 3 ;
                    aIndex( 2 , 18 ) = 0 ;
                    aIndex( 0 , 19 ) = 1 ;
                    aIndex( 1 , 19 ) = 3 ;
                    aIndex( 2 , 19 ) = 0 ;
                    aIndex( 0 , 20 ) = 3 ;
                    aIndex( 1 , 20 ) = 3 ;
                    aIndex( 2 , 20 ) = 1 ;
                    aIndex( 0 , 21 ) = 3 ;
                    aIndex( 1 , 21 ) = 3 ;
                    aIndex( 2 , 21 ) = 2 ;
                    aIndex( 0 , 22 ) = 0 ;
                    aIndex( 1 , 22 ) = 3 ;
                    aIndex( 2 , 22 ) = 1 ;
                    aIndex( 0 , 23 ) = 0 ;
                    aIndex( 1 , 23 ) = 3 ;
                    aIndex( 2 , 23 ) = 2 ;
                    aIndex( 0 , 24 ) = 1 ;
                    aIndex( 1 , 24 ) = 0 ;
                    aIndex( 2 , 24 ) = 3 ;
                    aIndex( 0 , 25 ) = 2 ;
                    aIndex( 1 , 25 ) = 0 ;
                    aIndex( 2 , 25 ) = 3 ;
                    aIndex( 0 , 26 ) = 0 ;
                    aIndex( 1 , 26 ) = 1 ;
                    aIndex( 2 , 26 ) = 3 ;
                    aIndex( 0 , 27 ) = 0 ;
                    aIndex( 1 , 27 ) = 2 ;
                    aIndex( 2 , 27 ) = 3 ;
                    aIndex( 0 , 28 ) = 3 ;
                    aIndex( 1 , 28 ) = 1 ;
                    aIndex( 2 , 28 ) = 3 ;
                    aIndex( 0 , 29 ) = 3 ;
                    aIndex( 1 , 29 ) = 2 ;
                    aIndex( 2 , 29 ) = 3 ;
                    aIndex( 0 , 30 ) = 2 ;
                    aIndex( 1 , 30 ) = 3 ;
                    aIndex( 2 , 30 ) = 3 ;
                    aIndex( 0 , 31 ) = 1 ;
                    aIndex( 1 , 31 ) = 3 ;
                    aIndex( 2 , 31 ) = 3 ;
                    aIndex( 0 , 32 ) = 1 ;
                    aIndex( 1 , 32 ) = 1 ;
                    aIndex( 2 , 32 ) = 0 ;
                    aIndex( 0 , 33 ) = 1 ;
                    aIndex( 1 , 33 ) = 2 ;
                    aIndex( 2 , 33 ) = 0 ;
                    aIndex( 0 , 34 ) = 2 ;
                    aIndex( 1 , 34 ) = 2 ;
                    aIndex( 2 , 34 ) = 0 ;
                    aIndex( 0 , 35 ) = 2 ;
                    aIndex( 1 , 35 ) = 1 ;
                    aIndex( 2 , 35 ) = 0 ;
                    aIndex( 0 , 36 ) = 1 ;
                    aIndex( 1 , 36 ) = 0 ;
                    aIndex( 2 , 36 ) = 1 ;
                    aIndex( 0 , 37 ) = 2 ;
                    aIndex( 1 , 37 ) = 0 ;
                    aIndex( 2 , 37 ) = 1 ;
                    aIndex( 0 , 38 ) = 2 ;
                    aIndex( 1 , 38 ) = 0 ;
                    aIndex( 2 , 38 ) = 2 ;
                    aIndex( 0 , 39 ) = 1 ;
                    aIndex( 1 , 39 ) = 0 ;
                    aIndex( 2 , 39 ) = 2 ;
                    aIndex( 0 , 40 ) = 0 ;
                    aIndex( 1 , 40 ) = 1 ;
                    aIndex( 2 , 40 ) = 1 ;
                    aIndex( 0 , 41 ) = 0 ;
                    aIndex( 1 , 41 ) = 1 ;
                    aIndex( 2 , 41 ) = 2 ;
                    aIndex( 0 , 42 ) = 0 ;
                    aIndex( 1 , 42 ) = 2 ;
                    aIndex( 2 , 42 ) = 2 ;
                    aIndex( 0 , 43 ) = 0 ;
                    aIndex( 1 , 43 ) = 2 ;
                    aIndex( 2 , 43 ) = 1 ;
                    aIndex( 0 , 44 ) = 3 ;
                    aIndex( 1 , 44 ) = 1 ;
                    aIndex( 2 , 44 ) = 1 ;
                    aIndex( 0 , 45 ) = 3 ;
                    aIndex( 1 , 45 ) = 2 ;
                    aIndex( 2 , 45 ) = 1 ;
                    aIndex( 0 , 46 ) = 3 ;
                    aIndex( 1 , 46 ) = 2 ;
                    aIndex( 2 , 46 ) = 2 ;
                    aIndex( 0 , 47 ) = 3 ;
                    aIndex( 1 , 47 ) = 1 ;
                    aIndex( 2 , 47 ) = 2 ;
                    aIndex( 0 , 48 ) = 2 ;
                    aIndex( 1 , 48 ) = 3 ;
                    aIndex( 2 , 48 ) = 1 ;
                    aIndex( 0 , 49 ) = 1 ;
                    aIndex( 1 , 49 ) = 3 ;
                    aIndex( 2 , 49 ) = 1 ;
                    aIndex( 0 , 50 ) = 1 ;
                    aIndex( 1 , 50 ) = 3 ;
                    aIndex( 2 , 50 ) = 2 ;
                    aIndex( 0 , 51 ) = 2 ;
                    aIndex( 1 , 51 ) = 3 ;
                    aIndex( 2 , 51 ) = 2 ;
                    aIndex( 0 , 52 ) = 1 ;
                    aIndex( 1 , 52 ) = 1 ;
                    aIndex( 2 , 52 ) = 3 ;
                    aIndex( 0 , 53 ) = 2 ;
                    aIndex( 1 , 53 ) = 1 ;
                    aIndex( 2 , 53 ) = 3 ;
                    aIndex( 0 , 54 ) = 2 ;
                    aIndex( 1 , 54 ) = 2 ;
                    aIndex( 2 , 54 ) = 3 ;
                    aIndex( 0 , 55 ) = 1 ;
                    aIndex( 1 , 55 ) = 2 ;
                    aIndex( 2 , 55 ) = 3 ;
                    aIndex( 0 , 56 ) = 1 ;
                    aIndex( 1 , 56 ) = 1 ;
                    aIndex( 2 , 56 ) = 1 ;
                    aIndex( 0 , 57 ) = 2 ;
                    aIndex( 1 , 57 ) = 1 ;
                    aIndex( 2 , 57 ) = 1 ;
                    aIndex( 0 , 58 ) = 2 ;
                    aIndex( 1 , 58 ) = 2 ;
                    aIndex( 2 , 58 ) = 1 ;
                    aIndex( 0 , 59 ) = 1 ;
                    aIndex( 1 , 59 ) = 2 ;
                    aIndex( 2 , 59 ) = 1 ;
                    aIndex( 0 , 60 ) = 1 ;
                    aIndex( 1 , 60 ) = 1 ;
                    aIndex( 2 , 60 ) = 2 ;
                    aIndex( 0 , 61 ) = 2 ;
                    aIndex( 1 , 61 ) = 1 ;
                    aIndex( 2 , 61 ) = 2 ;
                    aIndex( 0 , 62 ) = 2 ;
                    aIndex( 1 , 62 ) = 2 ;
                    aIndex( 2 , 62 ) = 2 ;
                    aIndex( 0 , 63 ) = 1 ;
                    aIndex( 1 , 63 ) = 2 ;
                    aIndex( 2 , 63 ) = 2 ;
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Invalid Element Type.");
                    break;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        TMatrix::populate_node_param_coords()
        {
            mNodeParamCoords.set_size( mNumberOfDimensions, mNumberOfNodes );

            // scaling factor
            real tScale = 1.0/( ( real ) mLagrangeOrder );

            for( uint k = 0; k < mNumberOfNodes; ++k )
            {

                // save coordinate into memory
                for( uint i = 0; i < mNumberOfDimensions; ++i )
                {
                    // fill in node ijk positions in element
                    mNodeParamCoords( i, k ) = 2. * tScale * mNodeIndex( i, k ) - 1.0;
                }
            }
        }

//------------------------------------------------------------------------------

        real
        TMatrix::shape_function_1d(
                            const uint & aK,
                           const real & aXi )
        {
            // max number of entries in lookup table
            uint tSteps = 2 * ( mBsplineOrder + 1 );

            Vector< real > tDeltaXi( tSteps, 0.0 );

            // temporary matrix that contains B-Spline segments
            for( uint i = 0; i < tSteps; ++i )
            {
                tDeltaXi( i ) = ( ( ( real ) i ) - ( ( real ) mBsplineOrder ) ) * 2.0 - 1.0;
            }

            // temporary matrix that contains evaluated values
            Vector< real > tN( mBsplineOrder + 1, 0.0 );

            // initialize zero order values
            for( uint i = 0; i <= mBsplineOrder; ++i )
            {
                if ( tDeltaXi( i + aK ) <= aXi && aXi < tDeltaXi( i + aK + 1 ) )
                {
                    tN( i ) = 1.0;
                }
            }

            // loop over all orders
            for( uint r = 1; r <= mBsplineOrder; ++r )
            {
                // copy values of tN into old matrix
                Vector< real > tNold( tN );

                // loop over all contributions
                for( uint i=0; i<=mBsplineOrder-r; ++i )
                {
                    // help values
                    real tA = aXi - tDeltaXi( i+aK );
                    real tB = tDeltaXi( i+aK + r + 1 ) - aXi;

                    tN( i ) = 0.5*( tA * tNold( i ) + tB * ( tNold( i + 1 ) ) ) / ( ( real ) r );
                }
            }

            // first value in entry is shape value
            return tN( 0 );
        }

//------------------------------------------------------------------------------

        void
        TMatrix::shape_function(
                const real        & aXi,
                Vector< real >    & aN )
        {
            // temporary value of function
            // evaluate contributions for xi
            Vector< real >  tNxi( mBsplineOrder+1 );

            for( uint i = 0; i <= mBsplineOrder; ++i )
            {
                tNxi( i ) = this->shape_function_1d( i,
                                                     aXi );
            }

            // resort indices
            for( uint k=0; k<mNumberOfBasis; ++k )
            {
                aN( k ) = tNxi( mBasisIndex( 0, k ) );
            }
        }

//------------------------------------------------------------------------------

        void
        TMatrix::shape_function(
                const real        & aXi,
                const real        & aEta,
                Vector< real >    & aN )
        {
            // temporary value of function
            // evaluate contributions for xi and eta

            Vector< real >  tNxi( mBsplineOrder+1 );
            Vector< real >  tNeta( mBsplineOrder+1 );

            for( uint i = 0; i <= mBsplineOrder; ++i )
            {
                tNxi( i ) = this->shape_function_1d( i,
                                                     aXi );
            }
            for( uint j = 0; j <= mBsplineOrder; ++j )
            {
                tNeta( j ) = this->shape_function_1d( j,
                                                      aEta );
            }

            // create shape vector in correct order
            for( uint k = 0; k < mNumberOfBasis; ++k )
            {
                aN( k ) = tNxi( mBasisIndex( 0, k ) )
                        * tNeta( mBasisIndex( 1, k ) );
            }
        }

//------------------------------------------------------------------------------

        void
        TMatrix::shape_function(
                        const real               & aXi,
                        const real               & aEta,
                        const real               & aZeta,
                        Vector< real >           & aN)
        {
            Vector< real >  tNxi( mBsplineOrder+1 );
            Vector< real > tNeta( mBsplineOrder+1 );
            Vector< real > tNzeta( mBsplineOrder+1 );

            for( uint i = 0; i <= mBsplineOrder; ++i )
            {
                tNxi( i ) = this->shape_function_1d(
                        i,
                        aXi );
            }
            for( uint j = 0; j <= mBsplineOrder; ++j )
            {
                tNeta( j ) = this->shape_function_1d(
                        j,
                        aEta );
            }
            for( uint k = 0; k <= mBsplineOrder; ++k )
            {
                tNzeta( k ) = this->shape_function_1d(
                        k,
                        aZeta );
            }

            // create shape vector in correct order
            for( uint k = 0; k < mNumberOfBasis; ++k )
            {
                aN( k ) =    tNxi(   mBasisIndex( 0, k ) )
                           * tNeta(  mBasisIndex( 1, k ) )
                           * tNzeta( mBasisIndex( 2, k ) );
            }
        }

//------------------------------------------------------------------------------

        void
        TMatrix::compute_lagrange_matrix()
        {
            // init the matrix
            mLagrangeMatrix.set_size(  mNumberOfNodes, mNumberOfBasis, 1.0 );

            // loop over all Lagrange nodes
            for( uint k = 0; k < mNumberOfNodes; ++k)
            {
                // loop over all B-Spline Basis
                for( uint j = 0; j < mNumberOfBasis; ++j )
                {
                    // loop over all dimensions
                    for( uint i = 0; i < mNumberOfDimensions; ++i )
                    {
                        mLagrangeMatrix( k, j ) *= this->shape_function_1d( mBasisIndex( i, j ),
                                                                             mNodeParamCoords( i, k ) );
                    }
                }
            }

        }

//------------------------------------------------------------------------------
    }
}