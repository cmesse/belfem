//
// Created by Christian Messe on 2019-01-17.
//
#include "assert.hpp"
#include "fn_intpoints.hpp"
#include "fn_intpoints_gauss_tri1.hpp"
#include "fn_intpoints_gauss_tri3.hpp"
#include "fn_intpoints_gauss_tri4.hpp"
#include "fn_intpoints_gauss_tri7.hpp"
#include "fn_intpoints_gauss_tri12.hpp"
#include "fn_intpoints_gauss_tri15.hpp"
#include "fn_intpoints_gauss_tri16.hpp"
#include "fn_intpoints_gauss_tri19.hpp"
#include "fn_intpoints_gauss_tri25.hpp"
#include "fn_intpoints_gauss_tri28.hpp"
#include "fn_intpoints_gauss_tri33.hpp"
#include "fn_intpoints_gauss_tri37.hpp"
#include "fn_intpoints_gauss_tri42.hpp"
#include "fn_intpoints_gauss_tri49.hpp"
#include "fn_intpoints_gauss_tri55.hpp"
#include "fn_intpoints_gauss_tri60.hpp"
#include "fn_intpoints_gauss_tri67.hpp"
#include "fn_intpoints_gauss_tri73.hpp"
#include "fn_intpoints_gauss_tri79.hpp"


#include "fn_intpoints_gauss_penta1.hpp"
#include "fn_intpoints_gauss_penta5.hpp"
#include "fn_intpoints_gauss_penta8.hpp"
#include "fn_intpoints_gauss_penta11.hpp"
#include "fn_intpoints_gauss_penta16.hpp"
#include "fn_intpoints_gauss_penta29.hpp"
#include "fn_intpoints_gauss_penta40.hpp"
#include "fn_intpoints_gauss_penta50.hpp"
#include "fn_intpoints_gauss_penta71.hpp"

// Hammer and Stroud
#include "fn_intpoints_gauss_hex6.hpp"
#include "fn_intpoints_gauss_hex14.hpp"
#include "fn_intpoints_gauss_hex34.hpp"

#include "fn_intpoints_gauss_tet1.hpp"
#include "fn_intpoints_gauss_tet4.hpp"

// classical Keast
#include "fn_intpoints_gauss_tet5.hpp"
#include "fn_intpoints_gauss_tet11.hpp"
#include "fn_intpoints_gauss_tet15.hpp"
#include "fn_intpoints_gauss_tet24.hpp"
#include "fn_intpoints_gauss_tet31.hpp"
#include "fn_intpoints_gauss_tet35.hpp"

// Shunn and Ham
#include "fn_intpoints_gauss_tet10.hpp"
#include "fn_intpoints_gauss_tet20.hpp"
#include "fn_intpoints_gauss_tet35.hpp"
#include "fn_intpoints_gauss_tet56.hpp"

// Witherden and Vincent
#include "fn_intpoints_gauss_tet46.hpp"
#include "fn_intpoints_gauss_tet81.hpp"

// Vioreanu and Rokhlin
#include "fn_intpoints_gauss_tet165.hpp"
#include "fn_intpoints_gauss_tet220.hpp"

// Zhang, Cui and Liu
#include "fn_intpoints_gauss_tet236.hpp"

#include "intpoints.hpp"

#include "fn_intpoints_gauss_quad8.hpp"
#include "fn_intpoints_gauss_quad12.hpp"
#include "fn_intpoints_gauss_quad20.hpp"
#include "fn_intpoints_gauss_quad28.hpp"
#include "fn_intpoints_gauss_quad37.hpp"

//------------------------------------------------------------------------------
namespace belfem
{
//------------------------------------------------------------------------------

    void
    intpoints( const enum IntegrationScheme aIntegrationScheme,
               const enum GeometryType      aGeometryType,
               const               uint     aOrder,
                          Vector< real >  & aWeights,
                          Matrix< real >  & aPoints )
    {
        // classic is only different from gauss with respect to
        // quad and hex elements
        IntegrationScheme tIntegrationScheme =
                aIntegrationScheme == IntegrationScheme::GAUSSCLASSIC ?
                IntegrationScheme::GAUSS : aIntegrationScheme;

        switch ( tIntegrationScheme )
        {
            case ( IntegrationScheme::GAUSS ):
            {
                switch ( aGeometryType )
                {
                    case ( GeometryType::LINE ) :
                    {
                        integration::gauss_line( aOrder, aWeights, aPoints );
                        break;
                    }
                    case ( GeometryType::QUAD ) :
                    {
                        if ( aIntegrationScheme == IntegrationScheme::GAUSSCLASSIC )
                        {
                            integration::gauss_quad( aOrder, aWeights, aPoints );
                            break;
                        }
                        else
                        {
                            switch ( aOrder )
                            {
                                case ( 4 ) :
                                case ( 5 ) :
                                {
                                    integration::gauss_quad8( aWeights, aPoints );
                                    break;
                                }
                                case ( 6 ) :
                                case ( 7 ) :
                                {
                                    integration::gauss_quad12( aWeights, aPoints );
                                    break;
                                }
                                case ( 8 ) :
                                case ( 9 ) :
                                {
                                    integration::gauss_quad20( aWeights, aPoints );
                                    break;
                                }
                                case ( 10 ) :
                                case ( 11 ) :
                                {
                                    integration::gauss_quad28( aWeights, aPoints );
                                    break;
                                }
                                case ( 12 ) :
                                case ( 13 ) :
                                {
                                    integration::gauss_quad37( aWeights, aPoints );
                                    break;
                                }
                                default :
                                {
                                    integration::gauss_quad( aOrder, aWeights, aPoints );
                                    break;
                                }
                            }
                            break;
                        }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        case ( GeometryType::HEX ) :
                        {
                            if ( aIntegrationScheme == IntegrationScheme::GAUSSCLASSIC )
                            {
                                integration::gauss_hex( aOrder, aWeights, aPoints );
                                break;
                            }
                            else
                            {
                                switch ( aOrder )
                                {
                                    case ( 3 ) :
                                    {
                                        integration::gauss_hex6( aWeights, aPoints );
                                        break;
                                    }
                                    case ( 4 ) :
                                    case ( 5 ) :
                                    {
                                        integration::gauss_hex14( aWeights, aPoints );
                                        break;
                                    }
                                    case ( 6 ) :
                                    case ( 7 ) :
                                    {
                                        integration::gauss_hex34( aWeights, aPoints );
                                        break;
                                    }
                                    default:
                                    {
                                        integration::gauss_hex( aOrder, aWeights, aPoints );
                                        break;
                                    }
                                }
                            }
                            break;
                        }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        case ( GeometryType::TRI ) :
                        {
                            switch ( aOrder )
                            {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                case ( 0 ) :
                                case ( 1 ) :
                                {
                                    integration::gauss_tri1( aWeights, aPoints );
                                    break;
                                }
                                case ( 2 ) :
                                {
                                    integration::gauss_tri3( aWeights, aPoints );
                                    break;
                                }
                                case ( 3 ) :
                                {
                                    integration::gauss_tri4( aWeights, aPoints );
                                    break;
                                }
                                case ( 4 ) :
                                case ( 5 ) :
                                {
                                    integration::gauss_tri7( aWeights, aPoints );
                                    break;
                                }
                                case ( 6 ) :
                                    integration::gauss_tri12( aWeights, aPoints );
                                    break;
                                case ( 7 ) :
                                {
                                    integration::gauss_tri15( aWeights, aPoints );
                                    break;
                                }
                                case ( 8 ) :
                                {
                                    integration::gauss_tri16( aWeights, aPoints );
                                    break;
                                }
                                case ( 9 ) :
                                {
                                    integration::gauss_tri19( aWeights, aPoints );
                                    break;
                                }
                                case ( 10 ) :
                                {
                                    integration::gauss_tri25( aWeights, aPoints );
                                    break;
                                }
                                case ( 11 ) :
                                {
                                    integration::gauss_tri28( aWeights, aPoints );
                                    break;
                                }
                                case ( 12 ) :
                                {
                                    integration::gauss_tri33( aWeights, aPoints );
                                    break;
                                }
                                case ( 13 ) :
                                {
                                    integration::gauss_tri37( aWeights, aPoints );
                                    break;
                                }
                                case ( 14 ) :
                                {
                                    integration::gauss_tri42( aWeights, aPoints );
                                    break;
                                }
                                case ( 15 ) :
                                {
                                    integration::gauss_tri49( aWeights, aPoints );
                                    break;
                                }
                                case ( 16 ) :
                                {
                                    integration::gauss_tri55( aWeights, aPoints );
                                    break;
                                }
                                case ( 17 ) :
                                {
                                    integration::gauss_tri60( aWeights, aPoints );
                                    break;
                                }
                                case ( 18 ) :
                                {
                                    integration::gauss_tri67( aWeights, aPoints );
                                    break;
                                }
                                case ( 19 ) :
                                {
                                    integration::gauss_tri73( aWeights, aPoints );
                                    break;
                                }
                                case ( 20 ) :
                                {
                                    integration::gauss_tri79( aWeights, aPoints );
                                    break;
                                }
                                default:
                                {
                                    BELFEM_ERROR( false, "Invalid Order for GeometryType::TRI" );
                                    break;
                                }
                            }
                            break;
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        } // end geometry type tri;
                        case ( GeometryType::PENTA ) :
                        {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                            switch ( aOrder )
                            {
                                case ( 0 ) :
                                case ( 1 ) :
                                {
                                    integration::gauss_penta1( aWeights, aPoints );
                                    break;
                                }
                                case ( 2 ) :
                                {
                                    integration::gauss_penta5( aWeights, aPoints );
                                    break;
                                }
                                case ( 3 ) :
                                {
                                    integration::gauss_penta8( aWeights, aPoints );
                                    break;
                                }
                                case ( 4 ) :
                                {
                                    integration::gauss_penta11( aWeights, aPoints );
                                    break;
                                }
                                case ( 5 ) :
                                {
                                    integration::gauss_penta16( aWeights, aPoints );
                                    break;
                                }
                                case ( 6 ) :
                                {
                                    integration::gauss_penta29( aWeights, aPoints );
                                    break;
                                }
                                case ( 7 ) :
                                {
                                    integration::gauss_penta40( aWeights, aPoints );
                                    break;
                                }
                                case ( 8 ) :
                                {
                                    integration::gauss_penta50( aWeights, aPoints );
                                    break;
                                }
                                case ( 9 ) :
                                {
                                    integration::gauss_penta71( aWeights, aPoints );
                                    break;
                                }
                                default:
                                {
                                    BELFEM_ERROR( false, "Invalid Order for GeometryType::TRI" );
                                    break;
                                }
                            }

                            // apparently, the weights need to be scaled like this
                            // so that the volume is correct.
                            aWeights /= 2.2;
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                            break;
                        }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        case ( GeometryType::TET ) :
                        {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                            switch ( aOrder )
                            {
                                case ( 0 ) :
                                case ( 1 ) :
                                {
                                    integration::gauss_tet1( aWeights, aPoints );
                                    break;
                                }
                                case ( 2 ) :
                                {
                                    integration::gauss_tet4( aWeights, aPoints );
                                    break;
                                }
                                case ( 3 ) :
                                {
                                    integration::gauss_tet5( aWeights, aPoints );
                                    break;
                                }
                                case ( 4 ) :
                                {
                                    integration::gauss_tet11( aWeights, aPoints );
                                    break;
                                }
                                case ( 5 ) :
                                {
                                    // Keast
                                    integration::gauss_tet15( aWeights, aPoints );
                                    break;
                                }
                                case ( 6 ) :
                                {
                                    // Keast
                                    integration::gauss_tet24( aWeights, aPoints );
                                    // Shunn and Ham
                                    // integration::gauss_tet24( aWeights, aPoints );
                                    break;
                                }
                                case ( 7 ) :
                                {
                                    // Keast
                                    //integration::gauss_tet31( aWeights, aPoints );
                                    // Shunn and Ham
                                    integration::gauss_tet35( aWeights, aPoints );
                                    break;
                                }
                                case ( 8 ) :
                                    // Witherden and Vincent
                                    integration::gauss_tet46( aWeights, aPoints );
                                    break;
                                case ( 9 ) :
                                {
                                    // Shunn and Ham
                                    integration::gauss_tet56( aWeights, aPoints );
                                    break;
                                }
                                case ( 10 ) :
                                {
                                    // Witherden and Vincent
                                    integration::gauss_tet81( aWeights, aPoints );
                                    break;
                                }
                                case ( 11 ) :
                                {
                                    // Vioreanu and Rokhlin :
                                    integration::gauss_tet165( aWeights, aPoints );
                                    break;
                                }
                                case( 12 ) :
                                case( 13 ) :
                                {
                                    // Vioreanu and Rokhlin :
                                    integration::gauss_tet220( aWeights, aPoints );
                                    break;
                                }
                                case( 14 ) :
                                {
                                    // Zhang, Cui and Liu :
                                    integration::gauss_tet236( aWeights, aPoints );
                                    break;
                                }
                                default:
                                {
                                    BELFEM_ERROR( false, "Invalid Order for GeometryType::TET" );
                                    break;
                                }
                            }
                            break;
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        } // end geometry type tet;
                        default :
                        {
                            BELFEM_ERROR( false, "Invalid geometry type" );
                            break;
                        }
                    }
                }
                break ;
            }
            // end integration type gauss
            case ( IntegrationScheme::LOBATTO ) :
            {
                switch ( aGeometryType )
                {
                    case ( GeometryType::LINE ) :
                    {
                        int tN = ceil( 0.5 * ( aOrder + 3 ));
                        aWeights.set_size( tN );
                        aPoints.set_size( 1, tN );
                        Vector< real > tPoints( tN );

                        intpoints_lobatto( &tN, aWeights.data(), tPoints.data() );
                        for( int k=0; k<tN; ++k )
                        {
                            aPoints( 0, k ) = tPoints( k );
                        }

                        break;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "Invalid geometry type" );
                        break;
                    }
                }
                break;
            }

            default :
            {
                BELFEM_ERROR( false, "Invalid integration type" );
                break;
            }
        }
    }
//------------------------------------------------------------------------------
    namespace integration
    {
//------------------------------------------------------------------------------

        void
        gauss_line(
                const int       aOrder,
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {
            int tN = ceil( 0.5 * ( aOrder ) + 1 );
            aWeights.set_size( tN );

            // points as vector
            Vector <real> tPoints( tN );
            intpoints_gauss( &tN, aWeights.data(), tPoints.data());

            // allocate points matrix
            aPoints.set_size( 1, tN );

            // write data into vector
            for ( int k = 0; k < tN; ++k )
            {
                aPoints( 0, k ) = tPoints( k );
            }
        }

//------------------------------------------------------------------------------

        void
        gauss_quad(
                const int       aOrder,
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {
            int tNi = ceil( 0.5 * ( aOrder ) + 1 );

            // 1D points
            Vector <real> tWeights( tNi );
            Vector <real> tPoints( tNi );
            intpoints_gauss( &tNi, tWeights.data(), tPoints.data() );

            // total number of points
            int tN = tNi * tNi;

            // allocate memory for vector
            aPoints.set_size( 2, tN );
            aWeights.set_size( tN );

            // write data into vector
            uint tCount = 0;
            for ( int j = 0; j < tNi; ++j )
            {
                for ( int i = 0; i < tNi; ++i )
                {
                    aPoints( 0, tCount ) = tPoints( i );
                    aPoints( 1, tCount ) = tPoints( j );
                    aWeights( tCount++ ) = tWeights( i ) * tWeights( j );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        gauss_hex(
                const int       aOrder,
                Vector <real> & aWeights,
                Matrix <real> & aPoints )
        {
            int tNi = ceil( 0.5 * ( aOrder ) + 1 );

            // 1D points
            Vector <real> tWeights( tNi );
            Vector <real> tPoints( tNi );
            intpoints_gauss( &tNi, tWeights.data(), tPoints.data());

            // total number of points
            int tN = tNi * tNi * tNi;

            // allocate memory for vector
            aPoints.set_size( 3, tN );
            aWeights.set_size( tN );

            // write data into vector
            uint tCount = 0;
            for ( int k = 0; k < tNi; ++k )
            {
                for ( int j = 0; j < tNi; ++j )
                {
                    for ( int i = 0; i < tNi; ++i )
                    {
                        aPoints( 0, tCount ) = tPoints( i );
                        aPoints( 1, tCount ) = tPoints( j );
                        aPoints( 2, tCount ) = tPoints( k );
                        aWeights( tCount++ ) = tWeights( i ) * tWeights( j ) * tWeights( k );
                    }
                }
            }
        }
//------------------------------------------------------------------------------
    }
}
