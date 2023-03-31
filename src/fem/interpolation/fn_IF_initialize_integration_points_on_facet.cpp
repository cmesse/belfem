//
// Created by Christian Messe on 16.06.20.
//

#include "fn_IF_initialize_integration_points_on_facet.hpp"

#include "fn_intpoints_auto_integration_order.hpp"
#include "fn_IF_initialize_integration_points.hpp"
#include "meshtools.hpp"
#include "fn_reverse.hpp"

namespace belfem
{
    namespace fem
    {

//------------------------------------------------------------------------------

        // by default we want twice the order of the element
        void
        initialize_integration_points_on_facet(
                const ElementType    & aElementType,
                const uint           & aSideIndex,
                Vector< real >       & aWeights,
                Matrix< real >       & aPoints,
                const uint             aIntegrationOrder,
                const IntegrationScheme  aIntegrationScheme )
        {

            // auto detect integration order if not set
            uint tIntegrationOrder = aIntegrationOrder == 0 ?
                                     auto_integration_order( aElementType ) :
                                     aIntegrationOrder ;

            // find out geometry type of this element
            GeometryType tGeometryType = mesh::geometry_type( aElementType );

            switch ( tGeometryType )
            {
                case( GeometryType::TRI  ) :
                {
                    facetintpoints::intpoints_tri( aSideIndex, aWeights, aPoints, tIntegrationOrder, aIntegrationScheme );
                    break;
                }
                case( GeometryType::QUAD ) :
                {
                    facetintpoints::intpoints_quad( aSideIndex, aWeights, aPoints, tIntegrationOrder, aIntegrationScheme  );
                    break ;
                }
                case( GeometryType::TET ) :
                {
                    facetintpoints::intpoints_tet( aSideIndex, aWeights, aPoints, tIntegrationOrder, aIntegrationScheme  );
                    break ;
                }
                case( GeometryType::PENTA ) :
                {
                    facetintpoints::intpoints_penta( aSideIndex, aWeights, aPoints, tIntegrationOrder, aIntegrationScheme  );
                    break ;
                }
                case( GeometryType::HEX ) :
                {
                    facetintpoints::intpoints_hex( aSideIndex, aWeights, aPoints, tIntegrationOrder, aIntegrationScheme  );
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Illegal Geometry Type" );
                    break;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        facetintpoints::intpoints_tri(
                const uint  aMasterIndex,
                Vector< real > & aWeights,
                Matrix< real > & aPoints,
                const uint  aIntegrationOrder,
                const IntegrationScheme  aIntegrationScheme )
        {
            // compute points for edge
            Matrix< real > tPoints;
            initialize_integration_points(
                    GeometryType::LINE,
                    aWeights,
                    tPoints,
                    aIntegrationOrder,
                    aIntegrationScheme );

            // scale and shift points, so that 0 <= xi <= 1
            tPoints *= 0.5;
            tPoints += 0.5;

            // allocate memory for output
            aPoints.set_size( 2, aWeights.length() , 0.0 );

            switch ( aMasterIndex )
            {
                case ( 0 ) :
                {
                    // eta
                    aPoints.set_row( 1, tPoints.row( 0 ));

                    // xi
                    tPoints *= -1.0;
                    tPoints += 1.0;
                    aPoints.set_row( 0, tPoints.row( 0 ));

                    break;
                }
                case ( 1 ) :
                {
                    // eta
                    tPoints *= -1.0;
                    tPoints += 1.0;
                    aPoints.set_row( 1, tPoints.row( 0 ));

                    // xi = 0

                    break;
                }
                case ( 2 ):
                {
                    // xi
                    aPoints.set_row( 0, tPoints.row( 0 ));

                    // eta = 0

                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Illegal side index for TRI : %u", ( unsigned int ) aMasterIndex );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        facetintpoints::intpoints_quad(
                const uint       aMasterIndex,
                Vector< real > & aWeights,
                Matrix< real > & aPoints,
                const uint       aIntegrationOrder,
                const IntegrationScheme  aIntegrationScheme )
        {
            // compute points for edge
            Matrix< real > tPoints;
            initialize_integration_points(
                    GeometryType::LINE,
                    aWeights,
                    tPoints,
                    aIntegrationOrder,
                    aIntegrationScheme );

            // allocate memory for output
            aPoints.set_size( 2, aWeights.length() );

            switch ( aMasterIndex )
            {
                case ( 0 ) :
                {
                    // xi
                    aPoints.set_row( 0, tPoints.row( 0 ));

                    // eta
                    tPoints.fill( -1.0 );
                    aPoints.set_row( 1, tPoints.row( 0 ));
                    break;
                }
                case ( 1 ) :
                {
                    // eta
                    aPoints.set_row( 1, tPoints.row( 0 ));

                    // xi
                    tPoints.fill( 1.0 );
                    aPoints.set_row( 0, tPoints.row( 0 ));

                    break;
                }
                case ( 2 ) :
                {
                    int n = aWeights.length() - 1 ;
                    uint c = 0 ;
                    for( int k = n ; k >= 0; k--)
                    {
                        // xi
                        aPoints( 0, c ) = tPoints( 0, k );

                        // eta
                        aPoints( 1, c++ ) = 1.0 ;
                    }
                    break;
                }
                case ( 3 ) :
                {

                    int n = aWeights.length() - 1 ;
                    uint c = 0 ;
                    for( int k = n ; k >= 0; k--)
                    {
                        // xi
                        aPoints( 0, c ) =  -1.0 ;

                        // eta
                        aPoints( 1, c++ ) = tPoints( 0, k );
                    }

                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Illegal side index for QUAD : %u", ( unsigned int ) aMasterIndex );
                }
            }

        }

//------------------------------------------------------------------------------

        void
        facetintpoints::intpoints_tet(
                const uint       aMasterIndex,
                Vector< real > & aWeights,
                Matrix< real > & aPoints,
                const uint       aIntegrationOrder,
                const IntegrationScheme  aIntegrationScheme )
        {
            Matrix< real > tPoints;
            initialize_integration_points(
                    GeometryType::TRI,
                    aWeights,
                    tPoints,
                    aIntegrationOrder,
                    aIntegrationScheme );

            // allocate memory
            uint tNumPoints = aWeights.length() ;
            aPoints.set_size( 4, tNumPoints, 0.0 );

            switch( aMasterIndex )
            {
                case( 0 ) :
                {
                    aPoints.set_row( 0, tPoints.row( 0 ) );
                    aPoints.set_row( 2, tPoints.row( 1 ) );
                    aPoints.set_row( 3, tPoints.row( 2 ) );
                    break ;
                }
                case( 1 ) :
                {
                    aPoints.set_row( 1, tPoints.row( 1 ) );
                    aPoints.set_row( 2, tPoints.row( 0 ) );
                    aPoints.set_row( 3, tPoints.row( 2 ) );
                    break ;
                }
                case( 2 ) :
                {
                    aPoints.set_row( 0, tPoints.row( 0 ) );
                    aPoints.set_row( 1, tPoints.row( 2 ) );
                    aPoints.set_row( 3, tPoints.row( 1 ) );
                    break ;
                }

                case( 3 ) :
                {
                    aPoints.set_row( 0, tPoints.row( 0 ) );
                    aPoints.set_row( 1, tPoints.row( 1 ) );
                    aPoints.set_row( 2, tPoints.row( 2 ) );
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Illegal side index for TET : %u", ( unsigned int ) aMasterIndex );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        facetintpoints::intpoints_penta(
                const uint       aMasterIndex,
                Vector< real > & aWeights,
                Matrix< real > & aPoints,
                const uint       aIntegrationOrder,
                const IntegrationScheme  aIntegrationScheme )
        {
            // values for Xi go from 0 to 1 if aSideIndex < 0
            Matrix< real > tXi ;

            // values for Eta go from 1 to 0 if aSideIndex < 0
            Matrix< real > tEta ;

            // values for zeta go from -1 to 1
            Matrix< real > tZeta ;

            // number of integration points
            uint tNumPoints ;

            if( aMasterIndex < 3 )
            {
                Matrix< real > tPoints ;

                // get points for quad
                initialize_integration_points(
                        GeometryType::QUAD,
                        aWeights,
                        tPoints,
                        aIntegrationOrder,
                        aIntegrationScheme );

                tNumPoints = aWeights.length() ;

                tXi.set_size( 1, tNumPoints );
                tEta.set_size( 1, tNumPoints );
                tZeta.set_size( 1, tNumPoints );

                tXi.set_row( 0, tPoints.row( 0 ) );
                tXi *= 0.5 ;
                tXi += 0.5 ;

                tEta.set_row( 0, tXi.row( 0 ) );
                tEta *= -1.0 ;
                tEta +=  1.0 ;

                tZeta.set_row( 0, tPoints.row( 1 ) );
            }
            else
            {
                Matrix< real > tPoints ;

                // get points for triangle
                initialize_integration_points(
                        GeometryType::TRI,
                        aWeights,
                        tPoints,
                        aIntegrationOrder,
                        aIntegrationScheme );

                tNumPoints = aWeights.length() ;

                tXi.set_size( 1, tNumPoints );
                tEta.set_size( 1, tNumPoints );
                tZeta.set_size( 1, tNumPoints );

                tXi.set_row( 0, tPoints.row( 0 ) );
                tEta.set_row( 0, tPoints.row( 0 ) );
            }

            aPoints.set_size( 3, tNumPoints, 0.0 );

            switch( aMasterIndex )
            {
                case ( 0 ) :
                {
                    // xi : 1 --> 0
                    aPoints.set_row( 0, tEta.row( 0 ) );

                    // eta: 0 --> 1
                    aPoints.set_row( 1, tXi.row( 0 ) );

                    // zeta: -1 --> 1

                    break ;
                }
                case( 1 ) :
                {
                    // xi : 0 == 0

                    // eta: 1 --> 0
                    aPoints.set_row( 1, tEta.row( 0 ) );

                    // zeta: -1 --> 1

                    break ;
                }
                case( 2 ) :
                {
                    // xi : 0 --> 1
                    aPoints.set_row( 0, tXi.row( 0 ) );

                    // eta: == 0

                    // zeta: -1 --> 1

                    break ;
                }
                case( 3 ) :
                {
                    // xi from triangle
                    aPoints.set_row( 0, tXi.row( 0 ) );

                    // eta from triangle
                    aPoints.set_row( 1, tEta.row( 0 ) );

                    // zeta is -1
                    tZeta.fill( -1.0 );

                    break ;
                }
                case( 4 ) :
                {
                    // xi from triangle
                    aPoints.set_row( 0, tXi.row( 0 ) );

                    // eta from triangle
                    aPoints.set_row( 1, tEta.row( 0 ) );

                    // zeta is +1
                    tZeta.fill( 1.0 );

                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Illegal side index for PENTA : %u", ( unsigned int ) aMasterIndex );
                }
            }

            // populate zeta value
            aPoints.set_row( 2, tZeta.row( 0 ) );
        }

//------------------------------------------------------------------------------

        void
        facetintpoints::intpoints_hex(
                const uint       aMasterIndex,
                Vector< real > & aWeights,
                Matrix< real > & aPoints,
                const uint       aIntegrationOrder,
                const IntegrationScheme  aIntegrationScheme )
        {
            Matrix< real > tPoints ;


            // get points for quad
            initialize_integration_points(
                    GeometryType::QUAD,
                    aWeights,
                    tPoints,
                    aIntegrationOrder,
                    aIntegrationScheme );

            aPoints.set_size( 3, aWeights.length() ) ;


            auto tXi  = tPoints.row( 0 ) ;
            auto tEta = tPoints.row( 1 );
            Matrix< real > tSide( 1, aWeights.length(), 1. );
            auto tZeta = tSide.row( 0 );

            switch( aMasterIndex )
            {
                case ( 0 ) :
                {
                    // xi
                    aPoints.set_row( 0, tXi );

                    // eta
                    tZeta *= -1 ;
                    aPoints.set_row( 1, tZeta );

                    // zeta
                    aPoints.set_row( 2, tEta );

                    break;
                }
                case ( 1 ) :
                {
                    // xi
                    aPoints.set_row( 0, tZeta);

                    // eta
                    aPoints.set_row( 1, tXi );

                    // zeta
                    aPoints.set_row( 2, tEta );

                    break;
                }
                case ( 2 ) :
                {
                    // xi
                    tXi *= -1. ;
                    aPoints.set_row( 0, tXi );

                    // eta
                    aPoints.set_row( 1, tZeta);

                    // zeta
                    aPoints.set_row( 2, tEta );

                    break;
                }
                case ( 3 ) :
                {
                    // xi
                    tZeta *= -1 ;
                    aPoints.set_row( 0, tZeta);

                    // zeta
                    aPoints.set_row( 2, tXi );

                    // eta
                    aPoints.set_row( 1, tEta );

                    break;
                }
                case ( 4 ) :
                {
                    // eta
                    aPoints.set_row( 0, tEta );

                    // xi
                    aPoints.set_row( 1, tXi );

                    // zeta
                    tZeta *= -1. ;
                    aPoints.set_row( 2, tZeta);

                    break;
                }
                case ( 5 ) :
                {
                    // xi
                    aPoints.set_row( 0, tXi );

                    // eta
                    aPoints.set_row( 1, tEta );

                    // zeta
                    aPoints.set_row( 2, tZeta);

                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Illegal side index for HEX : %u", ( unsigned int ) aMasterIndex );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        facetintpoints::intpoints_tet(
                const uint       aSlaveIndex,
                const uint       aOrientation,
                Vector< real > & aWeights,
                Matrix< real > & aPoints,
                const uint       aIntegrationOrder,
                const IntegrationScheme  aIntegrationScheme )
        {
            Matrix< real > tPoints;
            initialize_integration_points(
                    GeometryType::TRI,
                    aWeights,
                    tPoints,
                    aIntegrationOrder,
                    aIntegrationScheme );

            // allocate memory
            uint tNumPoints = aWeights.length() ;

            aPoints.set_size( 4, tNumPoints, 0.0 );

            // compute the case index
            uint tCase = aSlaveIndex * 3 + aOrientation ;

            // the magic orientation indices
            const Matrix< uint > tIndex{ { 0, 2, 3, 2, 1, 3, 0, 3, 1, 0, 1, 2 },
                                         { 3, 0, 2, 3, 2, 1, 1, 0, 3, 2, 0, 1 },
                                         { 2, 3, 0, 1, 3, 2, 3, 1, 0, 1, 2, 0 },
                                         { 1, 1, 1, 0, 0, 0, 2, 2, 2, 3, 3, 3 } };

            // set the integration points
            aPoints.set_row( tIndex( 0, tCase ), tPoints.row( 0 ) );
            aPoints.set_row( tIndex( 1, tCase ), tPoints.row( 1 ) );
            aPoints.set_row( tIndex( 2, tCase ), tPoints.row( 2 ) );

        }

//------------------------------------------------------------------------------

        void
        facetintpoints::intpoints_hex(
                const uint       aSlaveIndex,
                const uint       aOrientation,
                Vector< real > & aWeights,
                Matrix< real > & aPoints,
                const uint       aIntegrationOrder,
                const IntegrationScheme  aIntegrationScheme )
        {


            Matrix< real > tPoints;
            initialize_integration_points(
                    GeometryType::QUAD,
                    aWeights,
                    tPoints,
                    aIntegrationOrder,
                    aIntegrationScheme );

            Vector< real > tXi  = tPoints.row( 0 );
            Vector< real > tEta = tPoints.row( 1 );

            // allocate memory
            uint tNumPoints = aWeights.length() ;



            // compute the case index
            uint tCase = aSlaveIndex * 4 + aOrientation ;

            // the magic orientation indices
            //                              0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
            const Matrix< int >   tIndex{ { 3,-1,-3, 1, 3,-2,-3, 2, 3, 1,-3,-1, 2,-3,-2, 3, 1,-2,-1, 2, 2,-1,-2, 1},    // xi
                                          { 1, 3,-1,-3, 2, 3,-2,-3,-1, 3, 1,-3, 3, 2,-3,-2, 2, 1,-2,-1, 1, 2,-1,-2} };  // eta
            const Vector< real > tSign  = {-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1}  ;

            // set the integration points
            aPoints.set_size( 3, tNumPoints, tSign( tCase ) );
            if( tIndex( 0, tCase ) < 0 ) tXi  *= -1. ;
            if( tIndex( 1, tCase ) < 0 ) tEta *= -1. ;
            aPoints.set_row( std::abs(tIndex( 0, tCase ))-1, tXi );
            aPoints.set_row( std::abs(tIndex( 1, tCase ))-1, tEta );
        }

//------------------------------------------------------------------------------
    }
}