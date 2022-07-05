//
// Created by Christian Messe on 16.06.20.
//

#include "fn_FEM_initialize_integration_points_on_facet.hpp"

#include "fn_intpoints_auto_integration_order.hpp"
#include "fn_FEM_initialize_integration_points.hpp"
#include "meshtools.hpp"

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
                    facetintpoints::intpoints_tri( aSideIndex, aWeights, aPoints, tIntegrationOrder );
                    break;
                }
                case( GeometryType::QUAD ) :
                {
                    facetintpoints::intpoints_quad( aSideIndex, aWeights, aPoints, tIntegrationOrder );
                    break ;
                }
                case( GeometryType::TET ) :
                {
                    facetintpoints::intpoints_tet( aSideIndex, aWeights, aPoints, tIntegrationOrder );
                    break ;
                }
                case( GeometryType::PENTA ) :
                {
                    facetintpoints::intpoints_penta( aSideIndex, aWeights, aPoints, tIntegrationOrder );
                    break ;
                }
                case( GeometryType::HEX ) :
                {
                    facetintpoints::intpoints_hex( aSideIndex, aWeights, aPoints, tIntegrationOrder );
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
                const uint  aSideIndex,
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

            // scale weights, because edge length is now 1
            aWeights *= 0.5;

            // scale and shift points, so that 0 <= xi <= 1
            tPoints *= 0.5;
            tPoints += 0.5;

            // allocate memory for output
            aPoints.set_size( 2, aWeights.length() , 0.0 );

            switch ( aSideIndex )
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
                    BELFEM_ERROR( false, "Illegal side index for TRI : %u", ( unsigned int ) aSideIndex );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        facetintpoints::intpoints_quad(
                const uint       aSideIndex,
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

            switch ( aSideIndex )
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
                    // xi
                    aPoints.set_row( 0, tPoints.row( 0 ));

                    // eta
                    tPoints.fill( 1.0 );
                    aPoints.set_row( 1, tPoints.row( 0 ));

                    break;
                }
                case ( 3 ) :
                {
                    // eta
                    aPoints.set_row( 1, tPoints.row( 0 ));

                    // xi
                    tPoints.fill( -1.0 );
                    aPoints.set_row( 0, tPoints.row( 0 ));

                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Illegal side index for QUAD : %u", ( unsigned int ) aSideIndex );
                }
            }

        }

//------------------------------------------------------------------------------

        void
        facetintpoints::intpoints_tet(
                const uint       aSideIndex,
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
            aPoints.set_size( 3, tNumPoints, 0.0 );

            switch( aSideIndex )
            {
                case( 0 ) :
                {
                    // xi
                    aPoints.set_row( 0, tPoints.row( 0 ) );

                    // eta
                    aPoints.set_row( 1, tPoints.row( 1 ) );

                    // zeta = 0
                    break ;
                }
                case( 1 ) :
                {
                    // xi = 0

                    // eta
                    aPoints.set_row( 1, tPoints.row( 0 ) );

                    // zeta
                    aPoints.set_row( 2, tPoints.row( 1 ) );

                    break ;
                }
                case( 2 ) :
                {
                    // xi
                    aPoints.set_row( 0, tPoints.row( 1 ) );

                    // eta = 0

                    // zeta
                    aPoints.set_row( 2, tPoints.row( 0 ) );

                    break ;
                }

                case( 3 ) :
                {
                    // xi
                    aPoints.set_row( 0, tPoints.row( 0 ) );

                    // eta
                    for( uint k=0; k<tNumPoints; ++k )
                    {
                        aPoints( 1, k ) = 1.0 - tPoints( 0, k ) - tPoints( 1, k );
                    }

                    // zeta
                    aPoints.set_row( 2, tPoints.row( 1 ) );

                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Illegal side index for TET : %u", ( unsigned int ) aSideIndex );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        facetintpoints::intpoints_tet(
                const uint       aSideIndex,
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
            aPoints.set_size( 3, tNumPoints, 0.0 );

            // compute the case index
            uint tCase = aSideIndex* 3 + aOrientation ;

            // the magic orientation indices
            const Matrix< uint > tIndex{ { 1, 2, 3, 0, 0, 0, 3, 1, 2, 1, 2, 3 },
                                         { 3, 1, 2, 1, 2, 3, 0, 0, 0, 2, 3, 1 },
                                         { 0, 0, 0, 3, 1, 2, 1, 2, 3, 3, 1, 2 } };

            // loop over all coordinates
            for( uint i=0;i<3; ++i )
            {
                // checl case index
                if( tIndex( i, tCase ) == 1 )
                {
                    for( uint k=0; k<tNumPoints; ++k )
                    {
                        aPoints( i, k ) = tPoints( 0, k );
                    }
                }
                else if( tIndex( i, tCase ) == 2 )
                {
                    for( uint k=0; k<tNumPoints; ++k )
                    {
                        aPoints( i, k ) = tPoints( 1, k );
                    }
                }
                else if( tIndex( i, tCase ) == 3 )
                {
                    for( uint k=0; k<tNumPoints; ++k )
                    {
                        aPoints( i, k ) = 1.0 - tPoints( 0, k ) - tPoints( 1, k );
                    }
                }
            }

        }

//------------------------------------------------------------------------------

        void
        facetintpoints::intpoints_penta(
                const uint       aSideIndex,
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

            if( aSideIndex < 3 )
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

            switch( aSideIndex )
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
                    BELFEM_ERROR( false, "Illegal side index for PENTA : %u", ( unsigned int ) aSideIndex );
                }
            }

            // populate zeta value
            aPoints.set_row( 2, tZeta.row( 0 ) );


        }

//------------------------------------------------------------------------------

        void
        facetintpoints::intpoints_hex(
                const uint       aSideIndex,
                Vector< real > & aWeights,
                Matrix< real > & aPoints,
                const uint       aIntegrationOrder,
                const IntegrationScheme  aIntegrationScheme )
        {
            Matrix< real > tPoints ;
            Matrix< real > tSide ;

            // get points for quad
            initialize_integration_points(
                    GeometryType::QUAD,
                    aWeights,
                    tPoints,
                    aIntegrationOrder,
                    aIntegrationScheme );

            aPoints.set_size( 3, aWeights.length() ) ;
            tSide.set_size( 1, aWeights.length() );

            auto tXi = tPoints.row( 0 ) ;
            auto tEta = tPoints.row( 1 );

            switch( aSideIndex )
            {
                case( 0 ) :
                {
                    // xi
                    aPoints.set_row( 0, tPoints.row( 0 ) );

                    // eta
                    tSide.fill( -1.0 );
                    aPoints.set_row( 1, tSide.row( 0 ) );

                    // zeta
                    aPoints.set_row( 2, tPoints.row( 1 ) );

                    break ;
                }
                case( 1 ) :
                {
                    // xi
                    tSide.fill( 1.0 );
                    aPoints.set_row( 0, tSide.row( 0 ) );

                    // eta
                    aPoints.set_row( 1, tPoints.row( 0 ) );

                    // zeta
                    aPoints.set_row( 2, tPoints.row( 1 ) );

                    break ;
                }
                case( 2 ) :
                {
                    // xi
                    aPoints.set_row( 0, tPoints.row( 0 ) );

                    // eta
                    tSide.fill( 1.0 );
                    aPoints.set_row( 1, tSide.row( 0 ) );

                    // zeta
                    aPoints.set_row( 2, tPoints.row( 1 ) );

                    break ;
                }
                case( 3 ) :
                {
                    // xi
                    tSide.fill( -1.0 );
                    aPoints.set_row( 0, tSide.row( 0 ) );

                    // eta
                    aPoints.set_row( 1, tPoints.row( 0 ) );

                    // zeta
                    aPoints.set_row( 2, tPoints.row( 1 ) );

                    break ;
                }
                case( 4 ) :
                {
                    // xi
                    aPoints.set_row( 0, tPoints.row( 0 ) );

                    // eta
                    aPoints.set_row( 1, tPoints.row( 1 ) );

                    // zeta
                    tSide.fill( -1.0 );
                    aPoints.set_row( 2, tSide.row( 0 ) );

                    break ;
                }
                case( 5 ) :
                {
                    // xi
                    aPoints.set_row( 0, tPoints.row( 0 ) );

                    // eta
                    aPoints.set_row( 1, tPoints.row( 1 ) );

                    // zeta
                    tSide.fill( 1.0 );
                    aPoints.set_row( 2, tSide.row( 0 ) );

                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Illegal side index for HEX : %u", ( unsigned int ) aSideIndex );
                }
            }
        }

//------------------------------------------------------------------------------
    }
}