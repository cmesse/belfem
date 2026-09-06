//
// Created by Christian Messe on 16.06.20.
//

#include "fn_IF_initialize_integration_points_on_facet.hpp"

#include "fn_intpoints_auto_integration_order.hpp"
#include "fn_IF_initialize_integration_points.hpp"
#include "meshtools.hpp"

namespace belfem
{
    namespace fem
    {

//------------------------------------------------------------------------------

        // aIntegrationOrder == 0 selects auto_integration_order() ( 4 / 7 / 10 for linear / quadratic / cubic elements )
        void
        initialize_integration_points_on_facet(
                const ElementType    aElementType,
                const uint           aSideIndex,
                Vector< real >       & aWeights,
                Matrix< real >       & aPoints,
                const uint             aIntegrationOrder,
                const IntegrationScheme  aIntegrationScheme )
        {

            // auto detect integration order if not set
            uint tIntegrationOrder = aIntegrationOrder == 0 ?
                                     auto_integration_order( aElementType ) :
                                     aIntegrationOrder ;

            // Thin-shell elements now share their facet topology with the
            // corresponding volume element (see cl_Element_PENTA6TS.hpp etc.),
            // so we route through the geometry-type dispatch below.

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
                case( GeometryType::PYRA ) :
                {
                    facetintpoints::intpoints_pyra( aSideIndex, aWeights, aPoints, tIntegrationOrder, aIntegrationScheme  );
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Illegal Geometry Type" );
                    break;
                }
            }
        }

        void
        initialize_integration_points_on_facet(
            const ElementType    aElementType,
            const uint           aSideIndex,
            const uint           aOrientation,
            Vector< real > & aWeights,
            Matrix< real > & aPoints,
            const uint             aIntegrationOrder,
            const IntegrationScheme  aIntegrationScheme )
        {
            // auto detect integration order if not set
            uint tIntegrationOrder = aIntegrationOrder == 0 ?
                                     auto_integration_order( aElementType ) :
                                     aIntegrationOrder ;

            // Thin-shell elements share their facet topology with the volume
            // element of the same geometry, so route through the geometry-type
            // dispatch below.

            // find out geometry type of this element
            GeometryType tGeometryType = mesh::geometry_type( aElementType );

            switch ( tGeometryType )
            {
                case( GeometryType::QUAD ) :
                {
                    facetintpoints::intpoints_quad( aSideIndex, aOrientation, aWeights, aPoints, tIntegrationOrder, aIntegrationScheme  );
                    break ;
                }
                case( GeometryType::TET ) :
                {
                    facetintpoints::intpoints_tet( aSideIndex, aOrientation, aWeights, aPoints, tIntegrationOrder, aIntegrationScheme  );
                    break ;
                }
                case( GeometryType::PENTA ) :
                {
                    facetintpoints::intpoints_penta( aSideIndex, aOrientation , aWeights, aPoints, tIntegrationOrder, aIntegrationScheme  );
                    break ;
                }
                case( GeometryType::HEX ) :
                {
                    facetintpoints::intpoints_hex( aSideIndex, aOrientation,  aWeights, aPoints, tIntegrationOrder, aIntegrationScheme  );
                    break ;
                }
                case( GeometryType::PYRA ) :
                {
                    facetintpoints::intpoints_pyra( aSideIndex, aOrientation, aWeights, aPoints, tIntegrationOrder, aIntegrationScheme  );
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
            // values for Xi go from 0 to 1 on the quad faces ( aMasterIndex < 3 )
            Matrix< real > tXi ;

            // values for Eta go from 1 to 0 on the quad faces ( aMasterIndex < 3 )
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
                tEta.set_row( 0, tPoints.row( 1 ) );
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
                    // face [0,3,5,2]: ξ_q drives zeta, η_q drives xi
                    // xi = 0.5*(1 - η_q)
                    tEta.set_row( 0, tZeta.row( 0 ) );
                    tEta *= -0.5 ;
                    tEta +=  0.5 ;
                    aPoints.set_row( 0, tEta.row( 0 ) );

                    // eta: == 0

                    // zeta = ξ_q = 2*tXi - 1
                    tZeta.set_row( 0, tXi.row( 0 ) );
                    tZeta *= 2.0 ;
                    tZeta -= 1.0 ;

                    break ;
                }
                case( 3 ) :
                {
                    // face [0,2,1]: ξ→node0, η→node2, (1-ξ-η)→node1
                    // xi from triangle
                    aPoints.set_row( 0, tXi.row( 0 ) );

                    // eta = 1 - xi - eta  (third barycentric coordinate)
                    tEta *= -1.0 ;
                    tEta +=  1.0 ;
                    tEta -= tXi ;
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
        facetintpoints::intpoints_penta(
                const uint       aSlaveIndex,
                const uint       aOrientation,
                Vector< real > & aWeights,
                Matrix< real > & aPoints,
                const uint       aIntegrationOrder,
                const IntegrationScheme  aIntegrationScheme )
        {
            if( aSlaveIndex < 3 )
            {
                // quad faces (sides 0, 1, 2)
                Matrix< real > tPoints ;
                initialize_integration_points(
                        GeometryType::QUAD,
                        aWeights,
                        tPoints,
                        aIntegrationOrder,
                        aIntegrationScheme );

                uint tNumPoints = aWeights.length() ;
                aPoints.set_size( 3, tNumPoints, 0.0 );

                // map from [-1,1] to [0,1] and [1,0]
                Vector< real > tXi  = tPoints.row( 0 );
                Vector< real > tEta = tPoints.row( 1 );
                Vector< real > tA_xi( tNumPoints );
                Vector< real > tB_xi( tNumPoints );
                Vector< real > tA_eta( tNumPoints );
                Vector< real > tB_eta( tNumPoints );
                for( uint k=0; k<tNumPoints; ++k )
                {
                    tA_xi( k )  = 0.5 * ( 1.0 + tXi( k ) );
                    tB_xi( k )  = 0.5 * ( 1.0 - tXi( k ) );
                    tA_eta( k ) = 0.5 * ( 1.0 + tEta( k ) );
                    tB_eta( k ) = 0.5 * ( 1.0 - tEta( k ) );
                }

                uint tCase = aSlaveIndex * 4 + aOrientation ;

                switch( tCase )
                {
                    // face 0, orientation 0
                    case( 0 ) :
                    {
                        aPoints.set_row( 0, tB_eta );
                        aPoints.set_row( 1, tA_eta );
                        aPoints.set_row( 2, tXi );
                        break ;
                    }
                    // face 0, orientation 1
                    case( 1 ) :
                    {
                        aPoints.set_row( 0, tA_xi );
                        aPoints.set_row( 1, tB_xi );
                        aPoints.set_row( 2, tEta );
                        break ;
                    }
                    // face 0, orientation 2
                    case( 2 ) :
                    {
                        aPoints.set_row( 0, tA_eta );
                        aPoints.set_row( 1, tB_eta );
                        tXi *= -1.0 ;
                        aPoints.set_row( 2, tXi );
                        break ;
                    }
                    // face 0, orientation 3
                    case( 3 ) :
                    {
                        aPoints.set_row( 0, tB_xi );
                        aPoints.set_row( 1, tA_xi );
                        tEta *= -1.0 ;
                        aPoints.set_row( 2, tEta );
                        break ;
                    }
                    // face 1, orientation 0
                    case( 4 ) :
                    {
                        // eta(1) = 0 (already zero from set_size)
                        aPoints.set_row( 1, tB_eta );
                        aPoints.set_row( 2, tXi );
                        break ;
                    }
                    // face 1, orientation 1
                    case( 5 ) :
                    {
                        aPoints.set_row( 1, tA_xi );
                        aPoints.set_row( 2, tEta );
                        break ;
                    }
                    // face 1, orientation 2
                    case( 6 ) :
                    {
                        aPoints.set_row( 1, tA_eta );
                        tXi *= -1.0 ;
                        aPoints.set_row( 2, tXi );
                        break ;
                    }
                    // face 1, orientation 3
                    case( 7 ) :
                    {
                        aPoints.set_row( 1, tB_xi );
                        tEta *= -1.0 ;
                        aPoints.set_row( 2, tEta );
                        break ;
                    }
                    // face 2, orientation 0
                    case( 8 ) :
                    {
                        aPoints.set_row( 0, tB_xi );
                        // eta(2) = 0 (already zero)
                        aPoints.set_row( 2, tEta );
                        break ;
                    }
                    // face 2, orientation 1
                    case( 9 ) :
                    {
                        aPoints.set_row( 0, tB_xi );
                        tEta *= -1.0 ;
                        aPoints.set_row( 2, tEta );
                        break ;
                    }
                    // face 2, orientation 2
                    case( 10 ) :
                    {
                        aPoints.set_row( 0, tA_eta );
                        tXi *= -1.0 ;
                        aPoints.set_row( 2, tXi );
                        break ;
                    }
                    // face 2, orientation 3
                    case( 11 ) :
                    {
                        aPoints.set_row( 0, tA_eta );
                        aPoints.set_row( 2, tXi );
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false,
                            "Illegal case %u for PENTA quad slave (index %u, orientation %u)",
                            ( unsigned int ) tCase,
                            ( unsigned int ) aSlaveIndex,
                            ( unsigned int ) aOrientation );
                    }
                }
            }
            else
            {
                // tri faces (sides 3, 4)
                Matrix< real > tPoints ;
                initialize_integration_points(
                        GeometryType::TRI,
                        aWeights,
                        tPoints,
                        aIntegrationOrder,
                        aIntegrationScheme );

                uint tNumPoints = aWeights.length() ;
                aPoints.set_size( 3, tNumPoints );

                Vector< real > tXi  = tPoints.row( 0 );
                Vector< real > tEta = tPoints.row( 1 );

                // third barycentric coordinate
                Vector< real > tGamma( tNumPoints );
                for( uint k=0; k<tNumPoints; ++k )
                {
                    tGamma( k ) = 1.0 - tXi( k ) - tEta( k );
                }

                uint tCase = ( aSlaveIndex - 3 ) * 3 + aOrientation ;

                switch( tCase )
                {
                    // face 3, orientation 0
                    case( 0 ) :
                    {
                        aPoints.set_row( 0, tXi );
                        aPoints.set_row( 1, tEta );
                        aPoints.set_row( 2, Vector< real >( tNumPoints, -1.0 ) );
                        break ;
                    }
                    // face 3, orientation 1
                    case( 1 ) :
                    {
                        aPoints.set_row( 0, tGamma );
                        aPoints.set_row( 1, tXi );
                        aPoints.set_row( 2, Vector< real >( tNumPoints, -1.0 ) );
                        break ;
                    }
                    // face 3, orientation 2
                    case( 2 ) :
                    {
                        aPoints.set_row( 0, tEta );
                        aPoints.set_row( 1, tGamma );
                        aPoints.set_row( 2, Vector< real >( tNumPoints, -1.0 ) );
                        break ;
                    }
                    // face 4, orientation 0
                    case( 3 ) :
                    {
                        aPoints.set_row( 0, tXi );
                        aPoints.set_row( 1, tGamma );
                        aPoints.set_row( 2, Vector< real >( tNumPoints, 1.0 ) );
                        break ;
                    }
                    // face 4, orientation 1
                    case( 4 ) :
                    {
                        aPoints.set_row( 0, tEta );
                        aPoints.set_row( 1, tXi );
                        aPoints.set_row( 2, Vector< real >( tNumPoints, 1.0 ) );
                        break ;
                    }
                    // face 4, orientation 2
                    case( 5 ) :
                    {
                        aPoints.set_row( 0, tGamma );
                        aPoints.set_row( 1, tEta );
                        aPoints.set_row( 2, Vector< real >( tNumPoints, 1.0 ) );
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false,
                            "Illegal case %u for PENTA tri slave (index %u, orientation %u)",
                            ( unsigned int ) tCase,
                            ( unsigned int ) aSlaveIndex,
                            ( unsigned int ) aOrientation );
                    }
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

        void
        facetintpoints::intpoints_pyra(
                const uint       aMasterIndex,
                Vector< real > & aWeights,
                Matrix< real > & aPoints,
                const uint       aIntegrationOrder,
                const IntegrationScheme  aIntegrationScheme )
        {
            // PYRA5 reference node coordinates
            static const real sXi[]   = { -1.0,  1.0,  1.0, -1.0,  0.0 };
            static const real sEta[]  = { -1.0, -1.0,  1.0,  1.0,  0.0 };
            static const real sZeta[] = {  0.0,  0.0,  0.0,  0.0,  1.0 };

            if( aMasterIndex < 4 )
            {
                // TRI faces (facets 0-3)
                Matrix< real > tPoints ;
                initialize_integration_points(
                        GeometryType::TRI,
                        aWeights,
                        tPoints,
                        aIntegrationOrder,
                        aIntegrationScheme );

                uint tNumPoints = aWeights.length() ;
                aPoints.set_size( 3, tNumPoints );

                // canonical face node ordering from get_nodes_of_facet
                static const uint sFaces[ 4 ][ 3 ] = {
                    { 0, 1, 4 },
                    { 1, 2, 4 },
                    { 2, 3, 4 },
                    { 3, 0, 4 }
                };

                uint nA = sFaces[ aMasterIndex ][ 0 ];
                uint nB = sFaces[ aMasterIndex ][ 1 ];
                uint nC = sFaces[ aMasterIndex ][ 2 ];

                for( uint k = 0; k < tNumPoints; ++k )
                {
                    real tL1 = tPoints( 0, k );
                    real tL2 = tPoints( 1, k );
                    real tL3 = 1.0 - tL1 - tL2 ;

                    aPoints( 0, k ) = tL1 * sXi[ nA ]   + tL2 * sXi[ nB ]   + tL3 * sXi[ nC ];
                    aPoints( 1, k ) = tL1 * sEta[ nA ]  + tL2 * sEta[ nB ]  + tL3 * sEta[ nC ];
                    aPoints( 2, k ) = tL1 * sZeta[ nA ] + tL2 * sZeta[ nB ] + tL3 * sZeta[ nC ];
                }
            }
            else
            {
                // QUAD base (facet 4), nodes [0, 3, 2, 1]
                Matrix< real > tPoints ;
                initialize_integration_points(
                        GeometryType::QUAD,
                        aWeights,
                        tPoints,
                        aIntegrationOrder,
                        aIntegrationScheme );

                uint tNumPoints = aWeights.length() ;
                aPoints.set_size( 3, tNumPoints );

                for( uint k = 0; k < tNumPoints; ++k )
                {
                    real tXi  = tPoints( 0, k );
                    real tEta = tPoints( 1, k );

                    real tN1 = 0.25 * ( 1.0 - tXi ) * ( 1.0 - tEta );
                    real tN2 = 0.25 * ( 1.0 + tXi ) * ( 1.0 - tEta );
                    real tN3 = 0.25 * ( 1.0 + tXi ) * ( 1.0 + tEta );
                    real tN4 = 0.25 * ( 1.0 - tXi ) * ( 1.0 + tEta );

                    aPoints( 0, k ) = tN1 * sXi[ 0 ]  + tN2 * sXi[ 3 ]  + tN3 * sXi[ 2 ]  + tN4 * sXi[ 1 ];
                    aPoints( 1, k ) = tN1 * sEta[ 0 ] + tN2 * sEta[ 3 ] + tN3 * sEta[ 2 ] + tN4 * sEta[ 1 ];
                    aPoints( 2, k ) = 0.0 ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        facetintpoints::intpoints_pyra(
                const uint       aSlaveIndex,
                const uint       aOrientation,
                Vector< real > & aWeights,
                Matrix< real > & aPoints,
                const uint       aIntegrationOrder,
                const IntegrationScheme  aIntegrationScheme )
        {
            // PYRA5 reference node coordinates
            static const real sXi[]   = { -1.0,  1.0,  1.0, -1.0,  0.0 };
            static const real sEta[]  = { -1.0, -1.0,  1.0,  1.0,  0.0 };
            static const real sZeta[] = {  0.0,  0.0,  0.0,  0.0,  1.0 };

            if( aSlaveIndex < 4 )
            {
                // TRI faces (facets 0-3), 3 orientations each
                Matrix< real > tPoints ;
                initialize_integration_points(
                        GeometryType::TRI,
                        aWeights,
                        tPoints,
                        aIntegrationOrder,
                        aIntegrationScheme );

                uint tNumPoints = aWeights.length() ;
                aPoints.set_size( 3, tNumPoints );

                // node indices from orientation table (columns 0-11)
                static const uint sNodes[ 3 ][ 12 ] = {
                    { 0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0, 4 },
                    { 4, 0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0 },
                    { 1, 4, 0, 2, 4, 1, 3, 4, 2, 0, 4, 3 }
                };

                uint tCase = aSlaveIndex * 3 + aOrientation ;
                uint nA = sNodes[ 0 ][ tCase ];
                uint nB = sNodes[ 1 ][ tCase ];
                uint nC = sNodes[ 2 ][ tCase ];

                for( uint k = 0; k < tNumPoints; ++k )
                {
                    real tL1 = tPoints( 0, k );
                    real tL2 = tPoints( 1, k );
                    real tL3 = 1.0 - tL1 - tL2 ;

                    aPoints( 0, k ) = tL1 * sXi[ nA ]   + tL2 * sXi[ nB ]   + tL3 * sXi[ nC ];
                    aPoints( 1, k ) = tL1 * sEta[ nA ]  + tL2 * sEta[ nB ]  + tL3 * sEta[ nC ];
                    aPoints( 2, k ) = tL1 * sZeta[ nA ] + tL2 * sZeta[ nB ] + tL3 * sZeta[ nC ];
                }
            }
            else
            {
                // QUAD base (facet 4), 4 orientations
                Matrix< real > tPoints ;
                initialize_integration_points(
                        GeometryType::QUAD,
                        aWeights,
                        tPoints,
                        aIntegrationOrder,
                        aIntegrationScheme );

                uint tNumPoints = aWeights.length() ;
                aPoints.set_size( 3, tNumPoints );

                // node indices from orientation table (columns 12-15)
                static const uint sNodes[ 4 ][ 4 ] = {
                    { 0, 1, 2, 3 },
                    { 1, 2, 3, 0 },
                    { 2, 3, 0, 1 },
                    { 3, 0, 1, 2 }
                };

                uint nA = sNodes[ aOrientation ][ 0 ];
                uint nB = sNodes[ aOrientation ][ 1 ];
                uint nC = sNodes[ aOrientation ][ 2 ];
                uint nD = sNodes[ aOrientation ][ 3 ];

                for( uint k = 0; k < tNumPoints; ++k )
                {
                    real tXi  = tPoints( 0, k );
                    real tEta = tPoints( 1, k );

                    real tN1 = 0.25 * ( 1.0 - tXi ) * ( 1.0 - tEta );
                    real tN2 = 0.25 * ( 1.0 + tXi ) * ( 1.0 - tEta );
                    real tN3 = 0.25 * ( 1.0 + tXi ) * ( 1.0 + tEta );
                    real tN4 = 0.25 * ( 1.0 - tXi ) * ( 1.0 + tEta );

                    aPoints( 0, k ) = tN1 * sXi[ nA ]  + tN2 * sXi[ nB ]  + tN3 * sXi[ nC ]  + tN4 * sXi[ nD ];
                    aPoints( 1, k ) = tN1 * sEta[ nA ] + tN2 * sEta[ nB ] + tN3 * sEta[ nC ] + tN4 * sEta[ nD ];
                    aPoints( 2, k ) = 0.0 ;
                }
            }
        }

//------------------------------------------------------------------------------

        void
        facetintpoints::intpoints_quad(
                const uint       aSlaveIndex,
                const uint       aOrientation,
                Vector< real > & aWeights,
                Matrix< real > & aPoints,
                const uint       aIntegrationOrder,
                const IntegrationScheme  aIntegrationScheme )
        {
            // QUAD slave path: LINE facets with 2 orientations per facet.
            //   orientation 0 = forward along the facet-local parameter
            //   orientation 1 = reversed along the facet-local parameter
            //
            // Facet-to-reference mapping on the canonical QUAD [-1,+1]^2
            // (matches the master overload above):
            //   facet 0: bottom (eta=-1, xi varies)
            //   facet 1: right  (xi=+1,  eta varies)
            //   facet 2: top    (eta=+1, xi varies)
            //   facet 3: left   (xi=-1,  eta varies)
            Matrix< real > tPoints ;
            initialize_integration_points(
                    GeometryType::LINE,
                    aWeights,
                    tPoints,
                    aIntegrationOrder,
                    aIntegrationScheme );

            uint tNumPoints = aWeights.length() ;
            aPoints.set_size( 2, tNumPoints );

            // pick the varying-coord values (forward or reversed)
            Vector< real > tParam( tNumPoints );
            if ( aOrientation == 0 )
            {
                for( uint k = 0; k < tNumPoints; ++k )
                {
                    tParam( k ) = tPoints( 0, k );
                }
            }
            else
            {
                for( uint k = 0; k < tNumPoints; ++k )
                {
                    tParam( k ) = -tPoints( 0, k );
                }
            }

            switch ( aSlaveIndex )
            {
                case 0 :
                {
                    // bottom: eta = -1, xi = tParam
                    for( uint k = 0; k < tNumPoints; ++k )
                    {
                        aPoints( 0, k ) = tParam( k );
                        aPoints( 1, k ) = -1.0 ;
                    }
                    break ;
                }
                case 1 :
                {
                    // right: xi = +1, eta = tParam
                    for( uint k = 0; k < tNumPoints; ++k )
                    {
                        aPoints( 0, k ) = 1.0 ;
                        aPoints( 1, k ) = tParam( k );
                    }
                    break ;
                }
                case 2 :
                {
                    // top: eta = +1, xi = tParam
                    for( uint k = 0; k < tNumPoints; ++k )
                    {
                        aPoints( 0, k ) = tParam( k );
                        aPoints( 1, k ) = 1.0 ;
                    }
                    break ;
                }
                case 3 :
                {
                    // left: xi = -1, eta = tParam
                    for( uint k = 0; k < tNumPoints; ++k )
                    {
                        aPoints( 0, k ) = -1.0 ;
                        aPoints( 1, k ) = tParam( k );
                    }
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false,
                        "intpoints_quad (slave): invalid facet index %u (expected 0..3)",
                        ( unsigned int ) aSlaveIndex );
                }
            }
        }

//------------------------------------------------------------------------------
    }
}