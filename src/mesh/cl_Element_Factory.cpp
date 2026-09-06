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

#include "cl_Element_Factory.hpp"
#include "meshtools.hpp"

#include "cl_ElementTemplate.hpp"

#include "cl_Element_Vertex.hpp"
#include "cl_Element_LINE2.hpp"
#include "cl_Element_LINE3.hpp"
#include "cl_Element_LINE4.hpp"
#include "cl_Element_LINE5.hpp"
#include "cl_Element_TRI3.hpp"
#include "cl_Element_TRI6.hpp"
#include "cl_Element_TRI10.hpp"
#include "cl_Element_TRI15.hpp"
#include "cl_Element_QUAD4.hpp"
#include "cl_Element_QUAD8.hpp"
#include "cl_Element_QUAD9.hpp"
#include "cl_Element_QUAD16.hpp"
#include "cl_Element_TET4.hpp"
#include "cl_Element_TET10.hpp"
#include "cl_Element_TET20.hpp"
#include "cl_Element_TET35.hpp"
#include "cl_Element_PENTA6.hpp"
#include "cl_Element_PENTA15.hpp"
#include "cl_Element_PENTA18.hpp"
#include "cl_Element_PYRA5.hpp"
#include "cl_Element_PYRA13.hpp"
#include "cl_Element_PYRA14.hpp"
#include "cl_Element_HEX8.hpp"
#include "cl_Element_HEX20.hpp"
#include "cl_Element_HEX27.hpp"
#include "cl_Element_HEX64.hpp"
#include "cl_Element_QUAD4TS.hpp"
#include "cl_Element_QUAD9TS.hpp"
#include "cl_Element_PENTA6TS.hpp"
#include "cl_Element_PENTA18TS.hpp"
#include "cl_Element_HEX8TS.hpp"
#include "cl_Element_HEX8TB.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        Element *
        ElementFactory::create_element(
                const ElementType aType, const id_t aID ) const
        {
            switch ( aType )
            {
                case ( ElementType::VERTEX ) :
                {
                    return new ElementTemplate< 1, 1, 0, 0, 0 >( aID );
                }
                case ( ElementType::LINE2 ) :
                {
                    return new ElementTemplate< 2, 2, 1, 0, 0 >( aID );
                }
                case ( ElementType::LINE3 ) :
                {
                    return new ElementTemplate< 3, 2, 1, 0, 0 >( aID );
                }
                case ( ElementType::LINE4 ) :
                {
                    return new ElementTemplate< 4, 2, 1, 0, 0 >( aID );
                }
                case ( ElementType::LINE5 ) :
                {
                    return new ElementTemplate< 5, 2, 1, 0, 0 >( aID );
                }
                case ( ElementType::TRI3 ) :
                {
                    return new ElementTemplate< 3, 3, 3, 3, 1 >( aID );
                }
                case ( ElementType::TRI6 ) :
                {
                    return new ElementTemplate< 6, 3, 3, 3, 1 >( aID );
                }
                case ( ElementType::TRI10 ) :
                {
                    return new ElementTemplate< 10, 3, 3, 3, 1 >( aID );
                }
                case ( ElementType::TRI15 ) :
                {
                    return new ElementTemplate< 15, 3, 3, 3, 1 >( aID );
                }
                case ( ElementType::QUAD4 ) :
                {
                    return new ElementTemplate< 4, 4, 4, 4, 1 >( aID );
                }
                case ( ElementType::QUAD8 ) :
                {
                    return new ElementTemplate< 8, 4, 4, 4, 1 >( aID );
                }
                case ( ElementType::QUAD9 ) :
                {
                    return new ElementTemplate< 9, 4, 4, 4, 1 >( aID );
                }
                case ( ElementType::QUAD16 ) :
                {
                    return new ElementTemplate< 16, 4, 4, 4, 1 >( aID );
                }
                case ( ElementType::TET4 ) :
                {
                    return new ElementTemplate< 4, 4, 6, 4, 4 >( aID );
                }
                case ( ElementType::TET10 ) :
                {
                    return new ElementTemplate< 10, 4, 6, 4, 4 >( aID );
                }
                case ( ElementType::TET20 ) :
                {
                    return new ElementTemplate< 20, 4, 6, 4, 4 >( aID );
                }
                case ( ElementType::TET35 ) :
                {
                    return new ElementTemplate< 35, 4, 6, 4, 4 >( aID );
                }
                case ( ElementType::PENTA6 ) :
                {
                    return new ElementTemplate< 6, 6, 9, 5, 5 >( aID );
                }
                case ( ElementType::PENTA15 ) :
                {
                    return new ElementTemplate< 15, 6, 9, 5, 5 >( aID );
                }
                case ( ElementType::PYRA5 ) :
                {
                    return new ElementTemplate< 5, 5, 8, 5, 5 >( aID );
                }
                case ( ElementType::PYRA13 ) :
                {
                    return new ElementTemplate< 13, 5, 8, 5, 5 >( aID );
                }
                case ( ElementType::PYRA14 ) :
                {
                    return new ElementTemplate< 14, 5, 8, 5, 5 >( aID );
                }
                case ( ElementType::PENTA18 ) :
                {
                    return new ElementTemplate< 18, 6, 9, 5, 5 >( aID );
                }
                case ( ElementType::HEX8 ) :
                {
                    return new ElementTemplate< 8, 8, 12, 6, 6 >( aID );
                }
                case ( ElementType::HEX20 ) :
                {
                    return new ElementTemplate< 20, 8, 12, 6, 6 >( aID );
                }
                case ( ElementType::HEX27 ) :
                {
                    return new ElementTemplate< 27, 8, 12, 6, 6 >( aID );
                }
                case ( ElementType::HEX64 ) :
                {
                    return new ElementTemplate< 64, 8, 12, 6, 6 >( aID );
                }
                case ( ElementType::QUAD4TS ) :
                {
                    return new ElementTemplate< 4, 4, 2, 4, 1 >( aID );
                }
                case ( ElementType::QUAD9TS ) :
                {
                    return new ElementTemplate< 9, 4, 3, 4, 1 >( aID );
                }
                case ( ElementType::PENTA6TS ) :
                {
                    return new ElementTemplate< 6, 6, 6, 5, 1 >( aID );
                }
                case( ElementType::PENTA18TS ) :
                {
                    return new ElementTemplate< 18, 6, 9, 5, 3 >( aID );
                }
                case ( ElementType::HEX8TS ) :
                {
                    return new ElementTemplate< 8, 8, 8, 6, 1 >( aID );
                }
                case ( ElementType::HEX8TB ) :
                {
                    return new ElementTemplate< 8, 8, 4, 6, 1 >( aID );
                }
                default:
                {
                    BELFEM_ERROR( false, "Unknown element type" );
                    return nullptr ;
                }
            }
        }

//------------------------------------------------------------------------------

        ReferenceElement *
        ElementFactory::create_reference_element(
            const ElementType aType ) const
        {
            ReferenceElement * aElement = new ReferenceElement( this->create_element( aType, 1 ), { } );

            Matrix< real > tX ;
            this->create_unity_nodes( aType, tX ) ;

            uint n = tX.n_cols();
            Cell< Node * > & tNodes = aElement->nodes() ;
            tNodes.set_size( n, nullptr );



            for( uint k = 0 ; k < n ; ++k )
            {
                Node * tNode = new Node( k+1,
                       tX( 0, k ),
                    tX( 1, k ),
                    tX( 2, k ) ) ;
                tNodes( k ) = tNode ;
                aElement->element()->insert_node( tNode, k );
            }

            return aElement ;

        }

//------------------------------------------------------------------------------

        void
        ElementFactory::create_unity_nodes(
                const ElementType   aType,
                Matrix< real >    & aNodes ) const
        {
            Matrix< real > & X = aNodes;

            auto set_node = [&]( const uint aIndex,
                                 const real aX,
                                 const real aY = 0.0,
                                 const real aZ = 0.0 )
            {
                X( 0, aIndex ) = aX;
                X( 1, aIndex ) = aY;
                X( 2, aIndex ) = aZ;
            };

            const real tThird        = 1.0 / 3.0;
            const real tTwoThirds    = 2.0 / 3.0;

            // Equilateral triangle: side = 2/3^(1/4), area = 1
            const real tTriSide      = 2.0 / std::sqrt( std::sqrt( 3.0 ) );
            const real tTriHeight    = std::sqrt( std::sqrt( 3.0 ) );

            // Regular tetrahedron: side = (6*sqrt(2))^(1/3), volume = 1
            const real tTetSide      = std::cbrt( 6.0 * std::sqrt( 2.0 ) );
            const real tTetBaseH     = 0.5 * std::sqrt( 3.0 ) * tTetSide;
            const real tTetHeight    = std::sqrt( 2.0 / 3.0 ) * tTetSide;

            // Equilateral prism: side = (4/sqrt(3))^(1/3), volume = 1
            const real tPentaSide    = std::cbrt( 4.0 / std::sqrt( 3.0 ) );
            const real tPentaHeight  = 0.5 * std::sqrt( 3.0 ) * tPentaSide;

            // Square pyramid: scale = (3/4)^(1/3), volume = 1
            const real tPyraScale    = std::cbrt( 3.0 / 4.0 );
            const real tPyraZBase    = -0.25 * tPyraScale;

            // Corner coordinates (centroid = origin)
            const real tTriX[ 3 ] = { -0.5 * tTriSide,  0.5 * tTriSide, 0.0 };
            const real tTriY[ 3 ] = { -tTriHeight / 3.0, -tTriHeight / 3.0, 2.0 * tTriHeight / 3.0 };

            const real tTetX[ 4 ] = { -0.5 * tTetSide,  0.5 * tTetSide, 0.0, 0.0 };
            const real tTetY[ 4 ] = { -tTetBaseH / 3.0, -tTetBaseH / 3.0, 2.0 * tTetBaseH / 3.0, 0.0 };
            const real tTetZ[ 4 ] = { -0.25 * tTetHeight, -0.25 * tTetHeight,
                                      -0.25 * tTetHeight,  0.75 * tTetHeight };

            const real tPentaX[ 3 ] = { -0.5 * tPentaSide,  0.5 * tPentaSide, 0.0 };
            const real tPentaY[ 3 ] = { -tPentaHeight / 3.0, -tPentaHeight / 3.0,
                                         2.0 * tPentaHeight / 3.0 };

            // Barycentric-to-physical mapping for line elements
            auto load_line = [&]( const real * aXiHat, const uint aNumNodes )
            {
                X.set_size( 3, aNumNodes, 0.0 );
                for( uint k = 0; k < aNumNodes; ++k )
                {
                    set_node( k, 0.5 * aXiHat[ k ] );
                }
            };

            // Barycentric-to-physical mapping for triangular elements
            auto load_tri = [&]( const real ( * aXiHat )[ 2 ], const uint aNumNodes )
            {
                X.set_size( 3, aNumNodes, 0.0 );
                for( uint k = 0; k < aNumNodes; ++k )
                {
                    const real tXi  = aXiHat[ k ][ 0 ];
                    const real tEta = aXiHat[ k ][ 1 ];
                    const real tTau = 1.0 - tXi - tEta;

                    set_node( k,
                              tXi  * tTriX[ 0 ] + tEta * tTriX[ 1 ] + tTau * tTriX[ 2 ],
                              tXi  * tTriY[ 0 ] + tEta * tTriY[ 1 ] + tTau * tTriY[ 2 ] );
                }
            };

            // Tensor-product mapping for quad elements: [-1,1]^2 -> [-0.5,0.5]^2
            auto load_quad = [&]( const real ( * aXiHat )[ 2 ], const uint aNumNodes )
            {
                X.set_size( 3, aNumNodes, 0.0 );
                for( uint k = 0; k < aNumNodes; ++k )
                {
                    set_node( k, 0.5 * aXiHat[ k ][ 0 ], 0.5 * aXiHat[ k ][ 1 ] );
                }
            };

            // Barycentric-to-physical mapping for tet elements
            auto load_tet = [&]( const real ( * aXiHat )[ 3 ], const uint aNumNodes )
            {
                X.set_size( 3, aNumNodes, 0.0 );
                for( uint k = 0; k < aNumNodes; ++k )
                {
                    const real tXi   = aXiHat[ k ][ 0 ];
                    const real tEta  = aXiHat[ k ][ 1 ];
                    const real tZeta = aXiHat[ k ][ 2 ];
                    const real tTau  = 1.0 - tXi - tEta - tZeta;

                    set_node( k,
                              tXi   * tTetX[ 0 ] + tZeta * tTetX[ 1 ]
                            + tEta  * tTetX[ 2 ] + tTau  * tTetX[ 3 ],
                              tXi   * tTetY[ 0 ] + tZeta * tTetY[ 1 ]
                            + tEta  * tTetY[ 2 ] + tTau  * tTetY[ 3 ],
                              tXi   * tTetZ[ 0 ] + tZeta * tTetZ[ 1 ]
                            + tEta  * tTetZ[ 2 ] + tTau  * tTetZ[ 3 ] );
                }
            };

            // Prismatic mapping: triangle base x [-1,1] axial
            auto load_penta = [&]( const real ( * aXiHat )[ 3 ], const uint aNumNodes )
            {
                X.set_size( 3, aNumNodes, 0.0 );
                for( uint k = 0; k < aNumNodes; ++k )
                {
                    const real tXi   = aXiHat[ k ][ 0 ];
                    const real tEta  = aXiHat[ k ][ 1 ];
                    const real tZeta = aXiHat[ k ][ 2 ];
                    const real tTau  = 1.0 - tXi - tEta;

                    set_node( k,
                              tXi  * tPentaX[ 0 ] + tEta * tPentaX[ 1 ] + tTau * tPentaX[ 2 ],
                              tXi  * tPentaY[ 0 ] + tEta * tPentaY[ 1 ] + tTau * tPentaY[ 2 ],
                              0.5 * tPentaSide * tZeta );
                }
            };

            // Scaled pyramid mapping with centroid shift
            auto load_pyra = [&]( const real ( * aXiHat )[ 3 ], const uint aNumNodes )
            {
                X.set_size( 3, aNumNodes, 0.0 );
                for( uint k = 0; k < aNumNodes; ++k )
                {
                    set_node( k,
                              tPyraScale * aXiHat[ k ][ 0 ],
                              tPyraScale * aXiHat[ k ][ 1 ],
                              tPyraZBase + tPyraScale * aXiHat[ k ][ 2 ] );
                }
            };

            // Tensor-product mapping for hex elements: [-1,1]^3 -> [-0.5,0.5]^3
            auto load_hex = [&]( const real ( * aXiHat )[ 3 ], const uint aNumNodes )
            {
                X.set_size( 3, aNumNodes, 0.0 );
                for( uint k = 0; k < aNumNodes; ++k )
                {
                    set_node( k,
                              0.5 * aXiHat[ k ][ 0 ],
                              0.5 * aXiHat[ k ][ 1 ],
                              0.5 * aXiHat[ k ][ 2 ] );
                }
            };

            // Parametric coordinates below match the interpolation function headers
            // (cl_IF_*.hpp). Exodus II ordering for order <= 2, Gmsh for order > 2.
            switch ( aType )
            {
                case( ElementType::VERTEX ) :
                {
                    X.set_size( 3, 1, 0.0 );
                    break;
                }
                case( ElementType::LINE2 ) :
                {
                    const real tXiHat[] = { -1.0, 1.0 };
                    load_line( tXiHat, 2 );
                    break;
                }
                case( ElementType::LINE3 ) :
                {
                    const real tXiHat[] = { -1.0, 1.0, 0.0 };
                    load_line( tXiHat, 3 );
                    break;
                }
                case( ElementType::LINE4 ) :
                {
                    const real tXiHat[] = { -1.0, 1.0, -tThird, tThird };
                    load_line( tXiHat, 4 );
                    break;
                }
                case( ElementType::LINE5 ) :
                {
                    const real tXiHat[] = { -1.0, 1.0, -0.5, 0.0, 0.5 };
                    load_line( tXiHat, 5 );
                    break;
                }
                case( ElementType::TRI3 ) :
                {
                    const real tXiHat[][ 2 ] = {
                            { 1.0, 0.0 }, { 0.0, 1.0 }, { 0.0, 0.0 }
                    };
                    load_tri( tXiHat, 3 );
                    break;
                }
                case( ElementType::TRI6 ) :
                {
                    const real tXiHat[][ 2 ] = {
                            { 1.0, 0.0 }, { 0.0, 1.0 }, { 0.0, 0.0 },
                            { 0.5, 0.5 }, { 0.0, 0.5 }, { 0.5, 0.0 }
                    };
                    load_tri( tXiHat, 6 );
                    break;
                }
                case( ElementType::TRI10 ) :
                {
                    const real tXiHat[][ 2 ] = {
                            { 1.0, 0.0 }, { 0.0, 1.0 }, { 0.0, 0.0 },
                            { tTwoThirds, tThird }, { tThird, tTwoThirds },
                            { 0.0, tTwoThirds }, { 0.0, tThird },
                            { tThird, 0.0 }, { tTwoThirds, 0.0 },
                            { tThird, tThird }
                    };
                    load_tri( tXiHat, 10 );
                    break;
                }
                case( ElementType::TRI15 ) :
                {
                    const real tXiHat[][ 2 ] = {
                            { 1.0, 0.0 }, { 0.0, 1.0 }, { 0.0, 0.0 },
                            { 0.75, 0.25 }, { 0.50, 0.50 }, { 0.25, 0.75 },
                            { 0.00, 0.75 }, { 0.00, 0.50 }, { 0.00, 0.25 },
                            { 0.25, 0.00 }, { 0.50, 0.00 }, { 0.75, 0.00 },
                            { 0.50, 0.25 }, { 0.25, 0.50 }, { 0.25, 0.25 }
                    };
                    load_tri( tXiHat, 15 );
                    break;
                }
                case( ElementType::QUAD4 ) :
                case( ElementType::QUAD4TS ) :
                {
                    const real tXiHat[][ 2 ] = {
                            { -1.0, -1.0 }, {  1.0, -1.0 },
                            {  1.0,  1.0 }, { -1.0,  1.0 }
                    };
                    load_quad( tXiHat, 4 );
                    break;
                }
                case( ElementType::QUAD8 ) :
                {
                    const real tXiHat[][ 2 ] = {
                            { -1.0, -1.0 }, {  1.0, -1.0 }, {  1.0,  1.0 }, { -1.0,  1.0 },
                            {  0.0, -1.0 }, {  1.0,  0.0 }, {  0.0,  1.0 }, { -1.0,  0.0 }
                    };
                    load_quad( tXiHat, 8 );
                    break;
                }
                case( ElementType::QUAD9 ) :
                case( ElementType::QUAD9TS ) :
                {
                    const real tXiHat[][ 2 ] = {
                            { -1.0, -1.0 }, {  1.0, -1.0 }, {  1.0,  1.0 }, { -1.0,  1.0 },
                            {  0.0, -1.0 }, {  1.0,  0.0 }, {  0.0,  1.0 }, { -1.0,  0.0 },
                            {  0.0,  0.0 }
                    };
                    load_quad( tXiHat, 9 );
                    break;
                }
                case( ElementType::QUAD16 ) :
                {
                    const real tXiHat[][ 2 ] = {
                            { -1.0,   -1.0 }, {  1.0,   -1.0 }, {  1.0,    1.0 }, { -1.0,    1.0 },
                            { -tThird, -1.0 }, {  tThird, -1.0 }, {  1.0, -tThird }, {  1.0,  tThird },
                            {  tThird,  1.0 }, { -tThird,  1.0 }, { -1.0,  tThird }, { -1.0, -tThird },
                            { -tThird, -tThird }, {  tThird, -tThird }, {  tThird,  tThird }, { -tThird,  tThird }
                    };
                    load_quad( tXiHat, 16 );
                    break;
                }
                case( ElementType::TET4 ) :
                {
                    const real tXiHat[][ 3 ] = {
                            { 1.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 },
                            { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 0.0 }
                    };
                    load_tet( tXiHat, 4 );
                    break;
                }
                case( ElementType::TET10 ) :
                {
                    const real tXiHat[][ 3 ] = {
                            { 1.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 0.0 },
                            { 0.5, 0.0, 0.5 }, { 0.0, 0.5, 0.5 }, { 0.5, 0.5, 0.0 },
                            { 0.5, 0.0, 0.0 }, { 0.0, 0.0, 0.5 }, { 0.0, 0.5, 0.0 }
                    };
                    load_tet( tXiHat, 10 );
                    break;
                }
                case( ElementType::TET20 ) :
                {
                    const real tXiHat[][ 3 ] = {
                            { 1.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 0.0 },
                            { tTwoThirds, 0.0, tThird }, { tThird, 0.0, tTwoThirds },
                            { 0.0, tThird, tTwoThirds }, { 0.0, tTwoThirds, tThird },
                            { tThird, tTwoThirds, 0.0 }, { tTwoThirds, tThird, 0.0 },
                            { tThird, 0.0, 0.0 }, { tTwoThirds, 0.0, 0.0 },
                            { 0.0, tThird, 0.0 }, { 0.0, tTwoThirds, 0.0 },
                            { 0.0, 0.0, tThird }, { 0.0, 0.0, tTwoThirds },
                            { tThird, tThird, tThird }, { tThird, 0.0, tThird },
                            { tThird, tThird, 0.0 }, { 0.0, tThird, tThird }
                    };
                    load_tet( tXiHat, 20 );
                    break;
                }
                case( ElementType::TET35 ) :
                {
                    const real tXiHat[][ 3 ] = {
                            { 1.00, 0.00, 0.00 }, { 0.00, 0.00, 1.00 }, { 0.00, 1.00, 0.00 }, { 0.00, 0.00, 0.00 },
                            { 0.75, 0.00, 0.25 }, { 0.50, 0.00, 0.50 }, { 0.25, 0.00, 0.75 },
                            { 0.00, 0.25, 0.75 }, { 0.00, 0.50, 0.50 }, { 0.00, 0.75, 0.25 },
                            { 0.25, 0.75, 0.00 }, { 0.50, 0.50, 0.00 }, { 0.75, 0.25, 0.00 },
                            { 0.25, 0.00, 0.00 }, { 0.50, 0.00, 0.00 }, { 0.75, 0.00, 0.00 },
                            { 0.00, 0.25, 0.00 }, { 0.00, 0.50, 0.00 }, { 0.00, 0.75, 0.00 },
                            { 0.00, 0.00, 0.25 }, { 0.00, 0.00, 0.50 }, { 0.00, 0.00, 0.75 },
                            { 0.50, 0.25, 0.25 }, { 0.25, 0.50, 0.25 }, { 0.25, 0.25, 0.50 },
                            { 0.50, 0.00, 0.25 }, { 0.25, 0.00, 0.50 }, { 0.25, 0.00, 0.25 },
                            { 0.50, 0.25, 0.00 }, { 0.25, 0.25, 0.00 }, { 0.25, 0.50, 0.00 },
                            { 0.00, 0.25, 0.25 }, { 0.00, 0.25, 0.50 }, { 0.00, 0.50, 0.25 },
                            { 0.25, 0.25, 0.25 }
                    };
                    load_tet( tXiHat, 35 );
                    break;
                }
                case( ElementType::PENTA6 ) :
                case( ElementType::PENTA6TS ) :
                {
                    const real tXiHat[][ 3 ] = {
                            { 1.0, 0.0, -1.0 }, { 0.0, 1.0, -1.0 }, { 0.0, 0.0, -1.0 },
                            { 1.0, 0.0,  1.0 }, { 0.0, 1.0,  1.0 }, { 0.0, 0.0,  1.0 }
                    };
                    load_penta( tXiHat, 6 );
                    break;
                }
                case( ElementType::PENTA15 ) :
                {
                    const real tXiHat[][ 3 ] = {
                            { 1.0, 0.0, -1.0 }, { 0.0, 1.0, -1.0 }, { 0.0, 0.0, -1.0 },
                            { 1.0, 0.0,  1.0 }, { 0.0, 1.0,  1.0 }, { 0.0, 0.0,  1.0 },
                            { 0.5, 0.5, -1.0 }, { 0.0, 0.5, -1.0 }, { 0.5, 0.0, -1.0 },
                            { 1.0, 0.0,  0.0 }, { 0.0, 1.0,  0.0 }, { 0.0, 0.0,  0.0 },
                            { 0.5, 0.5,  1.0 }, { 0.0, 0.5,  1.0 }, { 0.5, 0.0,  1.0 }
                    };
                    load_penta( tXiHat, 15 );
                    break;
                }
                case( ElementType::PENTA18 ) :
                case( ElementType::PENTA18TS ) :
                {
                    const real tXiHat[][ 3 ] = {
                            { 1.0, 0.0, -1.0 }, { 0.0, 1.0, -1.0 }, { 0.0, 0.0, -1.0 },
                            { 1.0, 0.0,  1.0 }, { 0.0, 1.0,  1.0 }, { 0.0, 0.0,  1.0 },
                            { 0.5, 0.5, -1.0 }, { 0.0, 0.5, -1.0 }, { 0.5, 0.0, -1.0 },
                            { 1.0, 0.0,  0.0 }, { 0.0, 1.0,  0.0 }, { 0.0, 0.0,  0.0 },
                            { 0.5, 0.5,  1.0 }, { 0.0, 0.5,  1.0 }, { 0.5, 0.0,  1.0 },
                            { 0.5, 0.5,  0.0 }, { 0.0, 0.5,  0.0 }, { 0.5, 0.0,  0.0 }
                    };
                    load_penta( tXiHat, 18 );
                    break;
                }
                case( ElementType::PYRA5 ) :
                {
                    const real tXiHat[][ 3 ] = {
                            { -1.0, -1.0, 0.0 }, {  1.0, -1.0, 0.0 }, {  1.0,  1.0, 0.0 },
                            { -1.0,  1.0, 0.0 }, {  0.0,  0.0, 1.0 }
                    };
                    load_pyra( tXiHat, 5 );
                    break;
                }
                case( ElementType::PYRA13 ) :
                {
                    const real tXiHat[][ 3 ] = {
                            { -1.0, -1.0, 0.0 }, {  1.0, -1.0, 0.0 }, {  1.0,  1.0, 0.0 },
                            { -1.0,  1.0, 0.0 }, {  0.0,  0.0, 1.0 },
                            {  0.0, -1.0, 0.0 }, {  1.0,  0.0, 0.0 }, {  0.0,  1.0, 0.0 }, { -1.0,  0.0, 0.0 },
                            { -0.5, -0.5, 0.5 }, {  0.5, -0.5, 0.5 }, {  0.5,  0.5, 0.5 }, { -0.5,  0.5, 0.5 }
                    };
                    load_pyra( tXiHat, 13 );
                    break;
                }
                case( ElementType::PYRA14 ) :
                {
                    const real tXiHat[][ 3 ] = {
                            { -1.0, -1.0, 0.0 }, {  1.0, -1.0, 0.0 }, {  1.0,  1.0, 0.0 },
                            { -1.0,  1.0, 0.0 }, {  0.0,  0.0, 1.0 },
                            {  0.0, -1.0, 0.0 }, {  1.0,  0.0, 0.0 }, {  0.0,  1.0, 0.0 }, { -1.0,  0.0, 0.0 },
                            { -0.5, -0.5, 0.5 }, {  0.5, -0.5, 0.5 }, {  0.5,  0.5, 0.5 }, { -0.5,  0.5, 0.5 },
                            {  0.0,  0.0, 0.0 }
                    };
                    load_pyra( tXiHat, 14 );
                    break;
                }
                case( ElementType::HEX8 ) :
                {
                    const real tXiHat[][ 3 ] = {
                            { -1.0, -1.0, -1.0 }, {  1.0, -1.0, -1.0 }, {  1.0,  1.0, -1.0 }, { -1.0,  1.0, -1.0 },
                            { -1.0, -1.0,  1.0 }, {  1.0, -1.0,  1.0 }, {  1.0,  1.0,  1.0 }, { -1.0,  1.0,  1.0 }
                    };
                    load_hex( tXiHat, 8 );
                    break;
                }
                case( ElementType::HEX20 ) :
                {
                    const real tXiHat[][ 3 ] = {
                            { -1.0, -1.0, -1.0 }, {  1.0, -1.0, -1.0 }, {  1.0,  1.0, -1.0 }, { -1.0,  1.0, -1.0 },
                            { -1.0, -1.0,  1.0 }, {  1.0, -1.0,  1.0 }, {  1.0,  1.0,  1.0 }, { -1.0,  1.0,  1.0 },
                            {  0.0, -1.0, -1.0 }, {  1.0,  0.0, -1.0 }, {  0.0,  1.0, -1.0 }, { -1.0,  0.0, -1.0 },
                            { -1.0, -1.0,  0.0 }, {  1.0, -1.0,  0.0 }, {  1.0,  1.0,  0.0 }, { -1.0,  1.0,  0.0 },
                            {  0.0, -1.0,  1.0 }, {  1.0,  0.0,  1.0 }, {  0.0,  1.0,  1.0 }, { -1.0,  0.0,  1.0 }
                    };
                    load_hex( tXiHat, 20 );
                    break;
                }
                case( ElementType::HEX27 ) :
                {
                    const real tXiHat[][ 3 ] = {
                            { -1.0, -1.0, -1.0 }, {  1.0, -1.0, -1.0 }, {  1.0,  1.0, -1.0 }, { -1.0,  1.0, -1.0 },
                            { -1.0, -1.0,  1.0 }, {  1.0, -1.0,  1.0 }, {  1.0,  1.0,  1.0 }, { -1.0,  1.0,  1.0 },
                            {  0.0, -1.0, -1.0 }, {  1.0,  0.0, -1.0 }, {  0.0,  1.0, -1.0 }, { -1.0,  0.0, -1.0 },
                            { -1.0, -1.0,  0.0 }, {  1.0, -1.0,  0.0 }, {  1.0,  1.0,  0.0 }, { -1.0,  1.0,  0.0 },
                            {  0.0, -1.0,  1.0 }, {  1.0,  0.0,  1.0 }, {  0.0,  1.0,  1.0 }, { -1.0,  0.0,  1.0 },
                            {  0.0,  0.0,  0.0 }, {  0.0,  0.0, -1.0 }, {  0.0,  0.0,  1.0 },
                            { -1.0,  0.0,  0.0 }, {  1.0,  0.0,  0.0 }, {  0.0, -1.0,  0.0 }, {  0.0,  1.0,  0.0 }
                    };
                    load_hex( tXiHat, 27 );
                    break;
                }
                case( ElementType::HEX64 ) :
                {
                    const real tXiHat[][ 3 ] = {
                            { -1.0,   -1.0,   -1.0 }, {  1.0,   -1.0,   -1.0 }, {  1.0,    1.0,   -1.0 }, { -1.0,    1.0,   -1.0 },
                            { -1.0,   -1.0,    1.0 }, {  1.0,   -1.0,    1.0 }, {  1.0,    1.0,    1.0 }, { -1.0,    1.0,    1.0 },
                            { -tThird, -1.0,   -1.0 }, {  tThird, -1.0,   -1.0 }, { -1.0, -tThird, -1.0 }, { -1.0,  tThird, -1.0 },
                            { -1.0,   -1.0, -tThird }, { -1.0,   -1.0,  tThird }, {  1.0, -tThird, -1.0 }, {  1.0,  tThird, -1.0 },
                            {  1.0,   -1.0, -tThird }, {  1.0,   -1.0,  tThird }, {  tThird,  1.0, -1.0 }, { -tThird,  1.0, -1.0 },
                            {  1.0,    1.0, -tThird }, {  1.0,    1.0,  tThird }, { -1.0,    1.0, -tThird }, { -1.0,    1.0,  tThird },
                            { -tThird, -1.0,  1.0 }, {  tThird, -1.0,  1.0 }, { -1.0, -tThird,  1.0 }, { -1.0,  tThird,  1.0 },
                            {  1.0, -tThird,  1.0 }, {  1.0,  tThird,  1.0 }, {  tThird,  1.0,  1.0 }, { -tThird,  1.0,  1.0 },
                            { -tThird, -tThird, -1.0 }, { -tThird,  tThird, -1.0 }, {  tThird,  tThird, -1.0 }, {  tThird, -tThird, -1.0 },
                            { -tThird, -1.0, -tThird }, {  tThird, -1.0, -tThird }, {  tThird, -1.0,  tThird }, { -tThird, -1.0,  tThird },
                            { -1.0, -tThird, -tThird }, { -1.0, -tThird,  tThird }, { -1.0,  tThird,  tThird }, { -1.0,  tThird, -tThird },
                            {  1.0, -tThird, -tThird }, {  1.0,  tThird, -tThird }, {  1.0,  tThird,  tThird }, {  1.0, -tThird,  tThird },
                            {  tThird,  1.0, -tThird }, { -tThird,  1.0, -tThird }, { -tThird,  1.0,  tThird }, {  tThird,  1.0,  tThird },
                            { -tThird, -tThird,  1.0 }, {  tThird, -tThird,  1.0 }, {  tThird,  tThird,  1.0 }, { -tThird,  tThird,  1.0 },
                            { -tThird, -tThird, -tThird }, {  tThird, -tThird, -tThird }, {  tThird,  tThird, -tThird }, { -tThird,  tThird, -tThird },
                            { -tThird, -tThird,  tThird }, {  tThird, -tThird,  tThird }, {  tThird,  tThird,  tThird }, { -tThird,  tThird,  tThird }
                    };
                    load_hex( tXiHat, 64 );
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false,
                        "create_unity_nodes not implemented for this element type" );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        ElementFactory::create_orientation_table ( const ElementType aType, Matrix< uint > & aTable ) const
        {
            // placeholder for void matrix entries
            const uint x = BELFEM_UINT_MAX ;

            switch ( aType )
            {
                case ElementType::TET4  :
                {
                    aTable = {
                        { 0, 1, 3, 1, 2, 3, 0, 3, 2, 0, 2, 1 },
                               { 3, 0, 1, 3, 1, 2, 2, 0, 3, 1, 0, 2 },
                               { 1, 3, 0, 2, 3, 1, 3, 2, 0, 2, 1, 0 }
                    };
                    break;
                }
                case ElementType::TET10 :
                {
                    aTable = {
                         { 0, 1, 3, 1, 2, 3, 0, 3, 2, 0, 2, 1 },
                                { 3, 0, 1, 3, 1, 2, 2, 0, 3, 1, 0, 2 },
                                { 1, 3, 0, 2, 3, 1, 3, 2, 0, 2, 1, 0 },
                                { 7, 4, 8, 8, 5, 9, 6, 7, 9, 4, 6, 5 },
                                { 8, 7, 4, 9, 8, 5, 9, 6, 7, 5, 4, 6 },
                                { 4, 8, 7, 5, 9, 8, 7, 9, 6, 6, 5, 4 } };
                    break ;
                }
                case ElementType::TET20 :
                {
                    aTable = {
                        {  0,  1,  3,  1,  2,  3,  0,  3,  2,  0,  2,  1 },
                        {  3,  0,  1,  3,  1,  2,  2,  0,  3,  1,  0,  2 },
                        {  1,  3,  0,  2,  3,  1,  3,  2,  0,  2,  1,  0 },
                        { 11,  5, 14, 15,  7, 12,  9, 10, 13,  4,  8,  6 },
                        { 10,  4, 15, 14,  6, 13,  8, 11, 12,  5,  9,  7 },
                        { 14, 11,  5, 12, 15,  7, 13,  9, 10,  6,  4,  8 },
                        { 15, 10,  4, 13, 14,  6, 12,  8, 11,  7,  5,  9 },
                        {  5, 14, 11,  7, 12, 15, 10, 13,  9,  8,  6,  4 },
                        {  4, 15, 10,  6, 13, 14, 11, 12,  8,  9,  7,  5 },
                        { 17, 17, 17, 19, 19, 19, 18, 18, 18, 16, 16, 16 } };
                    break ;
                }
                case ElementType::TET35 :
                {
                    aTable = {
                        {  0,  1,  3,  1,  2,  3,  0,  3,  2,  0,  2,  1 },
                        {  3,  0,  1,  3,  1,  2,  2,  0,  3,  1,  0,  2 },
                        {  1,  3,  0,  2,  3,  1,  3,  2,  0,  2,  1,  0 },
                        { 15,  6, 19, 21,  9, 16, 12, 13, 18,  4, 10,  7 },
                        { 14,  5, 20, 20,  8, 17, 11, 14, 17,  5, 11,  8 },
                        { 13,  4, 21, 19,  7, 18, 10, 15, 16,  6, 12,  9 },
                        { 19, 15,  6, 16, 21,  9, 18, 12, 13,  7,  4, 10 },
                        { 20, 14,  5, 17, 20,  8, 17, 11, 14,  8,  5, 11 },
                        { 21, 13,  4, 18, 19,  7, 16, 10, 15,  9,  6, 12 },
                        {  6, 19, 15,  9, 16, 21, 13, 18, 12, 10,  7,  4 },
                        {  5, 20, 14,  8, 17, 20, 14, 17, 11, 11,  8,  5 },
                        {  4, 21, 13,  7, 18, 19, 15, 16, 10, 12,  9,  6 },
                        { 25, 26, 27, 32, 33, 31, 28, 29, 30, 22, 23, 24 },
                        { 27, 25, 26, 31, 32, 33, 30, 28, 29, 24, 22, 23 },
                        { 26, 27, 25, 33, 31, 32, 29, 30, 28, 23, 24, 22 } };
                    break ;
                }
                case ElementType::HEX8 :
                {
                    aTable = {
                        {  0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 0, 4, 7, 3, 0, 3, 2, 1, 4, 5, 6, 7 },
                         {  4, 0, 1, 5, 5, 1, 2, 6, 6, 2, 3, 7, 3, 0, 4, 7, 1, 0, 3, 2, 7, 4, 5, 6 },
                         {  5, 4, 0, 1, 6, 5, 1, 2, 7, 6, 2, 3, 7, 3, 0, 4, 2, 1, 0, 3, 6, 7, 4, 5 },
                         {  1, 5, 4, 0, 2, 6, 5, 1, 3, 7, 6, 2, 4, 7, 3, 0, 3, 2, 1, 0, 5, 6, 7, 4 } };
                    break ;
                }
                case ElementType::HEX20 :
                {
                    aTable = {
                        {  0,  1,  5,  4,  1,  2,  6,  5,  2,  3,  7,  6,  0,  4,  7,  3,  0,  3,  2,  1,  4,  5,  6,  7 },
                        {  4,  0,  1,  5,  5,  1,  2,  6,  6,  2,  3,  7,  3,  0,  4,  7,  1,  0,  3,  2,  7,  4,  5,  6 },
                        {  5,  4,  0,  1,  6,  5,  1,  2,  7,  6,  2,  3,  7,  3,  0,  4,  2,  1,  0,  3,  6,  7,  4,  5 },
                        {  1,  5,  4,  0,  2,  6,  5,  1,  3,  7,  6,  2,  4,  7,  3,  0,  3,  2,  1,  0,  5,  6,  7,  4 },
                        { 12,  8, 13, 16, 13,  9, 14, 17, 14, 10, 15, 18, 11, 12, 19, 15,  8, 11, 10,  9, 19, 16, 17, 18 },
                        { 16, 12,  8, 13, 17, 13,  9, 14, 18, 14, 10, 15, 15, 11, 12, 19,  9,  8, 11, 10, 18, 19, 16, 17 },
                        { 13, 16, 12,  8, 14, 17, 13,  9, 15, 18, 14, 10, 19, 15, 11, 12, 10,  9,  8, 11, 17, 18, 19, 16 },
                        {  8, 13, 16, 12,  9, 14, 17, 13, 10, 15, 18, 14, 12, 19, 15, 11, 11, 10,  9,  8, 16, 17, 18, 19 },
                        {  x,  x,  x,  x,  x,  x,  x,  x,  x,  x,  x,  x,  x,  x,  x,  x,  x,  x,  x,  x,  x,  x,  x,  x } };
                    break;
                }
                case ElementType::HEX27 :
                {
                    aTable = {
                        {  0,  1,  5,  4,  1,  2,  6,  5,  2,  3,  7,  6,  0,  4,  7,  3,  0,  3,  2,  1,  4,  5,  6,  7 },
                        {  4,  0,  1,  5,  5,  1,  2,  6,  6,  2,  3,  7,  3,  0,  4,  7,  1,  0,  3,  2,  7,  4,  5,  6 },
                        {  5,  4,  0,  1,  6,  5,  1,  2,  7,  6,  2,  3,  7,  3,  0,  4,  2,  1,  0,  3,  6,  7,  4,  5 },
                        {  1,  5,  4,  0,  2,  6,  5,  1,  3,  7,  6,  2,  4,  7,  3,  0,  3,  2,  1,  0,  5,  6,  7,  4 },
                        { 12,  8, 13, 16, 13,  9, 14, 17, 14, 10, 15, 18, 11, 12, 19, 15,  8, 11, 10,  9, 19, 16, 17, 18 },
                        { 16, 12,  8, 13, 17, 13,  9, 14, 18, 14, 10, 15, 15, 11, 12, 19,  9,  8, 11, 10, 18, 19, 16, 17 },
                        { 13, 16, 12,  8, 14, 17, 13,  9, 15, 18, 14, 10, 19, 15, 11, 12, 10,  9,  8, 11, 17, 18, 19, 16 },
                        {  8, 13, 16, 12,  9, 14, 17, 13, 10, 15, 18, 14, 12, 19, 15, 11, 11, 10,  9,  8, 16, 17, 18, 19 },
                        { 25, 25, 25, 25, 24, 24, 24, 24, 26, 26, 26, 26, 23, 23, 23, 23, 21, 21, 21, 21, 22, 22, 22, 22 } };
                    break;
                }
                case ElementType::HEX64 :
                {
                    aTable = {
                        {  0,  1,  5,  4,  1,  2,  6,  5,  2,  3,  7,  6,  0,  4,  7,  3,  0,  3,  2,  1,  4,  5,  6,  7 },
                        {  4,  0,  1,  5,  5,  1,  2,  6,  6,  2,  3,  7,  3,  0,  4,  7,  1,  0,  3,  2,  7,  4,  5,  6 },
                        {  5,  4,  0,  1,  6,  5,  1,  2,  7,  6,  2,  3,  7,  3,  0,  4,  2,  1,  0,  3,  6,  7,  4,  5 },
                        {  1,  5,  4,  0,  2,  6,  5,  1,  3,  7,  6,  2,  4,  7,  3,  0,  3,  2,  1,  0,  5,  6,  7,  4 },
                        { 12,  9, 17, 24, 16, 15, 21, 28, 20, 19, 23, 30, 10, 13, 27, 22,  8, 11, 18, 14, 26, 25, 29, 31 },
                        { 13,  8, 16, 25, 17, 14, 20, 29, 21, 18, 22, 31, 11, 12, 26, 23,  9, 10, 19, 15, 27, 24, 28, 30 },
                        { 24, 12,  9, 17, 28, 16, 15, 21, 30, 20, 19, 23, 22, 10, 13, 27, 14,  8, 11, 18, 31, 26, 25, 29 },
                        { 25, 13,  8, 16, 29, 17, 14, 20, 31, 21, 18, 22, 23, 11, 12, 26, 15,  9, 10, 19, 30, 27, 24, 28 },
                        { 17, 24, 12,  9, 21, 28, 16, 15, 23, 30, 20, 19, 27, 22, 10, 13, 18, 14,  8, 11, 29, 31, 26, 25 },
                        { 16, 25, 13,  8, 20, 29, 17, 14, 22, 31, 21, 18, 26, 23, 11, 12, 19, 15,  9, 10, 28, 30, 27, 24 },
                        {  9, 17, 24, 12, 15, 21, 28, 16, 19, 23, 30, 20, 13, 27, 22, 10, 11, 18, 14,  8, 25, 29, 31, 26 },
                        {  8, 16, 25, 13, 14, 20, 29, 17, 18, 22, 31, 21, 12, 26, 23, 11, 10, 19, 15,  9, 24, 28, 30, 27 },
                        { 36, 37, 38, 39, 44, 45, 46, 47, 48, 49, 50, 51, 40, 41, 42, 43, 32, 33, 34, 35, 52, 53, 54, 55 },
                        { 39, 36, 37, 38, 47, 44, 45, 46, 51, 48, 49, 50, 43, 40, 41, 42, 35, 32, 33, 34, 55, 52, 53, 54 },
                        { 38, 39, 36, 37, 46, 47, 44, 45, 50, 51, 48, 49, 42, 43, 40, 41, 34, 35, 32, 33, 54, 55, 52, 53 },
                        { 37, 38, 39, 36, 45, 46, 47, 44, 49, 50, 51, 48, 41, 42, 43, 40, 33, 34, 35, 32, 53, 54, 55, 52 } };
                    break;
                }
                case ElementType::PENTA6 :
                case ElementType::PENTA6TS :
                {
                    // Option-B canonicalization: PENTA6TS shares the volume
                    // PENTA6 facet topology and orientation layout
                    // (4+4+4+3+3 = 18 columns).
                    aTable = { { 0, 1, 4, 3, 1, 2, 5, 4, 0, 3, 5, 2, 0, 1, 2, 3, 4, 5 },
                               { 3, 0, 1, 4, 4, 1, 2, 5, 2, 5, 2, 5, 1, 2, 0, 5, 3, 4 },
                               { 4, 3, 0, 1, 5, 4, 1, 2, 5, 2, 0, 3, 2, 0, 1, 4, 5, 3 },
                               { 1, 4, 3, 0, 2, 5, 4, 1, 3, 0, 3, 0, x, x, x, x, x, x } };
                    break ;
                }
                case ElementType::PENTA15 :
                {
                    aTable = {
                        {  0,  1,  4,  3,  1,  2,  5,  4,  0,  3,  5,  2,  0,  1,  2,  3,  4,  5 },
                        {  3,  0,  1,  4,  4,  1,  2,  5,  2,  5,  2,  5,  1,  2,  0,  5,  3,  4 },
                        {  4,  3,  0,  1,  5,  4,  1,  2,  5,  2,  0,  3,  2,  0,  1,  4,  5,  3 },
                        {  1,  4,  3,  0,  2,  5,  4,  1,  3,  0,  3,  0,  x,  x,  x,  x,  x,  x },
                        {  9,  6, 10, 12, 10,  7, 11, 13,  8, 14, 11, 11,  6,  7,  8, 14, 12, 13 },
                        { 12,  9,  6, 10, 13, 10,  7, 11, 11, 11,  8, 14,  7,  8,  6, 13, 14, 12 },
                        { 10, 12,  9,  6, 11, 13, 10,  7, 14,  8,  9,  9,  8,  6,  7, 12, 13, 14 },
                        {  6, 10, 12,  9,  7, 11, 13, 10,  9,  9, 14,  8,  x,  x,  x,  x,  x,  x } };
                    break ;
                }
                case ElementType::PENTA18 :
                case ElementType::PENTA18TS :
                {
                    // Option-B canonicalization: PENTA18TS shares the volume
                    // PENTA18 facet topology and orientation layout
                    // (4+4+4+3+3 = 18 columns).
                    aTable = {
                        {  0,  1,  4,  3,  1,  2,  5,  4,  0,  3,  5,  2,  0,  1,  2,  3,  4,  5 },
                        {  3,  0,  1,  4,  4,  1,  2,  5,  2,  5,  2,  5,  1,  2,  0,  5,  3,  4 },
                        {  4,  3,  0,  1,  5,  4,  1,  2,  5,  2,  0,  3,  2,  0,  1,  4,  5,  3 },
                        {  1,  4,  3,  0,  2,  5,  4,  1,  3,  0,  3,  0,  x,  x,  x,  x,  x,  x },
                        {  9,  6, 10, 12, 10,  7, 11, 13,  8, 14, 11, 11,  6,  7,  8, 14, 12, 13 },
                        { 12,  9,  6, 10, 13, 10,  7, 11, 11, 11,  8, 14,  7,  8,  6, 13, 14, 12 },
                        { 10, 12,  9,  6, 11, 13, 10,  7, 14,  8,  9,  9,  8,  6,  7, 12, 13, 14 },
                        {  6, 10, 12,  9,  7, 11, 13, 10,  9,  9, 14,  8,  x,  x,  x,  x,  x,  x },
                        { 15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 17,  x,  x,  x,  x,  x,  x } };
                    break ;
                }
                case ElementType::PYRA5 :
                {
                    aTable =
                        {
                        { 0,1,4,1,2,4,2,3,4,3,0,4,0,1,2,3 },
                            { 4,0,1,4,1,2,4,2,3,4,3,0,1,2,3,0 },
                            { 1,4,0,2,4,1,3,4,2,0,4,3,2,3,0,1 },
                            { x,x,x,x,x,x,x,x,x,x,x,x,3,0,1,2 }
                                                    };
                    break;
                }
                case ElementType::PYRA13 :
                {
                    aTable = { { 0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0, 4, 0, 1, 2, 3 },
                                { 4, 0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0, 1, 2, 3, 0 },
                                { 1, 4, 0, 2, 4, 1, 3, 4, 2, 0, 4, 3, 2, 3, 0, 1 },
                                { 9, 5,10,10, 6,11,11, 7,12,12, 8, 9, 3, 0, 1, 2 },
                                {10, 9, 5,11,10, 6,12,11, 7, 9,12, 8, 5, 6, 7, 8 },
                                { 5,10, 9, 6,11,10, 7,12,11, 8, 9,12, 6, 7, 8, 5 },
                                { x, x, x, x, x, x, x, x, x, x, x, x, 7, 8, 5, 6 },
                                { x, x, x, x, x, x, x, x, x, x, x, x, 8, 5, 6, 7 }
                    };
                    break ;
                }
                case ElementType::PYRA14 :
                {
                    aTable = { { 0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0, 4, 0, 1, 2, 3 },
                                { 4, 0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0, 1, 2, 3, 0 },
                                { 1, 4, 0, 2, 4, 1, 3, 4, 2, 0, 4, 3, 2, 3, 0, 1 },
                                { 9, 5,10,10, 6,11,11, 7,12,12, 8, 9, 3, 0, 1, 2 },
                                {10, 9, 5,11,10, 6,12,11, 7, 9,12, 8, 5, 6, 7, 8 },
                                { 5,10, 9, 6,11,10, 7,12,11, 8, 9,12, 6, 7, 8, 5 },
                                { x, x, x, x, x, x, x, x, x, x, x, x, 7, 8, 5, 6 },
                                { x, x, x, x, x, x, x, x, x, x, x, x, 8, 5, 6, 7 },
                                { x, x, x, x, x, x, x, x, x, x, x, x,13,13,13,13 }
                    };
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false,
                        "create_orientation_table not implemented for this element type" );
                }
            }
        }


    }
 }
