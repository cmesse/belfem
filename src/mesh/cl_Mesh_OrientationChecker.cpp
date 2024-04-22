//
// Created by christian on 9/2/22.
//
#include "constants.hpp"
#include "cl_Mesh_OrientationChecker.hpp"
#include "meshtools.hpp"
#include "fn_sum.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"
#include "fn_rotation_matrix.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        OrientationChecker::OrientationChecker()
        {
            // normal vector of first plane
            mN.set_size( 3 );

            // rotation matrix
            mT.set_size( 3, 3 );

            // rotation axis
            mR.set_size( 3, 0.0 );

            // node coordinates
            mP.set_size( 3 );
            mQ.set_size( 3 );
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::set_element_type( const ElementType aType )
        {
            // remember the type
            mElementType = aType ;

            mNumNodesPerElement = number_of_nodes(aType  );

            mNumCornerNodesPerElement
                = number_of_nodes( linear_element_type(
                    aType ) );

            // allocate containers
            mX.set_size( mNumCornerNodesPerElement );
            mY.set_size( mNumCornerNodesPerElement );
            mZ.set_size( mNumCornerNodesPerElement );

            switch( aType )
            {
                case( ElementType::TRI3 ) :
                {
                    mFunCollect = & OrientationChecker::collect_node_coords_tri ;
                    mFunProcess = & OrientationChecker::process_2d_element ;
                    mFunFlip    = & OrientationChecker::flip_tri3 ;
                    break;
                }
                case( ElementType::TRI6 ) :
                {
                    mFunCollect = & OrientationChecker::collect_node_coords_tri ;
                    mFunProcess = & OrientationChecker::process_2d_element ;
                    mFunFlip    = & OrientationChecker::flip_tri6 ;
                    break;
                }
                case( ElementType::TRI10 ) :
                {
                    mFunCollect = & OrientationChecker::collect_node_coords_tri ;
                    mFunProcess = & OrientationChecker::process_2d_element ;
                    mFunFlip    = & OrientationChecker::flip_tri10 ;
                    break;
                }
                case( ElementType::TRI15 ) :
                {
                    mFunCollect = & OrientationChecker::collect_node_coords_tri ;
                    mFunProcess = & OrientationChecker::process_2d_element ;
                    mFunFlip    = & OrientationChecker::flip_tri15 ;
                    break;
                }
                case( ElementType::QUAD4 ) :
                {
                    mFunCollect = & OrientationChecker::collect_node_coords_quad ;
                    mFunProcess = & OrientationChecker::process_2d_element ;
                    mFunFlip    = & OrientationChecker::flip_quad4 ;
                    break;
                }
                case( ElementType::QUAD8 ) :
                case( ElementType::QUAD9 ) :
                {
                    mFunCollect = & OrientationChecker::collect_node_coords_quad ;
                    mFunProcess = & OrientationChecker::process_2d_element ;
                    mFunFlip    = & OrientationChecker::flip_quad9 ;
                    break;
                }
                case( ElementType::QUAD16 ) :
                {
                    mFunCollect = & OrientationChecker::collect_node_coords_quad ;
                    mFunProcess = & OrientationChecker::process_2d_element ;
                    mFunFlip    = & OrientationChecker::flip_quad16 ;
                    break;
                }
                case( ElementType::TET4 ) :
                {
                    mFunCollect = & OrientationChecker::collect_node_coords_tet ;
                    mFunProcess = & OrientationChecker::process_tet_element ;
                    mFunFlip    = & OrientationChecker::flip_tet4 ;
                    break ;
                }
                case( ElementType::TET10 ) :
                {
                    mFunCollect = & OrientationChecker::collect_node_coords_tet ;
                    mFunProcess = & OrientationChecker::process_tet_element ;
                    mFunFlip    = & OrientationChecker::flip_tet10 ;
                    break ;
                }
                case( ElementType::TET20 ) :
                {
                    mFunCollect = & OrientationChecker::collect_node_coords_tet ;
                    mFunProcess = & OrientationChecker::process_tet_element ;
                    mFunFlip    = & OrientationChecker::flip_tet20 ;
                    break ;
                }
                case( ElementType::TET35 ) :
                {
                    mFunCollect = & OrientationChecker::collect_node_coords_tet ;
                    mFunProcess = & OrientationChecker::process_tet_element ;
                    mFunFlip    = & OrientationChecker::flip_tet35 ;
                    break ;
                }
                case( ElementType::PENTA6 ) :
                {
                    mFunCollect = & OrientationChecker::collect_node_coords_penta ;
                    mFunProcess = & OrientationChecker::process_penta_element ;
                    mFunFlip    = & OrientationChecker::flip_penta6 ;
                    break ;
                }
                case( ElementType::PENTA15 ) :
                {
                    mFunCollect = & OrientationChecker::collect_node_coords_penta ;
                    mFunProcess = & OrientationChecker::process_penta_element ;
                    mFunFlip    = & OrientationChecker::flip_penta15 ;
                    break ;
                }
                case( ElementType::PENTA18 ) :
                {
                    mFunCollect = & OrientationChecker::collect_node_coords_penta ;
                    mFunProcess = & OrientationChecker::process_penta_element ;
                    mFunFlip    = & OrientationChecker::flip_penta18 ;
                    break ;
                }
                case( ElementType::HEX8 ) :
                {
                    mFunCollect = & OrientationChecker::collect_node_coords_hex ;
                    mFunProcess = & OrientationChecker::process_hex_element ;
                    mFunFlip    = & OrientationChecker::flip_hex8 ;
                    break ;
                }
                case( ElementType::HEX20 ) :
                {
                    mFunCollect = & OrientationChecker::collect_node_coords_hex ;
                    mFunProcess  = & OrientationChecker::process_hex_element ;
                    mFunFlip    = & OrientationChecker::flip_hex20 ;
                    break ;
                }
                case( ElementType::HEX27 ) :
                {
                    mFunCollect = & OrientationChecker::collect_node_coords_hex ;
                    mFunProcess = & OrientationChecker::process_hex_element ;
                    mFunFlip    = & OrientationChecker::flip_hex27 ;
                    break ;
                }
                case( ElementType::HEX64 ) :
                {
                    mFunCollect = & OrientationChecker::collect_node_coords_hex ;
                    mFunProcess = & OrientationChecker::process_hex_element ;
                    mFunFlip    = & OrientationChecker::flip_hex64 ;
                    break ;
                }
                default :
                {
                    BELFEM_ERROR( false, "Element Type not supported in orientation checker");
                }
            }
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::collect_node_coords_tri( Element * aElement )
        {
            mX( 0 ) = aElement->node( 0 )->x() ;
            mX( 1 ) = aElement->node( 1 )->x() ;
            mX( 2 ) = aElement->node( 2 )->x() ;

            mY( 0 ) = aElement->node( 0 )->y() ;
            mY( 1 ) = aElement->node( 1 )->y() ;
            mY( 2 ) = aElement->node( 2 )->y() ;
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::collect_node_coords_quad( Element * aElement )
        {
            mX( 0 ) = aElement->node( 0 )->x() ;
            mX( 1 ) = aElement->node( 1 )->x() ;
            mX( 2 ) = aElement->node( 2 )->x() ;
            mX( 3 ) = aElement->node( 3 )->x() ;

            mY( 0 ) = aElement->node( 0 )->y() ;
            mY( 1 ) = aElement->node( 1 )->y() ;
            mY( 2 ) = aElement->node( 2 )->y() ;
            mY( 3 ) = aElement->node( 3 )->y() ;
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::collect_node_coords_tet( Element * aElement )
        {
            mX( 0 ) = aElement->node( 0 )->x() ;
            mX( 1 ) = aElement->node( 1 )->x() ;
            mX( 2 ) = aElement->node( 2 )->x() ;
            mX( 3 ) = aElement->node( 3 )->x() ;

            mY( 0 ) = aElement->node( 0 )->y() ;
            mY( 1 ) = aElement->node( 1 )->y() ;
            mY( 2 ) = aElement->node( 2 )->y() ;
            mY( 3 ) = aElement->node( 3 )->y() ;

            mZ( 0 ) = aElement->node( 0 )->z() ;
            mZ( 1 ) = aElement->node( 1 )->z() ;
            mZ( 2 ) = aElement->node( 2 )->z() ;
            mZ( 3 ) = aElement->node( 3 )->z() ;
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::collect_node_coords_penta( Element * aElement )
        {
            mX( 0 ) = aElement->node( 0 )->x() ;
            mX( 1 ) = aElement->node( 1 )->x() ;
            mX( 2 ) = aElement->node( 2 )->x() ;
            mX( 3 ) = aElement->node( 3 )->x() ;
            mX( 4 ) = aElement->node( 4 )->x() ;
            mX( 5 ) = aElement->node( 5 )->x() ;


            mY( 0 ) = aElement->node( 0 )->y() ;
            mY( 1 ) = aElement->node( 1 )->y() ;
            mY( 2 ) = aElement->node( 2 )->y() ;
            mY( 3 ) = aElement->node( 3 )->y() ;
            mY( 4 ) = aElement->node( 4 )->y() ;
            mY( 5 ) = aElement->node( 5 )->y() ;

            mZ( 0 ) = aElement->node( 0 )->z() ;
            mZ( 1 ) = aElement->node( 1 )->z() ;
            mZ( 2 ) = aElement->node( 2 )->z() ;
            mZ( 3 ) = aElement->node( 3 )->z() ;
            mZ( 4 ) = aElement->node( 4 )->z() ;
            mZ( 5 ) = aElement->node( 5 )->z() ;
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::collect_node_coords_hex( Element * aElement )
        {
            mX( 0 ) = aElement->node( 0 )->x() ;
            mX( 1 ) = aElement->node( 1 )->x() ;
            mX( 2 ) = aElement->node( 2 )->x() ;
            mX( 3 ) = aElement->node( 3 )->x() ;
            mX( 4 ) = aElement->node( 4 )->x() ;
            mX( 5 ) = aElement->node( 5 )->x() ;
            mX( 6 ) = aElement->node( 6 )->x() ;
            mX( 7 ) = aElement->node( 7 )->x() ;


            mY( 0 ) = aElement->node( 0 )->y() ;
            mY( 1 ) = aElement->node( 1 )->y() ;
            mY( 2 ) = aElement->node( 2 )->y() ;
            mY( 3 ) = aElement->node( 3 )->y() ;
            mY( 4 ) = aElement->node( 4 )->y() ;
            mY( 5 ) = aElement->node( 5 )->y() ;
            mY( 6 ) = aElement->node( 6 )->y() ;
            mY( 7 ) = aElement->node( 7 )->y() ;

            mZ( 0 ) = aElement->node( 0 )->z() ;
            mZ( 1 ) = aElement->node( 1 )->z() ;
            mZ( 2 ) = aElement->node( 2 )->z() ;
            mZ( 3 ) = aElement->node( 3 )->z() ;
            mZ( 4 ) = aElement->node( 4 )->z() ;
            mZ( 5 ) = aElement->node( 5 )->z() ;
            mZ( 6 ) = aElement->node( 6 )->z() ;
            mZ( 7 ) = aElement->node( 7 )->z() ;
        }


//-----------------------------------------------------------------------------

        real
        OrientationChecker::compute_phi1( const uint aNumNodes )
        {
            real tXm = 0.0 ;
            real tYm = 0.0 ;
            for( uint k=0; k<aNumNodes; ++k )
            {
                tXm += mX( k );
                tYm += mY( k );
            }
            tXm /= aNumNodes ;
            tYm /= aNumNodes ;


            // polar coordinate of node 0
            real tPhi0 = std::atan2( mY( 0 )-tYm, mX( 0 )-tXm );

            // rotate node coordinates
            real tS = std::sin( tPhi0 );
            real tC = std::cos( tPhi0 );
            real tX1 =   tC * ( mX(1 ) - tXm ) + tS * ( mY( 1 ) - tYm ) ;
            real tY1 = - tS * ( mX(1 ) - tXm ) + tC * ( mY( 1 ) - tYm ) ;

            return std::atan2( tY1, tX1 ) ;
        }

//------------------------------------------------------------------------------

        index_t
        OrientationChecker::process_2d_element( Element * aElement )
        {
            // collect the node coordinates
            ( this->*mFunCollect )( aElement );
            // check if element needs to be flipped
            if( this->compute_phi1( mNumCornerNodesPerElement ) < 0 )
            {
                // flip the element
                ( this->*mFunFlip )( aElement );

                return 1;
            }
            else
            {
                return 0;
            }
        }

//------------------------------------------------------------------------------

        index_t
        OrientationChecker::process_tet_element( Element * aElement )
        {
            this->collect_node_coords_tet( aElement );

            // compute center of tet
            real tXm = sum( mX ) * 0.25 ;
            real tYm = sum( mY ) * 0.25 ;
            real tZm = sum( mZ ) * 0.25 ;

            // shift node coords
            mX -= tXm ;
            mY -= tYm ;
            mZ -= tZm ;

            // normal of first triangle plane
            mN( 0 ) = (mY(1)-mY(0))*(mZ(2)-mZ(0))
                            -(mY(2)-mY(0))*(mZ(1)-mZ(0));
            mN( 1 ) = (mX(2)-mX(0))*(mZ(1)-mZ(0))
                            -(mX(1)-mX(0))*(mZ(2)-mZ(0));
            mN( 2 ) = (mX(1)-mX(0))*(mY(2)-mY(0))
                            -(mX(2)-mX(0))*(mY(1)-mY(0));
            mN /= norm( mN );


            // rotation axis
            mR( 0 ) =  mN( 1 );
            mR( 1 ) = -mN( 0 );

            // rotation matrix to turn the first plane
            rotation_matrix( mR, std::atan2( 1.0, mN( 2 ) ), mT );

            // rotate the node coords so that the first triangle is in the x-y plane
            for( uint k=0; k<4; ++k )
            {
                mP( 0 ) = mX( k );
                mP( 1 ) = mY( k );
                mP( 2 ) = mZ( k );
                mQ = mT * mP ;
                mX( k ) = mQ( 0 );
                mY( k ) = mQ( 1 );
                mZ( k ) = mQ( 2 );
            }
            mZ -= mZ( 0 );

            // sanity check
            BELFEM_ASSERT( std::abs( mZ( 0 ) ) < BELFEM_MESH_EPSILON &&
                           std::abs( mZ( 1 ) ) < BELFEM_MESH_EPSILON &&
                           std::abs( mZ( 2 ) ) < BELFEM_MESH_EPSILON,
                          "error while trying to check tet element" );

            // check if element needs to be flipped
            if( this->compute_phi1( 3 ) * mZ( 3 ) < 0 )
            {
                // flip the element
                ( this->*mFunFlip )( aElement );

                return 1;
            }
            else
            {
                return 0;
            }
        }

//------------------------------------------------------------------------------

        index_t
        OrientationChecker::process_penta_element( Element * aElement )
        {
            this->collect_node_coords_penta( aElement );

            // compute center of penta
            real tXm = sum( mX ) / 6.0 ;
            real tYm = sum( mY ) / 6.0 ;
            real tZm = sum( mZ ) / 6.0 ;

            // shift node coords
            mX -= tXm ;
            mY -= tYm ;
            mZ -= tZm ;

            // normal of first triangle plane
            mN( 0 ) = (mY(1)-mY(0))*(mZ(2)-mZ(0))
                      -(mY(2)-mY(0))*(mZ(1)-mZ(0));
            mN( 1 ) = (mX(2)-mX(0))*(mZ(1)-mZ(0))
                      -(mX(1)-mX(0))*(mZ(2)-mZ(0));
            mN( 2 ) = (mX(1)-mX(0))*(mY(2)-mY(0))
                      -(mX(2)-mX(0))*(mY(1)-mY(0));
            mN /= norm( mN );


            // rotation axis
            mR( 0 ) =  mN( 1 );
            mR( 1 ) = -mN( 0 );

            // rotation matrix to turn the first plane
            rotation_matrix( mR, std::atan2( 1.0, mN( 2 ) ), mT );

            // rotate the node coords so that the first triangle is in the x-y plane
            for( uint k=0; k<6; ++k )
            {
                mP( 0 ) = mX( k );
                mP( 1 ) = mY( k );
                mP( 2 ) = mZ( k );
                mQ = mT * mP ;
                mX( k ) = mQ( 0 );
                mY( k ) = mQ( 1 );
                mZ( k ) = mQ( 2 );
            }
            mZ -= mZ( 0 );

            // sanity check
            BELFEM_ASSERT( std::abs( mZ( 0 ) ) < BELFEM_MESH_EPSILON &&
                           std::abs( mZ( 1 ) ) < BELFEM_MESH_EPSILON &&
                           std::abs( mZ( 2 ) ) < BELFEM_MESH_EPSILON,
                           "error while trying to check tet element" );


            // now we compute the center of the other triangle
            tZm = ( mZ( 0 ) + mZ( 1 ) + mZ( 2 ) ) / 3.0;

            // check if element needs to be flipped
            if( this->compute_phi1( 3 ) * tZm < 0 )
            {
                // flip the element
                ( this->*mFunFlip )( aElement );
                return 1;
            }
            else
            {
                return 0;
            }
        }
//------------------------------------------------------------------------------

        index_t
        OrientationChecker::process_hex_element( Element * aElement )
        {
            this->collect_node_coords_hex( aElement );

            // compute center of the hex
            real tXm = sum( mX ) * 0.125 ;
            real tYm = sum( mY ) * 0.125 ;
            real tZm = sum( mZ ) * 0.125 ;

            // shift node coords
            mX -= tXm ;
            mY -= tYm ;
            mZ -= tZm ;

            // compute the center of the first plane
            tXm = 0.25 * ( mX( 0 ) + mX( 1 ) + mX( 2 ) + mX( 3 ) );
            tYm = 0.25 * ( mY( 0 ) + mY( 1 ) + mY( 2 ) + mY( 3 ) );
            tZm = 0.25 * ( mZ( 0 ) + mZ( 1 ) + mZ( 2 ) + mZ( 3 ) );


            // compute an approximated normal for the first plane
            mN( 0 ) = (mY(2)-mY(0))*(mZ(3)-mZ(1))
                            -(mY(3)-mY(1))*(mZ(2)-mZ(0));

            mN( 1 ) = (mX(2)-mX(1))*(mZ(2)-mZ(0))
                            -(mX(2)-mX(0))*(mZ(3)-mZ(1));

            mN( 2 ) = (mX(2)-mX(0))*(mY(3)-mY(1))
                            -(mX(2)-mX(1))*(mY(2)-mY(0));
            mN /= norm( mN );

            // rotation axis
            mR( 0 ) =  mN( 1 );
            mR( 1 ) = -mN( 0 );

            // rotation matrix to turn the first plane
            rotation_matrix( mR, std::atan2( 1.0, mN( 2 ) ), mT );

            // rotate the node coords so that the first triangle is in the x-y plane
            for( uint k=0; k<8; ++k )
            {
                mP( 0 ) = mX( k );
                mP( 1 ) = mY( k );
                mP( 2 ) = mZ( k );
                mQ = mT * mP ;
                mX( k ) = mQ( 0 );
                mY( k ) = mQ( 1 );
                mZ( k ) = mQ( 2 );
            }

            // 4 * center of second surface
            tZm = mZ( 4 ) + mZ( 5 ) + mZ( 6 ) + mZ( 7 );

            // 4 * center of second surface
            tZm -= mZ( 0 ) + mZ( 1 ) + mZ( 2 ) + mZ( 3 );

            // check if element needs to be flipped
            if( this->compute_phi1( 4 ) * tZm < 0 )
            {
                // flip the element
                ( this->*mFunFlip )( aElement );

                return 1;
            }
            else
            {
                return 0;
            }
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::flip_tri3( Element * aElement )
        {
            // flip nodes
            this->swap( aElement, 1, 2 );
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::flip_tri6( Element * aElement )
        {
            // flip nodes
            this->swap( aElement, 1, 2 );
            this->swap( aElement, 3, 5 );
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::flip_tri10( Element * aElement )
        {
            // flip nodes
            this->swap( aElement, 1, 2 );
            this->swap( aElement, 3, 8 );
            this->swap( aElement, 4, 7 );
            this->swap( aElement, 5, 6 );
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::flip_tri15( Element * aElement )
        {
            // flip nodes
            this->swap( aElement,  1,  2 );
            this->swap( aElement,  3, 11 );
            this->swap( aElement,  4, 10 );
            this->swap( aElement,  5,  9 );
            this->swap( aElement,  6,  8 );
            this->swap( aElement, 13, 14 );
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::flip_quad4( Element * aElement )
        {
            // flip nodes
            this->swap( aElement,  1,  3 );
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::flip_quad9( Element * aElement )
        {
            // flip nodes
            this->swap( aElement,  1,  3 );
            this->swap( aElement,  4,  7 );
            this->swap( aElement,  5,  6 );
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::flip_quad16( Element * aElement )
        {
            this->swap( aElement,  1,  3 );
            this->swap( aElement,  4, 11 );
            this->swap( aElement,  5, 10 );
            this->swap( aElement,  6,  9 );
            this->swap( aElement,  7,  8 );
            this->swap( aElement, 13, 15 );
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::flip_tet4( Element * aElement )
        {
            this->swap( aElement, 1, 2 );
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::flip_tet10( Element * aElement )
        {
            this->swap( aElement, 1, 2 );
            this->swap( aElement, 4, 6 );
            this->swap( aElement, 8, 9 );
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::flip_tet20( Element * aElement )
        {
            this->swap( aElement,  1,  2 );
            this->swap( aElement,  4,  9 );
            this->swap( aElement,  5,  8 );
            this->swap( aElement,  6,  7 );
            this->swap( aElement, 12, 14 );
            this->swap( aElement, 13, 15 );
            this->swap( aElement, 17, 18 );
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::flip_tet35( Element * aElement )
        {
            this->swap( aElement,  1,  2 );
            this->swap( aElement,  4, 12 );
            this->swap( aElement,  5, 11 );
            this->swap( aElement,  6, 10 );
            this->swap( aElement,  7,  9 );
            this->swap( aElement, 16, 19 );
            this->swap( aElement, 17, 20 );
            this->swap( aElement, 18, 21 );
            this->swap( aElement, 23, 24 );
            this->swap( aElement, 25, 28 );
            this->swap( aElement, 26, 30 );
            this->swap( aElement, 27, 29 );
            this->swap( aElement, 32, 33 );
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::flip_penta6( Element * aElement )
        {
            this->swap( aElement,  1,  2 );
            this->swap( aElement,  4,  5 );
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::flip_penta15( Element * aElement )
        {
            this->swap( aElement,  1, 2 );
            this->swap( aElement,  4, 5 );
            this->swap( aElement,  6, 9 );
            this->swap( aElement, 10,11 );
            this->swap( aElement, 12,14 );
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::flip_penta18( Element * aElement )
        {
            this->swap( aElement,  1, 2 );
            this->swap( aElement,  4, 5 );
            this->swap( aElement,  6, 9 );
            this->swap( aElement, 10,11 );
            this->swap( aElement, 12,14 );
            this->swap( aElement, 15,17 );
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::flip_hex8( Element * aElement )
        {
            this->swap( aElement,  1, 3 );
            this->swap( aElement,  5, 7 );
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::flip_hex20( Element * aElement )
        {
            this->swap( aElement,  1, 3 );
            this->swap( aElement,  5, 7 );
            this->swap( aElement,  8,11 );
            this->swap( aElement,  9,12 );
            this->swap( aElement, 10,15 );
            this->swap( aElement, 16,18 );
            this->swap( aElement, 19,17 );
        }

//------------------------------------------------------------------------------

        void
        OrientationChecker::flip_hex27( Element * aElement )
        {
            this->swap( aElement,  1, 3 );
            this->swap( aElement,  5, 7 );
            this->swap( aElement,  8,11 );
            this->swap( aElement,  9,12 );
            this->swap( aElement, 10,15 );
            this->swap( aElement, 16,18 );
            this->swap( aElement, 19,17 );
            this->swap( aElement, 20,25 );
            this->swap( aElement, 22,23 );
        }
//------------------------------------------------------------------------------

        void
        OrientationChecker::flip_hex64( Element * aElement )
        {
            this->swap( aElement, 1, 3 );
            this->swap( aElement, 5, 7 );
            this->swap( aElement, 8, 10 );
            this->swap( aElement, 9, 11 );
            this->swap( aElement, 14, 19 );
            this->swap( aElement, 15, 18 );
            this->swap( aElement, 16, 22 );
            this->swap( aElement, 17, 23 );
            this->swap( aElement, 24, 26 );
            this->swap( aElement, 25, 27 );
            this->swap( aElement, 28, 31 );
            this->swap( aElement, 29, 30 );
            this->swap( aElement, 33, 35 );
            this->swap( aElement, 36, 40 );
            this->swap( aElement, 37, 43 );
            this->swap( aElement, 38, 42 );
            this->swap( aElement, 39, 41 );
            this->swap( aElement, 44, 49 );
            this->swap( aElement, 45, 48 );
            this->swap( aElement, 46, 51 );
            this->swap( aElement, 47, 50 );
            this->swap( aElement, 53, 55 );
            this->swap( aElement, 57, 59 );
            this->swap( aElement, 61, 63 );
        }

//------------------------------------------------------------------------------
    }
}