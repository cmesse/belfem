//
// Created by Christian Messe on 05.01.22.
//
#include "commtools.hpp"
#include "cl_Element.hpp"
#include "cl_Face.hpp"
#include "cl_Mesh_CurvedElementChecker.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        CurvedElementChecker::CurvedElementChecker(
                const uint aMeshDimension,
                Cell< Block * >   & aBlocks,
                Cell< SideSet * > & aSideSets ) :
                mMyRank( comm_rank() ),
                mMeshDimension( aMeshDimension ),
                mBlocks( aBlocks ),
                mSideSets( aSideSets )
        {
        }

//------------------------------------------------------------------------------

        index_t
        CurvedElementChecker::flag_curved_elements()
        {
            index_t aCount = 0 ;

            // loop over all blocks
            for( Block * tBlock : mBlocks )
            {
                // link the element function
                this->link_check_function( tBlock->element_type() );

                // get elements on block
                Cell< Element * > & tElements = tBlock->elements() ;

                // loop over all elements
                for( Element * tElement : tElements )
                {

                    // perform flagging
                    if( ( this->*mFunCheck )( tElement ) )
                    {
                        tElement->set_curved_flag() ;
                        ++aCount ;
                    }
                }
            }

            // loop over all sidesets and flag facets
            for( SideSet * tSideSet : mSideSets )
            {
                // get facets of sideset
                Cell< Facet * > & tFacets = tSideSet->facets() ;

                bool tMasterFlag ;
                bool tSlaveFlag ;

                // loop over all facets
                for( Facet * tFacet : tFacets )
                {
                    // check master
                    if( tFacet->has_master() )
                    {
                        tMasterFlag = tFacet->master()->is_curved() ;
                    }
                    else
                    {
                        tMasterFlag = false ;
                    }

                    // check slave
                    if ( tFacet->has_slave() )
                    {
                        tSlaveFlag = tFacet->slave()->is_curved() ;
                    }
                    else
                    {
                        tSlaveFlag = false ;
                    }

                    if( tMasterFlag || tSlaveFlag )
                    {
                        tFacet->element()->set_curved_flag() ;
                    }
                    else
                    {
                        tFacet->element()->unset_curved_flag() ;
                    }
                }
            }

            // return element counter
            return aCount ;
        }

//------------------------------------------------------------------------------

        void
        CurvedElementChecker::link_check_function( const ElementType aElementType )
        {
            if( mMeshDimension == 2 )
            {
                switch( aElementType )
                {
                    // linear elements
                    case( ElementType::LINE2 ) :
                    case( ElementType::TRI3 ) :
                    case( ElementType::QUAD4 ) :
                    {
                        mFunCheck = & CurvedElementChecker::check_linear ;
                        break ;
                    }
                    case( ElementType::TRI6 ) :
                    {
                        mFunCheck = & CurvedElementChecker::check_tri6_2d ;
                        break ;
                    }
                    case( ElementType::QUAD8 ) :
                    {
                        mFunCheck = & CurvedElementChecker::check_quad8_2d ;
                        break ;
                    }
                    case( ElementType::QUAD9 ) :
                    {
                        mFunCheck = & CurvedElementChecker::check_quad9_2d ;
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "No check function implemented for this element type" );
                    }
                }
            }
            else if ( mMeshDimension == 3 )
            {
                switch( aElementType )
                {
                    case( ElementType::LINE2 ) :
                    case( ElementType::TRI3 ) :
                    case( ElementType::TET4 ) :
                    {
                        mFunCheck = & CurvedElementChecker::check_linear ;
                        break ;
                    }
                    case( ElementType::TRI6 ) :
                    {
                        mFunCheck = & CurvedElementChecker::check_tri6 ;
                        break ;
                    }
                    case( ElementType::QUAD4 ) :
                    {
                        mFunCheck = & CurvedElementChecker::check_quad4 ;
                        break ;
                    }
                    case( ElementType::QUAD8 ) :
                    {
                        mFunCheck = & CurvedElementChecker::check_quad8 ;
                        break ;
                    }
                    case( ElementType::QUAD9 ) :
                    {
                        mFunCheck = & CurvedElementChecker::check_quad9 ;
                        break ;
                    }
                    case( ElementType::TET10 ) :
                    {
                        mFunCheck = & CurvedElementChecker::check_tet10 ;
                        break ;
                    }
                    case( ElementType::PENTA6 ) :
                    {
                        mFunCheck = & CurvedElementChecker::check_penta6 ;
                        break ;
                    }
                    case( ElementType::PENTA15 ) :
                    {
                        mFunCheck = & CurvedElementChecker::check_penta15 ;
                        break ;
                    }
                    case( ElementType::PENTA18 ) :
                    {
                        mFunCheck = & CurvedElementChecker::check_penta18 ;
                        break ;
                    }
                    case( ElementType::HEX8 ) :
                    {
                        mFunCheck = & CurvedElementChecker::check_hex8 ;
                        break ;
                    }
                    case( ElementType::HEX20 ) :
                    {
                        mFunCheck = & CurvedElementChecker::check_hex20 ;
                        break ;
                    }
                    case( ElementType::HEX27 ) :
                    {
                        mFunCheck = & CurvedElementChecker::check_hex27 ;
                        break ;
                    }
                    default :
                    {
                        BELFEM_ERROR( false, "No check function implemented for this element type" );
                    }
                }
            }
            else
            {
                BELFEM_ERROR( false, "Invalid mesh dimension: %u", ( unsigned int ) mMeshDimension );
            }
        }

//------------------------------------------------------------------------------

        bool
        CurvedElementChecker::check_tri6_2d( Element * aElement )
        {

            if( this->check_midpoint( aElement->node( 0 )->x(),
                                      aElement->node( 1 )->x(),
                                      aElement->node( 3 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 4 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 2 )->x(),
                                           aElement->node( 0 )->x(),
                                           aElement->node( 5 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->y(),
                                      aElement->node( 1 )->y(),
                                      aElement->node( 3 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 4 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 2 )->y(),
                                           aElement->node( 0 )->y(),
                                           aElement->node( 5 )->y() ) )
            {
                return true ;
            }

            return false ;

        }

//------------------------------------------------------------------------------

        bool
        CurvedElementChecker::check_quad8_2d( Element * aElement )
        {

            if( this->check_midpoint( aElement->node( 0 )->x(),
                                      aElement->node( 1 )->x(),
                                      aElement->node( 4 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 5 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 2 )->x(),
                                           aElement->node( 3 )->x(),
                                           aElement->node( 6 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 3 )->x(),
                                           aElement->node( 0 )->x(),
                                           aElement->node( 7 )->x()) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->y(),
                                      aElement->node( 1 )->y(),
                                      aElement->node( 4 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 5 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 2 )->y(),
                                           aElement->node( 3 )->y(),
                                           aElement->node( 6 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 3 )->y(),
                                           aElement->node( 0 )->y(),
                                           aElement->node( 7 )->y()) )
            {
                return true ;
            }

            return false ;
        }

//------------------------------------------------------------------------------

        bool
        CurvedElementChecker::check_quad9_2d( Element * aElement )
        {
            if( this->check_quad8_2d( aElement ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->x(),
                                           aElement->node( 1 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 3 )->x(),
                                           aElement->node( 8 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->y(),
                                           aElement->node( 1 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 3 )->y(),
                                           aElement->node( 8 )->y() ) )
            {
                return true ;
            }

            return false ;
        }

//------------------------------------------------------------------------------

        bool
        CurvedElementChecker::check_tri6( Element * aElement )
        {

            if( this->check_midpoint( aElement->node( 0 )->x(),
                                      aElement->node( 1 )->x(),
                                      aElement->node( 3 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 4 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 2 )->x(),
                                           aElement->node( 0 )->x(),
                                           aElement->node( 5 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->y(),
                                      aElement->node( 1 )->y(),
                                      aElement->node( 3 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 4 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 2 )->y(),
                                           aElement->node( 0 )->y(),
                                           aElement->node( 5 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->z(),
                                      aElement->node( 1 )->z(),
                                      aElement->node( 3 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->z(),
                                           aElement->node( 2 )->z(),
                                           aElement->node( 4 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 2 )->z(),
                                           aElement->node( 0 )->z(),
                                           aElement->node( 5 )->z() ) )
            {
                return true ;
            }

            return false ;
        }

//------------------------------------------------------------------------------

        bool
        CurvedElementChecker::check_quad4( Element * aElement )
        {
            // we create a plane from 0->1 and 0-> 2. Then we check
            // if point 3 is on this plane

            mX[ 0 ] = aElement->node( 0 )->x();
            mX[ 1 ] = aElement->node( 1 )->x();
            mX[ 2 ] = aElement->node( 2 )->x();
            mX[ 3 ] = aElement->node( 3 )->x();

            mY[ 0 ] = aElement->node( 0 )->y();
            mY[ 1 ] = aElement->node( 1 )->y();
            mY[ 2 ] = aElement->node( 2 )->y();
            mY[ 3 ] = aElement->node( 3 )->y();

            mZ[ 0 ] = aElement->node( 0 )->z();
            mZ[ 1 ] = aElement->node( 1 )->z();
            mZ[ 2 ] = aElement->node( 2 )->z();
            mZ[ 3 ] = aElement->node( 3 )->z();

            // let  a = p1-p0 and b = p2-p0
            // then n = cross( a, b )

            // distance of plane to origin = x * n
            real td = mZ[0] * ( mX[1] * mY[2] - mX[2]*mY[1] )
                      + mZ[1] * ( mX[2] * mY[0] - mX[0] * mY[2] )
                      + mZ[2] * ( mX[0] * mY[1] - mX[1] * mY[0] );

            // expression for point 3
            real tf = mX[0] * ( mZ[3] * ( mY[1] - mY[2] ) + mY[3] * ( mZ[2] - mZ[1]  ) )
                    + mX[1] * ( mZ[3] * ( mY[2] - mY[0] ) + mY[3] *( mZ[0] - mZ[2]  ) )
                    + mX[2] * ( mZ[3] * ( mY[0] - mY[1] ) + mY[3] *( mZ[1] - mZ[0] ) )
                    + mX[3] * ( mY[0] * ( mZ[1] - mZ[2] ) + mY[1] *( mZ[2] - mZ[0] )
                                + mY[2]*( mZ[0] - mZ[1] ) );
            // check if point does not sit on plane
            return std::abs( td - tf ) > mTwoMeshEpsilon ;
        }

//------------------------------------------------------------------------------

        bool
        CurvedElementChecker::check_quad4_face( Element * aElement, const uint aFaceIndex )
        {
            // we create a plane from 0->1 and 0-> 2. Then we check
            // if point 3 is on this plane

            aElement->get_corner_nodes_of_facet( aFaceIndex, mNodes );

            mX[ 0 ] = mNodes( 0 )->x();
            mX[ 1 ] = mNodes( 1 )->x();
            mX[ 2 ] = mNodes( 2 )->x();
            mX[ 3 ] = mNodes( 3 )->x();

            mY[ 0 ] = mNodes( 0 )->y();
            mY[ 1 ] = mNodes( 1 )->y();
            mY[ 2 ] = mNodes( 2 )->y();
            mY[ 3 ] = mNodes( 3 )->y();

            mZ[ 0 ] = mNodes( 0 )->z();
            mZ[ 1 ] = mNodes( 1 )->z();
            mZ[ 2 ] = mNodes( 2 )->z();
            mZ[ 3 ] = mNodes( 3 )->z();

            // let  a = p1-p0 and b = p2-p0
            // then n = cross( a, b )

            // distance of plane to origin = x * n
            real td = mZ[0] * ( mX[1] * mY[2] - mX[2]*mY[1] )
                      + mZ[1] * ( mX[2] * mY[0] - mX[0] * mY[2] )
                      + mZ[2] * ( mX[0] * mY[1] - mX[1] * mY[0] );

            // expression for point 3
            real tf = mX[0] * ( mZ[3] * ( mY[1] - mY[2] ) + mY[3] *( mZ[2] - mZ[1]  ) )
            + mX[1] * ( mZ[3] * ( mY[2] - mY[0] ) + mY[3] *( mZ[0] - mZ[2]  ) )
            + mX[2] * ( mZ[3] * ( mY[0] - mY[1] ) + mY[3] *( mZ[1] - mZ[0] ) )
            + mX[3] * ( mY[0] * ( mZ[1] - mZ[2] ) + mY[1] *( mZ[2] - mZ[0] )
                + mY[2]*( mZ[0] - mZ[1] ) );

            // check if point does not sit on plane
            return std::abs( td - tf ) > mTwoMeshEpsilon ;
        }

//------------------------------------------------------------------------------

        bool
        CurvedElementChecker::check_quad8( Element * aElement )
        {
            if( this->check_quad4( aElement ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->x(),
                                      aElement->node( 1 )->x(),
                                      aElement->node( 4 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 5 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 2 )->x(),
                                           aElement->node( 3 )->x(),
                                           aElement->node( 6 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 3 )->x(),
                                           aElement->node( 0 )->x(),
                                           aElement->node( 7 )->x()) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->y(),
                                      aElement->node( 1 )->y(),
                                      aElement->node( 4 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 5 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 2 )->y(),
                                           aElement->node( 3 )->y(),
                                           aElement->node( 6 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 3 )->y(),
                                           aElement->node( 0 )->y(),
                                           aElement->node( 7 )->y()) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->z(),
                                      aElement->node( 1 )->z(),
                                      aElement->node( 4 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->z(),
                                           aElement->node( 2 )->z(),
                                           aElement->node( 5 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 2 )->z(),
                                           aElement->node( 3 )->z(),
                                           aElement->node( 6 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 3 )->z(),
                                           aElement->node( 0 )->z(),
                                           aElement->node( 7 )->z()) )
            {
                return true ;
            }

            return false ;
        }

//------------------------------------------------------------------------------

        bool
        CurvedElementChecker::check_quad9( Element * aElement )
        {
            if( this->check_quad8( aElement ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->x(),
                                           aElement->node( 1 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 3 )->x(),
                                           aElement->node( 8 )->x() ) )
            {
                return true ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->y(),
                                           aElement->node( 1 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 3 )->y(),
                                           aElement->node( 8 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->z(),
                                           aElement->node( 1 )->z(),
                                           aElement->node( 2 )->z(),
                                           aElement->node( 3 )->z(),
                                           aElement->node( 8 )->z() ) )
            {
                return true ;
            }

            return false ;
        }

//------------------------------------------------------------------------------

        bool
        CurvedElementChecker::check_tet10( Element * aElement )
        {
            if( this->check_midpoint( aElement->node( 0 )->x(),
                                      aElement->node( 1 )->x(),
                                      aElement->node( 4 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->x(),
                                      aElement->node( 2 )->x(),
                                      aElement->node( 6 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->x(),
                                      aElement->node( 3 )->x(),
                                      aElement->node( 7 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->x(),
                                      aElement->node( 2 )->x(),
                                      aElement->node( 5 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->x(),
                                      aElement->node( 3 )->x(),
                                      aElement->node( 8 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 2 )->x(),
                                      aElement->node( 3 )->x(),
                                      aElement->node( 9 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->y(),
                                      aElement->node( 1 )->y(),
                                      aElement->node( 4 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->y(),
                                      aElement->node( 2 )->y(),
                                      aElement->node( 6 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->y(),
                                      aElement->node( 3 )->y(),
                                      aElement->node( 7 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->y(),
                                      aElement->node( 2 )->y(),
                                      aElement->node( 5 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->y(),
                                      aElement->node( 3 )->y(),
                                      aElement->node( 8 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 2 )->y(),
                                      aElement->node( 3 )->y(),
                                      aElement->node( 9 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->z(),
                                      aElement->node( 1 )->z(),
                                      aElement->node( 4 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->z(),
                                      aElement->node( 2 )->z(),
                                      aElement->node( 6 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->z(),
                                      aElement->node( 3 )->z(),
                                      aElement->node( 7 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->z(),
                                      aElement->node( 2 )->z(),
                                      aElement->node( 5 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->z(),
                                      aElement->node( 3 )->z(),
                                      aElement->node( 8 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 2 )->z(),
                                      aElement->node( 3 )->z(),
                                      aElement->node( 9 )->z() ) )
            {
                return true ;
            }

            return false ;
        }

//------------------------------------------------------------------------------

        bool
        CurvedElementChecker::check_penta6( Element * aElement )
        {
            for( uint k=0 ; k<3; ++k )
            {
                if( this->check_quad4_face( aElement, k ) )
                {
                    return true ;
                }
            }
            return false ;
        }

//------------------------------------------------------------------------------

        bool
        CurvedElementChecker::check_penta15( Element * aElement )
        {
            if( this->check_penta6( aElement ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->x(),
                                      aElement->node( 1 )->x(),
                                      aElement->node( 6 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 8 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->x(),
                                           aElement->node( 3 )->x(),
                                           aElement->node( 9 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 7 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->x(),
                                           aElement->node( 4 )->x(),
                                           aElement->node( 10 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 2 )->x(),
                                           aElement->node( 5 )->x(),
                                           aElement->node( 11 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 3 )->x(),
                                           aElement->node( 4 )->x(),
                                           aElement->node( 12 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 3 )->x(),
                                           aElement->node( 5 )->x(),
                                           aElement->node( 14 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 4 )->x(),
                                           aElement->node( 5 )->x(),
                                           aElement->node( 13 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->y(),
                                           aElement->node( 1 )->y(),
                                           aElement->node( 6 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 8 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->y(),
                                           aElement->node( 3 )->y(),
                                           aElement->node( 9 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 7 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->y(),
                                           aElement->node( 4 )->y(),
                                           aElement->node( 10 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 2 )->y(),
                                           aElement->node( 5 )->y(),
                                           aElement->node( 11 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 3 )->y(),
                                           aElement->node( 4 )->y(),
                                           aElement->node( 12 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 3 )->y(),
                                           aElement->node( 5 )->y(),
                                           aElement->node( 14 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 4 )->y(),
                                           aElement->node( 5 )->y(),
                                           aElement->node( 13 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->z(),
                                      aElement->node( 1 )->z(),
                                      aElement->node( 6 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->z(),
                                           aElement->node( 2 )->z(),
                                           aElement->node( 8 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->z(),
                                           aElement->node( 3 )->z(),
                                           aElement->node( 9 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->z(),
                                           aElement->node( 2 )->z(),
                                           aElement->node( 7 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->z(),
                                           aElement->node( 4 )->z(),
                                           aElement->node( 10 )->z() ) )
            {
                return true ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->z(),
                                           aElement->node( 5 )->z(),
                                           aElement->node( 11 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 3 )->z(),
                                           aElement->node( 4 )->z(),
                                           aElement->node( 12 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 3 )->z(),
                                           aElement->node( 5 )->z(),
                                           aElement->node( 14 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 4 )->z(),
                                           aElement->node( 5 )->z(),
                                           aElement->node( 13 )->z() ) )
            {
                return true ;
            }

            return false ;
        }

//------------------------------------------------------------------------------

        bool
        CurvedElementChecker::check_penta18( Element * aElement )
        {
            if( this->check_penta15( aElement ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 4 )->x(),
                                           aElement->node( 5 )->x(),
                                           aElement->node( 16 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 3 )->x(),
                                           aElement->node( 5 )->x(),
                                           aElement->node( 17 )->x() ) )
            {
                return true ;
            }

            if( this->check_midpoint( aElement->node( 0 )->y(),
                                           aElement->node( 1 )->y(),
                                           aElement->node( 3 )->y(),
                                           aElement->node( 4 )->y(),
                                           aElement->node( 15 )->y() ) )
            {
                return true ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 4 )->y(),
                                           aElement->node( 5 )->y(),
                                           aElement->node( 16 )->y() ) )
            {
                return true ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 3 )->y(),
                                           aElement->node( 5 )->y(),
                                           aElement->node( 17 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->z(),
                                           aElement->node( 1 )->z(),
                                           aElement->node( 3 )->z(),
                                           aElement->node( 4 )->z(),
                                           aElement->node( 15 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 1 )->z(),
                                           aElement->node( 2 )->z(),
                                           aElement->node( 4 )->z(),
                                           aElement->node( 5 )->z(),
                                           aElement->node( 16 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node( 0 )->z(),
                                           aElement->node( 2 )->z(),
                                           aElement->node( 3 )->z(),
                                           aElement->node( 5 )->z(),
                                           aElement->node( 17 )->z() ) )
            {
                return true ;
            }

            return false ;
        }

//------------------------------------------------------------------------------

        bool
        CurvedElementChecker::check_hex8( Element * aElement )
        {
            for( uint k=0 ; k<6; ++k )
            {
                if( this->check_quad4_face( aElement, k ) )
                {
                    return true ;
                }
            }
            return false ;
        }

//------------------------------------------------------------------------------

        bool
        CurvedElementChecker::check_hex20( Element * aElement )
        {
            if( this->check_hex8( aElement ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->x(),
                                      aElement->node(  1 )->x(),
                                      aElement->node(  8 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->x(),
                                      aElement->node(  3 )->x(),
                                      aElement->node( 11 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->x(),
                                      aElement->node(  4 )->x(),
                                      aElement->node( 12 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  1 )->x(),
                                      aElement->node(  2 )->x(),
                                      aElement->node(  9 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  1 )->x(),
                                      aElement->node(  5 )->x(),
                                      aElement->node( 13 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  2 )->x(),
                                      aElement->node(  3 )->x(),
                                      aElement->node( 10 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  2 )->x(),
                                      aElement->node(  6 )->x(),
                                      aElement->node( 14 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  3 )->x(),
                                      aElement->node(  7 )->x(),
                                      aElement->node( 15 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  4 )->x(),
                                      aElement->node(  5 )->x(),
                                      aElement->node( 16 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  4 )->x(),
                                      aElement->node(  7 )->x(),
                                      aElement->node( 19 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  5 )->x(),
                                      aElement->node(  6 )->x(),
                                      aElement->node( 17 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  6 )->x(),
                                      aElement->node(  7 )->x(),
                                      aElement->node( 18 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->y(),
                                      aElement->node(  1 )->y(),
                                      aElement->node(  8 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->y(),
                                      aElement->node(  3 )->y(),
                                      aElement->node( 11 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->y(),
                                      aElement->node(  4 )->y(),
                                      aElement->node( 12 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  1 )->y(),
                                      aElement->node(  2 )->y(),
                                      aElement->node(  9 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  1 )->y(),
                                      aElement->node(  5 )->y(),
                                      aElement->node( 13 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  2 )->y(),
                                      aElement->node(  3 )->y(),
                                      aElement->node( 10 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  2 )->y(),
                                      aElement->node(  6 )->y(),
                                      aElement->node( 14 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  3 )->y(),
                                      aElement->node(  7 )->y(),
                                      aElement->node( 15 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  4 )->y(),
                                      aElement->node(  5 )->y(),
                                      aElement->node( 16 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  4 )->y(),
                                      aElement->node(  7 )->y(),
                                      aElement->node( 19 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  5 )->y(),
                                      aElement->node(  6 )->y(),
                                      aElement->node( 17 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  6 )->y(),
                                      aElement->node(  7 )->y(),
                                      aElement->node( 18 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->y(),
                                      aElement->node(  1 )->y(),
                                      aElement->node(  8 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->y(),
                                      aElement->node(  3 )->y(),
                                      aElement->node( 11 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->y(),
                                      aElement->node(  4 )->y(),
                                      aElement->node( 12 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  1 )->y(),
                                      aElement->node(  2 )->y(),
                                      aElement->node(  9 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  1 )->y(),
                                      aElement->node(  5 )->y(),
                                      aElement->node( 13 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  2 )->y(),
                                      aElement->node(  3 )->y(),
                                      aElement->node( 10 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  2 )->y(),
                                      aElement->node(  6 )->y(),
                                      aElement->node( 14 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  3 )->y(),
                                      aElement->node(  7 )->y(),
                                      aElement->node( 15 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  4 )->y(),
                                      aElement->node(  5 )->y(),
                                      aElement->node( 16 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  4 )->y(),
                                      aElement->node(  7 )->y(),
                                      aElement->node( 19 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  5 )->y(),
                                      aElement->node(  6 )->y(),
                                      aElement->node( 17 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  6 )->y(),
                                      aElement->node(  7 )->y(),
                                      aElement->node( 18 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->z(),
                                      aElement->node(  1 )->z(),
                                      aElement->node(  8 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->z(),
                                      aElement->node(  3 )->z(),
                                      aElement->node( 11 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->z(),
                                      aElement->node(  4 )->z(),
                                      aElement->node( 12 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  1 )->z(),
                                      aElement->node(  2 )->z(),
                                      aElement->node(  9 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  1 )->z(),
                                      aElement->node(  5 )->z(),
                                      aElement->node( 13 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  2 )->z(),
                                      aElement->node(  3 )->z(),
                                      aElement->node( 10 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  2 )->z(),
                                      aElement->node(  6 )->z(),
                                      aElement->node( 14 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  3 )->z(),
                                      aElement->node(  7 )->z(),
                                      aElement->node( 15 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  4 )->z(),
                                      aElement->node(  5 )->z(),
                                      aElement->node( 16 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  4 )->z(),
                                      aElement->node(  7 )->z(),
                                      aElement->node( 19 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  5 )->z(),
                                      aElement->node(  6 )->z(),
                                      aElement->node( 17 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  6 )->z(),
                                      aElement->node(  7 )->z(),
                                      aElement->node( 18 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->z(),
                                      aElement->node(  1 )->z(),
                                      aElement->node(  8 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->z(),
                                      aElement->node(  3 )->z(),
                                      aElement->node( 11 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->z(),
                                      aElement->node(  4 )->z(),
                                      aElement->node( 12 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  1 )->z(),
                                      aElement->node(  2 )->z(),
                                      aElement->node(  9 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  1 )->z(),
                                      aElement->node(  5 )->z(),
                                      aElement->node( 13 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  2 )->z(),
                                      aElement->node(  3 )->z(),
                                      aElement->node( 10 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  2 )->z(),
                                      aElement->node(  6 )->z(),
                                      aElement->node( 14 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  3 )->z(),
                                      aElement->node(  7 )->z(),
                                      aElement->node( 15 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  4 )->z(),
                                      aElement->node(  5 )->z(),
                                      aElement->node( 16 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  4 )->z(),
                                      aElement->node(  7 )->z(),
                                      aElement->node( 19 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  5 )->z(),
                                      aElement->node(  6 )->z(),
                                      aElement->node( 17 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  6 )->z(),
                                      aElement->node(  7 )->z(),
                                      aElement->node( 18 )->z() ) )
            {
                return true ;
            }

            return false ;
        }

//------------------------------------------------------------------------------


        bool
        CurvedElementChecker::check_hex27( Element * aElement )
        {
            if( this->check_hex20( aElement ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->x(),
                                      aElement->node(  1 )->x(),
                                      aElement->node(  4 )->x(),
                                      aElement->node(  5 )->x(),
                                      aElement->node( 25 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  1 )->x(),
                                      aElement->node(  2 )->x(),
                                      aElement->node(  5 )->x(),
                                      aElement->node(  6 )->x(),
                                      aElement->node( 24 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  2 )->x(),
                                      aElement->node(  3 )->x(),
                                      aElement->node(  6 )->x(),
                                      aElement->node(  7 )->x(),
                                      aElement->node( 26 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->x(),
                                      aElement->node(  3 )->x(),
                                      aElement->node(  4 )->x(),
                                      aElement->node(  7 )->x(),
                                      aElement->node( 23 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->x(),
                                      aElement->node(  1 )->x(),
                                      aElement->node(  2 )->x(),
                                      aElement->node(  3 )->x(),
                                      aElement->node( 21 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  4 )->x(),
                                      aElement->node(  5 )->x(),
                                      aElement->node(  6 )->x(),
                                      aElement->node(  7 )->x(),
                                      aElement->node( 22 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->x(),
                                      aElement->node(  1 )->x(),
                                      aElement->node(  2 )->x(),
                                      aElement->node(  3 )->x(),
                                      aElement->node(  4 )->x(),
                                      aElement->node(  5 )->x(),
                                      aElement->node(  6 )->x(),
                                      aElement->node(  7 )->x(),
                                      aElement->node( 20 )->x() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->y(),
                                      aElement->node(  1 )->y(),
                                      aElement->node(  4 )->y(),
                                      aElement->node(  5 )->y(),
                                      aElement->node( 25 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  1 )->y(),
                                      aElement->node(  2 )->y(),
                                      aElement->node(  5 )->y(),
                                      aElement->node(  6 )->y(),
                                      aElement->node( 24 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  2 )->y(),
                                      aElement->node(  3 )->y(),
                                      aElement->node(  6 )->y(),
                                      aElement->node(  7 )->y(),
                                      aElement->node( 26 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->y(),
                                      aElement->node(  3 )->y(),
                                      aElement->node(  4 )->y(),
                                      aElement->node(  7 )->y(),
                                      aElement->node( 23 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->y(),
                                      aElement->node(  1 )->y(),
                                      aElement->node(  2 )->y(),
                                      aElement->node(  3 )->y(),
                                      aElement->node( 21 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  4 )->y(),
                                      aElement->node(  5 )->y(),
                                      aElement->node(  6 )->y(),
                                      aElement->node(  7 )->y(),
                                      aElement->node( 22 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->y(),
                                      aElement->node(  1 )->y(),
                                      aElement->node(  2 )->y(),
                                      aElement->node(  3 )->y(),
                                      aElement->node(  4 )->y(),
                                      aElement->node(  5 )->y(),
                                      aElement->node(  6 )->y(),
                                      aElement->node(  7 )->y(),
                                      aElement->node( 20 )->y() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->z(),
                                      aElement->node(  1 )->z(),
                                      aElement->node(  4 )->z(),
                                      aElement->node(  5 )->z(),
                                      aElement->node( 25 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  1 )->z(),
                                      aElement->node(  2 )->z(),
                                      aElement->node(  5 )->z(),
                                      aElement->node(  6 )->z(),
                                      aElement->node( 24 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  2 )->z(),
                                      aElement->node(  3 )->z(),
                                      aElement->node(  6 )->z(),
                                      aElement->node(  7 )->z(),
                                      aElement->node( 26 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->z(),
                                      aElement->node(  3 )->z(),
                                      aElement->node(  4 )->z(),
                                      aElement->node(  7 )->z(),
                                      aElement->node( 23 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->z(),
                                      aElement->node(  1 )->z(),
                                      aElement->node(  2 )->z(),
                                      aElement->node(  3 )->z(),
                                      aElement->node( 21 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  4 )->z(),
                                      aElement->node(  5 )->z(),
                                      aElement->node(  6 )->z(),
                                      aElement->node(  7 )->z(),
                                      aElement->node( 22 )->z() ) )
            {
                return true ;
            }
            if( this->check_midpoint( aElement->node(  0 )->z(),
                                      aElement->node(  1 )->z(),
                                      aElement->node(  2 )->z(),
                                      aElement->node(  3 )->z(),
                                      aElement->node(  4 )->z(),
                                      aElement->node(  5 )->z(),
                                      aElement->node(  6 )->z(),
                                      aElement->node(  7 )->z(),
                                      aElement->node( 20 )->z() ) )
            {
                return true ;
            }

            return false;
        }

//------------------------------------------------------------------------------
    }
}