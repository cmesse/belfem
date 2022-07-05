//
// Created by Christian Messe on 05.01.22.
//
#include "commtools.hpp"
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
                    ( this->*mFunCheck )( tElement );

                    // check if element is counted
                    if( tElement->is_curved() )
                    {
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
                    case( ElementType::TRI3 ) :
                    case( ElementType::QUAD4 ) :
                    case( ElementType::TET4 ) :
                    case( ElementType::PENTA6 ) :
                    case( ElementType::HEX8 ) :
                    {
                        mFunCheck = & CurvedElementChecker::check_linear ;
                        break ;
                    }
                    case( ElementType::TRI6 ) :
                    {
                        mFunCheck = & CurvedElementChecker::check_tri6 ;
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

        void
        CurvedElementChecker::check_tri6_2d( Element * aElement )
        {

            if( this->check_midpoint( aElement->node( 0 )->x(),
                                      aElement->node( 1 )->x(),
                                      aElement->node( 3 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 4 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->x(),
                                           aElement->node( 0 )->x(),
                                           aElement->node( 5 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->y(),
                                      aElement->node( 1 )->y(),
                                      aElement->node( 3 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 4 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->y(),
                                           aElement->node( 0 )->y(),
                                           aElement->node( 5 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else
            {
                aElement->unset_curved_flag() ;
            }
        }

//------------------------------------------------------------------------------

        void
        CurvedElementChecker::check_quad8_2d( Element * aElement )
        {

            if( this->check_midpoint( aElement->node( 0 )->x(),
                                      aElement->node( 1 )->x(),
                                      aElement->node( 4 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 5 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->x(),
                                           aElement->node( 3 )->x(),
                                           aElement->node( 6 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->x(),
                                           aElement->node( 0 )->x(),
                                           aElement->node( 7 )->x()) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->y(),
                                      aElement->node( 1 )->y(),
                                      aElement->node( 4 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 5 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->y(),
                                           aElement->node( 3 )->y(),
                                           aElement->node( 6 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->y(),
                                           aElement->node( 0 )->y(),
                                           aElement->node( 7 )->y()) )
            {
                aElement->set_curved_flag() ;
            }
            else
            {
                aElement->unset_curved_flag() ;
            }
        }

//------------------------------------------------------------------------------

        void
        CurvedElementChecker::check_quad9_2d( Element * aElement )
        {

            if( this->check_midpoint( aElement->node( 0 )->x(),
                                      aElement->node( 1 )->x(),
                                      aElement->node( 4 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 5 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->x(),
                                           aElement->node( 3 )->x(),
                                           aElement->node( 6 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->x(),
                                           aElement->node( 0 )->x(),
                                           aElement->node( 7 )->x()) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 3 )->x(),
                                           aElement->node( 4 )->x(),
                                           aElement->node( 8 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->y(),
                                      aElement->node( 1 )->y(),
                                      aElement->node( 4 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 5 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->y(),
                                           aElement->node( 3 )->y(),
                                           aElement->node( 6 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->y(),
                                           aElement->node( 0 )->y(),
                                           aElement->node( 7 )->y()) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 3 )->y(),
                                           aElement->node( 4 )->y(),
                                           aElement->node( 8 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else
            {
                aElement->unset_curved_flag() ;
            }
        }

//------------------------------------------------------------------------------

        void
        CurvedElementChecker::check_tri6( Element * aElement )
        {

            if( this->check_midpoint( aElement->node( 0 )->x(),
                                      aElement->node( 1 )->x(),
                                      aElement->node( 3 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 4 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->x(),
                                           aElement->node( 0 )->x(),
                                           aElement->node( 5 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->y(),
                                      aElement->node( 1 )->y(),
                                      aElement->node( 3 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 4 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->y(),
                                           aElement->node( 0 )->y(),
                                           aElement->node( 5 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->z(),
                                      aElement->node( 1 )->z(),
                                      aElement->node( 3 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->z(),
                                           aElement->node( 2 )->z(),
                                           aElement->node( 4 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->z(),
                                           aElement->node( 0 )->z(),
                                           aElement->node( 5 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else
            {
                aElement->unset_curved_flag() ;
            }
        }

//------------------------------------------------------------------------------

        void
        CurvedElementChecker::check_quad8( Element * aElement )
        {

            if( this->check_midpoint( aElement->node( 0 )->x(),
                                      aElement->node( 1 )->x(),
                                      aElement->node( 4 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 5 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->x(),
                                           aElement->node( 3 )->x(),
                                           aElement->node( 6 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->x(),
                                           aElement->node( 0 )->x(),
                                           aElement->node( 7 )->x()) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->y(),
                                      aElement->node( 1 )->y(),
                                      aElement->node( 4 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 5 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->y(),
                                           aElement->node( 3 )->y(),
                                           aElement->node( 6 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->y(),
                                           aElement->node( 0 )->y(),
                                           aElement->node( 7 )->y()) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->z(),
                                      aElement->node( 1 )->z(),
                                      aElement->node( 4 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->z(),
                                           aElement->node( 2 )->z(),
                                           aElement->node( 5 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->z(),
                                           aElement->node( 3 )->z(),
                                           aElement->node( 6 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->z(),
                                           aElement->node( 0 )->z(),
                                           aElement->node( 7 )->z()) )
            {
                aElement->set_curved_flag() ;
            }
            else
            {
                aElement->unset_curved_flag() ;
            }
        }

//------------------------------------------------------------------------------

        void
        CurvedElementChecker::check_quad9( Element * aElement )
        {

            if( this->check_midpoint( aElement->node( 0 )->x(),
                                      aElement->node( 1 )->x(),
                                      aElement->node( 4 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 5 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->x(),
                                           aElement->node( 3 )->x(),
                                           aElement->node( 6 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->x(),
                                           aElement->node( 0 )->x(),
                                           aElement->node( 7 )->x()) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 3 )->x(),
                                           aElement->node( 4 )->x(),
                                           aElement->node( 8 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->y(),
                                      aElement->node( 1 )->y(),
                                      aElement->node( 4 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 5 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->y(),
                                           aElement->node( 3 )->y(),
                                           aElement->node( 6 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->y(),
                                           aElement->node( 0 )->y(),
                                           aElement->node( 7 )->y()) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 3 )->y(),
                                           aElement->node( 4 )->y(),
                                           aElement->node( 8 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            if( this->check_midpoint( aElement->node( 0 )->z(),
                                      aElement->node( 1 )->z(),
                                      aElement->node( 4 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->z(),
                                           aElement->node( 2 )->z(),
                                           aElement->node( 5 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->z(),
                                           aElement->node( 3 )->z(),
                                           aElement->node( 6 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->z(),
                                           aElement->node( 0 )->z(),
                                           aElement->node( 7 )->z()) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->z(),
                                           aElement->node( 2 )->z(),
                                           aElement->node( 3 )->z(),
                                           aElement->node( 4 )->z(),
                                           aElement->node( 8 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else
            {
                aElement->unset_curved_flag() ;
            }
        }

//------------------------------------------------------------------------------

        void
        CurvedElementChecker::check_tet10( Element * aElement )
        {
            if( this->check_midpoint( aElement->node( 0 )->x(),
                                      aElement->node( 1 )->x(),
                                      aElement->node( 4 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->x(),
                                      aElement->node( 2 )->x(),
                                      aElement->node( 6 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->x(),
                                      aElement->node( 3 )->x(),
                                      aElement->node( 7 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->x(),
                                      aElement->node( 2 )->x(),
                                      aElement->node( 5 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->x(),
                                      aElement->node( 3 )->x(),
                                      aElement->node( 8 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->x(),
                                      aElement->node( 3 )->x(),
                                      aElement->node( 9 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->y(),
                                      aElement->node( 1 )->y(),
                                      aElement->node( 4 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->y(),
                                      aElement->node( 2 )->y(),
                                      aElement->node( 6 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->y(),
                                      aElement->node( 3 )->y(),
                                      aElement->node( 7 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->y(),
                                      aElement->node( 2 )->y(),
                                      aElement->node( 5 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->y(),
                                      aElement->node( 3 )->y(),
                                      aElement->node( 8 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->y(),
                                      aElement->node( 3 )->y(),
                                      aElement->node( 9 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->z(),
                                      aElement->node( 1 )->z(),
                                      aElement->node( 4 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->z(),
                                      aElement->node( 2 )->z(),
                                      aElement->node( 6 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->z(),
                                      aElement->node( 3 )->z(),
                                      aElement->node( 7 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->z(),
                                      aElement->node( 2 )->z(),
                                      aElement->node( 5 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->z(),
                                      aElement->node( 3 )->z(),
                                      aElement->node( 8 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->z(),
                                      aElement->node( 3 )->z(),
                                      aElement->node( 9 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else
            {
                aElement->unset_curved_flag();
            }
        }
//------------------------------------------------------------------------------

        void
        CurvedElementChecker::check_penta15( Element * aElement )
        {
            if( this->check_midpoint( aElement->node( 0 )->x(),
                                      aElement->node( 1 )->x(),
                                      aElement->node( 6 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 8 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->x(),
                                           aElement->node( 3 )->x(),
                                           aElement->node( 9 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 7 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->x(),
                                           aElement->node( 4 )->x(),
                                           aElement->node( 10 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->x(),
                                           aElement->node( 5 )->x(),
                                           aElement->node( 11 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->x(),
                                           aElement->node( 4 )->x(),
                                           aElement->node( 12 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->x(),
                                           aElement->node( 5 )->x(),
                                           aElement->node( 14 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 4 )->x(),
                                           aElement->node( 5 )->x(),
                                           aElement->node( 13 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->y(),
                                           aElement->node( 1 )->y(),
                                           aElement->node( 6 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 8 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->y(),
                                           aElement->node( 3 )->y(),
                                           aElement->node( 9 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 7 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->y(),
                                           aElement->node( 4 )->y(),
                                           aElement->node( 10 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->y(),
                                           aElement->node( 5 )->y(),
                                           aElement->node( 11 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->y(),
                                           aElement->node( 4 )->y(),
                                           aElement->node( 12 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->y(),
                                           aElement->node( 5 )->y(),
                                           aElement->node( 14 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 4 )->y(),
                                           aElement->node( 5 )->y(),
                                           aElement->node( 13 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->z(),
                                      aElement->node( 1 )->z(),
                                      aElement->node( 6 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->z(),
                                           aElement->node( 2 )->z(),
                                           aElement->node( 8 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->z(),
                                           aElement->node( 3 )->z(),
                                           aElement->node( 9 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->z(),
                                           aElement->node( 2 )->z(),
                                           aElement->node( 7 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->z(),
                                           aElement->node( 4 )->z(),
                                           aElement->node( 10 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->z(),
                                           aElement->node( 5 )->z(),
                                           aElement->node( 11 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->z(),
                                           aElement->node( 4 )->z(),
                                           aElement->node( 12 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->z(),
                                           aElement->node( 5 )->z(),
                                           aElement->node( 14 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 4 )->z(),
                                           aElement->node( 5 )->z(),
                                           aElement->node( 13 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else
            {
                aElement->unset_curved_flag() ;
            }
        }

//------------------------------------------------------------------------------

        void
        CurvedElementChecker::check_penta18( Element * aElement )
        {
            if( this->check_midpoint( aElement->node( 0 )->x(),
                                      aElement->node( 1 )->x(),
                                      aElement->node( 6 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 8 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->x(),
                                           aElement->node( 3 )->x(),
                                           aElement->node( 9 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 7 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->x(),
                                           aElement->node( 4 )->x(),
                                           aElement->node( 10 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->x(),
                                           aElement->node( 5 )->x(),
                                           aElement->node( 11 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->x(),
                                           aElement->node( 4 )->x(),
                                           aElement->node( 12 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->x(),
                                           aElement->node( 5 )->x(),
                                           aElement->node( 14 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 4 )->x(),
                                           aElement->node( 5 )->x(),
                                           aElement->node( 13 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->x(),
                                           aElement->node( 1 )->x(),
                                           aElement->node( 3 )->x(),
                                           aElement->node( 4 )->x(),
                                           aElement->node( 15 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 4 )->x(),
                                           aElement->node( 5 )->x(),
                                           aElement->node( 16 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->x(),
                                           aElement->node( 2 )->x(),
                                           aElement->node( 3 )->x(),
                                           aElement->node( 5 )->x(),
                                           aElement->node( 17 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->y(),
                                           aElement->node( 1 )->y(),
                                           aElement->node( 6 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 8 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->y(),
                                           aElement->node( 3 )->y(),
                                           aElement->node( 9 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 7 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->y(),
                                           aElement->node( 4 )->y(),
                                           aElement->node( 10 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->y(),
                                           aElement->node( 5 )->y(),
                                           aElement->node( 11 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->y(),
                                           aElement->node( 4 )->y(),
                                           aElement->node( 12 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->y(),
                                           aElement->node( 5 )->y(),
                                           aElement->node( 14 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 4 )->y(),
                                           aElement->node( 5 )->y(),
                                           aElement->node( 13 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->y(),
                                           aElement->node( 1 )->y(),
                                           aElement->node( 3 )->y(),
                                           aElement->node( 4 )->y(),
                                           aElement->node( 15 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 4 )->y(),
                                           aElement->node( 5 )->y(),
                                           aElement->node( 16 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->y(),
                                           aElement->node( 2 )->y(),
                                           aElement->node( 3 )->y(),
                                           aElement->node( 5 )->y(),
                                           aElement->node( 17 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->z(),
                                           aElement->node( 1 )->z(),
                                           aElement->node( 6 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->z(),
                                           aElement->node( 2 )->z(),
                                           aElement->node( 8 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->z(),
                                           aElement->node( 3 )->z(),
                                           aElement->node( 9 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->z(),
                                           aElement->node( 2 )->z(),
                                           aElement->node( 7 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->z(),
                                           aElement->node( 4 )->z(),
                                           aElement->node( 10 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 2 )->z(),
                                           aElement->node( 5 )->z(),
                                           aElement->node( 11 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->z(),
                                           aElement->node( 4 )->z(),
                                           aElement->node( 12 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 3 )->z(),
                                           aElement->node( 5 )->z(),
                                           aElement->node( 14 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 4 )->z(),
                                           aElement->node( 5 )->z(),
                                           aElement->node( 13 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->z(),
                                           aElement->node( 1 )->z(),
                                           aElement->node( 3 )->z(),
                                           aElement->node( 4 )->z(),
                                           aElement->node( 15 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 1 )->z(),
                                           aElement->node( 2 )->z(),
                                           aElement->node( 4 )->z(),
                                           aElement->node( 5 )->z(),
                                           aElement->node( 16 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node( 0 )->z(),
                                           aElement->node( 2 )->z(),
                                           aElement->node( 3 )->z(),
                                           aElement->node( 5 )->z(),
                                           aElement->node( 17 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else
            {
                aElement->unset_curved_flag() ;
            }
        }

//------------------------------------------------------------------------------

        void
        CurvedElementChecker::check_hex20( Element * aElement )
        {
            if( this->check_midpoint( aElement->node(  0 )->x(),
                                           aElement->node(  1 )->x(),
                                           aElement->node(  8 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->x(),
                                           aElement->node(  3 )->x(),
                                           aElement->node( 11 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->x(),
                                           aElement->node(  4 )->x(),
                                           aElement->node( 12 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  1 )->x(),
                                           aElement->node(  2 )->x(),
                                           aElement->node(  9 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  1 )->x(),
                                           aElement->node(  5 )->x(),
                                           aElement->node( 13 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  2 )->x(),
                                           aElement->node(  3 )->x(),
                                           aElement->node( 10 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  2 )->x(),
                                           aElement->node(  6 )->x(),
                                           aElement->node( 14 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  3 )->x(),
                                           aElement->node(  7 )->x(),
                                           aElement->node( 15 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  4 )->x(),
                                           aElement->node(  5 )->x(),
                                           aElement->node( 16 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  4 )->x(),
                                           aElement->node(  7 )->x(),
                                           aElement->node( 19 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  5 )->x(),
                                           aElement->node(  6 )->x(),
                                           aElement->node( 17 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  6 )->x(),
                                           aElement->node(  7 )->x(),
                                           aElement->node( 18 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->y(),
                                           aElement->node(  1 )->y(),
                                           aElement->node(  8 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->y(),
                                           aElement->node(  3 )->y(),
                                           aElement->node( 11 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->y(),
                                           aElement->node(  4 )->y(),
                                           aElement->node( 12 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  1 )->y(),
                                           aElement->node(  2 )->y(),
                                           aElement->node(  9 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  1 )->y(),
                                           aElement->node(  5 )->y(),
                                           aElement->node( 13 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  2 )->y(),
                                           aElement->node(  3 )->y(),
                                           aElement->node( 10 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  2 )->y(),
                                           aElement->node(  6 )->y(),
                                           aElement->node( 14 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  3 )->y(),
                                           aElement->node(  7 )->y(),
                                           aElement->node( 15 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  4 )->y(),
                                           aElement->node(  5 )->y(),
                                           aElement->node( 16 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  4 )->y(),
                                           aElement->node(  7 )->y(),
                                           aElement->node( 19 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  5 )->y(),
                                           aElement->node(  6 )->y(),
                                           aElement->node( 17 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  6 )->y(),
                                           aElement->node(  7 )->y(),
                                           aElement->node( 18 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->y(),
                                           aElement->node(  1 )->y(),
                                           aElement->node(  8 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->y(),
                                           aElement->node(  3 )->y(),
                                           aElement->node( 11 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->y(),
                                           aElement->node(  4 )->y(),
                                           aElement->node( 12 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  1 )->y(),
                                           aElement->node(  2 )->y(),
                                           aElement->node(  9 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  1 )->y(),
                                           aElement->node(  5 )->y(),
                                           aElement->node( 13 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  2 )->y(),
                                           aElement->node(  3 )->y(),
                                           aElement->node( 10 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  2 )->y(),
                                           aElement->node(  6 )->y(),
                                           aElement->node( 14 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  3 )->y(),
                                           aElement->node(  7 )->y(),
                                           aElement->node( 15 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  4 )->y(),
                                           aElement->node(  5 )->y(),
                                           aElement->node( 16 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  4 )->y(),
                                           aElement->node(  7 )->y(),
                                           aElement->node( 19 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  5 )->y(),
                                           aElement->node(  6 )->y(),
                                           aElement->node( 17 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  6 )->y(),
                                           aElement->node(  7 )->y(),
                                           aElement->node( 18 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else
            {
                aElement->unset_curved_flag();
            }
        }

//------------------------------------------------------------------------------


    void
    CurvedElementChecker::check_hex27( Element * aElement )
    {
            if( this->check_midpoint( aElement->node(  0 )->x(),
                                           aElement->node(  1 )->x(),
                                           aElement->node(  8 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->x(),
                                           aElement->node(  3 )->x(),
                                           aElement->node( 11 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->x(),
                                           aElement->node(  4 )->x(),
                                           aElement->node( 12 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  1 )->x(),
                                           aElement->node(  2 )->x(),
                                           aElement->node(  9 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  1 )->x(),
                                           aElement->node(  5 )->x(),
                                           aElement->node( 13 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  2 )->x(),
                                           aElement->node(  3 )->x(),
                                           aElement->node( 10 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  2 )->x(),
                                           aElement->node(  6 )->x(),
                                           aElement->node( 14 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  3 )->x(),
                                           aElement->node(  7 )->x(),
                                           aElement->node( 15 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  4 )->x(),
                                           aElement->node(  5 )->x(),
                                           aElement->node( 16 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  4 )->x(),
                                           aElement->node(  7 )->x(),
                                           aElement->node( 19 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  5 )->x(),
                                           aElement->node(  6 )->x(),
                                           aElement->node( 17 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  6 )->x(),
                                           aElement->node(  7 )->x(),
                                           aElement->node( 18 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->x(),
                                           aElement->node(  1 )->x(),
                                           aElement->node(  4 )->x(),
                                           aElement->node(  5 )->x(),
                                           aElement->node( 25 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  1 )->x(),
                                           aElement->node(  2 )->x(),
                                           aElement->node(  5 )->x(),
                                           aElement->node(  6 )->x(),
                                           aElement->node( 24 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  2 )->x(),
                                           aElement->node(  3 )->x(),
                                           aElement->node(  6 )->x(),
                                           aElement->node(  7 )->x(),
                                           aElement->node( 26 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->x(),
                                           aElement->node(  3 )->x(),
                                           aElement->node(  4 )->x(),
                                           aElement->node(  7 )->x(),
                                           aElement->node( 23 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->x(),
                                           aElement->node(  1 )->x(),
                                           aElement->node(  2 )->x(),
                                           aElement->node(  3 )->x(),
                                           aElement->node( 21 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  4 )->x(),
                                           aElement->node(  5 )->x(),
                                           aElement->node(  6 )->x(),
                                           aElement->node(  7 )->x(),
                                           aElement->node( 22 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->x(),
                                           aElement->node(  1 )->x(),
                                           aElement->node(  2 )->x(),
                                           aElement->node(  3 )->x(),
                                           aElement->node(  4 )->x(),
                                           aElement->node(  5 )->x(),
                                           aElement->node(  6 )->x(),
                                           aElement->node(  7 )->x(),
                                           aElement->node( 20 )->x() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->y(),
                                           aElement->node(  1 )->y(),
                                           aElement->node(  8 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->y(),
                                           aElement->node(  3 )->y(),
                                           aElement->node( 11 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->y(),
                                           aElement->node(  4 )->y(),
                                           aElement->node( 12 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  1 )->y(),
                                           aElement->node(  2 )->y(),
                                           aElement->node(  9 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  1 )->y(),
                                           aElement->node(  5 )->y(),
                                           aElement->node( 13 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  2 )->y(),
                                           aElement->node(  3 )->y(),
                                           aElement->node( 10 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  2 )->y(),
                                           aElement->node(  6 )->y(),
                                           aElement->node( 14 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  3 )->y(),
                                           aElement->node(  7 )->y(),
                                           aElement->node( 15 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  4 )->y(),
                                           aElement->node(  5 )->y(),
                                           aElement->node( 16 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  4 )->y(),
                                           aElement->node(  7 )->y(),
                                           aElement->node( 19 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  5 )->y(),
                                           aElement->node(  6 )->y(),
                                           aElement->node( 17 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  6 )->y(),
                                           aElement->node(  7 )->y(),
                                           aElement->node( 18 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->y(),
                                           aElement->node(  1 )->y(),
                                           aElement->node(  4 )->y(),
                                           aElement->node(  5 )->y(),
                                           aElement->node( 25 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  1 )->y(),
                                           aElement->node(  2 )->y(),
                                           aElement->node(  5 )->y(),
                                           aElement->node(  6 )->y(),
                                           aElement->node( 24 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  2 )->y(),
                                           aElement->node(  3 )->y(),
                                           aElement->node(  6 )->y(),
                                           aElement->node(  7 )->y(),
                                           aElement->node( 26 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->y(),
                                           aElement->node(  3 )->y(),
                                           aElement->node(  4 )->y(),
                                           aElement->node(  7 )->y(),
                                           aElement->node( 23 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->y(),
                                           aElement->node(  1 )->y(),
                                           aElement->node(  2 )->y(),
                                           aElement->node(  3 )->y(),
                                           aElement->node( 21 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  4 )->y(),
                                           aElement->node(  5 )->y(),
                                           aElement->node(  6 )->y(),
                                           aElement->node(  7 )->y(),
                                           aElement->node( 22 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->y(),
                                           aElement->node(  1 )->y(),
                                           aElement->node(  2 )->y(),
                                           aElement->node(  3 )->y(),
                                           aElement->node(  4 )->y(),
                                           aElement->node(  5 )->y(),
                                           aElement->node(  6 )->y(),
                                           aElement->node(  7 )->y(),
                                           aElement->node( 20 )->y() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->z(),
                                           aElement->node(  1 )->z(),
                                           aElement->node(  8 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->z(),
                                           aElement->node(  3 )->z(),
                                           aElement->node( 11 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->z(),
                                           aElement->node(  4 )->z(),
                                           aElement->node( 12 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  1 )->z(),
                                           aElement->node(  2 )->z(),
                                           aElement->node(  9 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  1 )->z(),
                                           aElement->node(  5 )->z(),
                                           aElement->node( 13 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  2 )->z(),
                                           aElement->node(  3 )->z(),
                                           aElement->node( 10 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  2 )->z(),
                                           aElement->node(  6 )->z(),
                                           aElement->node( 14 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  3 )->z(),
                                           aElement->node(  7 )->z(),
                                           aElement->node( 15 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  4 )->z(),
                                           aElement->node(  5 )->z(),
                                           aElement->node( 16 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  4 )->z(),
                                           aElement->node(  7 )->z(),
                                           aElement->node( 19 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  5 )->z(),
                                           aElement->node(  6 )->z(),
                                           aElement->node( 17 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  6 )->z(),
                                           aElement->node(  7 )->z(),
                                           aElement->node( 18 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->z(),
                                           aElement->node(  1 )->z(),
                                           aElement->node(  4 )->z(),
                                           aElement->node(  5 )->z(),
                                           aElement->node( 25 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  1 )->z(),
                                           aElement->node(  2 )->z(),
                                           aElement->node(  5 )->z(),
                                           aElement->node(  6 )->z(),
                                           aElement->node( 24 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  2 )->z(),
                                           aElement->node(  3 )->z(),
                                           aElement->node(  6 )->z(),
                                           aElement->node(  7 )->z(),
                                           aElement->node( 26 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->z(),
                                           aElement->node(  3 )->z(),
                                           aElement->node(  4 )->z(),
                                           aElement->node(  7 )->z(),
                                           aElement->node( 23 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->z(),
                                           aElement->node(  1 )->z(),
                                           aElement->node(  2 )->z(),
                                           aElement->node(  3 )->z(),
                                           aElement->node( 21 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  4 )->z(),
                                           aElement->node(  5 )->z(),
                                           aElement->node(  6 )->z(),
                                           aElement->node(  7 )->z(),
                                           aElement->node( 22 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else if( this->check_midpoint( aElement->node(  0 )->z(),
                                           aElement->node(  1 )->z(),
                                           aElement->node(  2 )->z(),
                                           aElement->node(  3 )->z(),
                                           aElement->node(  4 )->z(),
                                           aElement->node(  5 )->z(),
                                           aElement->node(  6 )->z(),
                                           aElement->node(  7 )->z(),
                                           aElement->node( 20 )->z() ) )
            {
                aElement->set_curved_flag() ;
            }
            else
            {
                aElement->unset_curved_flag() ;
            }
    }

//------------------------------------------------------------------------------
    }
}