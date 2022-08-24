//
// Created by Christian Messe on 26.10.19.
//



#include "cl_FEM_Block.hpp"
#include "cl_FEM_Field.hpp"


#include "meshtools.hpp"
#include "assert.hpp"

#include "fn_IF_initialize_integration_points.hpp"
#include "fn_IF_initialize_shape_function.hpp"
#include "cl_FEM_DofManager.hpp"

namespace belfem
{
    namespace fem
    {
//------------------------------------------------------------------------------

        // creates an empty block
        Block::Block(  DofManagerBase * aParent ) :
                Group( aParent,
                       GroupType::BLOCK,
                       ElementType::EMPTY,
                       0, 0 )
        {

        }

//------------------------------------------------------------------------------

        Block::Block( DofManagerBase * aParent, mesh::Block * aBlock,  const index_t aNumberOfElements ) :
            Group( aParent, GroupType::BLOCK, aParent->enforce_linear_interpolation() ?
                                              mesh::linear_element_type( aBlock->element_type() ) :
                                              aBlock->element_type(),
                                              aBlock->id(), aNumberOfElements ),
            mBlock( aBlock )
        {
            mIntegrationData = new IntegrationData( aBlock->element_type() );

            this->set_integration_order( aParent->block_integration_order() );
            this->initialize_elements();
            this->create_element_map();
        }

//------------------------------------------------------------------------------

        Block::~Block()
        {

            this->delete_pointers();
        }

//------------------------------------------------------------------------------

        void
        Block::set_integration_order( const uint aOrder )
        {
            // set the number of nodes per element
            mNumberOfNodesPerElement = mesh::number_of_nodes( mElementType ) ;

            mIntegrationData->populate( aOrder,
                            mParent->integration_scheme() );

            // note: sidesets must also be changed if assume_isogeometry is false
            this->assume_isogeometry();

        }

//------------------------------------------------------------------------------

        ElementType
        Block::element_type() const
        {
            return mElementType;
        }

//------------------------------------------------------------------------------

        void
        Block::initialize_elements()
        {
            // allocate element container
            mElements.set_size( mNumberOfElements, nullptr );

            // initialize counter
            index_t tCount = 0;

            // get pointer to element container on block
            Cell< mesh::Element * > & tElements = mBlock->elements();

            // legacy support
            switch ( mParent->type() )
            {
                case( DofManagerType::OLD ) :
                {
                    Field * tParent = reinterpret_cast< Field * >( mParent );

                    // create elements
                    for( mesh::Element * tElement : tElements )
                    {
                        // test if I own this element
                        if( tElement->owner() == mMyRank )
                        {
                            // create the new element
                            mElements( tCount++ ) = new fem::Element(
                                    this, tParent, tElement );
                        }
                    }
                    break;
                }
                case( DofManagerType::NEW ) :
                {
                    DofManager * tParent = reinterpret_cast< DofManager * >( mParent );
                    // create elements
                    for( mesh::Element * tElement : tElements )
                    {
                        // test if I own this element
                        if( tElement->owner() == mMyRank )
                        {
                            // create the new element
                            mElements( tCount++ ) = new fem::Element(
                                    this, tParent, tElement );
                        }
                    }
                    break;
                }
                default:
                {
                    BELFEM_ERROR( false, "Invalid DOF Manager type");
                }
                case DofManagerType::UNDEFINED:
                    break;
            }

            BELFEM_ASSERT( tCount == mNumberOfElements,
                   "Something went wrong while creating finite elements" );
        }

//------------------------------------------------------------------------------
    }
}