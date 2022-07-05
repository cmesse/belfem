//
// Created by Christian Messe on 08.10.19.
//

#include "cl_Mesh_Field.hpp"
#include "cl_Mesh.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        Field::Field(
                      Mesh            & aParent,
                const string          & aLabel,
                const index_t         & aIndex,
                const id_t            & aID,
                const EntityType        aEntityType,
                const FieldType         aFieldType ) :
                mParent( aParent ),
                mLabel( aLabel ),
                mIndex( aIndex ),
                mID( aID ),
                mEntityType( aEntityType ),
                mFieldType( aFieldType )
        {
            // todo: remove field type, it makes to sense
            BELFEM_ERROR( mFieldType == FieldType::SCALAR,
                    "only scalar fields are supported at this time" );

            switch( aEntityType )
            {
                case( EntityType::NODE ) :
                {
                    mData.set_size( aParent.number_of_nodes(), 0.0 );
                    break;
                }
                case( EntityType::ELEMENT ) :
                {
                    mData.set_size( aParent.number_of_elements(), 0.0 );
                    break;
                }
                case( EntityType::EDGE ) :
                case( EntityType::FACE ) :
                {
                    // size must be assigned later
                    break;
                }
                case( EntityType::FACET ) :
                {
                    mData.set_size( aParent.number_of_facets(), 0.0 );
                    break;
                }
                default :
                {
                    BELFEM_ERROR( false, "Unsupported field entity");
                    break;
                }
            }
        }

//------------------------------------------------------------------------------

        const string &
        Field::label() const
        {
            return mLabel;
        }

//------------------------------------------------------------------------------

        const EntityType &
        Field::entity_type()
        {
            return mEntityType;
        }

//------------------------------------------------------------------------------

        const FieldType &
        Field::field_type()
        {
            return mFieldType;
        }

//------------------------------------------------------------------------------

        Vector< real > &
        Field::data()
        {
            return mData;
        }

//------------------------------------------------------------------------------

        const id_t &
        Field::id() const
        {
            return mID;
        }

//------------------------------------------------------------------------------

        void
        Field::fill( const id_t aID, const real & aValue )
        {
            // unflag all elements on mesh
            mParent.unflag_all_elements() ;
            mParent.unflag_all_nodes() ;

            bool tFlag = false ;

            // check if this is a block ID
            if( mParent.block_exists( aID ) )
            {
                mParent.block( aID )->flag_nodes() ;
                tFlag = true ;
            }
            else if ( mParent.sideset_exists( aID ) )
            {
                mParent.sideset( aID )->flag_all_nodes() ;
                tFlag = true ;
            }

            // check if there is something to do
            if( tFlag )
            {
                // grab nodes on mesh
                Cell< Node * > & tNodes = mParent.nodes() ;

                // loop over all nodes and write value if node was flagged
                for( Node * tNode : tNodes )
                {
                    if( tNode->is_flagged() )
                    {
                        mData( tNode->index() ) = aValue ;
                    }
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Field::fill( const Vector< id_t > & aIDs, const real & aValue )
        {
            // unflag all elements on mesh
            mParent.unflag_all_elements() ;
            mParent.unflag_all_nodes() ;

            bool tFlag = false ;

            for( id_t tID : aIDs )
            {
                // check if this is a block ID
                if ( mParent.block_exists( tID ) )
                {
                    mParent.block( tID )->flag_nodes();
                    tFlag = true;
                }
                else if ( mParent.sideset_exists( tID ) )
                {
                    mParent.sideset( tID )->flag_all_nodes();
                    tFlag = true;
                }
            }

            // check if there is something to do
            if( tFlag )
            {
                // grab nodes on mesh
                Cell< Node * > & tNodes = mParent.nodes() ;

                // loop over all nodes and write value if node was flagged
                for( Node * tNode : tNodes )
                {
                    if( tNode->is_flagged() )
                    {
                        mData( tNode->index() ) = aValue ;
                    }
                }
            }
        }

//------------------------------------------------------------------------------
    }
}