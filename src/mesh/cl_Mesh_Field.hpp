//
// Created by Christian Messe on 08.10.19.
//

#ifndef BELFEM_CL_MESH_FIELD_HPP
#define BELFEM_CL_MESH_FIELD_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"

#include "Mesh_Enums.hpp"

namespace belfem
{
    class Mesh;

    namespace mesh
    {
        class Field
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

                   Mesh          & mParent;
            string mLabel;
            const index_t     mIndex ;
            const id_t        mID ;
            const EntityType  mEntityType;
            const FieldType   mFieldType;

            Vector< real > mData;

            // flag telling if field is written to output file
            bool mWriteToFile = true ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Field(       Mesh           & aParent,
                   const string         & aLabel,
                   const index_t        & aIndex,
                   const id_t           & aID,
                   const EntityType       aEntityType=EntityType::NODE,
                   const FieldType        aFieldType=FieldType::SCALAR );

//------------------------------------------------------------------------------

            virtual ~Field() = default;

//------------------------------------------------------------------------------

            const string &
            label() const;

//------------------------------------------------------------------------------

            const EntityType &
            entity_type();

//------------------------------------------------------------------------------

            const FieldType &
            field_type();

//------------------------------------------------------------------------------

            Vector< real > &
            data();

//------------------------------------------------------------------------------

            const id_t &
            id() const;

//------------------------------------------------------------------------------

            real &
            value( const index_t aIndex );

//------------------------------------------------------------------------------

            /**
             * returns the index of the field on the mesh
             * @return
             */
            const index_t &
            index() const ;

//------------------------------------------------------------------------------

            /**
             * reset a field based on a sideset or block id
             */
             void
             fill( const id_t aID, const real aValue );

//------------------------------------------------------------------------------

            /**
             * reset a field based on a sideset or block ids
             */
            void
            fill( const Vector< id_t > & aIDs, const real aValue );

//------------------------------------------------------------------------------

            bool
            write_field_to_file() const ;

//------------------------------------------------------------------------------

            void
            set_write_to_file_flag( const bool aFlag );

//------------------------------------------------------------------------------
        };


//------------------------------------------------------------------------------

        inline real &
        Field::value( const index_t aIndex )
        {
            return mData( aIndex );
        }

//------------------------------------------------------------------------------

        inline const index_t &
        Field::index() const
        {
            return mIndex ;
        }

//------------------------------------------------------------------------------

        inline bool
        Field::write_field_to_file() const
        {
            return mWriteToFile ;
        }

//------------------------------------------------------------------------------

        inline void
        Field::set_write_to_file_flag( const bool aFlag )
        {
            mWriteToFile = aFlag ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_MESH_FIELD_HPP
