//
// Created by Christian Messe on 09.10.19.
//

#ifndef BELFEM_CL_MESH_GLOBALVARIABLE_HPP
#define BELFEM_CL_MESH_GLOBALVARIABLE_HPP

#include "typedefs.hpp"

namespace belfem
{
    namespace mesh
    {
//------------------------------------------------------------------------------

        class GlobalVariable
        {
            string     mLabel;
            const id_t mID;
            real       mValue;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            GlobalVariable(
                    const string & aLabel,
                    const id_t   & aID,
                    const real     aValue = 0.0 );

//------------------------------------------------------------------------------

            ~GlobalVariable() = default;

//------------------------------------------------------------------------------

            inline const string &
            label() const;

//------------------------------------------------------------------------------

            inline const id_t &
            id() const;

//------------------------------------------------------------------------------

            inline real &
            value();

//------------------------------------------------------------------------------

            inline const real &
            value() const;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline const string &
        GlobalVariable::label() const
        {
            return mLabel;
        }

//------------------------------------------------------------------------------

        inline const id_t &
        GlobalVariable::id() const
        {
            return mID;
        }

//------------------------------------------------------------------------------

        inline real &
        GlobalVariable::value()
        {
            return mValue;
        }

//------------------------------------------------------------------------------

        inline const real &
        GlobalVariable::value() const
        {
            return mValue;
        }

//------------------------------------------------------------------------------

    }
}
#endif //BELFEM_CL_MESH_GLOBALVARIABLE_HPP
