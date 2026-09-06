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

            inline size_t
            memory() const;

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

        inline size_t
        GlobalVariable::memory() const
        {
            return sizeof( GlobalVariable ) + mLabel.capacity() * sizeof( char );
        }

//------------------------------------------------------------------------------

    }
}
#endif //BELFEM_CL_MESH_GLOBALVARIABLE_HPP
