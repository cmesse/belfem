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

#ifndef CL_FEM_DOMAIN_HPP
#define CL_FEM_DOMAIN_HPP

#include "typedefs.hpp"
#include "cl_Input_Section.hpp"
#include "cl_Vector.hpp"
#include "en_DomainType.hpp"

namespace belfem
{
    namespace fem
    {
        class Domain
        {
            DomainType mType ;

            Vector< id_t > mGroupIDs ;

            // if the group ids are used as terminal to
            // a thin shell, the next string is not empty. Otherwise it is
            string mThinShellLabel ;

            string mLabel ;

            string mMaterialLabel ;

            bool mIsBlock = true ;
            bool mIsSideSet = false ;

        public:

            Domain( const input::Section * aSection );

            ~Domain() = default;

            bool
            is_block() const ;

            bool
            is_sideset() const ;

            DomainType
            type() const ;

            const Vector< id_t > &
            group_ids() const ;

            const string &
            material() const ;

        private:

            //! aAllowSigns: thin-shell sidesets may carry gmsh-style negative
            //! signs ( sidesets : -5, -6, 7:20 ; ) that request an orientation
            //! flip. The signs are CONSUMED elsewhere
            //! ( MaxwellFactory::read_signed_sidesets ); here they are only
            //! stripped so the ids resolve — without this, "-5" wraps through
            //! the unsigned parse to 4294967291 and the sideset lookup aborts
            void
            read_groups( const input::Section * aSection,
                         bool aUseBlocks  = true,
                         bool aAllowSigns = false );

        };

        inline DomainType
        Domain::type() const
        {
            return mType ;
        }

        inline bool
        Domain::is_block() const
        {
            return mIsBlock ;
        }

        inline bool
        Domain::is_sideset() const
        {
            return mIsSideSet ;
        }

        inline const Vector< id_t > &
        Domain::group_ids() const
        {
            return mGroupIDs ;
        }

        inline const string &
        Domain::material() const
        {
            return mMaterialLabel ;
        }


    }
}
#endif //CL_FEM_DOMAIN_HPP
