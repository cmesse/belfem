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

#ifndef BELFEM_CL_IWGFACTORY_HPP
#define BELFEM_CL_IWGFACTORY_HPP

#include "cl_IWG.hpp"
#include "en_IWGs.hpp"

namespace belfem
{
    class Mesh ;
    namespace fem
    {
        /**
         * @brief Creates IWG instances by equation type.
         *
         * @ingroup grp_fem_iwg
         * @see @ref fem_iwg_iwg_usage_guide
         */
        class IwgFactory
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            Mesh * mMesh = nullptr ;
            const uint mNumberOfDimensions ;
            Vector< id_t > mAllBlockIDs ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            IwgFactory( Mesh & aMesh ) ;

            IwgFactory( Mesh * aMesh ) ;

//------------------------------------------------------------------------------

            virtual ~IwgFactory() = default ;

//------------------------------------------------------------------------------

            IWG *
            create_iwg( const IwgType aType, const ModelDimensionality aDimensionality=ModelDimensionality::TwoD ) const;

//------------------------------------------------------------------------------

            const Vector< id_t > &
            all_block_ids();

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            populate_block_ids();

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        inline const Vector< id_t > &
        IwgFactory::all_block_ids()
        {
            return mAllBlockIDs ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_IWGFACTORY_HPP
