//
// Created by christian on 7/6/21.
//

#ifndef BELFEM_CL_IWGFACTORY_HPP
#define BELFEM_CL_IWGFACTORY_HPP

#include "cl_IWG.hpp"
#include "en_IWGs.hpp"

namespace belfem
{
    class Mesh ;
    namespace fem
    {
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
